using Images, Colors, FixedPointNumbers, Plots
using FFTW
using JuMP
using NLopt
using Ipopt

include("DIC.jl")
include("DICPatternDistortion.jl")
include("XASGOSupport.jl")
plotly()

function CalcFDIC(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)

m = minimum(size(DefI))  # smallest dimension of image
ROIsize_px = 2*round(Int64,m*ROIsize/200) #ROI forced to be even and square
ccfilt = ccfilter(2,50, [ROIsize_px,ROIsize_px], 13) #filter masks computed only once
windowfunc = ccwindow(ROIsize_px)
qerror = zeros(2,size(ROIs)[1])
qs = zeros(3,size(ROIs)[1])
rs = zeros(3,size(ROIs)[1])

PD_p_px = PC_to_phosphor_frame(PD,m)
PR_p_px = PC_to_phosphor_frame(PR,m)
for i in 1:size(ROIs)[1]
  ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
  ROI_p_px = image_vec_to_phosphor_frame(ROI_px)
  r_p = ROI_p_px - PR_p_px

  shiftROI_px = offset_estimate_w_F(r_p,PD_p_px,PR_p_px,F_guess) + ROI_px

  trueq = shiftROI_px - ROI_px

  #thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)

  ppp = RunIcgn(DefI, RefI, ROI_px, trueq, ROIsize_px)
  thisq = ppp[1:2]

  qerror[:,i] = thisq - trueq
  #=
  if norm(thisq-trueq) > 4
    display(i)
    Plot_ROI(RefI, ROI_px, ROIsize_px)
    PlotCC(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
  end
  =#

  qstar_p = reverse_PC_offset(thisq, r_p, PD_p_px, PR_p_px)

  qs[:,i] = qstar_p/m #Why does this fail if I take out the /m?
  rs[:,i] = r_p/m
end

#=
println("mean X error: ", mean(qerror[1,:]))
println("mean Y error: ", mean(qerror[2,:]))
println("std X error: ", std(qerror[1,:]))
println("std Y error: ", std(qerror[2,:]))
=#

β = beta_calc_Ruggles(qs,rs)
β = VR_deviatoric(β+eye(3)) - eye(3)

#display(plot(qerror[1,:],qerror[2,:]))
#display(scatter3d(ROIs[:,1],ROIs[:,2],qerror[1,:]))
#display(scatter3d(ROIs[:,1],ROIs[:,2],qerror[2,:]))
return β + eye(3)

end

function CalcFNCorr(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)

m = minimum(size(DefI))  # smallest dimension of image
ROIsize_px = 2*round(Int64,m*ROIsize/200) #ROI forced to be even and square
ccfilt = ccfilter(2,50, [ROIsize_px,ROIsize_px], 13) #filter masks computed only once
windowfunc = ccwindow(ROIsize_px)
qerror = zeros(2,size(ROIs)[1])
qs = zeros(3,size(ROIs)[1])
rs = zeros(3,size(ROIs)[1])

PD_p_px = PC_to_phosphor_frame(PD,m)
PR_p_px = PC_to_phosphor_frame(PR,m)
for i in 1:size(ROIs)[1]
  ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
  ROI_p_px = image_vec_to_phosphor_frame(ROI_px)
  r_p = ROI_p_px - PR_p_px

  shiftROI_px = offset_estimate_w_F(r_p,PD_p_px,PR_p_px,F_guess) + ROI_px

  trueq = shiftROI_px - ROI_px

  thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)

  qerror[:,i] = thisq - trueq

  #=
  if norm(thisq-trueq) > 4
    display(i)
    Plot_ROI(RefI, ROI_px, ROIsize_px)
    PlotCC(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
  end
  =#

  qstar_p = reverse_PC_offset(thisq, r_p, PD_p_px, PR_p_px)

  qs[:,i] = qstar_p/m #Why does this fail if I take out the /m?
  rs[:,i] = r_p/m
end
β = beta_calc_Ruggles(qs,rs, false)
β = VR_deviatoric(β+eye(3)) - eye(3)

return β + eye(3)

end

function CalcFPD(ROI, DefI, PD, gD, RefI, PR, gR, F_guess)
  m = minimum(size(DefI))

  PD_p_px = PC_to_phosphor_frame(PD,m)
  PR_p_px = PC_to_phosphor_frame(PR,m)
  Δ_p_px =  PD_p_px - PR_p_px
  DD = PR_p_px[3]
  P_ivec = phosphor_frame_to_image_vec(PR_p_px)
  ΔDD = Δ_p_px[3]
  Δ_ivec = phosphor_frame_to_image_vec(Δ_p_px)

  pad_size = 50
  RefI_coeff = calc_B_coeffs(RefI, pad_size)
  df_ddp = GetSteepestDescentImagesPatternDistortion(RefI_coeff, ROI, pad_size, P_ivec, DD)
  f = EvaluateWarpedImagePatternDistortion(RefI_coeff, ROI, mattovec(eye(3)), pad_size, P_ivec, DD, [0,0], 0)
  f_m = mean(f)
  hess = ComputeHessian(f, f_m, df_ddp)

  Fgi = rotate_to_image_frame(F_guess)

  Fvec = RunIcgnPatternDistortion(DefI, RefI_coeff, df_ddp, pad_size, f, f_m, hess, ROI, mattovec(Fgi), P_ivec, DD, Δ_ivec, ΔDD)
  F = VR_deviatoric(vectomat(Fvec))
  F = rotate_to_phosframe_from_image(F)
end

function CalcFNCorrIterate(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)

m = minimum(size(DefI))  # smallest dimension of image
ROIsize_px = 2*round(Int64,m*ROIsize/200) #ROI forced to be even and square
ccfilt = ccfilter(2,50, [ROIsize_px,ROIsize_px], 13) #filter masks computed only once
windowfunc = ccwindow(ROIsize_px)
qerror = zeros(2,size(ROIs)[1])
qs = zeros(3,size(ROIs)[1])
rs = zeros(3,size(ROIs)[1])

PD_p_px = PC_to_phosphor_frame(PD,m)
PR_p_px = PC_to_phosphor_frame(PR,m)

converged = false
num_iterations = 0
F_new = eye(3)
while !converged && num_iterations < 20
  num_iterations += 1
  display(num_iterations)
  for i in 1:size(ROIs)[1]
    ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
    ROI_p_px = image_vec_to_phosphor_frame(ROI_px)
    r_p = ROI_p_px - PR_p_px

    shiftROI_px = offset_estimate_w_F(r_p,PD_p_px,PR_p_px,F_guess) + ROI_px

    trueq = shiftROI_px - ROI_px

    thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)

    qerror[:,i] = thisq - trueq

    #=
    if norm(thisq-trueq) > 4
      display(i)
      Plot_ROI(RefI, ROI_px, ROIsize_px)
      PlotCC(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
    end
    =#

    qstar_p = reverse_PC_offset(thisq, r_p, PD_p_px, PR_p_px)

    qs[:,i] = qstar_p/m #Why does this fail if I take out the /m?
    rs[:,i] = r_p/m
  end

  β = beta_calc_Ruggles(qs,rs, false)
  β = VR_deviatoric(β+eye(3)) - eye(3)

  F_new = β + eye(3)
  if norm(F_new-F_guess)<100e-6
    converged = true
  else
    F_guess = F_new
  end
end

return F_new

end

function CalcFNCorrRemap(ROIs, ROIsize, WarpROI, DefI, PD, gD, RefI, PR, gR, F_guess)

tic()
pad_size = 50
m = minimum(size(DefI))
RefI_coeffs = calc_B_coeffs(RefI, pad_size)


  # smallest dimension of image
ROIsize_px = 2*round(Int64,m*ROIsize/200) #ROI forced to be even and square
ccfilt = ccfilter(2,50, [ROIsize_px,ROIsize_px], 13) #filter masks computed only once
windowfunc = ccwindow(ROIsize_px)
qs = zeros(3,size(ROIs)[1])
rs = zeros(3,size(ROIs)[1])

PD_p_px = PC_to_phosphor_frame(PD,m)
PR_p_px = PC_to_phosphor_frame(PR,m)

#display(heatmap(RefI))
#display(heatmap(DefI))

converged = false
num_iterations = 0
F_guess = VR_deviatoric(F_guess)
F_new = eye(3)
toc()
while !converged && num_iterations < 20
  tic()
  num_iterations += 1
  display(num_iterations)
  RefI = patternremap(RefI_coeffs, pad_size, WarpROI, PR_p_px, PD_p_px-PR_p_px, F_guess, m)
  #display(heatmap(RefI))
  for i in 1:size(ROIs)[1]
    ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
    ROI_p_px = image_vec_to_phosphor_frame(ROI_px)
    r_p = ROI_p_px - PR_p_px

    thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, ROI_px, ROIsize_px, windowfunc, ccfilt)

    qs[:,i] = image_vec_to_phosphor_frame(thisq)/m
    rs[:,i] = r_p/m
  end

  β = beta_calc_Ruggles(qs,rs, false)
  β = VR_deviatoric(β+eye(3)) - eye(3)

  #display(β)

  F_new = (β + eye(3))*F_guess
  if norm(F_new-F_guess)<10e-6
    converged = true
  else
    F_guess = F_new
  end
  toc()
end

return F_new

end

function GetCC(ROICenter, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)
  m = minimum(size(DefI))
  ROIsize_px = 2*round(Int64,m*ROIsize/200)
  ROICenter_px = round.(m*ROICenter + [.5,.5]) -  [.5,.5]
  PD_p_px = PC_to_phosphor_frame(PD,m)
  PR_p_px = PC_to_phosphor_frame(PR,m)
  Δ_p_px =  PD_p_px - PR_p_px
  DD = PR_p_px[3]
  P_ivec = phosphor_frame_to_image_vec(PR_p_px)
  ΔDD = Δ_p_px[3]
  Δ_ivec = phosphor_frame_to_image_vec(Δ_p_px)
  CC = TestCCforPD(DefI, RefI, ROICenter_px, mattovec(rotate_to_image_frame(F_guess)), ROIsize_px, P_ivec, DD, Δ_ivec, ΔDD)
end

#=
RefI = prep_ebsp("ZeroCamElevation_x0y0.png")
PR = [.5;.5;.7]
gR = eye(3)

DefI = prep_ebsp("ZeroCamElevation_x500y500.png")
PD = [.48046875;.518353371499725;.706680080924329] #x500y500
#PD = [.5;.51835337;.70668008] #x0y500
#PD = [0.498046875000000;0.501835337149973;0.700668008092433]#x50y50
gD = eye(3)
=#

RefI = prep_ebsp(joinpath("al_ebsd","set1","ebsd_0.png"))
PR = [.5;.5;.625]
gR = euler_to_gmat(pi/12,pi/18,0)

DefI = prep_ebsp(joinpath("al_ebsd","set1","ebsd_2.png"))
PD = [.5;.5;.625]
gD = euler_to_gmat(pi/12,pi/18,8*pi/180)

m = minimum(size(DefI))

N = 48
ROIsize = 25
rad = .25

#ROIs = BullseyeROIs(N,rad)
ROIs = AnnularROIs(N,rad)
#ROIs = GridROIs(N,ROIsize,.1)

α = pi/2 - 7*pi/18 + pi/18
Qps=[0 cos(α) -sin(α);
        -1     0            0;
        0   -sin(α) -cos(α)]

F_exact = Qps'*gD'*gR*Qps

NN = rand(3,3) - rand(3,3)
#=
NN = [0.134621  0.86776   0.931237;
0.187134  0.555074  0.767708;
0.870081  0.973473  0.0432369]
=#
F_guess = F_exact #+ .01*NN/norm(NN)

#=
tic()
F_N = CalcFNCorrIterate(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)
toc()
display("Difference (iNcorr)")
display((F_N - F_exact)*1e6)
display(norm(F_N - F_exact)*1e6)
=#



F_N = CalcFNCorr(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)
display("Difference (Ncorr)")
display((F_N - F_exact)*1e6)
display(norm(F_N - F_exact)*1e6)



WarpROI = AnnularROI([m,m]/2,m*(rad-(ROIsize/100)/sqrt(2)), m*(rad+(ROIsize/100)/sqrt(2)), m) #m*(rad-(ROIsize/100)/sqrt(2))
#WarpROI = ROIUnion(ROIs,ROIsize,m)

F_N = CalcFNCorrRemap(ROIs, ROIsize, WarpROI, DefI, PD, gD, RefI, PR, gR, F_guess)
display("Difference (Ncorr w Remapping)")
display((F_N - F_exact)*1e6)
display(norm(F_N - F_exact)*1e6)







#=
tic()
F_DIC = CalcFDIC(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)
toc()
display("Difference (DIC)")
display((F_DIC - F_guess)*1e6)
=#

#=
ROIsize_px = 2*round(Int64,m*70/200)
ROICenter_px = round.(m*[.5, .5] + [.5,.5]) -  [.5,.5]
ROI = SquareRoiEvenSize(ROIsize_px).+ ROICenter_px'
=#

#=
ROI = AnnularROI(m*[.5,.5], m*.35, m*.375, m)#ROI = AnnularROI(m*[.5,.5], m*.125, m*.375, m)
#ROI = LymeROI(m*[.5,.5], m*.1, m*.125, m*.375)
#ROI = SquareROI(m*[.5,.5],.63*m)

F_PD = CalcFPD(ROI, DefI, PD, gD, RefI, PR, gR, F_guess)
display("Difference (PDDIC)")
display((F_PD - F_exact)*1e6)
display(norm(F_PD - F_exact)*1e6)
=#





#=
S = VR_poldec(F_N)
display((S[1] - eye(3))*1e6)
display(norm((S[1] - eye(3))*1e6))
display((S[2] - F_PD)*1e6)
=#
