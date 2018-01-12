using Images, Colors, FixedPointNumbers, Plots
using FFTW
using JuMP
using NLopt
using Ipopt

include("DIC.jl")
include("DICPatternDistortion.jl")
include("XASGOSupport.jl")
plotly()

function CalcF(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)

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
  #display(i)
  ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
  #ROI_px = round.(m*ROIs[i,:])

  ROI_p_px = image_vec_to_phosphor_frame(ROI_px)

  r_p = ROI_p_px - PR_p_px

  shiftROI_px = offset_estimate_w_F(r_p,PD_p_px,PR_p_px,F_guess) + ROI_px

  trueq = shiftROI_px - ROI_px

  #thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)

  ppp = RunIcgn(DefI, RefI, ROI_px, trueq, ROIsize_px)
  thisq = ppp[1:2]

  qerror[:,i] = thisq - trueq

  if norm(thisq-trueq) > 4
    display(i)
    Plot_ROI(RefI, ROI_px, ROIsize_px)
    PlotCC(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
  end

  qstar_p = reverse_PC_offset(thisq, r_p, PD_p_px, PR_p_px)

  qs[:,i] = qstar_p/m #Why does this fail if I take out the /m?
  rs[:,i] = r_p/m
end

println("mean X error: ", mean(qerror[1,:]))
println("mean Y error: ", mean(qerror[2,:]))
println("std X error: ", std(qerror[1,:]))
println("std Y error: ", std(qerror[2,:]))

β = beta_calc_Ruggles(qs,rs)
β = VR_deviatoric(β+eye(3)) - eye(3)

#display(plot(qerror[1,:],qerror[2,:]))
#display(scatter3d(ROIs[:,1],ROIs[:,2],qerror[1,:]))
#display(scatter3d(ROIs[:,1],ROIs[:,2],qerror[2,:]))
return β

end

function CalcFPD(ROICenter, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)
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
  Fvec = RunIcgnPatternDistortion(DefI, RefI, ROICenter_px, mattovec(F_guess), ROIsize_px, P_ivec, DD, Δ_ivec, ΔDD)
  β = vectomat(Fvec) - eye(3)
end

#RefI = prep_ebsp("ZeroCamElevation_x0y0.png")
RefI = prep_ebsp(joinpath("al_ebsd","set1","ebsd_0.png"))
#PR = [.5;.5;.7]
PR = [.5;.5;.625]
#gR = eye(3)
gR = euler_to_gmat(pi/12,pi/18,0)

#DefI = prep_ebsp("ZeroCamElevation_x500y500.png")
DefI = prep_ebsp(joinpath("al_ebsd","set1","ebsd_2.png"))
#PD = [.48046875;.518353371499725;.706680080924329] #x500y500
#PD = [.5;.51835337;.70668008] #x0y500
#PD = [0.498046875000000;0.501835337149973;0.700668008092433]#x50y50
PD = [.5;.5;.625]
#gD = eye(3)
gD = euler_to_gmat(pi/12,pi/18,8*pi/180)

N = 25
θs = 2*pi*collect(0:N-2)/(N-1)
ROIsize = 25
rad = .25
ROIs = Array{Float64}(N,2)
ROIs[1,:] = [.5 .5]
for i=2:N
  ROIs[i,:] = [.5 .5] + rad*[sin(θs[i-1]) cos(θs[i-1])]
end

F_guess = Qps'*gD'*gR*Qps - eye(3)

α = pi/2 - 7*pi/18 + pi/18
Qps=[0 cos(α) -sin(α);
        -1     0            0;
        0   -sin(α) -cos(α)]

#F = CalcF(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess+eye(3))

F = CalcFPD([.5,.5], ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess)

display(F)

display(F_guess)

display(1e6*(F_guess - F))

display(100*(F_guess - F)./F_guess)

display(norm(F_guess - F))

println("Done")
