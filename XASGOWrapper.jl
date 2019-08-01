using Images, Colors, FixedPointNumbers, Plots
using FFTW
#using JuMP
#using NLopt
#using Ipopt

include("XASGOCore.jl")
include("QuinticBSpline.jl")
#include("CubicBSpline.jl")
include("DIC.jl")
include("DICPatternDistortion.jl")
include("XASGOSupport.jl")
plotly()

function GrainFCalcCC(RefI, PR, gR, Patterns, PatternInfo, Qps; ROIs = BullseyeROIs(48,.25), ROIsize = 25, FAlgorithm = "ncorrshift", WarpROI = [], FG = [], numimax = 10)
  num_pats = size(Patterns, 1)
  m = minimum(size(RefI))
  if FAlgorithm[1:5] == "remap"
    pad_size = 50
    RefI_coeffs = calc_B_coeffs(RefI, pad_size)
    if WarpROI == []
      WarpROI = ROIUnion(ROIs,ROIsize,m)
    end
  end
  ROIsize_px = 2*round(Int64,m*ROIsize/200)
  ccfilt = ccfilter(2,50, [ROIsize_px,ROIsize_px], 13) #filter masks computed only once
  windowfunc = ccwindow(ROIsize_px)
  F_out = zeros(num_pats,3,3)

  for i=1:num_pats
    println("Pattern #: ", i)
    PD = PatternInfo[i,4:6]
    gD = euler_to_gmat([PatternInfo[i,1]], [PatternInfo[i,2]], [PatternInfo[i,3]])
    DefI = prep_ebsp(Patterns[i])
    if FG == []
        F_guess = Qps'*gD'*gR*Qps
    else
        F_guess = FG[i,:,:]
    end

    if FAlgorithm == "remap"
      F = F_NCorrRemap(ROIs, ROIsize, WarpROI, DefI, PD, gD, RefI_coeffs, pad_size, PR, gR, F_guess, ccfilt, windowfunc, numimax = numimax)
    elseif FAlgorithm == "ncorriterate"
      F = F_NCorrShiftIterate(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess, ccfilt, windowfunc, numimax = numimax)
    elseif FAlgorithm == "remapSingle"
      F_guess_new = F_NCorrShift(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess, ccfilt, windowfunc)
      F = F_NCorrRemap(ROIs, ROIsize, WarpROI, DefI, PD, gD, RefI_coeffs, pad_size, PR, gR, F_guess_new, ccfilt, windowfunc, numimax = 1)
    elseif FAlgorithm == "ncorrnoshift"
      F = F_NCorrShift(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, eye(3), ccfilt, windowfunc)
    else
        if angofF(F_guess)>12
            pa = angofF(F_guess)
            display(string("Problem angle: ",pa))
            F_guess = VR_poldec((F_guess - eye(3))*12/pa + eye(3))[2]
        end
      if FAlgorithm != "ncorrshift"
        display("Algorithm not recognized, using NCorrShift")
      end
      F = F_NCorrShift(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess, ccfilt, windowfunc)
    end
    #display(F)
    #display(F_guess)
    #display(((F - F_guess)*1e6))
    #display(norm((F - F_guess)*1e6))
    F_out[i,:,:] = F
  end
  return F_out
end

function GrainFCalcICGN(RefI, PR, gR, Patterns, PatternInfo, Qps; ROIinfo = [.5,.5,.325,.375], numimax = 10, FG = [], Hupdate = false)
  num_pats = size(Patterns, 1)
  m = minimum(size(RefI))
  if typeof(ROIinfo)==Array{Int,2}
    ROI = ROIinfo
    display("Using custom ROI")
  else
    ROI = AnnularROI(m*[ROIinfo[1],ROIinfo[2]], m*ROIinfo[3], m*ROIinfo[4], m)
  end

  PR_p_px = PC_to_phosphor_frame(PR,m)
  DD = PR_p_px[3]
  P_ivec = phosphor_frame_to_image_vec(PR_p_px)
  pad_size = 50

  RefI_coeff = calc_B_coeffs(RefI, pad_size)

  df_ddp = GetSteepestDescentImagesPatternDistortion(RefI_coeff, ROI, pad_size, P_ivec, DD)

  f = EvaluateWarpedImagePatternDistortion(RefI_coeff, ROI, mattovec(I), pad_size, P_ivec, DD, [0,0], 0)

  f_m = mean(f)
  hess = ComputeHessian(f, f_m, df_ddp)

  F_out = zeros(num_pats,3,3)

  for i=1:num_pats
    println("Pattern #: ", i)

    PD = PatternInfo[i,4:6]
    gD = euler_to_gmat([PatternInfo[i,1]], [PatternInfo[i,2]], [PatternInfo[i,3]])
    DefI = prep_ebsp(Patterns[i])

    if FG == []
        F_guess = Qps'*gD'*gR*Qps
    else
        F_guess = FG[i,:,:]
    end

    PD_p_px = PC_to_phosphor_frame(PD,m)
    Δ_p_px =  PD_p_px - PR_p_px
    ΔDD = Δ_p_px[3]
    Δ_ivec = phosphor_frame_to_image_vec(Δ_p_px)
    Fgi = rotate_to_image_frame(F_guess)


    #display(heatmap(DefI))

    Fvec = RunIcgnPatternDistortion(DefI, RefI_coeff, df_ddp, pad_size, f, f_m, hess, ROI, mattovec(Fgi), P_ivec, DD, Δ_ivec, ΔDD, numimax = numimax, Hupdate = Hupdate)
    F = VR_deviatoric(vectomat(Fvec))
    F = rotate_to_phosframe_from_image(F)

    F_out[i,:,:] = F
  end
  return F_out
end

#=
α = pi/2 - 7*pi/18 + pi/18
Qps=[0 cos(α) -sin(α);
        -1     0            0;
        0   -sin(α) -cos(α)]

RefI = prep_ebsp(joinpath("al_ebsd","set1","ebsd_0.png"))
PR = [.5;.5;.625]
gR = euler_to_gmat(pi/12,pi/18,0)
m = minimum(size(RefI))

#dx = -2/960
#dy = -1.8/960
#dz = .6/960

dx = -20/960
dy = -18/960
dz = 6/960
PatternInfo = [pi/12 pi/18 0 .5+dx .5+dy .625+dz;
              pi/12 pi/18 2*pi/180 .5+dx .5+dy .625+dz;
              pi/12 pi/18 8*pi/180 .5+dx .5+dy .625+dz;
              pi/12 pi/18 10*pi/180 .5+dx .5+dy .625+dz;
              pi/12 pi/18 12*pi/180 .5+dx .5+dy .625+dz;
              pi/12 pi/18 14*pi/180 .5+dx .5+dy .625+dz;]

Patterns = String[]
for i=0:5
  push!(Patterns, joinpath("al_ebsd","set3",string("ebsd_",i,".png")))
end

#display(heatmap(RefI))

ROIinfo = [.5,.5,0,m*(rad+(ROIsize/100)/sqrt(2))j]
tic()
F = GrainFCalcICGN(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIinfo = ROIinfo)
toc()
#display("Done")

#=
N = 200
ROIsize = 25
rad = .25
ROIs = BullseyeROIs(N,rad)
WarpROI = AnnularROI([m,m]/2,0, m*(rad+(ROIsize/100)/sqrt(2)), m)
tic()
F = GrainFCalcCC(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIs = ROIs, ROIsize = ROIsize, FAlgorithm = "remap", WarpROI = WarpROI)
toc()
=#

#=
xmap = zeros(100,100)
ymap = zeros(100,100)
for i=1:100
  for j=1:100
    ymap[i,j] = j
    xmap[i,j] = i
  end
end
display(heatmap(xmap))
display(heatmap(ymap))
=#

=#
