using Images, Colors, FixedPointNumbers, Plots
using FFTW
using JuMP
using NLopt
using Ipopt

include("XASGOCore.jl")
include("DIC.jl")
include("DICPatternDistortion.jl")
include("XASGOSupport.jl")
plotly()

function GrainFCalcCC(RefI, PR, gR, Patterns, PatternInfo, Qps; ROIs = BullseyeROIs(48,.25), ROIsize = 25, FAlgorithm = " ", WarpROI = [])
  num_pats = size(PatternInfo, 1)
  m = minimum(size(RefI))
  if FAlgorithm == "remap"
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
    PD = PatternInfo[i,4:6]
    gD = euler_to_gmat(PatternInfo[i,1], PatternInfo[i,2], PatternInfo[i,3])
    DefI = prep_ebsp(Patterns[i])
    F_guess = Qps'*gD'*gR*Qps
    if FAlgorithm == "remap"
      F = F_NCorrRemap(ROIs, ROIsize, WarpROI, DefI, PD, gD, RefI_coeffs, pad_size, PR, gR, F_guess, ccfilt, windowfunc)
    elseif FAlgorithm == "ncorriterate"
      F = F_NCorrShiftIterate(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess, ccfilt, windowfunc)
    else
      F = F_NCorrShift(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess, ccfilt, windowfunc)
    end
    #display((F - F_guess)*1e6)
    display(norm((F - F_guess)*1e6))
    F_out[i,:,:] = F
  end
  return F_out
end

function GrainFCalcICGN(RefI, PR, gR, Patterns, PatternInfo, Qps; ROIinfo = [.5,.5,.325,.375])
  num_pats = size(PatternInfo, 1)
  m = minimum(size(RefI))
  ROI = AnnularROI(m*[ROIinfo[1],ROIinfo[2]], m*ROIinfo[3], m*ROIinfo[4], m)

  PR_p_px = PC_to_phosphor_frame(PR,m)
  DD = PR_p_px[3]
  P_ivec = phosphor_frame_to_image_vec(PR_p_px)
  pad_size = 50
  RefI_coeff = calc_B_coeffs(RefI, pad_size)
  df_ddp = GetSteepestDescentImagesPatternDistortion(RefI_coeff, ROI, pad_size, P_ivec, DD)
  f = EvaluateWarpedImagePatternDistortion(RefI_coeff, ROI, mattovec(eye(3)), pad_size, P_ivec, DD, [0,0], 0)
  f_m = mean(f)
  hess = ComputeHessian(f, f_m, df_ddp)

  F_out = zeros(num_pats,3,3)
  for i=1:num_pats
    PD = PatternInfo[i,4:6]
    gD = euler_to_gmat(PatternInfo[i,1], PatternInfo[i,2], PatternInfo[i,3])
    DefI = prep_ebsp(Patterns[i])

    F_guess = Qps'*gD'*gR*Qps


    PD_p_px = PC_to_phosphor_frame(PD,m)
    Δ_p_px =  PD_p_px - PR_p_px
    ΔDD = Δ_p_px[3]
    Δ_ivec = phosphor_frame_to_image_vec(Δ_p_px)
    Fgi = rotate_to_image_frame(F_guess)

    Fvec = RunIcgnPatternDistortion(DefI, RefI_coeff, df_ddp, pad_size, f, f_m, hess, ROI, mattovec(Fgi), P_ivec, DD, Δ_ivec, ΔDD)
    F = VR_deviatoric(vectomat(Fvec))
    F = rotate_to_phosframe_from_image(F)

    #F = F_ICGN(ROI, DefI, PD, gD, RefI_coeff, df_ddp, pad_size, f, f_m, hess, PR, gR, F_guess)
    #display((F - F_guess)*1e6)
    #display(norm((F - F_guess)*1e6))

    F_out[i,:,:] = F
  end
  return F_out
end

α = pi/2 - 7*pi/18 + pi/18
Qps=[0 cos(α) -sin(α);
        -1     0            0;
        0   -sin(α) -cos(α)]

RefI = prep_ebsp(joinpath("al_ebsd","set1","ebsd_0.png"))
PR = [.5;.5;.625]
gR = euler_to_gmat(pi/12,pi/18,0)
PatternInfo = [pi/12 pi/18 2*pi/180 .5 .5 .625;
              pi/12 pi/18 8*pi/180 .5 .5 .625;
              pi/12 pi/18 10*pi/180 .5 .5 .625;
              pi/12 pi/18 12*pi/180 .5 .5 .625;
              pi/12 pi/18 14*pi/180 .5 .5 .625;]
Patterns = String[]
for i=1:5
  push!(Patterns, joinpath("al_ebsd","set1",string("ebsd_",i,".png")))
end

ROIinfo = [.5,.5,.35,.375]
F = GrainFCalcICGN(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIinfo = ROIinfo)
#display("Done")

#=
N = 49
ROIsize = 25
rad = .25
ROIs = BullseyeROIs(N,rad)
WarpROI = AnnularROI([m,m]/2,0, m*(rad+(ROIsize/100)/sqrt(2)), m)
F = GrainFCalcCC(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIs = ROIs, ROIsize = ROIsize, FAlgorithm = "remap", WarpROI = WarpROI)
=#
