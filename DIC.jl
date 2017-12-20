#using Plots
include("QuinticBSpline.jl")


function RunIcgn(G, F, ref_c, initial_guess_u, ROIsize)
  # notation:
  # F: intensity of entire undeformed image
  # G: intensity of entire deformed image
  # F_coeff: biquintic spline coefficients of entire undeformed image
  # G_coeff: biquintic spline coefficients of entire deformed image
  # ROI: region of interest (nx2 array of x,y locations) 
  # ROIrelative: region of interest described relative to center of ROI
  # ref_c: [x, y] location of center of ROI in reference image
  # f: intensity of ROI in undeformed image
  # g: intensity of ROI in deformed image
  # p: vector of deformations to be solved for
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  G_coeff = calc_B_coeffs(G, pad_size)
  
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = GetSteepestDescentImages(F_coeff, ref_c, ROIrelative, pad_size)
  f = EvaluateWarpedImage(F_coeff, ref_c, ROIrelative, zeros(6), pad_size)
  f_m = mean(f)
  hess = ComputeHessian(f, f_m, df_ddp)
  p_old = TranslateInitialGuess(initial_guess_u) 
  converged = false
  num_iterations = 0
  while !converged && num_iterations < 20
    num_iterations += 1
    g = EvaluateWarpedImage(G_coeff, ref_c, ROIrelative, p_old, pad_size)
    g_m = mean(g) 
    grad = ComputeGradient(f, f_m, g, g_m, df_ddp)
    dp = (hess\(-grad))
    p_old = UpdateP(p_old, dp)
    converged = (norm(dp) < 1.0e-6)
  end

  return p_old[[1, 2]]
end


function TranslateInitialGuess(initial_guess)
  p = zeros(6)
  for i in 1 : min(6, length(initial_guess))
    p[i] = initial_guess[i]
  end
  return p
end


function SquareRoiEvenSize(ROIsize)
  ROIrelative = Array{Float64}(ROIsize*ROIsize, 2)
  count = 0
  for i in -ROIsize/2 : ROIsize/2 - 1  
    for j in -ROIsize/2 : ROIsize/2 - 1 
      count += 1
      ROIrelative[count, 1] = i+0.5
      ROIrelative[count, 2] = j+0.5
    end
  end
  return ROIrelative 
end


function GetSteepestDescentImages(F_coeff, ref_c, ROIrelative, pad_size)
  # this function assumes x_ref_tilde = x_ref
  # ddp = [ddu, ddv, ddu_dx, ddu_dy, ddv_dx, ddv_dy]
  df_ddp = zeros(size(ROIrelative, 1), 6)
  
  ROI = ROIrelative .+ ref_c' 
  #x = ROIrange + ref_c[1] + p[1] + p[3]*ROIrange + p[4]*ROIrange'
  #y = ROIrange' + ref_c[2] + p[2] + p[5]*ROIrange + p[6]*ROIrange'
  
  df_dxy = SplineDerivative(F_coeff, ROI, pad_size)
  df_ddp[:, 1] = copy(df_dxy[:, 1])
  df_ddp[:, 2] = copy(df_dxy[:, 2])
  df_ddp[:, 3] = df_dxy[:, 1] .* ROIrelative[:, 1]
  df_ddp[:, 4] = df_dxy[:, 1] .* ROIrelative[:, 2]
  df_ddp[:, 5] = df_dxy[:, 2] .* ROIrelative[:, 1]
  df_ddp[:, 6] = df_dxy[:, 2] .* ROIrelative[:, 2]
  
  return df_ddp
end


function ComputeHessian(f, f_m, df_ddp)
  hessian = zeros(6,6)
  for i in 1:length(f)
    hessian += df_ddp[i, :] * df_ddp[i, :]'
  end
  hessian *= 2.0 / sum((f-f_m).^2)
  return factorize(hessian)
end


function ComputeGradient(f, f_m, g, g_m, df_ddp)
  gradient = zeros(6)
  norm_f = sum((f-f_m).^2)^0.5 
  norm_g = sum((g-g_m).^2)^0.5 
  for i in 1:length(f)
    gradient += ((f[i]-f_m)/norm_f - (g[i]-g_m)/norm_g)* df_ddp[i, :]
  end
  gradient *= 2/norm_f 
  return gradient
end


function ComputeCorrelationCriteria(f, f_m, g, g_m)
  norm_f = sum((f-f_m).^2)^0.5 
  norm_g = sum((g-g_m).^2)^0.5 
  C = 0
  for i in 1:length(f)
    C += ((f[i]-f_m)/norm_f - (g[i]-g_m)/norm_g)^2
  end
  return C
end


function UpdateP(p_old, dp)
  w_old = [1+p_old[3]  p_old[4]  p_old[1]; 
           p_old[5]  1+p_old[6]  p_old[2]; 
           0  0  1]
  w_dp = [1+dp[3]  dp[4]  dp[1]; 
          dp[5]  1+dp[6]  dp[2];
          0  0  1]
  w_new = w_old*inv(w_dp)
  p_new = [w_new[1, 3], w_new[2, 3], w_new[1, 1] - 1, 
           w_new[1, 2], w_new[2, 1], w_new[2, 2] - 1]
           
  return p_new
end


function EvaluateWarpedImage(F, ref_c, ROIrelative, p, pad_size)
  ROI = ROIrelative .+ ref_c' .+ p[1:2]' + 
        ROIrelative .* p[[3, 6]]' +
        ROIrelative[:, [2, 1]] .* p[[4, 5]]'
  return SplineEvaluate(F, ROI, pad_size)
end


function TestGradient(G, F, ref_c, p_old, ROIsize)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  G_coeff = calc_B_coeffs(G, pad_size)
  
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = GetSteepestDescentImages(F_coeff, ref_c, ROIrelative, pad_size)
  f = EvaluateWarpedImage(F_coeff, ref_c, ROIrelative, zeros(6), pad_size)
  f_m = mean(f)
  g = EvaluateWarpedImage(G_coeff, ref_c, ROIrelative, p_old, pad_size)
  g_m = mean(g)
  grad = ComputeGradient(f, f_m, g, g_m, df_ddp)
end


function TestCorrelation(G, F, ref_c, p_old, ROIsize)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  G_coeff = calc_B_coeffs(G, pad_size)
  
  ROIrelative = SquareRoiEvenSize(ROIsize)
  f = EvaluateWarpedImage(F_coeff, ref_c, ROIrelative, zeros(6), pad_size)
  f_m = mean(f)
  g = EvaluateWarpedImage(G_coeff, ref_c, ROIrelative, p_old, pad_size)
  g_m = mean(g)
  C_ls = ComputeCorrelationCriteria(f, f_m, g, g_m)
end


function FiniteDiffGradient(G, F, ref_c, p_0, ROIsize, delta)
  grad = zeros(6)
  C_0 = TestCorrelation(G, F, ref_c, p_0, ROIsize)
  for i=1:6
    p_1 = copy(p_0)
    p_1[i] += delta
    C_1 = TestCorrelation(G, F, ref_c, p_1, ROIsize)
    grad[i] = (C_1 - C_0)/delta
  end
  return grad
end


function Testdfddp(F, ref_c, ROIsize)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = zeros(size(ROIrelative, 1), 6)
  ROI = ROIrelative .+ ref_c'
  df_ddp = GetSteepestDescentImages(F_coeff, ref_c, ROIrelative, pad_size)
end


function FiniteDiffdfddp(F, ref_c, ROIsize, delta)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = zeros(size(ROIrelative, 1), 6)
  p_0 = zeros(6)
  f_0 = EvaluateWarpedImage(F_coeff, ref_c, ROIrelative, p_0, pad_size)
  for i=1:6
    p_1 = copy(p_0)
    p_1[i] += delta
    f_1 = EvaluateWarpedImage(F_coeff, ref_c, ROIrelative, p_1, pad_size)
    df_ddp[:, i] = (f_1 - f_0)/delta
  end
  return df_ddp
end




function DICtest()

    plotly()

    FF = ones(50, 50)
    FF[11:16,11:16] = [1 1 1 1 1 1;
                      1 .9 .7 .7 .9 1;
                      1 .7 .2 .2 .7 1;
                      1 .7 .2 .2 .7 1;
                      1 .9 .7 .7 .9 1;
                      1 1 1 1 1 1]
    GG = ones(50, 50)
    tmp = EvaluateWarpedImage(calc_B_coeffs(FF, 50), [15.5, 15.5], 
                              SquareRoiEvenSize(12), 
                              [0.00,0.0,0.,0,0.5,0], 50)
    GG[10:21,10:21] = reshape(tmp,12, 12)'

    #display(heatmap(FF'))
    #display(heatmap(GG'))

    RunIcgn(FF, GG, [15.5, 15.5], [0.0, 0.0], 12)

    p=zeros(6)
    println("---Gradient calc:")
    gcalc = TestGradient(FF, GG, [15.5, 15.5], p, 12)
    display(gcalc)
    display(gcalc./gcalc[1])
    display(gcalc./gcalc[2])
    println("---Gradient finitediff:")
    gdiff = FiniteDiffGradient(FF, GG, [15.5, 15.5], p, 12, 1e-10)
    display(gdiff)
    display(gdiff./gdiff[1])
    display(gdiff./gdiff[2])


    #=
    Gcalc = zeros(10, 6)
    Gdiff = zeros(10, 6)
    for i = 1:10
      Gdiff[i, :] = FiniteDiffGradient(FF, GG, [15.5, 15.5], p, 12, 1.0*10.0^(-5-i))
      Gcalc[i, :] = TestGradient(FF, GG, [15.5, 15.5], p, 12)
    end
    display(plot(Gcalc, title="calc"))
    display(plot(Gdiff, title="diff"))
    =#

    #=
    println("Gradient calc:")
    display(heatmap(Testdfddp(F, [15.5, 15.5], 2)))
    println("Gradient finitediff:")
    display(heatmap(FiniteDiffdfddp(F, [15.5, 15.5], 2, 1e-8)))
    =#
end


