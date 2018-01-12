function RunIcgnPatternDistortion(G, F, ref_c, initial_guess_M, ROIsize,P_ivec, DD, Δ_ivec, ΔDD)
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
  df_ddp = GetSteepestDescentImagesPatternDistortion(F_coeff, ROIrelative .+ ref_c', pad_size, P_ivec, DD)
  f = EvaluateWarpedImagePatternDistortion(F_coeff, ROIrelative .+ ref_c', zeros(9), pad_size, P_ivec, DD, [0,0], 0)
  f_m = mean(f)
  hess = ComputeHessian(f, f_m, df_ddp)
  p_old = TranslateInitialGuess(initial_guess_M)
  converged = false
  num_iterations = 0
  while !converged && num_iterations < 20
    num_iterations += 1
    g = EvaluateWarpedImagePatternDistortion(G_coeff, ROIrelative .+ ref_c', p, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
    g_m = mean(g)
    grad = ComputeGradient(f, f_m, g, g_m, df_ddp)
    dp = (hess\(-grad))
    p_old = UpdatePPatternDistortion(p_old, dp)
    converged = (norm(dp) < 1.0e-6)
  end

  return p_old
end

function GetSteepestDescentImagesPatternDistortion(F_coeff, ROIabsolute, pad_size, P_ivec, DD)
  # this function assumes x_ref_tilde = x_ref
  # ddp = [ddβ_ij]
  df_ddp = zeros(size(ROIabsolute, 1), 9)
  x = [ROIabsolute[:,1]-P_ivec[1] ROIabsolute[:,2]-P_ivec[2]]

  df_dxy = SplineDerivative(F_coeff, ROIabsolute, pad_size)
  f1 = copy(df_dxy[:, 1])
  f2 = copy(df_dxy[:, 2])
  df_ddp[:, 1] = f1.*x[:,1]
  df_ddp[:, 2] = f1.*x[:,2]
  df_ddp[:, 3] = -DD*f1
  df_ddp[:, 4] = f2.*x[:,1]
  df_ddp[:, 5] = f2.*x[:,2]
  df_ddp[:, 6] = -DD*f2
  df_ddp[:, 7] = (f1.*x[:,1] + f2.*x[:,2]).*x[:,1]/DD
  df_ddp[:, 8] = (f1.*x[:,1] + f2.*x[:,2]).*x[:,2]/DD
  df_ddp[:, 9] = -f1.*x[:,1] - f2.*x[:,2]
  return df_ddp
end

function EvaluateWarpedImagePatternDistortion(F, ROIabsolute, p, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
  M = [p[1] p[2] p[3];p[4] p[5] p[6];p[7] p[8] p[9]]
  x = [ROIabsolute[:,1]-P_ivec[1],ROIabsolute[:,2]-P_ivec[2]]
  Mx = M*[x[1];x[2];-DD] #this is a mess, I need to just loop through ROI absolute
  ROI = P_ivec + D_ivec - Mx*(DD + ΔDD)/Mx[3]
  return SplineEvaluate(F, ROI, pad_size)
end

function UpdatePPatternDistortion(p_old, dp)
  M = vectomat(p_old)
  dM = vectomat(dp)
  p_new = mattovec(M*inv(dM))
end

function vectomat(v)
  M = [v[1] v[2] v[3];v[4] v[5] v[6];v[7] v[8] v[9]]
end

function mattovec(M)
  v = [M[1,1],M[1,2],M[1,3],M[2,1],M[2,2],M[2,3],M[3,1],M[3,2],M[3,3]]
end
