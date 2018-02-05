function RunIcgnPatternDistortion(G, F_coeff, df_ddp, pad_size, f, f_m, hess, ROI, initial_guess_M, P_ivec, DD, Δ_ivec, ΔDD)
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

  p_old = initial_guess_M#initial_guess_M
  G_coeff = calc_B_coeffs(G, pad_size)

  converged = false
  num_iterations = 0
  while !converged && num_iterations < 20
    num_iterations += 1
    display(num_iterations)
    g = EvaluateWarpedImagePatternDistortion(G_coeff, ROI, p_old, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
    g_m = mean(g)
    grad = ComputeGradient(f, f_m, g, g_m, df_ddp)
    dp = (hess\(-grad))
    p_new = UpdatePPatternDistortion(p_old, dp)
    converged = (norm(p_new - p_old) < 10.0e-6)
    p_old = p_new
  end
  display("fin")
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
  M = vectomat(p)
  xs = [ROIabsolute[:,1]-P_ivec[1] ROIabsolute[:,2]-P_ivec[2]]
  ROI = zeros(size(ROIabsolute))
  for i=1:size(ROIabsolute,1)
    x = xs[i,:]
    Mr = M*[x[1];x[2];-DD]
    ROInew = P_ivec + Δ_ivec - Mr[1:2]*(DD + ΔDD)/Mr[3]
    ROI[i,:] = ROInew
  end
  return SplineEvaluate(F, ROI, pad_size)
end

function UpdatePPatternDistortion(p_old, dp)
  M = vectomat(p_old)
  dM = vectomat(dp)+eye(3)
  p_new = mattovec(VR_deviatoric(M*inv(dM)))
end

function vectomat(v)
  M = [v[1] v[2] v[3];v[4] v[5] v[6];v[7] v[8] v[9]]
  #M = [v[1] v[2] v[3];v[4] v[5] v[6];v[7] v[8] 1]
end

function mattovec(M)
  v = [M[1,1],M[1,2],M[1,3],M[2,1],M[2,2],M[2,3],M[3,1],M[3,2],M[3,3]]
  #v = [M[1,1],M[1,2],M[1,3],M[2,1],M[2,2],M[2,3],M[3,1],M[3,2]]
end



function DICtestPD()

    plotly()

    FF = ones(50, 50)
    FF[11:16,11:16] = [1 1 1 1 1 1;
                      1 .9 .7 .7 .9 1;
                      1 .7 .2 .2 .7 1;
                      1 .7 .2 .2 .7 1;
                      1 .9 .7 .7 .9 1;
                      1 1 1 1 1 1]
    GG = ones(50, 50)
    tmp = EvaluateWarpedImagePatternDistortion(calc_B_coeffs(FF, 50), SquareRoiEvenSize(12) .+ [15.5, 15.5]',
                              [1.,0.0,0.,0,1.,0.,0.,0.,1.], 50, [25.3,25.8], 30, [3,3], 3)
    GG[10:21,10:21] = reshape(tmp,12, 12)'

    p = mattovec(eye(3))
    println("---Gradient calc:")
    gcalc = TestGradientPatternDistortion(FF, GG, [15.5, 15.5], p, 12, [25.3,25.8], 30, [0,0], 0)
    display(gcalc)
    println("---Gradient finitediff:")
    #gdiff = FiniteDiffGradientPatternDistortion(FF, GG, [15.5, 15.5], p, 12, 1e-10, [25.3,25.8], 30, [3,3], 3)
    #display(gdiff)
    ds = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9]
    display(size(ds))
    vec1=zeros(size(ds,1),9)
    for i=1:size(ds,1)
      delta = ds[i]
      gdiff = FiniteDiffGradientPatternDistortion(FF, GG, [15.5, 15.5], p, 12, delta, [25.3,25.8], 30, [0,0], 0)
      display(gdiff)
      vec1[i,:] = gdiff
    end
    display(scatter(log10(ds),log10(vec1)))

    #display(heatmap(TestdfddpPatternDistortion(F, [15.5, 15.5], 2,[25.3,25.8], 30)))
    #display(heatmap(FiniteDiffdfddpPatternDistortion(F, [15.5, 15.5], 2, 1e-8, [25.3,25.8], 30)))

end

function TestdfddpPatternDistortion(F, ref_c, ROIsize,  P_ivec, DD)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = zeros(size(ROIrelative, 1), 9)
  ROI = ROIrelative .+ ref_c'
  df_ddp = GetSteepestDescentImagesPatternDistortion(F_coeff, ROI, pad_size, P_ivec, DD)
end

function FiniteDiffdfddpPatternDistortion(F, ref_c, ROIsize, delta, P_ivec, DD)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = zeros(size(ROIrelative, 1), 9)
  p_0 = mattovec(eye(3))
  f_0 = EvaluateWarpedImagePatternDistortion(F_coeff, ref_c' .+ ROIrelative, p_0, pad_size, P_ivec, DD, [0,0],0)
  for i=1:9
    p_1 = copy(p_0)
    p_1[i] += delta
    f_1 = EvaluateWarpedImagePatternDistortion(F_coeff, ref_c' .+ ROIrelative, p_1, pad_size, P_ivec, DD, [0,0],0)
    df_ddp[:, i] = (f_1 - f_0)/delta
  end
  return df_ddp
end

function TestGradientPatternDistortion(G, F, ref_c, p_old, ROIsize, P_ivec, DD, Δ_ivec, ΔDD)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  G_coeff = calc_B_coeffs(G, pad_size)

  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = GetSteepestDescentImagesPatternDistortion(F_coeff, ref_c' .+ ROIrelative, pad_size, P_ivec, DD)
  f = EvaluateWarpedImagePatternDistortion(F_coeff, ref_c' .+ ROIrelative, mattovec(eye(3)), pad_size, P_ivec, DD, [0.0,0.0], 0)
  f_m = mean(f)
  g = EvaluateWarpedImagePatternDistortion(G_coeff, ref_c' .+ ROIrelative, p_old, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
  g_m = mean(g)
  grad = ComputeGradient(f, f_m, g, g_m, df_ddp)
end

function FiniteDiffGradientPatternDistortion(G, F, ref_c, p_0, ROIsize, delta, P_ivec, DD, Δ_ivec, ΔDD)
  grad = zeros(9)
  C_0 = TestCorrelationPatternDistortion(G, F, ref_c, p_0, ROIsize, P_ivec, DD, Δ_ivec, ΔDD)
  for i=1:9
    p_1 = copy(p_0)
    p_1[i] += delta
    C_1 = TestCorrelationPatternDistortion(G, F, ref_c, p_1, ROIsize, P_ivec, DD, Δ_ivec, ΔDD)
    grad[i] = (C_1 - C_0)/delta
  end
  return grad
end

function TestCorrelationPatternDistortion(G, F, ref_c, p_old, ROIsize, P_ivec, DD, Δ_ivec, ΔDD)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  G_coeff = calc_B_coeffs(G, pad_size)

  ROIrelative = SquareRoiEvenSize(ROIsize)
  f = EvaluateWarpedImagePatternDistortion(F_coeff, ref_c' .+ ROIrelative, mattovec(eye(3)), pad_size, P_ivec, DD, [0.0,0.0], 0)
  f_m = mean(f)
  g = EvaluateWarpedImagePatternDistortion(G_coeff, ref_c' .+ ROIrelative, p_old, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
  g_m = mean(g)
  C_ls = ComputeCorrelationCriteria(f, f_m, g, g_m)
end

function TestCCforPD(G, F, ref_c, initial_guess_M, ROIsize,P_ivec, DD, Δ_ivec, ΔDD)
  pad_size = 50
  F_coeff = calc_B_coeffs(F, pad_size)
  G_coeff = calc_B_coeffs(G, pad_size)

  ROIrelative = SquareRoiEvenSize(ROIsize)
  f = EvaluateWarpedImagePatternDistortion(F_coeff, ROIrelative .+ ref_c', mattovec(eye(3)), pad_size, P_ivec, DD, [0,0], 0)
  f_m = mean(f)

  p_old = initial_guess_M
  g = EvaluateWarpedImagePatternDistortion(G_coeff, ROIrelative .+ ref_c', p_old, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
  g_m = mean(g)

  CC = ComputeCorrelationCriteria(f, f_m, g, g_m)
end
