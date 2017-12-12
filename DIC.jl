
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
  
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = GetSteepestDescentImages(F_coeff, ref_c, ROIrelative)
  
  f = EvaluateWarpedImage(F_coeff, ref_c, ROIrelative, zeros(6))
  f_m = mean(f)
  hess = ComputeHessian(f, f_m, df_ddp)
  p_old = [initial_guess_u[1], initial_guess_u[2], 0, 0, 0, 0]  
  converged = false
  num_iterations = 0
  
  while !converged && num_iterations < 100
    println(num_iterations)
    num_iterations += 1
    g = EvaluateWarpedImage(G_coeff, ref_c, ROIrelative, p_old)
    g_m = mean(g)
    grad = ComputeGradient(f, f_m, g, g_m, df_ddp)
    dp = hess\(-grad)
    p_old = UpdateP(p_old, dp)
    converged = (norm(dp) < 1.0e-6)
  end

  println(p_old)
  return p_old[[1, 2]]
end


function SquareRoiEvenSize(ROIsize)
  ROIrelative = Array{Int}(ROIsize*ROIsize, 2)
  count = 0
  for i in ROIsize : ROIsize/2 - 1  
    for i in ROIsize : ROIsize/2 - 1 
      count += 1
      ROIrelative[count, :] = [i, j]
    end
  end
  return ROIrelative 
end


function GetSteepestDescentImages(F, ref_c, ROIrelative)
  # this function assumes x_ref_tilde = x_ref
  # ddp = [ddu, ddv, ddu_dx, ddu_dy, ddv_dx, ddv_dy]
  df_ddp = zeros(size(ROIrelative, 1), 6)
  
  ROI = ROIrelative .+ ref_c' 
  #x = ROIrange + ref_c[1] + p[1] + p[3]*ROIrange + p[4]*ROIrange'
  #y = ROIrange' + ref_c[2] + p[2] + p[5]*ROIrange + p[6]*ROIrange'
  
  df_dxy = SplineDerivative(F, ROI) # need to implement this function
  df_ddp[:, 1] = df_dxy[:, 1]
  df_ddp[:, 2] = df_dxy[:, 2]
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
  for i in length(f)
    gradient += ((f[i]-f_m)/norm_f - (g[i]-g_m)/norm_g)* df_ddp[i, :]
  end
  gradient *= 2/norm_f 
  return gradient
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
end


function EvaluateWarpedImage(F, ref_c, ROIrelative, p)
  ROI = ROIrelative .+ ref_c' .+ p[1:2]' + 
        ROIrelative .* p[[3, 5]]' +
        ROIrelative[:, [2, 1]] .* p[[4, 6]]'
  return SplineEvaluate(F, ROI)
end


function SplineDerivative(F_coeff, ROI)
  pad_size = 2
  QK = [1/120  13/60  11/20  13/60  1/120  0;
        -1/24  -5/12  0  5/12  1/24  0;
        1/12  1/6  -1/2  1/6  1/12  0;
        -1/12  1/6  0  -1/6  1/12  0;
        1/24  -1/6 1/4  -1/6  1/24  0;
        -1/120  1/24  -1/12  1/12  -1/24  1/120]
  df = Array{float}(size(ROI))
  vec1 = zeros(6)
  vec1[1] = 1
  vec2 = zeros(6)
  vec2[2] = 1
  for i in 1:size(ROI,1)
    x_floor = convert(Array{Int,1}, floor(ROI[i,:]))
    
    c = F_coeff[pad_size + x_floor[1]-2:pad_size + x_floor[1] + 3,
                pad_size + x_floor[2]-2:pad_size + x_floor[2] + 3]
    df[i,1] = vec1' * QK * c * QK' * vec2
    df[i,2] = vec2' * QK * c * QK' * vec1
  end
  return f
end


function SplineEvaluate(F_coeff, ROI)
  pad_size = 2
  QK = [1/120  13/60  11/20  13/60  1/120  0;
        -1/24  -5/12  0  5/12  1/24  0;
        1/12  1/6  -1/2  1/6  1/12  0;
        -1/12  1/6  0  -1/6  1/12  0;
        1/24  -1/6 1/4  -1/6  1/24  0;
        -1/120  1/24  -1/12  1/12  -1/24  1/120]
  f = Array{float}(size(ROI, 1))
  for i in 1:size(ROI,1)
    x_floor = convert(Array{Int,1}, floor(ROI[i,:]))
    dx, dy = ROI[i, :] - x_floor
    dx_vec = [dx^n for i in 0:5]
    dy_vec = [dy^n for i in 0:5]
    
    c = F_coeff[pad_size + x_floor[1]-2:pad_size + x_floor[1] + 3,
                pad_size + x_floor[2]-2:pad_size + x_floor[2] + 3]
    f[i] = dy_vec' * QK * c * QK' * dx_vec
  end
  return f
end
