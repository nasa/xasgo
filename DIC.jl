using Plots
plotly()

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
  pad_size = 500
  F_coeff = calc_B_coeffs(F, pad_size)
  G_coeff = calc_B_coeffs(G, pad_size)
  
  ROIrelative = SquareRoiEvenSize(ROIsize)
  df_ddp = GetSteepestDescentImages(F_coeff, ref_c, ROIrelative, pad_size)
  f = EvaluateWarpedImage(F_coeff, ref_c, ROIrelative, zeros(6), pad_size)
  display(heatmap(reshape(f, ROIsize, ROIsize)))
  f_m = mean(f)
  hess = ComputeHessian(f, f_m, df_ddp)
  p_old = [initial_guess_u[1], initial_guess_u[2], .00, 0, 0, 0]  
  converged = false
  num_iterations = 0
  println("0: ", p_old)
  while !converged && num_iterations < 100
    num_iterations += 1
    g = EvaluateWarpedImage(G_coeff, ref_c, ROIrelative, p_old, pad_size)
    #if num_iterations == 10
      #display(heatmap(reshape(g, ROIsize, ROIsize)))
    #end
    g_m = mean(g)
    C_ls = ComputeCorrelationCriteria(f, f_m, g, g_m)
    grad = ComputeGradient(f, f_m, g, g_m, df_ddp)
    dp = (hess\(-grad))
    p_old = UpdateP(p_old, dp)
    #println("   dp: ", norm(dp), dp)
    println(num_iterations, ": ", C_ls, "    ", norm(dp), "    ", p_old)
    println("grad: ", grad)
    converged = (norm(dp) < 1.0e-3)
  end

  println(p_old)
  return p_old[[1, 2]]
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
  df_ddp[:, 1] = df_dxy[:, 1]
  df_ddp[:, 2] = df_dxy[:, 2]
  df_ddp[:, 3] = df_dxy[:, 1] .* ROIrelative[:, 1]
  df_ddp[:, 4] = df_dxy[:, 1] .* ROIrelative[:, 2]
  df_ddp[:, 5] = df_dxy[:, 2] .* ROIrelative[:, 1]
  df_ddp[:, 6] = df_dxy[:, 2] .* ROIrelative[:, 2]
  
  #display(heatmap(reshape(ROIrelative[:,1], 12,12), title="roix"))
  #display(heatmap(reshape(ROIrelative[:,2], 12,12), title="roiy"))
  #
  #display(heatmap(reshape(df_ddp[:,1], 12,12), title="dfddu"))
  #display(heatmap(reshape(df_ddp[:,2], 12,12), title="dfddv"))
  #display(heatmap(reshape(df_ddp[:,3], 12,12), title="dfddudx"))
  #display(heatmap(reshape(df_ddp[:,4], 12,12), title="dfddudy"))
  #display(heatmap(reshape(df_ddp[:,5], 12,12), title="dfddvdx"))
  #display(heatmap(reshape(df_ddp[:,6], 12,12), title="dfddvdy"))
  
  return df_ddp
end


function ComputeHessian(f, f_m, df_ddp)
  hessian = zeros(6,6)
  for i in 1:length(f)
    hessian += df_ddp[i, :] * df_ddp[i, :]'
  end
  hessian *= 2.0 / sum((f-f_m).^2)
  display(hessian)
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
  #w_old = [1+p_old[6]  p_old[5]  p_old[2]; 
  #         p_old[4]  1+p_old[3]  p_old[1]; 
  #         0  0  1]
  #w_dp = [1+dp[6]  dp[5]  dp[2]; 
  #        dp[4]  1+dp[3]  dp[1];
  #        0  0  1]
  w_new = w_old*inv(w_dp)
  p_new = [w_new[1, 3], w_new[2, 3], w_new[1, 1] - 1, 
           w_new[1, 2], w_new[2, 1], w_new[2, 2] - 1]
  #p_new = [w_new[2, 3], w_new[1, 3], w_new[2, 2] - 1, 
  #         w_new[2, 1], w_new[1, 2], w_new[1, 1] - 1]
end


function EvaluateWarpedImage(F, ref_c, ROIrelative, p, pad_size)
  ROI = ROIrelative .+ ref_c' .+ p[1:2]' + 
        ROIrelative .* p[[3, 5]]' +
        ROIrelative[:, [2, 1]] .* p[[4, 6]]'
  return SplineEvaluate(F, ROI, pad_size)
end


"""
r = rand(49,49)
for i = 1:49
  for j = 1:49
    r[i, j] = (i/50*j/50)^2 
  end
end
"""

F = ones(50, 50)
F[11:16,11:16] = [1 1 1 1 1 1;1 .95 .35 .02 .24 .85;1 .49 0 0 0 .26;1 .41 0 0 0 .18;1 .84 .06 0 .01 .64;1 1 .92 .71 .87 1]
F[11:16,11:16] = [1 1 1 1 1 1;
                  1 .95 .85 .85 .95 1;
                  1 .85 .6 .6 .85 1;
                  1 .85 .6 .6 .85 1;
                  1 .95 .85 .85 .95 1;
                  1 1 1 1 1 1]


#F[1:49,1:49] = r

G = ones(50, 50)
#G[13:18,11:16] = [1 1 1 1 1 1;1 .95 .35 .02 .24 .85;1 .49 0 0 0 .26;1 .41 0 0 0 .18;1 .84 .06 0 .01 .64;1 1 .92 .71 .87 1]
#G[2:50,1:49] = r
tmp = EvaluateWarpedImage(calc_B_coeffs(F[8:19,8:19], 500), [6.5, 6.5], SquareRoiEvenSize(12), 
                          [0,-0.5,0,0,0,0], 500)
G[8:19,8:19] = reshape(tmp,12, 12)'

display(heatmap(F'))
display(heatmap(G'))

RunIcgn(G, F, [15.5, 15.5], [0.0, 0.5], 12)



