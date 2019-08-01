function HyperCCGrad(F_coeff, G_coeff, ROI, p_old, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
    g = EvaluateWarpedImagePatternDistortion(G_coeff, ROI, p_old, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
    g_m = mean(g)
    dfdp = zeros(9)
    h = 1e-12
    for i = 1:9
        hvec = complex(zeros(9))
        hvec[i] = h*im
        f = EvaluateWarpedImagePatternDistortionComplex(F_coeff, ROI, mattovec(I) + hvec, pad_size, P_ivec, DD, [0,0], 0)
        f_m = mean(f)
        dfdp[i] = imag(ComputeCorrelationCriteria(f, f_m, g, g_m))/h
    end
    return dfdp
end

function EvaluateWarpedImagePatternDistortionComplex(F, ROIabsolute, p, pad_size, P_ivec, DD, Δ_ivec, ΔDD)
  M = VR_deviatoric(vectomat(p))

  xs = [ROIabsolute[:,1].-P_ivec[1] ROIabsolute[:,2].-P_ivec[2]]
  ROI = complex(zeros(size(ROIabsolute)))
  mn = size(F)
  boundx = mn[1] - 2*pad_size
  boundy = mn[2] - 2*pad_size
  outofbounds = false
  for i=1:size(ROIabsolute,1)
    x = xs[i,:]
    Mr = M*[x[1];x[2];-DD]
    ROInew = P_ivec + Δ_ivec - Mr[1:2]*(DD + ΔDD)/Mr[3]
    ROInew[2] = conj(ROInew[2])#WHY???!!!!????!!!!
    ROI[i,:] = ROInew
    if ~(boundx > real(ROInew[1]) > 1) || ~(boundy > real(ROInew[2]) > 1)
        outofbounds = true
    end
  end
  if ~outofbounds
      f = SplineEvaluateComplex(F, ROI, pad_size)
      return f
  else
      display("Out of bounds")
      return []
  end
end

function SplineEvaluateComplex(B_coeff, ROI, pad_size)
  QK = get_QK()
  f = complex(Array{Float64}(undef,size(ROI, 1)))
  for i in 1:size(ROI,1)
    x_floor = convert(Array{Int,1}, floor.(real(ROI[i,:])))
    dx, dy = ROI[i, :] - x_floor
    dx_vec = [dx^n for n in 0:5]
    dy_vec = [dy^n for n in 0:5]

    c = B_coeff[pad_size + x_floor[1]-2:pad_size + x_floor[1] + 3,
                pad_size + x_floor[2]-2:pad_size + x_floor[2] + 3]'
    f[i] = dy_vec' * QK * c * QK' * dx_vec
  end
  return f
end
