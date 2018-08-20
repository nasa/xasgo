function F_NCorrShift(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess, ccfilt, windowfunc)
  m = minimum(size(DefI))  # smallest dimension of image
  ROIsize_px = 2*round(Int64,m*ROIsize/200) #ROI forced to be even and square
  qs = zeros(3,size(ROIs)[1])
  rs = zeros(3,size(ROIs)[1])
  PD_p_px = PC_to_phosphor_frame(PD,m)
  PR_p_px = PC_to_phosphor_frame(PR,m)
  for i in 1:size(ROIs)[1]
    ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
    ROI_p_px = image_vec_to_phosphor_frame(ROI_px)
    r_p = ROI_p_px - PR_p_px
    shiftROI_px = offset_estimate_w_F(r_p,PD_p_px,PR_p_px,F_guess) + ROI_px
    thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
    qstar_p = reverse_PC_offset(thisq, r_p, PD_p_px, PR_p_px)
    qs[:,i] = qstar_p/m #Why does this fail if I take out the /m?
    rs[:,i] = r_p/m
  end
  β = beta_calc_Ruggles(qs,rs, false)
  F = VR_deviatoric(β+eye(3))
end

function F_NCorrShiftIterate(ROIs, ROIsize, DefI, PD, gD, RefI, PR, gR, F_guess, ccfilt, windowfunc; numimax = 10)
  m = minimum(size(DefI))  # smallest dimension of image
  ROIsize_px = 2*round(Int64,m*ROIsize/200) #ROI forced to be even and square
  qs = zeros(3,size(ROIs)[1])
  rs = zeros(3,size(ROIs)[1])
  PD_p_px = PC_to_phosphor_frame(PD,m)
  PR_p_px = PC_to_phosphor_frame(PR,m)
  converged = false
  num_iterations = 0
  F_new = eye(3)
  while !converged && num_iterations < numimax
    num_iterations += 1
    for i in 1:size(ROIs)[1]
      ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
      ROI_p_px = image_vec_to_phosphor_frame(ROI_px)
      r_p = ROI_p_px - PR_p_px
      shiftROI_px = offset_estimate_w_F(r_p,PD_p_px,PR_p_px,F_guess) + ROI_px
      thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
      qstar_p = reverse_PC_offset(thisq, r_p, PD_p_px, PR_p_px)
      qs[:,i] = qstar_p/m #Why does this fail if I take out the /m?
      rs[:,i] = r_p/m
    end
    β = beta_calc_Ruggles(qs,rs, false)
    F_new = VR_deviatoric(β+eye(3))
    if norm(F_new-F_guess)<100e-6
      converged = true
    else
      F_guess = F_new
    end
  end
  if num_iterations > 1
    println("Number of iterations: ", num_iterations)
  end
  return F_new
end

function F_NCorrRemap(ROIs, ROIsize, WarpROI, DefI, PD, gD, RefI_coeffs, pad_size, PR, gR, F_guess, ccfilt, windowfunc; numimax = 10)
  m = minimum(size(DefI))
  ROIsize_px = 2*round(Int64,m*ROIsize/200) #ROI forced to be even and square
  qs = zeros(3,size(ROIs)[1])
  rs = zeros(3,size(ROIs)[1])
  PD_p_px = PC_to_phosphor_frame(PD,m)
  PR_p_px = PC_to_phosphor_frame(PR,m)
  converged = false
  num_iterations = 0
  F_guess = VR_deviatoric(F_guess)
  F_new = eye(3)

  while !converged && num_iterations < numimax
    num_iterations += 1
    #display(num_iterations)
    F_remap = F_guess#VR_poldec(F_guess)[2]
    RefI = patternremap(RefI_coeffs, pad_size, WarpROI, PR_p_px, PD_p_px-PR_p_px, F_remap, m)
    if RefI == []
        F_new = eye(3) + ones(3,3)
        converged = true
    else
        for i in 1:size(ROIs)[1]
          ROI_px = round.(m*ROIs[i,:] + [.5,.5]) -  [.5,.5]
          ROI_p_px = image_vec_to_phosphor_frame(ROI_px)
          r_p = ROI_p_px - PR_p_px
          thisq = MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, ROI_px, ROIsize_px, windowfunc, ccfilt)
          qs[:,i] = image_vec_to_phosphor_frame(thisq)/m
          rs[:,i] = r_p/m
        end
        if num_iterations>100
            qnorms = m*sqrt.(sum(qs.*qs,1))
            #display(histogram(squeeze(qnorms,1)))
            display(mean(qnorms))
        end
        β = beta_calc_Ruggles(qs,rs, false)
        F_new = VR_deviatoric(β+eye(3))*F_remap
        #display(norm(F_new-F_guess))
        if norm(F_new-F_guess)<10e-6
          converged = true
        else
          F_guess = F_new
        end
    end
  end
  if (num_iterations >= 20) && (~converged)
      F_new = eye(3) + ones(3,3)
  end
  if num_iterations > 1
    println("Number of iterations: ", num_iterations)
  end
  return F_new
end
