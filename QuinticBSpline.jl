#using Plots
using FFTW

function sample_quintic_kernel(N)
  return fft([11./20  13./60 1./120 0 zeros(1,ceil(Int,(N-6))) 1./120 13./60 ])
end


function pad_image_w_border(I,padding)
  mn = size(I)
  PI = zeros(mn[1]+2*padding, mn[2]+2*padding)
  PI[padding+1:padding+mn[1],padding+1:padding+mn[2]] = I;
  for i=1:padding
    PI[padding+1:padding+mn[1],i] = I[1:end,1]
    PI[padding+1:padding+mn[1],padding + mn[2] + i] = I[1:end,end]
  end
  for i=1:padding
    PI[i,:] = PI[padding+1,1:end]
    PI[padding + mn[1] + i,:] = PI[padding+mn[1],1:end]
  end
  return PI
end


function convolve_for_spline_coeffs(I, FK)
  mn = size(I)
  B = zeros(mn[1],mn[2])
  thisrow = Array{Complex{Float64}}(1,mn[2])
  for i=1:mn[1]
    thisrow[1,:] = fft(I[i,:])
    B[i,:] = real(ifft(thisrow./FK))
  end
  return B
end


function calc_B_coeffs(I,padding)
  PI = pad_image_w_border(I,padding)
  mn = size(PI)
  FKr = sample_quintic_kernel(mn[2])
  rB = convolve_for_spline_coeffs(PI, FKr)
  FKc = sample_quintic_kernel(mn[1])
  B = convolve_for_spline_coeffs(rB',FKc)'
end


function get_QK()
  QK = [1/120 13/60 11/20 13/60 1/120 0;
        -1/24 -5/12 0 5/12 1/24 0;
        1/12 1/6 -1/2 1/6 1/12 0;
        -1/12 1/6 0 -1/6 1/12 0
        1/24 -1/6 1/4 -1/6 1/24 0;
        -1/120 1/24 -1/12 1/12 -1/24 1/120]
end


function SplineEvaluate(B_coeff, ROI, pad_size)
  QK = get_QK()
  f = Array{Float64}(size(ROI, 1))
  for i in 1:size(ROI,1)
    x_floor = convert(Array{Int,1}, floor.(ROI[i,:]))
    dx, dy = ROI[i, :] - x_floor
    dx_vec = [dx^n for n in 0:5]
    dy_vec = [dy^n for n in 0:5]
    
    c = B_coeff[pad_size + x_floor[1]-2:pad_size + x_floor[1] + 3,
                pad_size + x_floor[2]-2:pad_size + x_floor[2] + 3]'
    f[i] = dy_vec' * QK * c * QK' * dx_vec
  end
  return f
end


function SplineDerivative(B_coeff, ROI, pad_size)
  QK = get_QK()
  df = Array{Float64}(size(ROI))
  vec1 = zeros(6)
  vec1[1] = 1
  vec2 = zeros(6)
  vec2[2] = 1
  for i in 1:size(ROI,1)
    x_floor = convert(Array{Int,1}, floor.(ROI[i,:]))
    
    c = B_coeff[pad_size + x_floor[1]-2:pad_size + x_floor[1] + 3,
                pad_size + x_floor[2]-2:pad_size + x_floor[2] + 3]'
    df[i,1] = vec1' * QK * c * QK' * vec2
    df[i,2] = vec2' * QK * c * QK' * vec1
  end
  return df
end


function TestSplineIntepolation()
    plotly()

    A = [1 1 1 1 1 1;
         1 .95 .35 .02 .24 .85;
         1 .49 0 0 0 .26;
         1 .41 0 0 0 .18;
         1 .84 .06 0 .01 .64;
         1 1 .92 .71 .87 1]
    #A = Array(124:-1:1).*Array(1:124)'
    #A = rand(124,124)
    A = ones(50, 50)
    A[20:30, 30:40]=0

    pad_size=2

    B = calc_B_coeffs(A,pad_size)

    #display(heatmap(B))

    Abig = zeros(size(A,1)*10, size(A,2)*10)
    AshiftY = zeros(size(A))
    AshiftX = zeros(size(A))
    AbigX = zeros(size(A))
    AbigY = zeros(size(A))
    for i=1:size(A, 1) - 1
      for j=1:size(A, 2) - 1 
        for k=0:9
          for l=0:9
            Abig[i*10+k,j*10+l] = SplineEvaluate(B, [i+k/10 j+l/10], pad_size)[1]
          end
        end
        AshiftX[i,j] = SplineEvaluate(B, [i+0.5 j], pad_size)[1]
        AshiftY[i,j] = SplineEvaluate(B, [i j+0.5], pad_size)[1]
        AbigX[i,j] = SplineDerivative(B, [i*1.0 j*1.0], pad_size)[1, 1]
        AbigY[i,j] = SplineDerivative(B, [i*1.0 j*1.0], pad_size)[1, 2]
      end
    end

    #display(heatmap(pad_image_w_border(A,2)))

    display(heatmap(A[1:end-1, 1:end-1]', title="original"))
    display(heatmap(AshiftX[1:end-1, 1:end-1]', title="shift X"))
    display(heatmap(AshiftY[1:end-1, 1:end-1]', title="shift Y"))
    display(heatmap(Abig[1:end-10, 1:end-10]', title="interpolation"))
    display(heatmap(AbigX[1:end-1, 1:end-1]', title="dx"))
    display(heatmap(AbigY[1:end-1, 1:end-1]', title="dy"))
end

