using Plots, FFTW
plotly()

function sample_quintic_kernel(N)
  #return fft([zeros(1,floor(Int,(N-6)/2)-2) 1./120 13./60 11./20 13./60 1./120 0 zeros(1,ceil(Int,(N-6)/2)+2)])
  return fft([11./20 13./60 1./120 0 zeros(1,ceil(Int,(N-6))) 1./120 13./60])
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
  thisrow = Array(Complex{Float64},1,mn[2])
  for i=1:mn[1]
    thisrow[1,:] = fft(I[i,:])
    B[i,:] = (ifft(thisrow./FK)) #Do things need to get real?
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

function get_val(B, x, y)
  xpix = floor(Int,x)
  ypix = floor(Int,y)
  dx = x - xpix;
  dy = y - ypix
  C = B[xpix-2:xpix+3,ypix-2:ypix+3]
  QK = get_QK()
  ypowvec = ones(1,6)
  xpowvec = ones(6,1)
  for i=1:5
    display(i)
    for j=i+1:6
      display(j)
      xpowvec[j,1] = xpowvec[j,1]*dx
      ypowvec[1,j] = ypowvec[1,j]*dy
    end
  end
  val = ypowvec*QK*C*QK'*xpowvec
  val[1,1]
end

A = [1 1 1 1 1 1;1 .95 .35 .02 .24 .85;1 .49 0 0 0 .26;1 .41 0 0 0 .18;1 .84 .06 0 .01 .64;1 1 .92 .71 .87 1]

B = calc_B_coeffs(A,2)

display(heatmap(B))

Abig = zeros(50,50)
for i=1:50
  for j=1:50
    Abig[i,j] = get_val(B,(i*4.0)/50+3,(j*4.0)/50+3)
  end
end

display(heatmap(pad_image_w_border(A,2)))
display(heatmap(Abig))
