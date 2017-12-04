using Images, Colors, FixedPointNumbers, Plots
using FFTW
plotly()

function CalcF(ROIs, ROIsize, DefI, PD, RefI, PR)

m = minimum(size(DefI))
qs = zeros(size(ROIs))
Qvp=[-1 0 0; 0 -1 0; 0 0 1]
ROIpix = 2*convert(Int64,m*ROIsize/200) #ROI forced to be even and square
ccfilt = ccfilter(2,50, [ROIpix,ROIpix], 13)
windowfunc = ccwindow(ROIpix)
qerror = zeros(2,size(ROIs)[1])
for i in 1:size(ROIs)[1]

  shiftROI = shift_estimate_simple(ROIs[i,:],PD,PR, m, Qvp)

  trueq = shiftROI - round.(m*[ROIs[i,2];ROIs[i,1]])
  truepeak = [trueq[2];trueq[1]] + [ROIpix+1,ROIpix+1]
  #println([trueq[2], trueq[1]])

  thisq = MeasureShiftClassic_w_Shift(DefI, RefI, m*ROIs[i,:], shiftROI, ROIpix, windowfunc, ccfilt)

  #println(thisq)
  #println(" ")
  qerror[:,i] = thisq - [trueq[2], trueq[1]]

end
#display(scatter(qerror[1,:],qerror[2,:]))
println(mean(qerror[1,:]))
println(mean(qerror[2,:]))
println(std(qerror[1,:]))
println(std(qerror[2,:]))
end

function MeasureShiftClassic_w_Shift(DefI, RefI, ROI, shiftROI, ROIpix, windowfunc, ccfilt)
  rrange=round(Int,ROI[1]-ROIpix/2):round(Int,ROI[1]-ROIpix/2)+ROIpix-1
  crange=round(Int,ROI[2]-ROIpix/2):round(Int,ROI[2]-ROIpix/2)+ROIpix-1

  srrange=round(Int,shiftROI[2] - ROIpix/2):round(Int,shiftROI[2] - ROIpix/2)+ROIpix-1
  scrange=round(Int,shiftROI[1] - ROIpix/2):round(Int,shiftROI[1] - ROIpix/2)+ROIpix-1

  DROI = DefI[srrange,scrange]
  RROI = RefI[rrange,crange]

  FD = rfft(windowfunc.*DROI)
  FR = rfft(windowfunc.*RROI)#windowfunc.*

  CC = fftshift(brfft(ccfilt[1:round(Int,ROIpix/2)+1,:].*FD.*conj(FR),ROIpix))#ccfilt[1:round(Int,ROIpix/2)+1,:].*

  q = findpeak2(CC) - [rrange[1] - srrange[1],crange[1] - scrange[1]] - ROIpix/2 - 1
end

function shift_estimate_simple(ROI, PD, PR, m, Qvp)
  r = Qvp*round.(m*[ROI[2];ROI[1];0]) + m*[PR[1];1-PR[2];-PR[3]]
  δP = PD - PR
  δP = m*[-δP[1];δP[2];δP[3]]
  shiftROI3D = Qvp'*(-m*[PR[1];1-PR[2];-PR[3]] + δP + r*PD[3]/PR[3])
  shiftROI = [shiftROI3D[1];shiftROI3D[2]]
end

function findpeak(A)
  mn = size(A)
  maxdat = findmax(A)
  xyind = [(maxdat[2]-1)%mn[1] + 1,convert(Int,ceil(maxdat[2]/mn[1]))]
  X = [1 -1 1;1 0 0;1 1 1]
  xs = [A[xyind[1]-1,xyind[2]],A[xyind[1],xyind[2]],A[xyind[1]+1,xyind[2]]]
  ys = [A[xyind[1],xyind[2]-1],A[xyind[1],xyind[2]],A[xyind[1],xyind[2]+1]]
  Bx = X\xs
  By = X\ys
  fittedpeak = xyind -.5*[Bx[2]/Bx[3],By[2]/By[3]]
  return fittedpeak
end

function findpeak2(A)
  mn = size(A)
  maxdat = findmax(A)
  xyind = [(maxdat[2]-1)%mn[1] + 1,convert(Int,ceil(maxdat[2]/mn[1]))]
  X = [1 -1 -1 1 1 1;
        1 -1 0 0 1 0;
        1 -1 1 -1 1 1;
        1 0 -1 0 0 1;
        1 0 0 0 0 0;
        1 0 1 0 0 1;
        1 1 -1 -1 1 1;
        1 1 0 0 1 0;
        1 1 1 1 1 1]
  count = 0;
  V = zeros(9,1)
  for i in -1:1
    for j in -1:1
      count = count + 1
      V[count] = A[xyind[1] + i, xyind[2] + j]
    end
  end
  C = X\V
  denom = 4*C[6]*C[5] - C[4]*C[4]
  fittedpeak = xyind + [C[4]*C[3] - 2.*C[6]*C[2],C[4]*C[2] - 2.*C[5]*C[3]]/denom
end

function ccfilter(lb,ub, mn, smooth)
  ccfilt = zeros(mn[1],mn[2])
  mid = round(Int,(mn[1]+1.001)/2)
  for i in 1:mn[1]
    for j in 1:mn[2]
      dist = sqrt((i-mid)^2 + (j-mid)^2)
      ccfilt[i,j] = 0.0*(dist>(ub+smooth)) + erfc(pi*(dist-ub)/smooth)*((ub+smooth)>=dist>ub) + 1.0*(ub>=dist>lb) + erfc(-pi*(dist-lb)/smooth)*(lb>=dist>lb-smooth)
    end
  end
  ccfilt[mid,mid] = 0.0
  ccfilt = fftshift(ccfilt)
end

function filterimage(I,lb,ub,smooth)
  mn = size(I)
  filt = ccfilter(lb,ub,mn,smooth)
  I = irfft(filt[1:round(Int,mn[1]/2)+1,:].*rfft(I),mn[1])
  return I
end

function ccwindow(L)
  D = zeros(L,L)
  for i=1:L
    for j=1:L
      D[i,j] = i;
    end
  end
  mid = (L+1)/2
  windowfunc = cos.((D-mid)*pi/L).*cos.((D'-mid)*pi/L)
  return windowfunc
end

cd("D:\\XASGO")
RefI = load("ZeroCamElevation_x0y0.png")
#display(heatmap(RefI))
RefI = convert(Array{ColorTypes.Gray{FixedPointNumbers.Normed{UInt8,8}},2},RefI)
RefI = convert(Array{Float64},RefI)
RefI = filterimage(RefI,9,90,25)
PR = [.5;.5;.7]

DefI = load("ZeroCamElevation_x500y500.png")
DefI = convert(Array{ColorTypes.Gray{FixedPointNumbers.Normed{UInt8,8}},2},DefI)
DefI = convert(Array{Float64},DefI)
DefI = filterimage(DefI,9,90,25)
PD = [.48046875;.518353371499725;.706680080924329] #x500y500
#PD = [.5;.51835337;.70668008] #x0y500
#PD = [0.498046875000000;0.501835337149973;0.700668008092433]#x50y50

N = 48

θs = 2*pi*collect(0:N-2)/(N-1)

ROIsize = 25
rad = .25

ROIs = Array{Float64}(N,2)
ROIs[1,:] = [.5 .5]
for i=2:N
  ROIs[i,:] = [.5 .5] + rad*[sin(θs[i-1]) cos(θs[i-1])]
end


CalcF(ROIs, ROIsize, DefI, PD, RefI, PR)

println("Done")