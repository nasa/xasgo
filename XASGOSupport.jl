function finite_rotation_and_small_strain_to_F(e11,e12,e13,e22,e23,e33,a,b,c,d)
  e = [e11 e12 e13;e12 e22 e23;e12 e23 e33]
  R = versor_to_R([a,b,c,d])
  return F = (eye(3)+e)*R
end

function versor_to_R(v)
  a = v[1]
  b = v[2]
  c = v[3]
  d = v[4]

  R =[1-2*(c*c+d*d) 2*b*c+2*a*d 2*b*d-2*a*c;
    2*b*c-2*a*d 1-2*(b*b+d*d) 2*c*d+2*a*b;
    2*b*d+2*a*c 2*c*d-2*a*b 1-2*(b*b+c*c)]
end

function euler_to_gmat(phi1, PHI, phi2)
  #Input angles may be vectors
  g = zeros(3,3,length(phi1));
  cp1 = cos(phi1);
  sp1 = sin(phi1);
  cp2 = cos(phi2);
  sp2 = sin(phi2);
  cP = cos(PHI);
  sP = sin(PHI);

  g[1,1,:]= cp1.*cp2-sp1.*sp2.*cP;
  g[1,2,:] = sp1.*cp2+cp1.*sp2.*cP;
  g[1,3,:] = sp2.*sP;
  g[2,1,:]= -cp1.*sp2-sp1.*cp2.*cP;
  g[2,2,:]= -sp1.*sp2+cp1.*cp2.*cP;
  g[2,3,:]= cp2.*sP;
  g[3,1,:]=  sp1.*sP;
  g[3,2,:]= -cp1.*sP;
  g[3,3,:]=  cP;

  if size(g,3) == 1
    g = squeeze(g,3)
  end

  return g
end

function LSfit(qs,rs,F)
  sse = 0;
  for i in 1:size(qs)[2]
    r = rs[:,i]
    q = [qs[1,i],qs[2,i],0]
    Fr = F*r
    rp = Fr*r[3]/Fr[3]
    errorvec = rp - q - r
    sse += sum(errorvec.*errorvec)
  end
  return sse
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

function squareimage(I)
  mn = size(I)
  buffer = round(Int64,(mn[2] - mn[1])/2)
  I_square = I[:,buffer+1:buffer+mn[1]]
end

function VR_poldec(A)
  #Assume square matrix
  S = svd(A)
  R = S[1]*(S[3]')
  V = S[1]*diagm(S[2])*(S[1]')
  VR = [V, R]
end

function VR_deviatoric(A)
  VR = VR_poldec(A)
  return (eye(3) + VR[1] - trace(VR[1])*eye(3)/3.)*VR[2]
end

function beta_calc_Ruggles(qs,rs)
  r1 = rs[1,:]
  r2 = rs[2,:]
  r3 = rs[3,:] #r3 is just a constant times a vector of ones...that's a wee bit silly
  q1 = qs[1,:]
  q2 = qs[2,:]
  zerovec = zeros(size(r1))

  A1 = [r1.*r3 r2.*r3 r3.*r3 zerovec zerovec zerovec -(r1.*r1 + q1.*r1) -(r1.*r2 + q1.*r2) -(r1.*r3 + q1.*r3)]
  A2 = [zerovec zerovec zerovec r1.*r3 r2.*r3 r3.*r3 -(r2.*r1 + q2.*r1) -(r2.*r2 + q2.*r2) -(r2.*r3 + q2.*r3)]
  b1 = q1.*r3;
  b2 = q2.*r3;

  b3 = 0
  A3 = [1 0 0 0 1 0 0 0 1]

  A = [A1;A2;A3]
  b = [b1;b2;b3]

  X = A\b
  β = [X[1] X[2] X[3];X[4] X[5] X[6];X[7] X[8] X[9]]
end

function beta_calc_FiniteR(qs,rs, v0)
  #m = Model(solver=NLoptSolver(algorithm=:LD_MMA))
  m = Model(solver=IpoptSolver(print_level=0))

  @variable(m, e11)
  @variable(m, e12)
  @variable(m, e13)
  @variable(m, e22)
  @variable(m, e23)
  @variable(m, e33)
  @variable(m, 1 >= a >= -1)
  @variable(m, 1 >= b >= -1)
  @variable(m, 1 >= c >= -1)
  @variable(m, 1 >= d >= -1)

  setvalue(e11, 0.0)
  setvalue(e12, 0.0)
  setvalue(e13, 0.0)
  setvalue(e22, 0.0)
  setvalue(e23, 0.0)
  setvalue(e33, 0.0)
  setvalue(a, v0[1])
  setvalue(b, v0[2])
  setvalue(c, v0[3])
  setvalue(d, v0[4])

  CMfitfun(e11,e12,e13,e22,e23,e33,a,b,c,d) = LSfit(qs,rs,finite_rotation_and_small_strain_to_F(e11,e12,e13,e22,e23,e33,a,b,c,d))
  JuMP.register(m,:CMfitfun,10,CMfitfun,autodiff=true)
  @NLobjective(m,Min,CMfitfun(e11,e12,e13,e22,e23,e33,a,b,c,d))

  @NLconstraint(m, a*a + b*b + c*c + d*d == 1.0)
  @NLconstraint(m, e11 + e22 + e33 == 0.0) #Not actually non-linear...

  status = solve(m)

  display(status)
  display([getvalue(e11),getvalue(e12),getvalue(e13),getvalue(e22),getvalue(e23),getvalue(e33),getvalue(a),getvalue(b),getvalue(c),getvalue(d)])
  F = finite_rotation_and_small_strain_to_F(getvalue(e11),getvalue(e12),getvalue(e13),getvalue(e22),getvalue(e23),getvalue(e33),getvalue(a),getvalue(b),getvalue(c),getvalue(d))

  β = F - eye(3)
end

function find_change_in_frame(A,B)
  #Assume A and B are the same tensor in different frames, find g such that A = g*B*inv(g)
  DA = eig(A)
  DB = eig(B)
  AA = DA[2]
  BB = DB[2]
  g = zeros(3,3)
  for i=1:3
    for j=1:3
      g[i,j] = real(dot(AA[i,:],BB[j,:]))
    end
  end
  return g
end

function MeasureShiftClassic_w_Shift(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
  RS2 = round(Int, ROIsize_px/2)

  ROIC = round.(Int,ROI_px + [0.5;0.5])
  sROIC = round.(Int,shiftROI_px + [0.5;0.5])

  rrange=ROIC[1] - RS2:ROIC[1] - RS2 + ROIsize_px - 1
  crange=ROIC[2] - RS2:ROIC[2] - RS2 + ROIsize_px - 1

  srrange=sROIC[1] - RS2:sROIC[1] - RS2 + ROIsize_px - 1
  scrange=sROIC[2] - RS2:sROIC[2] - RS2 + ROIsize_px - 1

  DROI = DefI[srrange,scrange]
  RROI = RefI[rrange,crange]

  FD = rfft(windowfunc.*DROI)
  FR = rfft(windowfunc.*RROI)#windowfunc.*

  CC = fftshift(brfft(ccfilt[1:round(Int,ROIsize_px/2)+1,:].*FD.*conj(FR),ROIsize_px))#ccfilt[1:round(Int,ROIsize_px/2)+1,:].*

  #println(findpeak2(CC) - ROIsize_px/2 - 1)

  q = findpeak2(CC) - [rrange[1] - srrange[1],crange[1] - scrange[1]] - ROIsize_px/2 - 1
end

function Plot_ROI(I, ROI_px, ROIsize_px)
  RS2 = round(Int, ROIsize_px/2)
  ROIC = round.(Int,ROI_px + [0.5;0.5])
  rrange=ROIC[1] - RS2:ROIC[1] - RS2 + ROIsize_px - 1
  crange=ROIC[2] - RS2:ROIC[2] - RS2 + ROIsize_px - 1
  ROI = I[rrange,crange]
  display(heatmap(ROI))
end

function PlotCC(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px, windowfunc, ccfilt)
  RS2 = round(Int, ROIsize_px/2)

  ROIC = round.(Int,ROI_px + [0.5;0.5])
  sROIC = round.(Int,shiftROI_px + [0.5;0.5])

  rrange=ROIC[1] - RS2:ROIC[1] - RS2 + ROIsize_px - 1
  crange=ROIC[2] - RS2:ROIC[2] - RS2 + ROIsize_px - 1

  srrange=sROIC[1] - RS2:sROIC[1] - RS2 + ROIsize_px - 1
  scrange=sROIC[2] - RS2:sROIC[2] - RS2 + ROIsize_px - 1

  DROI = DefI[srrange,scrange]
  RROI = RefI[rrange,crange]

  FD = rfft(DROI)
  FR = rfft(RROI)#windowfunc.*

  CC = fftshift(brfft(FD.*conj(FR),ROIsize_px))#ccfilt[1:round(Int,ROIsize_px/2)+1,:].*

  display(heatmap(CC))
end

function image_vec_to_phosphor_frame(Iv)
  Pv = [-Iv[2],-Iv[1],0]
end

function phosphor_frame_to_image_vec(Pv)
  Iv = [-Pv[2],-Pv[1]]
end

function rotate_to_image_frame(Mp)
  Qfix = [0 -1 0;-1 0 0;0 0 1]
  Mi = Qfix*Mp*(Qfix')
end

function rotate_to_phosframe_from_image(Mi)
  Qfix = [0 -1 0;-1 0 0;0 0 1]
  Mp = (Qfix')*Mi*Qfix
end

function PC_to_phosphor_frame(P,m)
  P_p = -m*[P[1]; 1-P[2]; -P[3]]
end

function offset_estimate_simple(r, PD, PR)
  δP = PD - PR
  q_p = δP + r*δP[3]/PR[3]
  q_px = phosphor_frame_to_image_vec(q_p)
end

function offset_estimate_w_F(r, PD, PR, F)
  δP = PD - PR
  Fr = F*r
  q_p = δP - r + Fr*(-PD[3])/Fr[3]
  q_px = phosphor_frame_to_image_vec(q_p)
end

function reverse_PC_offset(q, r_p,PD,PR)
  q_p = image_vec_to_phosphor_frame(q)
  δP = PD - PR
  qstar_p = (q_p + r_p - δP)*PR[3]/PD[3] - r_p
end

function prep_ebsp(filepath)
  I = load(filepath)
  I = squareimage(I)
  I = convert(Array{ColorTypes.Gray{FixedPointNumbers.Normed{UInt8,8}},2},I)
  I = convert(Array{Float64},I)
  I = filterimage(I,9,90,25)
end
