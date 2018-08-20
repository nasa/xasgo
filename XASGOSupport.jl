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

function angofF(F)
  return 180*acos((min(trace(F),3) - 1)/2)/pi
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

function gmat_to_euler(g)
  P = acos(sign(g[3,3])*min(abs(g[3,3]),1))
  p1 = atan2(g[3,1],-g[3,2])
	p2 = atan2(g[1,3], g[2,3])
  return [p1,P,p2]
end

function axisang_to_gmat(n,θ)
  K = [0 -n[3] n[2];n[3] 0 -n[1];-n[2] n[1] 0]
  g = eye(3) + sin(θ)*K + (1-cos(θ))*K*K
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

function findpeak3(A; n=1)
  mn = size(A)
  maxdat = findmax(A)
  xyind = [(maxdat[2]-1)%mn[1] + 1,convert(Int,ceil(maxdat[2]/mn[1]))]
  count = 0;
  V = zeros((2*n+1)^2,1)
  X = zeros((2*n+1)^2,6)
  if (mn[1]-n)>xyind[1]>n && (mn[2]-n)>xyind[2]>n
    for i in -n:n
      for j in -n:n
        count = count + 1
        V[count] = A[xyind[1] + i, xyind[2] + j]
        X[count,:] = [1 i j i*j i*i j*j]
      end
    end
    C = X\V
    denom = 4*C[6]*C[5] - C[4]*C[4]
    fittedpeak = xyind + [C[4]*C[3] - 2.*C[6]*C[2],C[4]*C[2] - 2.*C[5]*C[3]]/denom
  else
    fittedpeak = xyind
  end
  return fittedpeak
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

function RU_poldec(A)
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

function beta_calc_Ruggles(qs,rs, robust_e)
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
  if robust_e
    X = RobustLS(A,b)
  else
    X = A\b
  end

  β = [X[1] X[2] X[3];X[4] X[5] X[6];X[7] X[8] X[9]]
end

function beta_calc_Landon(qs,rs)
  rps = zeros(size(rs))
  for i=1:size(qs,2)
    #display(size([qs[:,i], 0]))
    #display(size(rs[:,i]))
    thisrp = [qs[1,i], qs[2,i],0] + rs[:,i]
    thisrp = thisrp/norm(thisrp)
    rps[:,i] = thisrp
  end
  r1 = rs[1,:]
  r2 = rs[2,:]
  r3 = rs[3,:] #r3 is just a constant times a vector of ones...that's a wee bit silly
  q1 = qs[1,:]
  q2 = qs[2,:]
  q3 = zeros(size(q1))
  rp1 = rps[1,:]
  rp2 = rps[2,:]
  rp3 = rps[3,:]

  A1 = [r1.*rp1.*rp1-r1 r2.*rp1.*rp1-r2 r3.*rp1.*rp1-r3 r1.*rp2.*rp1 r2.*rp2.*rp1 r3.*rp2.*rp1 r1.*rp3.*rp1 r2.*rp3.*rp1 r3.*rp3.*rp1]
  A2 = [r1.*rp1.*rp2 r2.*rp1.*rp2 r3.*rp1.*rp2 r1.*rp2.*rp2-r1 r2.*rp2.*rp2-r2 r3.*rp2.*rp2-r3 r1.*rp3.*rp2 r2.*rp3.*rp2 r3.*rp3.*rp2]
  A3 = [r1.*rp1.*rp3 r2.*rp1.*rp3 r3.*rp1.*rp3 r1.*rp2.*rp3 r2.*rp2.*rp3 r3.*rp2.*rp3 r1.*rp3.*rp3-r1 r2.*rp3.*rp3-r2 r3.*rp3.*rp3-r3]
  b1 = -q1+(q1.*rp1+q2.*rp2+q3.*rp3).*rp1
  b2 = -q2+(q1.*rp1+q2.*rp2+q3.*rp3).*rp2
  b3 = -q3+(q1.*rp1+q2.*rp2+q3.*rp3).*rp3

  b4 = 0
  A4 = [1 0 0 0 1 0 0 0 1]

  A = [A1;A2;A3;A4]
  b = [b1;b2;b3;b4]

  X = A\b
  #X = RobustLS(A,b)
  β = [X[1] X[2] X[3];X[4] X[5] X[6];X[7] X[8] X[9]]
end

#=function beta_calc_FiniteR(qs,rs, v0)
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
end=#

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

  q = findpeak3(CC) - [rrange[1] - srrange[1],crange[1] - scrange[1]] - ROIsize_px/2 - 1

  return q
end

function MeasureShiftNCC(DefI, RefI, ROI_px, shiftROI_px, ROIsize_px; ccsp = 20)
  RS2 = round(Int, ROIsize_px/2)

  ROIC = round.(Int,ROI_px + [0.5;0.5])
  sROIC = round.(Int,shiftROI_px + [0.5;0.5])

  rrange=ROIC[1] - RS2:ROIC[1] - RS2 + ROIsize_px - 1
  crange=ROIC[2] - RS2:ROIC[2] - RS2 + ROIsize_px - 1

  RROI = RefI[rrange,crange]

  q_guess = sROIC - ROIC

  CC = NCC(sROIC[1]-ccsp-RS2:sROIC[1]+ccsp-RS2,sROIC[2]-ccsp-RS2:sROIC[2]+ccsp-RS2,RROI,DefI)

  q = findpeak2(CC) + q_guess - ccsp

  return q
end

function NCC(us,vs,t,f)
  CC = zeros(size(us)[1],size(vs)[1])
  for i=1:size(us)[1]
    for j=1:size(vs)[1]
      CC[i,j] = PointNCC(us[i], vs[j],t, f)
    end
  end
  return CC
end

function PointNCC(u,v,t,f)
  mn = size(t)
  fmuv = mean(f[u+1:u+mn[1],v+1:v+mn[2]])
  tm = mean(t)
  n_gamuv = 0
  df = 0
  dt = 0
  for x=u+1:u+mn[1]
    for y = v+1:v+mn[2]
      n_gamuv += (f[x,y] - fmuv)*(t[x-u,y-v] - tm)
      df += (f[x,y] - fmuv)^2
      dt += (t[x-u,y-v] - tm)^2
    end
  end
  return gamuv = n_gamuv/sqrt(df*dt)
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

  FD = rfft(windowfunc.*DROI)
  FR = rfft(windowfunc.*RROI)#windowfunc.*

  CC = fftshift(brfft(ccfilt[1:round(Int,ROIsize_px/2)+1,:].*FD.*conj(FR),ROIsize_px))#ccfilt[1:round(Int,ROIsize_px/2)+1,:].*

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

function xstar_to_EMsoft(P,m,n)
  pcx = m*(P[1] - .5)
  pcy = m*P[2] - n/2
  pcz = m*P[3]
  P_EM = [pcx,pcy,pcz]
end

function EMsoft_to_image(P, m, n)
  P_i = [P[2]+n/2,P[1]+m/2, P[3]]
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

function get_ideal_qs(rs,PD,PR,F)
    qs = zeros(size(rs))
    for i=1:size(rs)[2]
        r = rs[:,i]
        qi = offset_estimate_w_F(r, PD, PR, F)
        q = image_vec_to_phosphor_frame(qi)
        qs[:,i] = q
    end
    return qs
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
  #I = flattopfilter(I,3.5)
  I = filterimage(I,9,90,25)
end

function flattopfilter(I,c)
  μ = mean(I)
  σ = std(I)
  #I = sqrt(pi)*(c*σ)*erf.((I-μ)/(c*σ))/2 + μ
  I[I.>(μ+c*σ)] = μ+c*σ
  I[I.<(μ-c*σ)] = μ-c*σ
  return I
end

function addnoiseUINT8(I,σ)
  I = convert(Array{Float64},I)
  I += σ*randn(size(I))
  I[I.<0] = 0
  I[I.>1] = 1
  I = round.(I*256)/256
end

function dodgy_prep_ebsp(filepath,σ)
  I = load(filepath)
  I = squareimage(I)
  I = convert(Array{ColorTypes.Gray{FixedPointNumbers.Normed{UInt8,8}},2},I)
  I = convert(Array{Float64},I)
  I = addnoiseUINT8(I,σ)
  I = filterimage(I,9,90,25)
end

function AnnularROI(C, Ri, Ro, m)
  C = round.(Int,C)
  Ro2 = Ro*Ro
  Ri2 = Ri*Ri
  count = 0
  ROI = round.(Int,zeros(ceil(Int,1.1*pi*(Ro2-Ri2))+1,2))
  for i = floor(Int,C[1] - Ro):ceil(Int,C[1] + Ro)
    i2 = (i-C[1])*(i-C[1])
    for j = floor(Int,C[2] - Ro):ceil(Int,C[2] + Ro)
      r2 = (i2 + (j-C[2])*(j-C[2]))
      if Ro2 > r2 > Ri2
        count +=1
        ROI[count,1] = i
        ROI[count,2] = j
      end
    end
  end
  return ROI[1:count,:]
end

function SquareROI(C,S)
  H = S/2
  count = 0
  ROI = round.(Int,zeros(ceil(Int,1.1*S*S)+1,2))
  for i = floor(Int,C[1] - H):ceil(Int,C[1] + H)
    for j = floor(Int,C[2] - H):ceil(Int,C[2] + H)
      if i > floor(Int,C[1] - H) && i < ceil(Int,C[1] + H) && j > floor(Int,C[1] - H) && j < ceil(Int,C[1] + H)
        count +=1
        ROI[count,1] = i
        ROI[count,2] = j
      end
    end
  end
  return ROI[1:count,:]
end

function patternremap(I_coeffs, pad_size, ROI, P, Δ, F, m)
  P_ivec = phosphor_frame_to_image_vec(P+Δ)
  DD = P[3] + Δ[3]
  Δ_ivec = phosphor_frame_to_image_vec(-Δ)
  ΔDD = -Δ[3]
  Fgi = rotate_to_image_frame(inv(F))
  f = EvaluateWarpedImagePatternDistortion(I_coeffs, ROI, mattovec(Fgi), pad_size, P_ivec, DD, Δ_ivec, ΔDD)
  rmI = zeros(m,m)
  if f==[]
      return []
  else
      for i=1:size(ROI,1)
        rmI[ROI[i,1],ROI[i,2]] = f[i]
      end
      return rmI
  end
end

function RobustLS(X,y)
  #tune = 4.685
  tune = 3.
  n = size(y)[1]
  w = ones(n)
  H = X*inv(X'*X)*X'
  h = zeros(n)
  for i=1:n
    h[i] = H[i,i]
  end
  converged = false
  num_iterations = 0
  β_old = zeros(size(X,2))
  β_new = []
  while !converged && num_iterations < 20
    num_iterations += 1
    Xp = diagm(w)*X
    yp = diagm(w)*y
    β_new = Xp\yp
    if norm(β_new-β_old) < 1e-6
      converged = true
    else
      resid = X*β_new - y
      s = median(abs.(resid - median(resid)))/.6745
      r = resid./(tune*s*sqrt.(1-h))
      w = (abs.(r).<1) .* (1 - r.^2).^2
      w = w.*(w.>=.01) + .01*(w.<.01)
      β_old = β_new
    end
  end
  display(num_iterations)
  return β_new
end

function BullseyeROIs(N, rad)
  θs = 2*pi*collect(0:N-2)/(N-1)
  ROIs = Array{Float64}(N,2)
  ROIs[1,:] = [.5 .5]
  for i=2:N
    ROIs[i,:] = [.5 .5] + rad*[sin(θs[i-1]) cos(θs[i-1])]
  end
  return ROIs
end

function AnnularROIs(N, rad; centerpoint = [.5 .5])
  θs = 2*pi*collect(0:N-1)/N
  ROIs = Array{Float64}(N,2)
  for i=1:N
    ROIs[i,:] = centerpoint + rad*[sin(θs[i]) cos(θs[i])]
  end
  return ROIs
end

function GridROIs(N, ROIsize, border)
  S = floor(Int,sqrt(N))
  ROIs = zeros(S*S,2)
  count = 0
  for i=1:S
    for j=1:S
      count += 1
      ROIs[count,1] = border + ROIsize/200 + (i-1)*(1 - 2*border - ROIsize/100)/(S-1)
      ROIs[count,2] = border + ROIsize/200 + (j-1)*(1 - 2*border - ROIsize/100)/(S-1)
    end
  end
  return ROIs
end

function ROIUnion(ROIs,ROIsize,m)
  ROIsize_px = 2*round(Int64,m*ROIsize/200)
  RS2 = round(Int, ROIsize_px/2)
  Union = zeros(Int,m*m,2)
  N = size(ROIs)[1]
  ks = zeros(Int,N)
  count = 0
  skip = 8
  for j=1:skip
    for i=1:floor(Int,N/skip)
      count += 1
      ks[count] = skip*(i-1) + j
    end
  end
  ks[floor(Int,N/skip)*skip+1:N] = floor(Int,N/skip)*skip+1:N

  count = 0
  for i=round(Int,m*minimum(ROIs[:,1])) - RS2 - 1:round(Int,m*maximum(ROIs[:,1])) + RS2 + 1
    for j=round(Int,m*minimum(ROIs[:,2])) - RS2 - 1:round(Int,m*maximum(ROIs[:,2])) + RS2 + 1
      JoinUnion = false
      for q in 1:N
        k = ks[q]
        ROI_px = round.(m*ROIs[k,:] + [.5,.5]) -  [.5,.5]
        ROIC = round.(Int,ROI_px + [0.5;0.5])

        rrange=ROIC[1] - RS2:ROIC[1] - RS2 + ROIsize_px - 1
        crange=ROIC[2] - RS2:ROIC[2] - RS2 + ROIsize_px - 1
        if (i in rrange)&&(j in crange)
          JoinUnion = true
          break
        end
      end
      if JoinUnion
        count+=1
        Union[count,1] = i
        Union[count,2] = j
      end
    end
  end
  return Union[1:count,:]
end

function rotateNx3x3(F,Q; tpose = false)
  rotF = zeros(size(F))
  for i=1:size(F)[1]
    thisF = F[i,:,:]
    if typeof(Q) == Array{Float64,3}
      thisQ = Q[i,:,:]
    elseif typeof(Q) == Array{Float64,2} && size(Q) == (3,3)
      thisQ = Q
    elseif typeof(Q) == Array{Float64,2} && size(Q) == (size(F)[1],3)
      thisQ = euler_to_gmat(Q[i,1],Q[i,2],Q[i,3])
    else
      display("Error in rotateNx3x3, Q in unrecognized format")
    end
    if tpose
      thisQ = thisQ'
    end
    newF = thisQ*thisF*(thisQ')
    rotF[i,:,:] = newF
  end
  return rotF
end

function rotPtoC(P,Qps, eulers)
    S = rotateNx3x3(P,Qps)
    C = rotateNx3x3(S,eulers)
end

function stresscalc(F, C)
    σ = zeros(size(F))
    for i=1:size(F)[1]
        thisF = F[i,:,:]
        E = VR_poldec(thisF)[1] - eye(3)
        Evec = [E[1,1], E[2,2], E[3,3], 2*E[2,3], 2*E[1,3], 2*E[1,2]]
        s = C*Evec
        σ[i,:,:] = [s[1] s[6] s[5];s[6] s[2] s[4];s[5] s[4] s[3]]
    end
    return σ
end

function normNx3x3(F; normtype = "operator2")
  n = zeros(size(F)[1])
  for i=1:size(F)[1]
    thisF = F[i,:,:]
    if normtype == "maxnorm"
      n[i] = maximum(abs.(thisF))
    elseif normtype == "operator2"
      n[i] = norm(thisF)
  elseif normtype == "mises"
        ssum = (thisF[1,1] - thisF[2,2])^2
        ssum += (thisF[3,3] - thisF[2,2])^2
        ssum += (thisF[1,1] - thisF[3,3])^2
        ssum += 6*(thisF[1,2]^2 + thisF[2,3]^2 + thisF[3,1]^2)
      n[i] = sqrt(ssum/2)
    else
      display("Error in norm calculation, normtype not valid")
      n[i] = norm(thisF)
    end
  end
  return n
end

function plotNx3x3Fcomponents(x,F;microstrain = true, minusI = true)
  toplot = []
  if minusI
    toplot = [F[:,1,1]-1 F[:,1,2] F[:,1,3] F[:,2,1] F[:,2,2]-1 F[:,2,3] F[:,3,1] F[:,3,2] F[:,3,3]-1]
  else
    toplot = [F[:,1,1] F[:,1,2] F[:,1,3] F[:,2,1] F[:,2,2] F[:,2,3] F[:,3,1] F[:,3,2] F[:,3,3]]
  end
  if microstrain
    toplot = toplot*1e6
  end
  display(plot(x,toplot))
end

function heatmapNx3x3(F,mn;microstrain = true, minusI = true)
  if minusI
    F = F - eyeNx3x3(size(F)[1])
  end
  if microstrain
    F *= 1e6
  end
  for i=1:3
    for j=1:3
      display(heatmap(reshape(F[:,i,j],mn[1],mn[2])))
    end
  end
end

function convertFtostrain(F)
  E = zeros(size(F))
  for i=1:size(F)[1]
    thisF = F[i,:,:]
    VR = VR_poldec(thisF)
    E[i,:,:] = VR[1] - eye(3)
  end
  return E
end

function misangNx3x3(F)
  misang = zeros(size(F)[1])
  for i=1:size(F)[1]
    thisF = F[i,:,:]
    if thisF != -eye(3)
        VR = VR_poldec(thisF)
        misang[i] = angofF(VR[2])
    end
  end
  return misang
end

function tetragNx3x3(F)
    tetrag = zeros(size(F)[1])
    for i=1:size(F)[1]
      thisF = F[i,:,:]
      thisE = VR_poldec(thisF)[1]
      thisD = eig(thisE)[1]
      tetrag[i] = maximum(thisD) - (sum(thisD) - maximum(thisD))/2
    end
    return tetrag
end

function ffinvNx3x3(F1,F2)
  Fout = zeros(size(F1))
  for i=1:size(F1)[1]
    thisF1 = F1[i,:,:]
    thisF2 = F2[i,:,:]
    Fout[i,:,:] = thisF1*inv(thisF2)
  end
  return Fout
end

function eyeNx3x3(N)
  I = zeros(N,3,3)
  for i=1:N
    I[i,:,:] = eye(3)
  end
  return I
end

function perturbeulers(p1,P,p2,ang,n)
  g0 = euler_to_gmat(p1,P,p2)
  gp = axisang_to_gmat(n,ang)
  g1 = gp*g0
  angs = gmat_to_euler(g1)
end

function BFGSupdate(Bk, gradkp1, gradk, xkp1, xk)
  yk = gradkp1 - gradk
  sk = xkp1 - xk
  Bkp1 = Bk + (yk*yk')/(yk'*sk) - Bk*(sk*sk')*Bk/(sk'*Bk*sk)
  return Bkp1
end
