#Notices:
#Copyright 2018 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. All Rights Reserved.
#
#Disclaimers
#No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS." 
# 
#Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.



using Images, Colors, FixedPointNumbers, Plots
using FFTW
#using JuMP
#using NLopt
#using Ipopt

include("XASGOWrapper.jl")
plotly()


α = pi/2 - 7*pi/18 + pi/18
Qps=[0 cos(α) -sin(α);
        -1     0            0;
        0   -sin(α) -cos(α)]


#=RefI = prep_ebsp(joinpath("i02-H1 8-30_03_Images","i02-H1 8-30_03_001.jpg"))
PR = [0.522339027595;0.92900;0.6057818659658]
gR = euler_to_gmat(5.461484295+pi/2,0.043380159,0.007648033)
m = minimum(size(RefI))=#

RefI = prep_ebsp(joinpath("al_ebsd","set4",string("ebsd_",0,".png")))
PR = [.5 .5 .625]
#gR = euler_to_gmat(pi/12,	pi/18,	0)
gR = euler_to_gmat(0,0,0)
m = minimum(size(RefI))

NP = 41
PatternInfo = zeros(NP,6)
include("LongEulers.jl")
#include("ShortEulers.jl")
Patterns = String[]
for i=1:NP
  push!(Patterns, joinpath("al_ebsd","set5",string("ebsd_",i-1,".png")))
  p1 = Pangles[i,1]
  P = Pangles[i,2]
  p2 = Pangles[i,3]

  #n = rand(3)
  #n = n/norm(n)
  #angs = perturbeulers(p1,P,p2,pi/900,n)
  #p1 = angs[1]
  #P = angs[2]
  #p2 = angs[3]



  PatternInfo[i,:] = [p1 P p2 .5 .5 .625]
end

#display(heatmap(RefI))

#=
N = 48
ROIsize = 15
rad = .2689
ROIinfo = [.5,.5,rad-(ROIsize/100)/sqrt(2),rad+(ROIsize/100)/sqrt(2)]#rad+(ROIsize/100)/sqrt(2)
=#


N = 100
ROIsize = 15
rad = .264
ROIinfo = [.5,.5,.14287,.355]


ROIs = AnnularROIs(N,rad)
WarpROI = AnnularROI([m,m]/2,m*(rad-(ROIsize/100)/sqrt(2))-3, m*(rad+(ROIsize/100)/sqrt(2))+1, m)

tic()
F = GrainFCalcICGN(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIinfo = ROIinfo, Hupdate = false)
toc()
tic()
FR = GrainFCalcCC(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIs = ROIs, ROIsize = ROIsize, FAlgorithm = "remap", WarpROI = WarpROI)
toc()
tic()
FN = GrainFCalcCC(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIs = ROIs, ROIsize = ROIsize, FAlgorithm = "ncorrshift")
toc()
tic()
FR1 = GrainFCalcCC(RefI, PR, gR, Patterns, PatternInfo, Qps, ROIs = ROIs, ROIsize = ROIsize, FAlgorithm = "remapSingle", WarpROI = WarpROI)
toc()




display("Done")
FT = zeros(size(F))
FG = zeros(size(F))
Pseudo = eye(3) + [-2310 0 0;0 2640 0;0 0 -330]/1e6
for i=1:size(F)[1]
  gD = euler_to_gmat(Pangles[i,1], Pangles[i,2], Pangles[i,3])
  thisFT = Qps'*gD'*Pseudo*gR*Qps
  FT[i,:,:] = thisFT

  gG = euler_to_gmat(PatternInfo[i,1], PatternInfo[i,2], PatternInfo[i,3])
  thisFG = Qps'*gG'*gR*Qps
  FG[i,:,:] = thisFG

end


#plotNx3x3Fcomponents(0:.3:(NP-1)*.3,convertFtostrain(rotateNx3x3(F - FT,eye(3))+eyeNx3x3(NP)),microstrain = true, minusI = false)

#plotNx3x3Fcomponents(0:.3:(NP-1)*.3,F - FT,microstrain = true, minusI = false)
#plotNx3x3Fcomponents(0:.3:(NP-1)*.3,FR - FT,microstrain = true, minusI = false)


#display(plot(0:.3:(NP-1)*.3, 1e6*[normNx3x3(convertFtostrain(F)) normNx3x3(convertFtostrain(FR))]))
#display(plot(0:.3:(NP-1)*.3, 1e6*[normNx3x3(F-FT) normNx3x3(FR-FT)]))
#display(plot(0:.3:(NP-1)*.3, [misangNx3x3(ffinvNx3x3(F,FT)) misangNx3x3(ffinvNx3x3(FR,FT))]))



#display(plot(misangNx3x3(ffinvNx3x3(FN,FT))))

#EF = normNx3x3(convertFtostrain(F))
#ER = normNx3x3(convertFtostrain(FR))
#EN = normNx3x3(convertFtostrain(FN))
#E1 = normNx3x3(convertFtostrain(FR1))
EF = normNx3x3(convertFtostrain(ffinvNx3x3(F,FT)))
ER = normNx3x3(convertFtostrain(ffinvNx3x3(FR,FT)))
EN = normNx3x3(convertFtostrain(ffinvNx3x3(FN,FT)))
E1 = normNx3x3(convertFtostrain(ffinvNx3x3(FR1,FT)))
display(plot(0:.3:(NP-1)*.3, 1e6*[EN E1 ER EF]))

MF = misangNx3x3(ffinvNx3x3(F,FT))
MR = misangNx3x3(ffinvNx3x3(FR,FT))
MN = misangNx3x3(ffinvNx3x3(FN,FT))
M1 = misangNx3x3(ffinvNx3x3(FR1,FT))
display(plot(0:.3:(NP-1)*.3, [MN M1 MR MF]))


plotNx3x3Fcomponents(0:.3:(NP-1)*.3,ffinvNx3x3(F,FT),microstrain = true, minusI = true)
plotNx3x3Fcomponents(0:.3:(NP-1)*.3,ffinvNx3x3(FR,FT),microstrain = true, minusI = true)

relF_p = ffinvNx3x3(F,FT)
relF_s = rotateNx3x3(relF_p,Qps,tpose = false)
relF_c = rotateNx3x3(relF_s,Pangles[1:NP,:])
plotNx3x3Fcomponents(0:.3:(NP-1)*.3,relF_c,microstrain = true, minusI = true)

relFR_p = ffinvNx3x3(FR,FT)
relFR_s = rotateNx3x3(relFR_p,Qps,tpose = false)
relFR_c = rotateNx3x3(relFR_s,Pangles[1:NP,:])
plotNx3x3Fcomponents(0:.3:(NP-1)*.3,relFR_c,microstrain = true, minusI = true)

relFR1_p = ffinvNx3x3(FR1,FT)
relFR1_s = rotateNx3x3(relFR1_p,Qps,tpose = false)
relFR1_c = rotateNx3x3(relFR1_s,Pangles[1:NP,:])
plotNx3x3Fcomponents(0:.3:(NP-1)*.3,relFR1_c,microstrain = true, minusI = true)

relFN_p = ffinvNx3x3(FN,FT)
relFN_s = rotateNx3x3(relFN_p,Qps,tpose = false)
relFN_c = rotateNx3x3(relFN_s,Pangles[1:NP,:])
plotNx3x3Fcomponents(0:.3:(NP-1)*.3,relFN_c,microstrain = true, minusI = true)

@save "Set5.jld"
