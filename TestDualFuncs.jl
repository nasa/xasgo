#using LinearAlgebra
#include("DualNumbers.jl")
#using Main.Dual

#A = rand(3,3) + imd*rand(3,3)


function eigenD(A)
    Ar = real(A)
    Ai = imag(A)

    realeig = eigen(Ar,permute=false,scale=false)
    λr = realeig.values
    vr = realeig.vectors

    λi = zeros(size(λr))
    vi = zeros(size(vr))
    N = length(λr)

    if true#length(unique(λr))==length(λr)
        M = zeros(N+1,N+1)
        rhs = zeros(N+1)
        for i=1:N
            M[1:N,1:N] = Ar - I*λr[i]
            M[1:N,N+1] = -vr[:,i]
            M[N+1,1:N] = transpose(vr[:,i])
            rhs[1:N] = -Ai*vr[:,i]
            viandλi = M\rhs
            vi[:,i] = viandλi[1:N]
            λi[i] = viandλi[N+1]
        end
    else
        Anew = copy(A)
        for i=1:20
            qrout = qr(Anew)
            Anew = qrout.R*qrout.Q
        end
        idiag = imag(diag(Anew))
        λi = idiag[sort_to_match(real(diag(Anew)),λr)]
        M = zeros(N+1,N)
        rhs = zeros(N+1)
        for i=1:N
            M = [Ar - λr[i]*I;transpose(vr[:,i])]
            rhs[1:N] = (λi[i]*I-Ai)*vr[:,i]
            thisvi = M\rhs
            vi[:,i] = thisvi
        end

    end
    v = vr + vi*imd
    λ = λr + λi*imd
    return Eigen(λ,v)
end

function sort_to_match(A,B; tol = .001)
    N = length(A)
    permvec = convert(Array{Int64,1},1:N)
    for i=1:N
        for j=1:N
            if abs(B[j] - A[i])/(abs(A[i])+.0000001) < tol
                toreplaceindex = findnext(permvec.==j,1)
                permvec[toreplaceindex] = permvec[i]
                permvec[i] = j
                break
            end
        end
    end
    return permvec
end
