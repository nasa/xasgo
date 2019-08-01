module Dual

#Constructors
struct DualNumber{T<:Real} <: Number
  rpart::T
  ipart::T
end
DualNumber(x::Real,y::Real) = DualNumber(promote(x,y)...)
DualNumber(x::Real) = DualNumber(x,zero(x))

const imd = DualNumber(0,1)


#Type conversion
DualNumber{T}(x::Real) where {T<:Real} = DualNumber{T}(x,0)
DualNumber{T}(z::Complex) where {T<:Real} = DualNumber{T}(real(z),imag(z))
DualNumber{T}(z::DualNumber) where {T<:Real} = DualNumber{T}(real(z),imag(z))
Complex{T}(z::DualNumber) where {T<:Real} = Complex{T}(real(z),imag(z))

DualNumber(z::DualNumber) = z

(::Type{T})(z::DualNumber) where {T<:Real} =
isreal(z) ? T(real(z))::T : throw(InexactError(Symbol(string(T)), T, z))

import Base.promote_rule
import Base.widen
import Base.float
promote_rule(::Type{DualNumber{T}}, ::Type{S}) where {T<:Real,S<:Real} =
    DualNumber{promote_type(T,S)}
promote_rule(::Type{DualNumber{T}}, ::Type{DualNumber{S}}) where {T<:Real,S<:Real} =
    DualNumber{promote_type(T,S)}

widen(::Type{DualNumber{T}}) where {T} = DualNumber{widen(T)}

float(::Type{DualNumber{T}}) where {T<:AbstractFloat} = DualNumber{T}
float(::Type{DualNumber{T}}) where {T} = DualNumber{float(T)}

#Add methods to various common complex operators
import Base.real
import Base.imag
import Base.isreal
import Base.isinteger
import Base.isfinite
import Base.isnan
import Base.isinf
import Base.iszero
import Base.isone
import Base.flipsign
import Base.show
import Base.read
import Base.write
import Base.bswap
import Base.in
import Base.==
import Base.isequal
import Base.+
import Base.-
import Base.*
import Base./
import Base.^
import Base.conj
import Base.abs
import Base.abs2
import Base.inv
import Base.sin
import Base.cos
import Base.hash
import Base.sqrt


real(x::DualNumber) = x.rpart
real(::Type{DualNumber{T}}) where {T<:Real} = T
imag(x::DualNumber) = x.ipart

isreal(z::DualNumber{<:Real}) = iszero(imag(z))
isinteger(z::DualNumber) = isreal(z) & isinteger(real(z))
isfinite(z::DualNumber) = isfinite(real(z)) & isfinite(imag(z))
isnan(z::DualNumber) = isnan(real(z)) | isnan(imag(z))
isinf(z::DualNumber) = isinf(real(z)) | isinf(imag(z))
iszero(z::DualNumber) = iszero(real(z)) & iszero(imag(z))
isone(z::DualNumber) = isone(real(z)) & iszero(imag(z))

dualify(z::DualNumber) = z
dualify(x::Real) = DualNumber(x)
dualify(x::Real,y::Real) = DualNumber(x,y)

dualify(A::AbstractArray{<:DualNumber}) = A

function dualify(A::AbstractArray{T}) where T
    if !isconcretetype(T)
        error("`dualify` not defined on abstractly-typed arrays; please convert to a more specific type")
    end
    convert(AbstractArray{typeof(dualify(zero(T)))}, A)
end

dualify(::Type{T}) where {T<:Real} = DualNumber{T}
dualify(::Type{DualNumber{T}}) where {T<:Real} = DualNumber{T}

flipsign(x::DualNumber, y::Real) = ifelse(signbit(y), -x, x)

function show(io::IO, z::DualNumber)
    r, i = reim(z)
    compact = get(io, :compact, false)
    show(io, r)
    if signbit(i) && !isnan(i)
        i = -i
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    show(io, i)
    if !(isa(i,Integer) && !isa(i,Bool) || isa(i,AbstractFloat) && isfinite(i))
        print(io, "*")
    end
    print(io, "ϵ")
end
show(io::IO, z::DualNumber{Bool}) =
    print(io, z == im ? "ϵ" : "DualNumber($(z.re),$(z.im))")

function show_unquoted(io::IO, z::DualNumber, ::Int, prec::Int)
    if operator_precedence(:+) <= prec
        print(io, "(")
        show(io, z)
        print(io, ")")
    else
        show(io, z)
    end
end

function read(s::IO, ::Type{DualNumber{T}}) where T<:Real
    r = read(s,T)
    i = read(s,T)
    DualNumber{T}(r,i)
end
function write(s::IO, z::DualNumber)
    write(s,real(z),imag(z))
end

bswap(z::DualNumber) = DualNumber(bswap(real(z)), bswap(imag(z)))

==(z::DualNumber, w::DualNumber) = (real(z) == real(w)) & (imag(z) == imag(w))
==(z::DualNumber, x::Real) = isreal(z) && real(z) == x
==(x::Real, z::DualNumber) = isreal(z) && real(z) == x

isequal(z::DualNumber, w::DualNumber) = isequal(real(z),real(w)) & isequal(imag(z),imag(w))

in(x::DualNumber, r::AbstractRange{<:Real}) = isreal(x) && real(x) in r

if UInt === UInt64
    const h_imag = 0x32a7a07f3e7cd1f9
else
    const h_imag = 0x3e7cd1f9
end
const hash_0_imag = hash(0, h_imag)

function hash(z::DualNumber, h::UInt)
    # TODO: with default argument specialization, this would be better:
    # hash(real(z), h ⊻ hash(imag(z), h ⊻ h_imag) ⊻ hash(0, h ⊻ h_imag))
    hash(real(z), h ⊻ hash(imag(z), h_imag) ⊻ hash_0_imag)
end

#Basic functions
conj(z::DualNumber) = DualNumber(real(z),-imag(z))
abs(z::DualNumber)  = abs(z.rpart)
abs2(z::DualNumber) = real(z)*real(z)
inv(z::DualNumber)  = DualNumber(inv(z.rpart),z.ipart/(z.rpart*z.rpart))
#inv(z::Complex{<:Integer}) = inv(float(z))

+(x::DualNumber,y::DualNumber) = DualNumber(x.rpart+y.rpart,x.ipart+y.ipart)
+(x::DualNumber,y::Real) = DualNumber(x.rpart+y,x.ipart)
+(x::Real,y::DualNumber)=DualNumber(x+y.rpart,y.ipart)

-(x::DualNumber)=DualNumber(-x.rpart,-x.ipart)
-(x::DualNumber,y::DualNumber)=DualNumber(x.rpart-y.rpart,x.ipart-y.ipart)
-(x::DualNumber,y::Real)=DualNumber(x.rpart-y,x.ipart)
-(x::Real,y::DualNumber)=DualNumber(x-y.rpart,-y.ipart)

function *(x::DualNumber,y::DualNumber)
    rpart = x.rpart*y.rpart
    ipart = x.rpart*y.ipart + x.ipart*y.rpart
    DualNumber(rpart,ipart)
end
function *(x::DualNumber,y::Real)
    rpart = x.rpart*y
    ipart = x.ipart*y
    DualNumber(rpart,ipart)
end
function *(x::Real,y::DualNumber)
    rpart = x*y.rpart
    ipart = x*y.ipart
    DualNumber(rpart,ipart)
end
#Their division is way better than mine
function /(x::DualNumber,y::DualNumber)
    rpart = x.rpart/y.rpart
    ipart = (x.ipart*y.rpart - x.rpart*y.ipart)/(y.rpart*y.rpart)
    DualNumber(rpart,ipart)
end
function /(x::DualNumber,y::Real)
    rpart = x.rpart/y
    ipart = x.ipart/y
    DualNumber(rpart,ipart)
end
function /(x::Real,y::DualNumber)
    rpart = x/y.rpart
    ipart = (-x*y.ipart)/(y.rpart*y.rpart)
    DualNumber(rpart,ipart)
end

#Did not include all the Boolean stuff from complex.jl

#Need to add in muladd

sin(x::DualNumber)=DualNumber(sin(x.rpart),x.ipart*cos(x.rpart))
cos(x::DualNumber)=DualNumber(cos(x.rpart),-x.ipart*sin(x.rpart))

function sqrt(x::DualNumber)
    sqrtreal = sqrt(real(x))
    return DualNumber(sqrtreal,imag(x)/(2*sqrtreal))
end

function ^(x::DualNumber, p::Real)
    xpr = real(x)^p
    return DualNumber(xpr,p*imag(x)*xpr/real(x))
end

function ^(x::DualNumber, p::Integer)
    xpr = real(x)^p
    return DualNumber(xpr,p*imag(x)*xpr/real(x))
end

using LinearAlgebra
import LinearAlgebra.eigen!
import LinearAlgebra.eigen
import LinearAlgebra.svd

function eigen(A::Array{DualNumber{T},2}) where T <: Real
    #Only works for square matrices
    #No idea what's going on when Ar has complex eigenvalues
    Ar = real(A)
    Ai = imag(A)

    realeig = eigen(Ar,permute=false,scale=false)
    λr = realeig.values
    vr = realeig.vectors
    λi = zeros(size(λr))
    vi = zeros(size(vr))
    N = length(λr)
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
    v = vr + vi*imd
    λ = λr + λi*imd
    return Eigen(λ,v)
end

function svd(A::Array{DualNumber{T},2}) where T <: Real
    m,n = size(A)
    p = min(m,n)
    AATeig = eigen(A*transpose(A))
    u = AATeig.vectors
    sv = sqrt.(AATeig.values)
    v = dualify(zeros(n,n))
    for i=1:p
        thisv = transpose(A)*u[:,i]/sv[i]
        v[:,i] = thisv
    end
    return SVD(u,sv[1:p],v)
end


export DualNumber, dualify, imd

end
