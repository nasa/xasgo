module BiDual
#Constructors
struct BiDualNumber{T<:Real} <: Number
  rpart::T
  i1part::T
  i2part::T
  i12part::T
end
BiDualNumber(x::Real,y::Real,z::Real,w::Real) = BiDualNumber(promote(x,y,z,w)...)
BiDualNumber(x::Real) = BiDualNumber(x,zero(x),zero(x),zero(x))

const im1 = BiDualNumber(0,1,0,0)
const im2 = BiDualNumber(0,0,1,0)
const im12 = BiDualNumber(0,0,0,1)


#Type conversion
BiDualNumber{T}(x::Real) where {T<:Real} = BiDualNumber{T}(x,0,0,0)
BiDualNumber{T}(z::BiDualNumber) where {T<:Real} = BiDualNumber{T}(z.rpart,z.i1part,z.i2part,z.i12part)

BiDualNumber(z::BiDualNumber) = z

(::Type{T})(z::BiDualNumber) where {T<:Real} =
isreal(z) ? T(real(z))::T : throw(InexactError(Symbol(string(T)), T, z))

import Base: promote_rule, widen, float, real
import Base: isreal, isinteger, isfinite, isnan, isinf, iszero, isone
import Base: flipsign, show, read, write, bswap, in
import Base: ==, isequal, +, -, *, /, conj, sin, cos, hash

promote_rule(::Type{BiDualNumber{T}}, ::Type{S}) where {T<:Real,S<:Real} =
    BiDualNumber{promote_type(T,S)}
promote_rule(::Type{BiDualNumber{T}}, ::Type{BiDualNumber{S}}) where {T<:Real,S<:Real} =
    BiDualNumber{promote_type(T,S)}

widen(::Type{BiDualNumber{T}}) where {T} = BiDualNumber{widen(T)}

float(::Type{BiDualNumber{T}}) where {T<:AbstractFloat} = BiDualNumber{T}
float(::Type{BiDualNumber{T}}) where {T} = BiDualNumber{float(T)}

#Add methods to various common complex operators


real(x::BiDualNumber) = x.rpart
real(::Type{BiDualNumber{T}}) where {T<:Real} = T
imag1(z::BiDualNumber) = z.i1part
imag2(z::BiDualNumber) = z.i2part
imag12(z::BiDualNumber) = z.i12part

isreal(z::BiDualNumber{<:Real}) = iszero(z.i1part)&iszero(z.i2part)&iszero(z.i12part)
isinteger(z::BiDualNumber) = isreal(z) & isinteger(real(z))
isfinite(z::BiDualNumber) = isfinite(real(z)) & isfinite(z.i1part)& isfinite(z.i2part)& isfinite(z.i12part)
isnan(z::BiDualNumber) = isnan(real(z)) | isnan(z.i1part)| isnan(z.i2part)| isnan(z.i12part)
isinf(z::BiDualNumber) = isinf(real(z)) | isinf(z.i1part)| isinf(z.i2part)| isinf(z.i12part)
iszero(z::BiDualNumber) = iszero(real(z)) & isreal(z)
isone(z::BiDualNumber) = isone(real(z)) & isreal(z)

bidualify(z::BiDualNumber) = z
bidualify(x::Real) = BiDualNumber(x)
bidualify(x::Real,y::Real,z::Real,w::Real) = BiDualNumber(x,y,z,w)

bidualify(A::AbstractArray{<:Complex}) = A

function bidualify(A::AbstractArray{T}) where T
    if !isconcretetype(T)
        error("`bidualify` not defined on abstractly-typed arrays; please convert to a more specific type")
    end
    convert(AbstractArray{typeof(bidualify(zero(T)))}, A)
end

bidualify(::Type{T}) where {T<:Real} = BiDualNumber{T}
bidualify(::Type{BiDualNumber{T}}) where {T<:Real} = BiDualNumber{T}

flipsign(x::BiDualNumber, y::Real) = ifelse(signbit(y), -x, x)

function show(io::IO, z::BiDualNumber)
    r = real(z)
    i1 = imag1(z)
    i2 = imag2(z)
    i12 = imag12(z)
    compact = get(io, :compact, false)
    show(io, r)
    if signbit(i1) && !isnan(i1)
        i1 = -i1
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    show(io, i1)
    print(io, "*")
    print(io, "i1")
    if signbit(i2) && !isnan(i2)
        i2 = -i2
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    show(io, i2)
    print(io, "*")
    print(io, "i2")
    if signbit(i12) && !isnan(i12)
        i12 = -i12
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    show(io, i12)
    print(io, "*")
    print(io, "i1*i2")
end

function show_unquoted(io::IO, z::BiDualNumber, ::Int, prec::Int)
    if operator_precedence(:+) <= prec
        print(io, "(")
        show(io, z)
        print(io, ")")
    else
        show(io, z)
    end
end

==(z::BiDualNumber, w::BiDualNumber) = (real(z) == real(w)) & (imag1(z) == imag1(w)) & (imag2(z) == imag2(w)) & (imag12(z) == imag12(w))
==(z::BiDualNumber, x::Real) = isreal(z) && real(z) == x
==(x::Real, z::BiDualNumber) = isreal(z) && real(z) == x


isequal(z::BiDualNumber, w::BiDualNumber) = isequal(real(z),real(w)) & isequal(imag1(z),imag1(w)) & isequal(imag2(z),imag2(w)) & isequal(imag12(z),imag12(w))

in(x::BiDualNumber, r::AbstractRange{<:Real}) = isreal(x) && real(x) in r

#Basic functions
conj(z::BiDualNumber) = BiDualNumber(real(z),-imag1(z),-imag2(z),-imag12(z))#???????
abs(z::BiDualNumber)  = abs(z.rpart)

+(x::BiDualNumber,y::BiDualNumber) = BiDualNumber(x.rpart+y.rpart,x.i1part+y.i1part,x.i2part+y.i2part,x.i12part+y.i12part)
+(x::BiDualNumber,y::Real) = BiDualNumber(x.rpart+y,x.i1part,x.i2part,x.i12part)
+(x::Real,y::BiDualNumber)=BiDualNumber(x+y.rpart,y.i1part,y.i2part,y.i12part)

-(x::BiDualNumber)=BiDualNumber(-x.rpart,-x.i1part,-x.i2part,-x.i12part)
-(x::BiDualNumber,y::BiDualNumber)=BiDualNumber(x.rpart-y.rpart,x.i1part-y.i1part,x.i2part-y.i2part,x.i12part-y.i12part)
-(x::BiDualNumber,y::Real)=BiDualNumber(x.rpart-y,x.i1part,x.i2part,x.i12part)
-(x::Real,y::BiDualNumber)=BiDualNumber(x-y.rpart,-y.i1part,-y.i2part,-y.i12part)

function *(x::BiDualNumber,y::BiDualNumber)
    a = x.rpart
    b = x.i1part
    c = x.i2part
    d = x.i12part
    e = y.rpart
    f = y.i1part
    g = y.i2part
    h = y.i12part
    rpart = a*e
    i1part = a*f+e*b
    i2part = a*g+e*c
    i12part = a*h + b*g + c*f + d*e
    BiDualNumber(rpart,i1part,i2part,i12part)
end
function *(x::BiDualNumber,y::Real)
    a = x.rpart
    b = x.i1part
    c = x.i2part
    d = x.i12part
    BiDualNumber(a*y,b*y,c*y,d*y)
end
function *(y::Real,x::BiDualNumber)
    a = x.rpart
    b = x.i1part
    c = x.i2part
    d = x.i12part
    BiDualNumber(a*y,b*y,c*y,d*y)
end
#Their division is way better than mine
function /(x::BiDualNumber,y::BiDualNumber)
    a = x.rpart
    b = x.i1part
    c = x.i2part
    d = x.i12part
    e = y.rpart
    f = y.i1part
    g = y.i2part
    h = y.i12part
    e2 = e*e
    rpart = a/e
    i1part = (e*b-a*f)/e2
    i2part = (e*c - a*g)/e2
    i12part = (d*e - a*h - b*g - c*f)/e2 + 2*a*g*f/(e2*e)
    BiDualNumber(rpart,i1part,i2part,i12part)
end
function /(x::BiDualNumber,y::Real)
    a = x.rpart
    b = x.i1part
    c = x.i2part
    d = x.i12part
    rpart = a/y
    i1part = b/y
    i2part = c/y
    i12part = d/y
    BiDualNumber(rpart,i1part,i2part,i12part)
end
function /(x::Real,y::BiDualNumber)
    a = x
    e = y.rpart
    f = y.i1part
    g = y.i2part
    h = y.i12part
    e2 = e*e
    rpart = a/e
    i1part = -a*f/e2
    i2part = -a*g/e2
    i12part = -a*h/e2 + 2*a*g*f/(e2*e)
    BiDualNumber(rpart,i1part,i2part,i12part)
end

export BiDualNumber, bidualify, imag1, imag2, imag12, im1, im2, im12

end
