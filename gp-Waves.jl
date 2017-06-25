function aniGP(d::Int64, n::Int64)
 # Description :

 # Returns a matrix of size d x n, representing a grand circle 
 # on the unit-sphere in n steps, starting at random location
 # Given a kernel (covariance matrix of domain location),
 # this can be turnedinto a tour through the sample space, 
 # simply by calling chol(K)' * K.

 # Phillip Hennig, September 2012
 
 x = randn(d)                         # starting sample
 r = sqrt(sum(x.^2))
 x ./= r                              # projection onto sphere
 
 t = randn(d)                         # sample tangent direction
 t .-= (t' * x)[] * x                 # orthogonalise by Gram-Schmidt
 t ./= sqrt(sum(t.^2))                # standardise

 s = Vector(linspace(0, 2*π, n + 1)); # space to span
 s = s[1:end - 1]';                   # span linspace in direction of t
 t = broadcast(*, s, t)               # projection onto sphere, re-scale

 return (r .* ExpMap(x, t))
end

function ExpMap(µ::Vector{Float64}, E::Matrix{Float64})
 # Description :

 # Computes the exponential-map on the sphere
 # Acknowledgement to Soren Hauberg
 # Phillip Hennig, September 2012
 
 D = size(E, 1)
 θ = sqrt(sum(E.^2, 1))
 M = µ * cos(θ) + E .* repmat(sin(θ)./θ, D, 1)

 # M[:, (abs(θ) .<= 1e-7)[:]] = µ ...
 if any(abs(θ) .<= 1e-7)
    for a in find(abs(θ) .<= 1e-7)
        M[:, a] = µ
    end
 end
 
 return(M)
end
