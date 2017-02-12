# ~ implementation in julia-language

# · packages and my configurations
using Cairo, Gadfly       
set_default_plot_format(:png);


# · prior on θ (parameters or weights)
F    = 5;
ø(a) = (broadcast(.>, a', 0:F-1) .* broadcast(.<=, a', 1:F))'     
# ø(a) = broadcast(^, a', 0:F-1)';     # phi function ø(a) = [1 a]^T
µ    = zeros(F, 1);                    # prior mean
Σ    = 5.*eye(F);                         # prior variance, π(θ) = N(θ|µ, Σ)

# · prior on f(x; θ)
np    = 100; xp = Vector(linspace(-6, 6, np));       # domain points to predict function values
phixp = ø(xp);                                       # ø(·) function
µpr   = phixp * µ;
kxpp  = phixp * Σ * phixp';                                               # unconditional distribution w.r.t θ, π(f(xp)) = N(f(xp)|µpr, kxpp))
s     = broadcast(+, µpr, chol(kxpp + 1e-08 * eye(np))' * rand(np, 3));   # sampling from the distribution above
stdpr = sqrt(diag(kxpp));                                                 # standart-deviation from the prior

# · data
# x  = [-5.0, -4.0, -1.0,  0.0, 2.0, 3.0, 4.0];
x  = [ 0.1,  0.5,  1.5,  1.7, 2.5, 4.1, 4.9];
y  = [-5.5, -5.4, -4.0, -1.0, 0.0, 1.0, 3.4]; 
n  = length(y);
σ² = 1;

# · the posterior
phix = ø(x);                            # ø(·) function
M    = phix * µ;
kxx  = phix * Σ * phix';                # unconditional distribution w.r.t θ, π(f(x)|M, kxx)

G = kxx + σ² * eye(n);                  # unconditional distribution w.r.t θ, π(y|M, kxx + Iσ²)
R = chol(G);                            # most expensive step: O(n³)

kxpx = phixp * Σ * phix';               # cov(f(x), f(xp)) = Kxxp
A    = kxpx / R;                        # pre-compute for re-use

µpost = µpr + A * (R' \ (y - M));       # π(f(x)|y) = N(f(x)|µpr + kxpx (kxx + Iσ²)-¹ (y - M)
Σpost = kxpp - A * A';                  # kxpxp - kxpx (kxx + Iσ²)-¹ kxxp)

spost = broadcast(+, µpost, chol(Σpost + 1e-08*eye(np))' * rand(np, 3));    # samples
stdpo = sqrt(diag(Σpost));                                                  # standart-deviation from the posterior