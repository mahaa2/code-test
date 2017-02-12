# The bivariate Ricker function
function ricker(N::Vector{Float64}, B::Matrix{Float64}, α::Vector{Float64})
    a = B * N;
    Nt = α .* N .* exp(-a)

 return(Nt)
end

# carrying capacity and interspecific parameters
B = [0.001 0.0005; 0.0008 0.001];

# intrinsic growth rate 
α = [2.0; 2.0];

# equilibrium population (is it stable ?)
Nhat = inv(B) * log(α)

# function handle
R(N::Vector{Float64}) = ricker(N, B, α);

# test population growth
Npred = [1000.0 1.0]

for i in 1:50
    Npred = [Npred; R(Npred[end, :])']
end

# visualize
set_default_plot_format(:png)

plot(layer(x = collect(1:size(Npred, 1)), y = Npred[:, 1], Geom.point, Geom.line, Theme(default_color = color("orange"))),
     layer(x = collect(1:size(Npred, 1)), y = Npred[:, 2], Geom.point, Geom.line), Theme(background_color = color("White")))
