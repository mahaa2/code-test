nF = 30;      # number of frames
ns = 1;       # number of sample-frames
# np = 50;     # np-dimensional Gaussian

s = Array{Array{Float64, 2}}(ns);   # Array of arrays 
for i in 1:ns                      # create animations
    s[i] = aniGP(np, nF);
end

# 2-dimensional domain example
# kSE(a, b) = exp(-sum((repeat(reshape(a, size(a, 1), 1, size(a, 2)), outer = [1, size(b, 1), 1]) - 
#                       repeat(reshape(b, 1, size(b, 1), size(b, 2)), outer = [size(a, 1), 1, 1])
#                       ).^2, 3)[:, :]);

# x1 = x2 = Vector(linspace(-5, 5, 20));
# Δ = [repeat(x1, inner= [size(x2, 1)]) repeat(x2, outer = [size(x1, 1)])];

# K = kSE(Δ, Δ);
# L = chol(K + eye(np) .* 1e-08)';

# initialize a 3D plot with 1 empty series
# plt = path3d(1, xlim=(x1[1], x1[end]), ylim=(x2[1], x2[end]), zlim=(-2, 2),
#                 xlab = "x", ylab = "y", zlab = "z",
#                 title = "Gaussian process waves", marker = 1)

# # build an animated gif, saving every 10th frame
# @gif for i in 1:nF
#     z = L * s[1][:, i]
#     push!(plt, Δ[:, 1], Δ[:, 2], z)
# end every 1

cols = ["orange", "purple", "red"] # colors

# µpost, Σpost
L = chol(Σpost + eye(np) .* 1e-08)';

# Gadfly.plot([layer(x = xp, y = µpost + L * s1[:, f], Geom.line) for f in 1:2] ..., Theme(background_color = color("White")))

for j in 1:nF
    p = plot([layer(x = xp, y = µpost + L * s[f][:, j], Geom.line, Theme(default_color = color(cols[f]))) for f in 1:ns] ...,
              layer(x = x, y = y, Geom.point), Coord.Cartesian(ymin = minimum(y)-2, ymax = maximum(y)+2),
              Theme(background_color = color("White")));

    draw(PNG(string(j, ".png"), 100, 100, dpi = 200) , p);
end

# creat gif in gnu-linux terminal
# convert -delay 10 -loop 0 *.png animation.gif