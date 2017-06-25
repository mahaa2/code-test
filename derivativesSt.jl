function f1(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

    -1/σ^3 * (1 + 1/ν) * (2/ν*z/(zB^2)) * (4/zB - 1)
end

f1(0.0, 1.0, 2.0, 3.0)

function f2(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

    2/(σ^3) * (1 + 1/ν)/(zB^2) * (1 - (1/zB)*(4/ν)*zA)
end

function f3(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

    -1/σ^2 * (-1/ν^2 * (2./(zB).^2 - 1./zB) + (1 + 1/ν)*(4/ν^2*zA/zB^3 - 1/ν^2*zA/zB^2))

end

function f4(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

   -2/σ^3 + 2/σ^3*(1 + 1/ν)/zB*(2*zA + 4*zA/zB - zA^2/(ν*zB) - 4*zA^2/(ν*zB^2))

end

function f5(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

    (1 + 1/ν)*z/(zB^2)*(1/σ^3)*(6 - 8*zA/(ν*zB))
end

function f6(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

   2/(σ^2 * zB^2) * (z/ν^2 - 2*(1 + 1/ν)*z^3/(ν^2 * zB))
end

function f7(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

    1/σ * (zA-1)*zA*(2*zA/(ν^4*zB^3) - 2/(ν^3 * zB^2))

end

function f8(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

    1/σ * (zA-1)*z*(2*zA/(ν^4*zB^3) - 2/(ν^3 * zB^2))
end

function f9(y, µ, σ, ν)
    z  = (y - µ)/σ
    zA = z^2
    zB = 1 + 1/ν * zA

    -1/σ^2 * (-zA/ν^2 * (2/zB^2 + 1/zB) + zA^2*(1 + 1/ν)/ν^2 * (4/zB^3 + 1/zB^2))
end

function f10(y, m, s, v)
    z  = (y - m)/s
    zA = z^2
    zB = 1 + 1/v * zA

    (polygamma(2, (v + 1)/2)/8 - polygamma(2, v/2)/8 - 1/v^3 - 2*zA/(v^3 * zB) + zA^2/(v^4 * zB^2) - 
       zA/2 * ( (-2/v^3 - 3/v^4) * (2 + 1/v*zA)/(zB^2) + (v+1)/v^3 * (2*zA/(v^2*zB^3) + zA/(v^2*zB^2)) ) )
end