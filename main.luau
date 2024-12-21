local function map(t, a, b, x, y)
    local s = (t - a)/(b - a)
    local i = (b - t)/(b - a)
    if s < 0.5 then
        return x + s*(y - x)
    elseif i < 0.5 then
        return y + i*(x - y)
    else
        return (x + y)/2
    end
end


local function chebyshevCoefficients(f, a, b, n)
    local c = {}
    for k = 0, n-1 do
        local s = 0
        for j = 1, n do
            local t = math.pi*(j - 0.5)/n
            local x_j = math.cos(t)
            --TODO: Cache all of this.
            --Addendum: x_(j+1) = x_j cos(pi/(2n)) + sqrt(1 - x_j^2) sin(pi/(2n)).
            s += f(map(x_j, -1, 1, a, b))*math.cos(k*t)
        end
        table.insert(c, 2/n*s)
    end
    return c
end

local function clenshawSumChebyshev(u, c, n)
    --u \in [-1, 1]
    local b_n1, b_n2 = 0, 0
    local b_n
    for k = n, 2, -1 do
        b_n = c[k] + 2*u*b_n1 - b_n2
        b_n2 = b_n1
        b_n1 = b_n
    end
    return c[1]/2 + u*b_n1 - b_n2
end


local function constructChebyshevApproximation(f, a, b, n)
    local c = chebyshevCoefficients(f, a, b, n)
    return c, function (x)
        return clenshawSumChebyshev(map(x, a, b, -1, 1), c, n)
    end
end

local f = function(x)
    return math.sqrt(x)
end
local a, b = 1, 4
local n = 6
local _, approximant = constructChebyshevApproximation(f, a, b, n)

for i=a,b,0.1 do
    print(string.format("x = %.8f, f(x) = %.8f, P(x) = %.8f", i, f(i), approximant(i)))
end
