using Plots, Distributions

# include main code
include("../src/functions.jl")

# definig distribution functions 
function FR12(u, b1, b2)
    lambda1 = 0.25
    lambda2 = 0.5
    return (-b1 * lambda2 * exp(-u * lambda1 / b1) + b2 * lambda1 * exp(-lambda2 * u / b2) + b1 * lambda2 - b2 * lambda1) / (b1 * lambda2 - b2 * lambda1)
end

F1 = SampleDistFun(x -> cdf(Gamma(3, 0.8), x), 40.0, 0.005)
F2(x) = 1 - exp(-0.7071067812 * sqrt(x))
FR1(x, a1) = F1(x / a1)
FR2(x, a2) = F2(x / a2)

FR(r, x) =
    0.3333333333333333 * FR1(x, r[1]) +
    0.4 * FR2(x, r[2]) +
    0.26666666666666666 * FR12(x, r[1], r[2])


ER1(a1) = 2.4 * a1
ER2(a2) = 4 * a2
ER12(a1, a2) = 4 * a1 + 2 * a2
ER(r) =
    0.3333333333333333 * ER1(r[1]) +
    0.4 * ER2(r[2]) +
    0.26666666666666666 * ER12(r[1], r[2])


res = bisection_hjbsolver(
    a=166.0,
    b=333.0,
    max_it=20,
    NRe=2,
    β=5.0 + 6.0 + 4.0,
    FR=FR,
    ER=ER,
    EnoR=ER([1.0, 1.0]),
    partition=(0.0, 19.0, 200.0, 500.0),
    Δx=0.1,
    SearchRanges=([(0.0, 1.0), (0.0, 1.0)], [(0.0, 1.0), (1.0, 1.0)]),
    η=0.25,
    η1=0.3,
    δ=0.15,
    MaxEvals=1500,
    p=120,
    L=300
)

# reinsurance strategy line 1
plot(res[end][1], res[end][4][:, 1], size=(600, 300), xlims=(0, 70))
# line 2
plot(res[end][1], res[end][4][:, 2], size=(600, 300), xlims=(0, 5), ylims=(0, 1))

# value function
plot(res[1][1], res[1][2], size = (600, 300))
for i = 2:length(res)
    plot!(res[i][1], res[i][2], size = (600, 300))
end
δ = 0.15
p = 120
plot!(x -> x / δ, xlims = (0, 400))
plot!(x -> x / δ + p / δ^2, xlims = (0, 400))