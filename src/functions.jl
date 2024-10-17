using BlackBoxOptim

function hjbsolver(NRe, β, FR, ER, EnoR, partition, Δx,
    SearchRanges, η, η1, δ, v0, MaxEvals, p, L)
    End1 = partition[end-1]
    End = partition[end]
    xv = 0:Δx:End
    n = length(xv)
    v = zeros(n + 1)
    vhat = zeros(n)
    v[1] = v0
    Re = zeros(n, NRe)

    pR(r) = (1 + η) * β * EnoR - ((1 + η1) * β * (EnoR - ER(r)))
    V_inf_res = 0
    function H0(r)
        pRR = pR(r)
        if pRR > 0
            s = v[1] * FR(r, Δx)
            return ((δ + β) * v[1] - β * s) / pRR
        else
            return Inf
        end
    end
    MIN = bboptimize(H0; SearchRange=SearchRanges[1], TraceMode=:silent)
    if isnan(best_fitness(MIN))
        MIN = bboptimize(H0; SearchRange=SearchRanges[1], TraceMode=:silent)
    end
    vhat[1] = best_fitness(MIN)
    B_C = best_candidate(MIN)
    v[2] = v[1] + Δx * vhat[1]
    Re[1, :] = B_C
    k = length(SearchRanges)

    TIMER = time()
    l = n
    for i in 2:n
        function H(r)
            pRR = pR(r)
            if pRR > 0
                FRΔ = FR(r, 0.0)
                FRΔ1 = FR(r, Δx)
                s = 0.0
                for j in 1:(i-1)
                    s += (v[i-j+1] + v[i-j]) * (FRΔ1 - FRΔ)
                    FRΔ = FRΔ1
                    FRΔ1 = FR(r, (j + 1) * Δx)
                end
                return ((δ + β) * v[i] - β * 0.5 * s - (i - 1) * Δx) / pRR
            else
                return Inf
            end
        end
        rein_inf = Re[i-1, :]
        for m in 1:k
            if partition[m] <= i * Δx < partition[m+1]
                MIN = bboptimize(H; SearchRange=SearchRanges[m], MaxFuncEvals=MaxEvals, TraceMode=:silent)
                if isnan(best_fitness(MIN))
                    MIN = bboptimize(H; SearchRange=SearchRanges[m], MaxFuncEvals=MaxEvals, TraceMode=:silent)
                end
                vhat[i] = best_fitness(MIN)
                B_C = best_candidate(MIN)

                break
            end

            # Hrein_inf = H(rein_inf)
            # if Hrein_inf < vhat[i]
            #     vhat[i] =  Hrein_inf
            #     B_C = rein_inf
            # end
        end
        if i * Δx >= partition[k+1]
            #rein_inf = Re[i-1,:]
            vhat[i] = H(rein_inf)
            B_C = rein_inf
        end

        Re[i, :] = B_C
        # Hi1 = H(Re[i-1,:])
        # if Hi1 < vhat[i]
        #     Re[i,:] = Re[i-1,:]
        #     vhat[i] = Hi1
        # end
        v[i+1] = v[i] + Δx * vhat[i]

        if (v[i] > i * Δx / δ + p / δ^2 + L)
            l = i
            V_inf_res = 1
            break
        end
        if (v[i] < i * Δx / δ - L)
            l = i
            V_inf_res = -1
            break
        end
        if (i % (div(n, 10)) == 0) | (i == 1)
            # IJulia.clear_output(true)
            print(string(round(Int, 100 * i / n)) * "%, ")
            println("V(" * string(round(i * Δx, sigdigits=4)) * ") = " *
                    string(round(v[i]; sigdigits=5)) * "   Time = " * string(round(time() - TIMER, sigdigits=3)) * " seconds")
        end
    end
    return xv[1:l], v[1:l], vhat[1:l], Re[1:l, :], V_inf_res
end


# function RAND_GEN(rand_dists)

#     R1() = rand(Gamma(3.0, 0.8))
#     mu1 = 2.4

#     R2() = rand(Weibull(2.0, 2.3))
#     mu2 = 2.038321928541343

#     R3() = rand(Exponential(4)) + rand(Exponential(2))
#     mu3 = 6.0
#     RAND_FUNS = [R1, R2, R3]
#     MEANS = [mu1, mu2, mu3]
#     return RAND_FUNS, MEANS
# end

# function simCS(η, η1, δ, α, x, βs, RAND_GEN)
#     RAND_MEAN = RAND_GEN()
#     M = RAND_MEAN[2]
#     RAND_FUNS = RAND_MEAN[1]
#     μsum = sum(M)
#     p = ((1 + η) - (1 - α) * (1 + η1)) * sum(βs .* M)
#     β = sum(βs)
#     k = length(βs)
#     δi2 = 1 / δ^2

#     intX(t1, t2, y) = -(-δ * p * t1 - δ * x + δ * y - p) * exp(-δ * t1) * δi2 + (-δ * p * t2 - δ * x + δ * y - p) * exp(-δ * t2) * δi2
#     S = 0.0
#     t1 = 0.0
#     t2 = rand(Exponential(1 / β))
#     X0 = x
#     y = 0.0
#     INT = 11.0
#     while (X0 >= 0) & (INT > 0.0000001)
#         INT = intX(t1, t2, y)
#         S += INT
#         j = wsample(1:k, βs)
#         y += α * RAND_FUNS[j]()
#         X0 = x + p * t2 - y
#         t1 = t2
#         t2 += rand(Exponential(1 / β))
#     end
#     return S
# end

# function SimulateVconst(η, η1, δ, α, x, n)
#     S = 0.0
#     for i in 1:n
#         S += simCS(η, η1, δ, α, x)
#     end
#     return S / n
# end

function SampleDistFun(F::Function, End::Float64, Δx::Float64)
    VF = F.(0:Δx:End)
    Δxi = 1 / Δx
    function FSample(x)
        if x >= End
            return 1.0
        end
        idx = floor(Int, Δxi * x) + 1
        Fa = VF[idx]
        Fb = VF[idx+1]
        return Fa + Δxi * (Fb - Fa) * (x - Δx * (idx - 1))
    end
    return FSample
end


function bisection_hjbsolver(; a, b, max_it, NRe, β, FR, ER, EnoR, partition, Δx,
    SearchRanges, η, η1, δ, MaxEvals=false, p, L)
    RES_ALL = []
    for i = 1:max_it
        V0 = (a + b) / 2
        println((a, V0, b))
        RES = hjbsolver(NRe, β, FR, ER, EnoR, partition, Δx,
            SearchRanges, η, η1, δ, V0, MaxEvals, p, L)
        push!(RES_ALL, RES)
        if RES[5] == 0
            break
        end
        if RES[5] == -1
            a = V0
        end
        if RES[5] == 1
            b = V0
        end
    end
    return RES_ALL
end