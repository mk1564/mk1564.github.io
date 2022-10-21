## HW 1 (ECON605, Prof. Carlos Esquivel)
# Min Kim (min.kim@rutgers.edu)
# January 23, 2021

using Plots

## Model
struct Model
    β::Float64
    σ::Float64
    α::Float64
    δ::Float64
    A::Float64
    N::Int64
    kgrid::StepRangeLen
    vguess::Array{Float64,1}
    u:: Function 
    f:: Function 
end

function construct_model(;σ=1.5, δ=.2)
    β = .6
    σ = σ # RRA
    α = .3 # capital share
    δ = δ # depreciation rate
    A = 1.0 # TFP

    kss = (α*A/(1/β  - (1-δ)))^(1/(1 - α)) # steady state capital level
    kbar = (A/δ)^(1/(1-α)) # maximum sutatiable capital level

    N = 300 # # of grids in capital

    kgrid = range(0.001, 3*kss, length=N) # equal-spaced grid for capital

    
    a0 = 1/(1-β) * (log(1-α*β)+α*β/(1-α*β)*log(α*β))
    a1 = α/(1-α*β)
    
    vguess = a0 .+ a1*log.(kgrid)

    u(c) = (c.^(1.0-σ)-1.0)/(1.0-σ)
    f(k) = A*k.^α
    return Model(β, σ, α, δ, A, N, kgrid, vguess, u, f)
end

function VFI(m::Model)
    v0 = zeros(m.N)
    v1 = zeros(m.N)
    rhs = zeros(m.N)
    kpolicy = collect(similar(m.kgrid))
    diff = 1.0
    maxiter = 100
    tol = 1e-16 # tolerance for convergence

    for iter in 1:maxiter
        for (i, k) in enumerate(m.kgrid)
            for (j, kp) in enumerate(m.kgrid)
                c = m.f(k)+(1-m.δ)*k - kp # consumption
                if c > 0.0 # check negativity
                    rhs[j] = m.u(c) + m.β * v0[j]
                else 
                    rhs[j] = -Inf # penalty if negative
                end
            end
            v1[i], idx = findmax(rhs)
            kpolicy[i] = m.kgrid[idx]
        end
        diff = maximum(abs.(v1-v0))
        println("Iteration $iter: difference $diff")
        if diff < tol
            break
        end
        v0 = copy(v1)
    end
        return v1, kpolicy
end

## Question 3
m = construct_model()
@time v1, kpolicy = VFI(m)

# figures
p = plot(m.kgrid,v1, label = "Value function")
plot!(p, legend =:bottomright)
display(p)
savefig(p,"v_q3.svg")

p1 = plot(m.kgrid,kpolicy, label = "Policy function")
plot!(p1, m.kgrid,m.kgrid, color = :black, linestyle=:dot, label ="45 degree line")
plot!(p1, legend =:bottomright)
display(p1)
savefig(p1,"kp_q3.svg")

## Question 4
m2 = construct_model(;σ = 1.00001, δ = 1)
@time v1_2, kpolicy_2 = VFI(m2)

p_2 = plot(m2.kgrid,v1_2, label = "Value function")
plot!(p_2,m2.kgrid,m2.vguess, label = "guess")
plot!(p_2, legend =:bottomright)
display(p_2)
savefig(p_2,"v_q4.svg")

p_21 = plot(m2.kgrid,kpolicy_2, label = "Policy function")
plot!(p_21, m2.kgrid,m2.kgrid, color = :black, linestyle=:dot, label ="45 degree line")
plot!(p_21, legend =:bottomright)
display(p_21)
savefig(p_21,"kp_q4.svg")