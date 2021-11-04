include("./Distribution.jl")
using LinearAlgebra

mutable struct Optimisation_problem
    #### This small structure is only to make Newton's method mor readable ####
    δ::Array{Float64, 1} # describes the step
    cost::Float64 # value of the cost function
end

function make_a_step!(dist::Distribution, prb::Optimisation_problem) 
    #### modifies β, γ, M -- assumes the step δ is correct ####
    # getting the constants we will need
    n_pop = size(dist.J)[1]
    lmax_half = size(dist.J)[3]

    # let's perform the step
    dist.β -= prb.δ[1]
    dist.γ -= prb.δ[2]
    for l_half = 1:lmax_half
        for a = 1:n_pop
            line = (l_half - 1) * n_pop + a + 2
            dist.M[line - 2] -= prb.δ[line]
        end
    end
end

function update_Newton!(dist::Distribution, prb::Optimisation_problem) 
    #### Compute the step of Newton's method, then makes it                                    ####
    #### Assumes that cost and jacobian have already been computed when we enter this function ####

    ldiv!(prb.δ, lu(dist.jacobian), dist.cost) # this computes the step to performs

    make_a_step!(dist, prb) # this makes the step and modifies β, γ, M

    Cost_and_Jacobian!(dist) # this is makes sure the cost and jacobian have been computed when we next call update_Newton!
    prb.cost = norm(dist.cost) # computes how close to the root we are AFTER the step
end

function Root_finder_Newton!(dist::Distribution) 
    #### Performs the entire Newton method                                                                    ####
    #### However, it doesn't choose the initial point itself, this should be done by a wrapper to this script ####
    # getting the constants we will need
    n_pop = size(dist.J)[1]
    lmax_half = size(dist.J)[3]

    # preallocating an array we will need
    δ = zeros(lmax_half * n_pop + 2)

    # Let's choose a precision goal and an maximal number of iterations
    ϵ = 10^(-6) # maximum acceptable cost
    n_max = 100 # maximumm acceptable number of iterations

    # Computing the initial cost because the Newton step assumes that when the step begins the cost has already been computed
    Cost_and_Jacobian!(dist)
    cost = norm(dist.cost)

    # Initializing our small Optimization problem with this cost
    prb = Optimisation_problem(δ, cost)

    # Let's perform the Newton method now
    n = 0
    while ((n < n_max) && (prb.cost > ϵ))
        n += 1
        update_Newton!(dist, prb)
    end
end
