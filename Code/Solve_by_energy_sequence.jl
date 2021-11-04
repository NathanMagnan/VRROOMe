#### This script is the highest wrapper, that's why I overcommented it ####
include("./Distribution.jl")
include("./Newton.jl")
include("./Spherical_harmonics.jl")
include("./No_energy_constraint_case.jl")
using LinearAlgebra
using FastGaussQuadrature

function solve_by_energy_sequence(J::Array{Float64, 3}, L::Array{Float64, 1}, N::Array{Float64, 1}, Etot::Float64, Ltot::Float64) 
    #### Returns a Distribution. It is the solution to the microcanonical problem Etot, Ltot IF it exists. ####
    #### One should always check that the cost of that dist is less than 10⁻⁶ !                                           ####
    #### This function solves the problem by starting at the β = 0 point which is easy to find,                           ####
    #### then moving through the sequence of equilibrium to the point of interest .                                       ####

    # let's get all the constants we need
    n_pop = size(J)[1] # number of populations
    lmax_half = size(J)[3] # number of interesting modes
    lmax = 2 * lmax_half # number of modes
    n_quad = max(100, 2 * lmax) # number of quadratures
    size_mag = lmax_half * n_pop # size of the magnitudes vector
    size_jac = lmax_half * n_pop + 2 # size of the jacobian matrix

    # let's preallocate the arrays we will need
    List_x , List_w = gausslegendre(n_quad) # will hold the position and weight (respectively) of the quadrature nodes
    List_Y = zeros(lmax_half, n_quad) # Will hold the spherical harmonics

    # let's compute the spherical harmonics
    Pₗ = zeros(lmax_half) # This is temporary (allocates as early as possible an array used by Spherical_harmonics_a!)
    Yₗ = zeros(lmax_half) # this is temporary (allocates as early as possible an array used by Spherical_harmonics_a!)
    for i = 1:n_quad
        Spherical_harmonics_a!(List_x[i], lmax, Yₗ = Yₗ, Pₗ = Pₗ) # Yₗ now holds the spherical harmonics at this point x
        List_Y[:, i] = Yₗ
    end

    # let's get the initial point (β, γ, M)
    β = 0.0
    γ = Solution_no_Etot_dichotomy(L, N, Ltot)
    M = Magnetisations_no_Etot(L, N, γ, lmax, List_x, List_w, List_Y)

    # let's create the initial distribution
    Etot_first_guess = 1.0 # at first, we create a distribution with a mock energy
    dist = Distribution(J, L, N, Etot_first_guess, Ltot, β = β, γ = γ, M = M)

    # let's compute the initial energy and use it to correct the energy of the initial distribution
    Cost_and_Jacobian!(dist)
    Etot_start = Etot_first_guess * (1.0 - dist.cost[1])
    dist.Etot = Etot_start

    # let's choose our steps
    n_steps = 50 # this number is quite arbitrary
    List_Etot = - 10 .^ (range(log10(- Etot_start), stop = log10(- Etot), length = n_steps))

    # let's perform the steps
    for n = 2:n_steps
        dist.Etot = List_Etot[n] # We do not create a new distribution at each step, it would cost too much memory
        Root_finder_Newton!(dist)
    end

    # Let's return the final distribution
    return(dist)
end
