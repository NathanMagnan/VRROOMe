function cost_no_Etot(L::Array{Float64, 1}, N::Array{Float64, 1}, Ltot::Float64, γ::Float64) 
    #### Returns a Float64. It is the cost function evaluated in γ (in the case where there is no constrain on the energy) ####
    # getting the constants we will need
    n_pop = length(L)

    # Let's start the integration
    int = 0.0
    for a = 1:n_pop
        Lₐ = L[a]
        nₐ = N[a]

        x = γ * Lₐ

        int += nₐ * (Lₐ * coth(x) - 1 / γ)
    end

    return(Ltot - int)
end

function Dichotomy_step(L::Array{Float64, 1}, N::Array{Float64, 1}, Ltot::Float64, γ_min::Float64, γ_max::Float64) 
    #### Returns a tuple of 2 floats. They are the new upper and lower bounds at the end of the dichotomy step ####
    # Let's get the current middle of the segment and the value of the cost function there
    γ_mid = (γ_max + γ_min) / 2
    c_mid = cost_no_Etot(L, N, Ltot, γ_mid)

    # Let's get the new bounds
    if (c_mid > 0) # We use the fact that the cost function is monotonous and decreasing
        return (γ_mid, γ_max)
    elseif (c_mid < 0)
        return (γ_min, γ_mid)
    else
        return (γ_mid, γ_mid)
    end
end

function Solution_no_Etot_dichotomy(L::Array{Float64, 1}, N::Array{Float64, 1}, Ltot::Float64) 
    #### Returns a Float64. It is the optimal γ found by the dichotomy method ####
    #### This function performs the whole dichotomy                           ####

    # first let's make define the starting bounds of the dichotomy
    min_L, useless = findmin(L)
    max_L, useless = findmax(L)

    γ_min = 10.0^(-7) / max_L # these bounds come from looking at the shape of coth(x) - 1/x
    γ_max = 10.0^(2) / min_L
    γ_mid = (γ_max + γ_min) / 2

    # then let's define an acceptable relative error and maximum number of steps
    ϵ = 10^(-6)
    n_step = 100

    # then let's do the dichotomy method itself
    i = 0
    cost = Inf
    while ((cost > ϵ) && (i <= n_step))
        γ_min, γ_max = Dichotomy_step(L, N, Ltot, γ_min, γ_max) # update the bounds

        γ_mid = (γ_max + γ_min) / 2
        cost = abs(cost_no_Etot(L, N, Ltot, γ_mid)) # update the cost and test wether it is acceptable

        i += 1
    end

    return(γ_mid)
end

function Magnetisations_no_Etot(L::Array{Float64, 1}, N::Array{Float64, 1}, γ::Float64, 
                                lmax::Int64, 
                                List_x::Array{Float64,1}, List_w::Array{Float64,1}, List_Y::Array{Float64, 2}) 
    #### Returns an array of size n_pop $ lmax_half. Each cell contains one (axisymmetric, even) magnetisation                               ####
    #### This function computes the magnetization from a solution given by the dichotomy method, but doesn't run the dichotomy method itself ####
    # getting the constants we will need
    n_pop = length(L)
    lmax_half = div(lmax, 2)
    n_quad = length(List_x)

    # Preallocating the arrays we will need
    M = zeros(n_pop * lmax_half) # Will hold the magnetisations

    # Let's start to compute the magnetisations
    for star = 1:n_pop
        Lₛ = L[star]
        nₛ = N[star]

        # first let's compute the renrmalizing factor (lower integral)
        int_down = 0.0
        for i = 1:n_quad
            weight = List_w[i]
            x = List_x[i]

            int_down += weight * exp(γ * Lₛ * x)
        end

        # then let's compute the upper integral
        for l_half = 1:lmax_half
            index = (star - 1) * lmax_half + l_half

            for i = 1:n_quad
                spherical_harmonic = List_Y[l_half, i]
                weight = List_w[i]
                x = List_x[i]

                M[index] += weight * nₛ * spherical_harmonic * exp(γ * Lₛ * x)
            end

            M[index] /= int_down
        end
    end

    return(M)
end