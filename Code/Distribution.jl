include("./Spherical_harmonics.jl")
using FastGaussQuadrature
using Plots

mutable struct Distribution
    #### This structure will hold many arrays, it reduce the memory cost and time cost (through allocations)        ####
    #### But the user doesn't really have to deal with most of them. I think only Etot, Ltot, J, L, N are important ####
    #### V could be usefull to compute an entropy quickly, but one can easily do without                            ####

    # Orbital parameters
    J::Array{Float64, 3} # coupling coefficients
    L::Array{Float64, 1} # angular momentums

    # Quadrature parameters
    List_x::Array{Float64,1} # nodes of the quadratures
    List_w::Array{Float64, 1} # weigth of the quadratures
    List_Y::Array{Float64, 2} # spherical harmonics at the nodes

    # Microcanonical parameters
    N::Array{Float64, 1} # number of particles in each population
    Etot::Float64 # Total energy
    Ltot::Float64 # Total angular momentum

    # Order parameters
    β::Float64 # Lagrange multiplier for the energy
    γ::Float64 # Lagrange multiplier for the momentum
    M::Array{Float64, 1} # magnetisations -- one dimension but of size n_pop * lmax_half
    #### About the magnetization vector M : we store lexicographically first by the number of the population, then by the number of the mode ####
    #### i.e. M_{i * n_pop + l} corresponds to star i and mode l (not exactly, but that's the idea).                                         ####

    # pre-allocations
    U::Array{Float64, 1} # Will hold an array that is usefull to compute the cost function (upper integrals)
    V::Array{Float64, 1} # Will hold an array that is usefull to compute the cost function (lower integrals)
    dU::Array{Float64, 2} # Will hold an array that is usefull to compute the Jacobian
    dV::Array{Float64, 2} # Will hold an array that is usefull to compute the Jacobian

    cost::Array{Float64, 1} # Will hold the cost vector
    #### About the cost vector : first line will hold the consistency cost on energy, ####
    #### second line the cost on momentum, and next lines on magnetisations           ####
    jacobian::Array{Float64, 2} # Will hold the Jacobian
    #### About the Jacobian : Jᵢ,ⱼ = ∂Costᵢ / ∂xⱼ                                                                                     ####
    #### i.e First line will contains the derivative of Etot, second line the derivatives of Ltot and next lines the derivatives of M ####
    #### and first column will contain the derivatives with respect to β, second column w.r.t. γ and next comumns w.r.t. M            ####

    # normalisations
    Norms::Array{Float64, 1} # to normalize the cost function
end

function Distribution(J::Array{Float64, 3}, L::Array{Float64, 1}, N::Array{Float64, 1}, Etot::Float64, Ltot::Float64; 
                                        β = nothing, γ = nothing, M = nothing) 
    #### Helps to create a distribution object when some parameters are unknown ####

    # getting the constants we will need
    n_pop = size(J)[1] # number of populations
    lmax_half = size(J)[3] # number of interesting modes
    lmax = 2 * lmax_half # number of modes
    size_mag = lmax_half * n_pop # size of the magnitudes vector
    size_jac = lmax_half * n_pop + 2 # size of the jacobian matrix
    n_quad = max(100, 2 * lmax) # number of quadratures
    N_stars = sum(N) # total number of stars

    # let's preallocate the arrays we will need
    List_x , List_w = gausslegendre(n_quad) # will hold the position and weight (respectively) of the quadrature nodes
    List_Y = zeros(lmax_half, n_quad) # Will hold the spherical harmonics
    U = zeros(size_jac) # Will hold an array that is usefull to compute the cost function
    V = zeros(n_pop) # Will hold an array that is usefull to compute the cost function
    dU = zeros(size_jac, size_jac) # Will hold an array that is usefull to compute the Jacobian
    dV = zeros(size_jac, size_jac) # Will hold an array that is usefull to compute the Jacobian
    cost = zeros(size_jac) # Will hold the cost vector
    jacobian = zeros(size_jac, size_jac) # Will hold the Jacobian of the cost function
    Norms = zeros(size_mag) # Will hold the vector used to renormalize every line of the cost function to 1


    # let's compute the spherical harmonics
    Pₗ = zeros(lmax_half) # this is temporary (allocates as early as possible an array used by Spherical_harmonics_a!)
    Yₗ = zeros(lmax_half) # this is temporary (allocates as early as possible an array used by Spherical_harmonics_a!)
    for i = 1:n_quad
        Spherical_harmonics_a!(List_x[i], lmax, Yₗ = Yₗ, Pₗ = Pₗ) # Yₗ now holds the spherical harmonics at this point x
        List_Y[:, i] = Yₗ
    end

    # let's compute the norms
    for star = 1:n_pop
        for l_half = 1:lmax_half
            l = 2 * l_half
            index = (star - 1) * lmax_half + l_half

            Norms[index] = (1 / N[star]) * sqrt((4 * π) / (2 * l + 1)) # || Mₗ[K] || <= N[K] * sqrt((2l + 1) / 4 π)
        end
    end

    # Let's get the initial point β, γ, M, by choice at random
    if isnothing(β)
        r = 2 * rand()
        β = r * (10 / (N_stars * J[1, 1, 1])) # Choice around an order of magnitude, from Takacs (fig. 4 left)
    end

    if isnothing(γ)
        r = 2 * rand()
        γ = r * (1.0 / L[1]) # Choice around an order of magnitude, from Roupas (fig. 7 right)
    end

    if isnothing(M)
        M = zeros(size_mag)

        for star = 1:n_pop # I will try to use a and b as indexes for all loops on populations
            for l_half = 1:lmax_half
                l = 2 * l_half
                index = (star - 1) * lmax_half + l_half

                r = 2 * rand() - 1 # random value between -1 and 1
                M[index] = N[star] * sqrt((2 * l + 1) / (4 * π)) * r
            end
        end
    else
        if (length(M) != lmax_half * n_pop)
            println("Error : lenght(M) = $(length(M)) whereas lmax = $(lmax) and n_pop = $(n_pop)")
        end
    end

    # Let's call the standard initializer for Distribution in order to create our distribution
    dist = Distribution(J, L, List_x, List_w, List_Y, N, Etot, Ltot, β, γ, M, U, V, dU, dV, cost, jacobian, Norms)

    return(dist)
end

function reset!(dist::Distribution) 
    #### Resets V, cost, jacobian to the value they have when all integrals are 0.0 ####
    # getting the constants we will need
    n_pop = size(dist.J)[1] # number of populations
    lmax_half = size(dist.J)[3] # number of interesting modes
    lmax = 2 * lmax_half # number of modes
    size_mag = lmax_half * n_pop # size of the magnitudes vector
    size_jac = lmax_half * n_pop + 2 # size of the jacobian matrix

    # Access the dist's array we will need, using passage by reference
    Etot, Ltot, M, V = dist.Etot, dist.Ltot, dist.M, dist.V
    cost, jacobian = dist.cost, dist.jacobian

    # reset V to 0
    fill!(V, 0.0)

    # resest cost to [0, 0, M]
    cost[1] = Etot
    cost[2] = Ltot
    for a = 1:n_pop
        for l_half = 1:lmax_half
            line = (a - 1) * lmax_half + l_half + 2
            cost[line] = M[line - 2]
        end
    end

    # reset jacobian to [0, 0, ....; 0, 0, ....; 0, 0, Identity]
    for column= 1:(size_jac)
        for line = 1:(size_jac)
            jacobian[line, column] = 0.0

            if ((column == line) && (line > 2))
                jacobian[line, column] += 1.0
            end
        end
    end
end

function reset_v2!(dist::Distribution, star_1::Int64) 
    #### Resets the lines of U, dU, dV we will use in this iteration over star_1 to 0.0 ####
    # getting the constants we will need
    n_pop = size(dist.J)[1] # number of stellar populations
    lmax_half = size(dist.J)[3] # number of interesting modes
    size_jac = lmax_half * n_pop + 2 # size of the jacobian matrix

    # Access the dist's array we will need, using passage by reference
    U, dU, dV = dist.U, dist.dU, dist.dV

    # Reset the lines of U
    U[1] = 0.0
    U[2] = 0.0
    for l_half = 1:lmax_half
        line = (star_1 - 1) * lmax_half + l_half + 2
        U[line] = 0.0
    end

    # Reset the lines of dU
    for column = 1:size_jac
        dU[1, column] = 0.0
        dU[2, column] = 0.0
        for l_half = 1:lmax_half
            line = (star_1 - 1) * lmax_half + l_half + 2
            dU[line, column] = 0.0
        end
    end

    # Reset the lines of dV
    for column = 1:size_jac
        dV[column] = 0.0
    end
end

function Cost_and_Jacobian!(dist::Distribution) 
    ####  Update the cost function and its jacobian, which are 2 of the dist's attributes  ####
    #### NOTE : in all our integrals we forget the factor 2π coming from the integral on ϕ ####
    ####        because it is both on the upper and lower integral.                        ####
    # getting the constants we will need
    n_pop = size(dist.J)[1] # number of stellar populations
    lmax_half = size(dist.J)[3] # number of interesting modes
    n_quad = length(dist.List_x) # number of quadrature nodes
    size_mag = lmax_half * n_pop # size of the magnitudes vector
    size_jac = lmax_half * n_pop + 2 # size of the jacobian matrix

    # Access the dist's array we will need, using passage by reference
    J, L, List_x, List_w, List_Y, N, Etot, Ltot, β, γ, M, U, V, dU, dV, Norms = dist.J, dist.L, dist.List_x, dist.List_w, dist.List_Y, dist.N, dist.Etot, dist.Ltot, dist.β, dist.γ, dist.M, dist.U, dist.V, dist.dU, dist.dV, dist.Norms
    cost, jacobian = dist.cost, dist.jacobian

    # Reseting dist.V, dist.cost, and dist.jacobian -- It does modify V, cost and jacobian here because they were passed by reference
    reset!(dist)

    # Now we can start to comput the integrals
    for star_1 = 1:n_pop
        # reset U[K], dU[K], dV[K] to 0
        # V[K] = dist.V[star_1] has aready been reset to 0
        reset_v2!(dist, star_1)

        # First let's compute U, V, dU, dV for this population
        #
        # U is the upper integral of the cost function : ∫dL (N[K]  * [ϵ; L; Y] * exp)
        # V is the lower integral of the cost function : ∫dL (exp)
        # the cost function looks something like C = ∫dK (U / V)
        #
        # dU is the jacobian of U
        # dV is the gradient of V
        # the Jacobian looks something like J = ∫dK (dU / V - U dVᵀ / V²)

        for i = 1:n_quad # We use a Gauss-Legendre quadrature to to integrate over L

            ##### computing the local parameters (wᵢ, ϵ, momentum)
            # local node of the quadrature
            cos_θ = List_x[i] # position
            weight = List_w[i] # weight

            # ϵ the local one-particle energy
            ϵ = 0.0
            for l_half = 1:lmax_half
                spherical_harmonic = List_Y[l_half, i] # local spherical harmonics

                for star_2 = 1:n_pop
                    index = (star_2 - 1) * lmax_half + l_half
                    ϵ -= spherical_harmonic * M[index] * J[star_2, star_1, l_half]
                end
            end

            # local one-particle momentum
            momentum = L[star_1] * cos_θ

            # local unnormalized distribution function
            exponential = exp(- β * ϵ + γ * momentum)


            #### precomputing some prefactors
            Ne_dθ = N[star_1] * weight * exponential # this prefactor is f * dL
            e_dθ = weight * exponential # this prefactor is f * dL / N[k]


            ###### updating all the upper integrals (U), lower integrals (V), and derivatives (dU, dV)

            ### update U[K]
            U[1] += 0.5 * Ne_dθ * ϵ
            U[2] += Ne_dθ * momentum

            for l_half = 1:lmax_half
                # here, notice that Mₗ[K] = ∫dL (N[K] Yₗ exp) does not need an integration over K.
                # that why each iteration on star_1 will fill diferrent cells of U, but there is no need for a second loop on star_2
                line = (star_1 - 1) * lmax_half + l_half + 2
                spherical_harmonic = List_Y[l_half, i]

                U[line] += Ne_dθ * spherical_harmonic
            end

            ### update V[K]
            V[star_1] += e_dθ

            ### update dU[K]
            # some prefactors we will need
            specific_pref_dM_dβ = -1.0 * Ne_dθ * ϵ
            specific_pref_dM_dγ = Ne_dθ * momentum
            specific_pref_dE_dM = -0.5 * Ne_dθ * (1 - β * ϵ)
            specific_pref_dL_dM = Ne_dθ * β * momentum
            specific_pref_dM_dM = Ne_dθ * β

            # upper left square of the matrix
            dU[1, 1] += -0.5 * Ne_dθ * ϵ^(2) # derivative of C_E with respect to β
            dU[1, 2] += 0.5 * Ne_dθ * momentum * ϵ # derivative of C_E with respect to γ
            dU[2, 1] += -1.0 * Ne_dθ * ϵ * momentum # derivative of C_L with respect to β
            dU[2, 2] += Ne_dθ * momentum^(2) # derivative of C_L with respect to γ

            # upper lines of the matrix
            for l_half = 1:lmax_half
                spherical_harmonic = List_Y[l_half, i]

                for star_2 = 1:n_pop
                    # here, E[K] and L[K] can be derived with respect to Mₗ[K'], that's why we need a second loop on star_2
                    column = (star_2 - 1) * lmax_half + l_half + 2
                    
                    coupling_coefficient = J[star_2, star_1, l_half]
                    dU[1, column] += specific_pref_dE_dM * coupling_coefficient * spherical_harmonic # derivative of C_E with respect to Mₗ[K]
                    dU[2, column] += specific_pref_dL_dM * coupling_coefficient * spherical_harmonic # derivative of C_L with respect to Mₗ[K]
                end
            end

            # left column of the matrix
            for l_half = 1:lmax_half
                # here, the derivatives of Mₗ[K] WITH RESPECT TO β AND γ do not need a second integration over star_2
                line = (star_1 - 1) * lmax_half + l_half + 2
                spherical_harmonic = List_Y[l_half, i]

                dU[line, 1] += specific_pref_dM_dβ * spherical_harmonic # derivative of C_Mₗ[K] with respect to β
                dU[line, 2] += specific_pref_dM_dγ * spherical_harmonic # derivative of C_Mₗ[K] with respect to γ
            end

            # lower right square of the matrix
            for l_half = 1:lmax_half
                spherical_harmonic_column = List_Y[l_half, i]

                for star_2 = 1:n_pop
                    # here we have to derivate Mₗ[K] with respect to Mₗ[K'] so we need a second loop on star_2
                    # but the derivative of Mₗ[K] with respect to Mₗ[K'] does not need a third integration over star_3
                    column = (star_2 - 1) * lmax_half + l_half + 2
                    coupling_coefficient = J[star_2, star_1, l_half]

                    @inbounds @simd for ll_half = 1:lmax_half
                        line = (star_1 - 1) * lmax_half + ll_half + 2
                        spherical_harmonic_line = List_Y[ll_half, i]

                        dU[line, column] += specific_pref_dM_dM * coupling_coefficient * spherical_harmonic_column * spherical_harmonic_line
                    end
                end
            end

            ### update dV[K]
            # some prefactors we will need
            specific_pref_dM = e_dθ * β

            # leftmost 2 coefficients
            dV[1] += -1.0 * e_dθ * ϵ
            dV[2] += e_dθ * momentum

            # rightmost lmax_half * n_pop coefficients
            for l_half = 1:lmax_half
                spherical_harmonic = List_Y[l_half, i]

                for star_2 = 1:n_pop
                    # here, V[K] can be derived with respect to any Mₗ[K'] so we need a second loop on star_2
                    column = (star_2 - 1) * lmax_half + l_half +2

                    dV[column] += specific_pref_dM * J[star_2, star_1, l_half] * spherical_harmonic
                end
            end

        end ##### end of the integrals over L

        # Now that U[K], V[K], dU[K] and dV[K] are known, we can integrate them over K to update cost and jacobian
        # of course, we are inside a loop on star_1, which represents the integration over K

        V_of_K = V[star_1]

        ### update cost 
        # Precisely, Cost = [Etot, Ltot, M] - ∫dK (U[K] / V[K]), and remember that cost has been reinitialized to [Etot; Ltot; M]
        cost[1] -= U[1] / V_of_K
        cost[2] -= U[2] / V_of_K
        for l_half = 1:lmax_half
            line = (star_1 - 1) * lmax_half + l_half + 2

            cost[line] -= U[line] / V_of_K
        end

        #### update Jacobian
        # Precisely, Jacobian = [0, 0, ...; 0, 0, ...; 0, 0, I] - ∫dK (dU'[K] / V[K] - U[K] dV[K]ᵀ / V[K]²),
        # and remember that Jacobian has already been reinitialized to [0, 0, ...; 0, 0, ...; 0, 0, I]
        
        ## part U' / V

        # upper left square of the martrix
        jacobian[1, 1] -= (dU[1, 1] / V_of_K)
        jacobian[1, 2] -= (dU[1, 2] / V_of_K)
        jacobian[2, 1] -= (dU[2, 1] / V_of_K)
        jacobian[2, 2] -= (dU[2, 2] / V_of_K)

        # upper lines of the matrix
        for star_2 = 1:n_pop
            for l_half = 1:lmax_half
                column = (star_2 - 1) * lmax_half + l_half + 2

                jacobian[1, column] -= (dU[1, column] / V_of_K)
                jacobian[2, column] -= (dU[2, column] / V_of_K)
            end
        end

        # left columns of the matrix
        for l_half = 1:lmax_half
            line = (star_1 - 1) * lmax_half + l_half + 2

            jacobian[line, 1] -= (dU[line, 1] / V_of_K)
            jacobian[line, 2] -= (dU[line, 2] / V_of_K)
        end

        # lower right square or the matrix
        for star_2 = 1:n_pop
            for l_half = 1:lmax_half
                column = (star_2 - 1) * lmax_half + l_half + 2

                @inbounds @simd for ll_half = 1:lmax_half
                    line = (star_1 - 1) * lmax_half + ll_half + 2

                    jacobian[line, column] -= (dU[line, column] / V_of_K)
                end
            end
        end

        ### part U dVᵀ / V²

        # upper two lines of the matrix
        for column = 1:(size_jac)
            jacobian[1, column] += (U[1] * dV[column] / (V_of_K^(2)))
            jacobian[2, column] += (U[2] * dV[column] / (V_of_K^(2)))
        end

        # lower n_pop * lmax_half lines of the matrix
        for column = 1:(size_jac)
            for l_half = 1:lmax_half
                line = (star_1 - 1) * lmax_half + l_half + 2

                jacobian[line, column] += (U[line] * dV[column] / (V_of_K^(2)))
            end
        end

        ##### End of the treatment of population a
    end

    ##### Normalisations 
    # We want all dimensions of the cost vector to be normalized to 1
    # To do that, we will make the renormalizations {Etot -> 1; Ltot -> 1; ∀ l, ∀ K,  N[K] sqrt((2l + 1) / 4 π) -> 1}

    
    # for the cost
    cost[1] /= Etot
    cost[2] /= Ltot

    for star_1 = 1:n_pop
        for l_half = 1:lmax_half
            line = (star_1 - 1) * lmax_half + l_half + 2

            if (isinf(Norms[line - 2]) && (cost[line] == 0.0))
                cost[line] = 0.0
            else
                cost[line] *= Norms[line - 2]
            end
        end
    end

    # for the jacobian
    for column = 1:(size_jac)
        jacobian[1, column] /= Etot
        jacobian[2, column] /= Ltot

        for star_1 = 1:n_pop
            for l_half = 1:lmax_half
                line = (star_1 - 1) * lmax_half + l_half + 2

                if (isinf(Norms[line - 2]) && (jacobian[line, column] == 0.0))
                    jacobian[line, column] = 0.0
                else
                    jacobian[line, column] *= Norms[line - 2]
                end
            end
        end
    end

    # I don't need to update dist.U, dist.V, dist.dU, dist.dV, dist.cost, dist.jacobian because they were passed by reference
end

function get_f(dist::Distribution, star::Int64, i::Int64)
    #### Returns a Float. It is the distribution function for this distribution, for a given population (star), at a given quadrature (index i) ####
    #### This shouldn't be used on an unoptimized distribution                                                                                  ####
    # Passing by reference the dist's fields we will need
    J, L, List_x, List_Y, N, β, γ, M, V = dist.J, dist.L, dist.List_x, dist.List_Y, dist.N, dist.β, dist.γ, dist.M, dist.V
    
    # getting the constants we will need
    n_pop = size(J)[1]
    lmax_half = size(J)[3]

    # getting ϵ : the local one-particle energy
    ϵ = 0.0
    for l_half = 1:lmax_half
        spherical_harmonic = List_Y[l_half, i] # local spherical harmonics

        for star_2 = 1:n_pop
            index = (star_2 - 1) * lmax_half + l_half
            ϵ -= spherical_harmonic * M[index] * J[star, star_2, l_half]
        end
    end

    # getting the local one-particle momentum
    momentum = L[star] * List_x[i] # local momentum

    # getting the local exponential
    exponential = exp(- β * ϵ + γ * momentum)

    # return
    f = N[star] * exponential / V[star]
    return(f) # f(K, L) = N[K] * exp / ∫dL exp
end

function plot_that_dist(x_axis::String, 
                        Etot_norm::Float64, s::Float64, lmax::Int64, 
                        dist::Distribution, 
                        Bins_m::Array{Tuple{Float64,Float64}, 1}, n_pop_m::Int64,
                        Bins_a::Array{Tuple{Float64,Float64}, 1}, n_pop_a::Int64,
                        Bins_e::Array{Tuple{Float64,Float64}, 1}, n_pop_e::Int64)
    ####  Returns a plot. It is a heatmap of the distribution function with y = cos(θ) and x = m, a or e ####
    #### Options for input x_axis are "m", "a" and "e"
    # Getting the constants we will need
    n_pop = n_pop_m * n_pop_a * n_pop_e

    # Let's get our x axis
    if (x_axis == "m")
        List_m = zeros(n_pop_m)

        for k_m = 1:n_pop_m
            m_down, m_up = Bins_m[k_m]
            List_m[k_m] = 10^((log10(m_up) + log10(m_down)) / 2)
        end
    elseif (x_axis == "a")
        List_a = zeros(n_pop_a)

        for k_a = 1:n_pop_a
            a_down, a_up = Bins_a[k_a]
            List_a[k_a] = 10^((log10(a_up) + log10(a_down)) / 2)
        end
    elseif (x_axis == "e")
        List_e = zeros(n_pop_e)
        
        for k_e = 1:n_pop_e
            e_down, e_up = Bins_e[k_e]
            List_e[k_e] = (e_up + e_down) / 2
        end
    else
        return "Error : x_axis should be m, a or e"
    end

    # Let's get our y axis : the positions of the quadratures
    n_quad = max(100, 2 * lmax)
    List_x, List_w = gausslegendre(n_quad)

    # Create the array that will become a heatmap
    if (x_axis == "m")
        List_f = zeros(n_quad, n_pop_m)

        for k_m = 1:n_pop_m

            #first let's get the raw vector of densities
            i_true = 1
            for i = 1:n_quad
                x = List_x[i]

                # first we find the index of the distribution's quadrature closest to our x
                while (x > dist.List_x[i_true]) 
                    i_true +=1
                end
        
                # then we read the distribution function there
                for k_a = 1:n_pop_a
                    for k_e = 1:n_pop_e
                        index = k_e + n_pop_e * (k_a - 1) + n_pop_e * n_pop_a * (k_m - 1)
                        f = get_f(dist, index, i_true)
        
                        List_f[i, k_m] += f # To get the distribution in mass, we have to integrate over a and e
                    end
                end
            end
        
            # now let's normalize with max
            max_f = -Inf
            for i = 1:n_quad
                f = List_f[i, k_m]
        
                if (f > max_f)
                    max_f = f
                end
            end
        
            for i = 1:n_quad
                List_f[i, k_m] /= max_f
            end
        end

    elseif (x_axis == "a")
        List_f = zeros(n_quad, n_pop_a)
        
        for k_a = 1:n_pop_a

            #first let's get the raw vector of densities
            i_true = 1
            for i = 1:n_quad
                x = List_x[i]

                # first we find the index of the distribution's quadrature closest to our x
                while (x > dist.List_x[i_true])
                    i_true +=1
                end
        
                # then we read the distribution function there
                for k_m = 1:n_pop_m
                    for k_e = 1:n_pop_e
                        index = k_e + n_pop_e * (k_a - 1) + n_pop_e * n_pop_a * (k_m - 1)
                        f = get_f(dist, index, i_true)
        
                        List_f[i, k_a] += f # To get the distribution in sma, we have to integrate over m and e
                    end
                end
            end
        
            # now let's normalize with max
            max_f = -Inf
            for i = 1:n_quad
                f = List_f[i, k_a]
        
                if (f > max_f)
                    max_f = f
                end
            end
        
            for i = 1:n_quad
                List_f[i, k_a] /= max_f
            end
        end

    elseif (x_axis == "e")
        List_f = zeros(n_quad, n_pop_e)
        
        for k_e = 1:n_pop_e

            #first let's get the raw vector of densities
            i_true = 1
            for i = 1:n_quad
                x = List_x[i]

                # first we find the index of the distribution's quadrature closest to our x
                while (x > dist.List_x[i_true])
                    i_true +=1
                end
        
                # then we read the distribution function there
                for k_m = 1:n_pop_m
                    for k_a = 1:n_pop_a
                        index = k_e + n_pop_e * (k_a - 1) + n_pop_e * n_pop_a * (k_m - 1)
                        f = get_f(dist, index, i_true)
        
                        List_f[i, k_e] += f # To get the distribution in eccentricity, we have to integrate over m and a
                    end
                end
            end
        
            # now let's normalize with max
            max_f = -Inf
            for i = 1:n_quad
                f = List_f[i, k_e]
        
                if (f > max_f)
                    max_f = f
                end
            end
        
            for i = 1:n_quad
                List_f[i, k_e] /= max_f
            end
        end

    else
        return "Error : x_axis should be m, a or e"
    end

    # Let's make the plot
    plotlyjs()

    if (x_axis == "m")
        my_plot = heatmap(List_m, List_x, List_f,
                        c = :thermal, colorbar_title = "f_eq",
                        xscale = :log10,
                        xlabel = "m", ylabel = "cos(θ)",
                        title = "E = $(round(Etot_norm, digits = 3)), s = $(round(s, digits = 2))")

    elseif (x_axis == "a")
        my_plot = heatmap(List_a, List_x, List_f,
                        c = :thermal, colorbar_title = "f_eq",
                        xscale = :log10,
                        xlabel = "a", ylabel = "cos(θ)",
                        title = "E = $(round(Etot_norm, digits = 3)), s = $(round(s, digits = 2))")
    
    elseif (x_axis == "e")
        my_plot = heatmap(List_e, List_x, List_f,
                        c = :thermal, colorbar_title = "f_eq",
                        xlabel = "e", ylabel = "cos(θ)",
                        title = "E = $(round(Etot_norm, digits = 3)), s = $(round(s, digits = 2))")
    
    else
        return "Error : x_axis should be m, a or e"
    end

    # Return
    return(my_plot)
end