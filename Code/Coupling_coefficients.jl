include("./Legendre_polynomials.jl")

function Coupling_coefficients(K::Array{Tuple{Float64,Float64,Float64},1}, lmax::Int64; G::Float64 = 1.0) 
    #### Returns an array of size n_pop * n_pop * lmax_half. Each cell holds one (even) coupling coefficient.    ####
    #### The difference with Coupling_coefficients_v2 is that here we take K the list of all orbital populations ####
    #### And return the list of all coupling coefficients. _v2 treats only 2 populations.                        ####

    # getting all the constants we will need
    n_pop = length(K)
    lmax_half = div(lmax, 2)
    n_bin = 10^5 # number of Riemann nodes we will use

    # preallocating all the arrays we will need
    J = zeros(n_pop, n_pop, lmax_half) # will store all the coupling coefficients
    P_l = zeros(lmax_half) # these will be the legendre polynomials evaluted in 0.0
    Coef_l = zeros(lmax_half) # these will be the multiplicative coefficients in front of the integrals (modified at each iteration)
    I_l = zeros(lmax_half) # these will store the integrals
    R₁ = zeros(n_bin) # will store the list of radii at Riemann nodes for the first annuli
    R₂ = zeros(n_bin) # will store the list of radii at Riemann nodes for the second annuli

    # Computing the prefactors
    Legendre_polynomials!(0.0, lmax, Pₗ = P_l) # P_l now holds the list of even Legendre polynomials at point 0.0

    # Computing the Riemann nodes
    cos_E = zeros(n_bin)
    for i = 1:n_bin
        cos_E[i] = cos(π * (2 * i - 1) / (2 * n_bin)) # These are the nodes of the Riemann sums (mid-point rule for faster convergence)
    end

    # Starting the integrals
    for star_1 = 1:n_pop
        m₁, a₁, e₁ = K[star_1]

        # getting the radii
        for i = 1:n_bin
            R₁[i] = a₁ * (1.0 - e₁ * cos_E[i])
        end

        for star_2 = 1:star_1 # We use the symmetry Jₗ[K₁, K₂] = Jₗ[K₂, K₁]
            m₂, a₂, e₂ = K[star_2]

            # getting the radii
            for i = 1:n_bin
                R₂[i] = a₂ * (1.0 - e₂ * cos_E[i])
            end

            # resetting the array of integrals
            for l_half = 1:lmax_half
                I_l[l_half] = 0.0
            end

            # Computing the prefactors
            for l_half = 1:lmax_half
                l = 2 * l_half
                Coef_l[l_half] = G * 4 * π * m₁ * m₂ / (2 * l + 1)
            end

            # starting the integrals -- We use the multipole algorithm
            for l_half = 1:lmax_half

                # let's start with the pⱼ
                wⱼ = 0
                pⱼ = 0.0
                for j = 1:n_bin
                    # first part of the recurrence : multiplication
                    if j > 1 
                        pⱼ = pⱼ * (R₂[j - 1] / R₂[j])^((2 * l_half) + 1)
                    end

                    # second part of the recurrence : addition of δp
                    i = wⱼ
                    while ((i + 1 <= n_bin) && (R₁[i + 1] <= R₂[j]))
                        i += 1
                        μ₁ᵢ = R₁[i] / a₁
                        pⱼ += μ₁ᵢ * (R₁[i] / R₂[j])^(2 * l_half) / R₂[j]
                    end
                    wⱼ = i

                    μ₂ⱼ = R₂[j] / a₂
                    I_l[l_half] += μ₂ⱼ * pⱼ
                end

                # now continue with the qⱼ
                wⱼ = n_bin + 1
                qⱼ = 0.0
                for j = n_bin:-1:1
                    # first part of the recurrence : multiplication
                    if j < n_bin
                        qⱼ = qⱼ * (R₂[j] / R₂[j + 1])^(2 * l_half)
                    end

                    # second part of the recurrence : addition of δq
                    i = wⱼ
                    while ((i - 1 >= 1) && (R₁[i - 1] > R₂[j]))
                        i -= 1
                        μ₁ᵢ = R₁[i] / a₁
                        qⱼ += μ₁ᵢ * (R₂[j] / R₁[i])^(2 * l_half) / R₁[i]
                    end
                    wⱼ = i

                    μ₂ⱼ = R₂[j] / a₂
                    I_l[l_half] += μ₂ⱼ * qⱼ
                end
            end

            # Let's store the results of the integral in J
            for l_half = 1:lmax_half
                J_K₁_K₂_l = Coef_l[l_half] * I_l[l_half] * (P_l[l_half]^(2)) / (n_bin^2) # Jₗ[K₁, K₂] = G m₁ m₂ 4 π Pₗ(0)² / (2 * l + 1) * int et int = double Riemann sum / n_bin^2
                J[star_1, star_2, l_half] = J_K₁_K₂_l # We use the symmetry once again there
                J[star_2, star_1, l_half] = J_K₁_K₂_l
            end
        end
    end

    return(J)
end

function Coupling_coefficients_v2!(K₁::Tuple{Float64,Float64,Float64}, K₂::Tuple{Float64,Float64,Float64}, lmax::Int64; 
                    G::Float64 = 1.0,
                    Pₗ = nothing, Cₗ = nothing, R₁ = nothing, R₂ = nothing, Iₗ = nothing, Jₗ = nothing) 
    #### returns an array of size lmax_half. Each cell contains one (even) coupling coupling coefficient between K₁ and K₂. ####
    #### The difference with Coupling_coefficients is that here we take 2 populations K₁ and K₂ and return the coupling     ####
    #### coefficients between them, whereas Coupling_coefficients works with any number n_pop of populations.               ####

    # getting the constants we will need
    lmax_half = div(lmax, 2)
    n_bin = 10^3 # number of Riemann nodes we will use

    # let's preallocate the arrays that have not been given as argument
    if isnothing(Pₗ)
        Pₗ = zeros(lmax_half) # Will store the (even) Legendre polynomials evaluated in 0.0
    end
    if isnothing(Cₗ)
        Cₗ = zeros(lmax_half) # Will store the prefactors G m₁ m₂ 4 π / (2l + 1)
    end
    if isnothing(R₁)
        R₁ = zeros(n_bin) # Will store the radii at Riemann nodes for the first annuli
    end
    if isnothing(R₂)
        R₂ = zeros(n_bin) # Will store the radii at Riemann nodes for the second annuli
    end
    if isnothing(Iₗ)
        Iₗ = zeros(lmax_half) # Will store the Riemann integrals
    end
    if isnothing(Jₗ)
        Jₗ = zeros(lmax_half) # Will store the coupling coefficients
    end

    # let's compute the prefactors
    m₁, a₁, e₁ = K₁
    m₂, a₂, e₂ = K₂
    for l_half = 1:lmax_half
        Cₗ[l_half] = G * 4 * π * m₁ * m₂ / (2 * (2 * l_half) + 1)
    end

    Legendre_polynomials!(0.0, lmax, Pₗ = Pₗ) # Pₗ now holds the even Legendre polynomials at point 0.0

    # let's compute the radii at the nodes of the Riemann sum
    for j = 1:n_bin
        Ecc = π * (2 * j - 1) / (2 * n_bin)
        R₁[j] = a₁ * (1.0 - e₁ * cos(Ecc))
        R₂[j] = a₂ * (1.0 - e₂ * cos(Ecc))
    end

    # let's reset the array containing the integrals
    for l_half = 1:lmax_half
        Iₗ[l_half] = 0.0
    end

    # Now let's compute the integrals -- We use the multipole algorithm
    for l_half = 1:lmax_half
        # let's start with the pⱼ
        wⱼ = 0
        pⱼ = 0.0
        for j = 1:n_bin
            # first part of the recurrence : multiplication
            if j > 1 
                pⱼ = pⱼ * (R₂[j - 1] / R₂[j])^((2 * l_half) + 1)
            end

            # second part of the recurrence : addition of δp
            i = wⱼ
            while ((i + 1 <= n_bin) && (R₁[i + 1] <= R₂[j]))
                i += 1
                μ₁ᵢ = R₁[i] / a₁
                pⱼ += μ₁ᵢ * (R₁[i] / R₂[j])^(2 * l_half) / R₂[j]
            end
            wⱼ = i

            μ₂ⱼ = R₂[j] / a₂
            Iₗ[l_half] += μ₂ⱼ * pⱼ / (n_bin^2)
        end

        # now continue with the qⱼ
        wⱼ = n_bin + 1
        qⱼ = 0.0
        for j = n_bin:-1:1
            # first part of the recurrence : multiplication
            if j < n_bin
                qⱼ = qⱼ * (R₂[j] / R₂[j + 1])^(2 * l_half)
            end

            # second part of the recurrence : addition of δq
            i = wⱼ
            while ((i - 1 >= 1) && (R₁[i - 1] > R₂[j]))
                i -= 1
                μ₁ᵢ = R₁[i] / a₁
                qⱼ += μ₁ᵢ * (R₂[j] / R₁[i])^(2 * l_half) / R₁[i]
            end
            wⱼ = i

            μ₂ⱼ = R₂[j] / a₂
            Iₗ[l_half] += μ₂ⱼ * qⱼ / (n_bin^2)
        end
    end

    # Let's save the integrals in J
    for l_half = 1:lmax_half
        Jₗ[l_half] = Cₗ[l_half] * Iₗ[l_half] * (Pₗ[l_half]^(2)) # Jₗ[K₁, K₂] = G m₁ m₂ 4 π Pₗ(0)² / (2 l + 1) * int
    end
end