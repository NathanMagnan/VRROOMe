include("./Coupling_coefficients.jl")
include("./Coupling_coefficients_interpolated.jl")
include("./Legendre_polynomials.jl")
include("./Spherical_harmonics.jl")
using Distributions
using Plots

mutable struct Cluster_initial
    #### This structure will hold all the data necessary to represent an astrophysical cluster ####
    #### But I will only use it to represent initial clusters i.e. clusters of 16 disks        ####
    # physical parameters
    G::Float64 # constant of graviation, typically 1.0
    M::Float64 # mass of the central black hole, typically 1.0

    # size parameters of the problem
    n_disks::Int64 # number of disks,typically 16
    N::Int64 # number of stars per disk, typically 512
    lmax::Int64 # number of modes considered in the multipole expansion

    # parameters of the distribution
    m_min::Float64 # minimal mass, typically 1.0
    m_max::Float64 # maximal mass, typically 100.0
    γₘ::Float64 # index of the power law for the distribution in mass, typically -2.0

    a_min::Float64 # minimal semi-major axis, typically 1.0
    a_max::Float64 # maximal semi-major axis, typically 100.0
    γₐ::Float64 # index of the power law for the distribution in semi-major axis, typically 0.0

    e_min::Float64 # minimal eccentricity, typically 0.0
    e_max::Float64 # maximal eccentricity, typically 0.3
    γₑ::Float64 # index of the power law for the distribution in eccentricity, typically 1.0

    σ::Float64 # thickness of the disks i.e. for any star in any disk, we have L_star . L_disk >= (1 - σ). σ is typically worth 0.006
    κ::Float64 # von Mises - Fisher parameter for the distribution of the disks, typically worth 0.0

    # description of the cluster
    K::Array{Tuple{Float64,Float64,Float64}, 1} # List of the orbital parameters of all stars in the cluster -- to be chosen to match the distribution
    L::Array{Tuple{Float64,Float64}, 1} # list of the initial orientation of the stellar annulii -- to be chosen to match the distribution

    # invariants of motion
    Etot::Float64 # total energy of the cluster -- invariant of motion -- to be computed once the orientations are chosen
    Ltot::Float64 # total energy of the cluster -- invariant of motion -- to be computed once the orientations are chosen
end

function Cluster_initial(n_disks::Int64, N::Int64, lmax::Int64; 
                        G::Float64 = 1.0, M::Float64 = 1.0, 
                        m_min::Float64 = 1.0, m_max::Float64 = 100.0, γₘ::Float64 = -2.0,
                        a_min::Float64 = 1.0, a_max::Float64 = 100.0, γₐ::Float64 = 0.0,
                        e_min::Float64 = 0.0, e_max::Float64 = 0.3, γₑ::Float64 = 1.0,
                        σ::Float64 = 0.006, κ::Float64 = 0.0) 
    ############ 
    # This function creates a cluster structure demanding only the size parameters n_disks, and lmax. 
    # To do that, it provides typical values for the the parameters of the distribution, from Szölgyen 2018 page 2 right column
    # Then it samples at random n_disks of N stars from the distribution
    # And finally it calls 2 functions that compute the total energy and total angular momentum
    ############

    # chose the orbital parameters
    K = Array{Tuple{Float64,Float64,Float64}, 1}(undef, n_disks * N)

    for i = 1:n_disks # i will always denote the number of the disk
        for j = 1:N # j will always denote the number of the star, within the disk i
            # draw m, a, e
            m = draw_from_power_law(m_min, m_max, γₘ)
            a = draw_from_power_law(a_min, a_max, γₐ)
            e = draw_from_power_law(e_min, e_max, γₑ)

            index = (i - 1) * N + j # index will always denote the absolute number of the star
            K[index] = (m, a, e) # I hope this is a tuple
        end
    end

    # chose the orbital orientations
    L = Array{Tuple{Float64,Float64}, 1}(undef, n_disks * N)

    # first let's choose the position of the disks' centers
    List_centers = Array{Tuple{Float64,Float64}, 1}(undef, n_disks)

    if (κ == 0)
        for i = 1:n_disks
            θ_disk = acos(2 * (rand() - 0.5)) # I want a uniform distribution on the sphere, not on θ
            ϕ_disk = 2 * π * rand()

            List_centers[i] = (θ_disk, ϕ_disk)
        end

    elseif (κ > 0)
        L_z = [0.0, 0.0, 1.0]
        vMF = VonMisesFisher(L_z, κ) # Defines the distribution we will use to draw the disks' centers. 
                                     # With κ = 0, it is the uniform distrbution on the sphere.
        for i = 1:n_disks
            l_drawn = rand(vMF, 1) # this gives the cartesian coordinates

            θ_disk = acos(l_drawn[3]) # WHich I need to translate back to angular coordinates
            if (l_drawn[2] > 0)
                ϕ_disk = acos(l_drawn[1] / sin(θ_disk))
            elseif (l_drawn[2] < 0)
                ϕ_disk = - acos(l_drawn[1] / sin(θ_disk))
            elseif ((l_drawn[2] == 0) && (l_drawn[1] > 0))
                ϕ_disk = 0
            elseif ((l_drawn[2] == 0) && (l_drawn[1] < 0))
                ϕ_disk = π
            elseif ((l_drawn[2] == 0) && (l_drawn[2] == 0))
                ϕ_disk = 0 # I prefer to explicitely break the pole issue, but it should never happen anyways
            else
                println("There is something wrong : $(l_drawn[1]), $(l_drawn[2]), $(l_drawn[3])")
            end

            List_centers[i] = (θ_disk, ϕ_disk) # now I can save it
        end

    else
        println("κ should be positive : κ = $(κ)")
    end

    for i = 1:n_disks
        θ_disk, ϕ_disk = List_centers[i]

        for j = 1:N
            # first I choose point uniformly on the northern polar cap of the standard sphere
            θ_1 = acos(1 - σ * rand())
            ϕ_1 = 2 * π * rand()

            # then I describe this point in standard cartesian coordinates
            x_1 = sin(θ_1) * cos(ϕ_1)
            y_1 = sin(θ_1) * sin(ϕ_1)
            z_1 = cos(θ_1)

            # then I perfom the first rotation θ = 0 -> θ = θ_disk i.e. the rotation around axis y
            x_2 = x_1 * cos(θ_disk) + z_1 * sin(θ_disk)
            y_2 = y_1
            z_2 = - x_1 * sin(θ_disk) + z_1 * cos(θ_disk)

            # then I perform the second rotation ϕ = 0 -> ϕ = ϕ_disk i.e. the rotation around the new axis z
            x_3 = x_2 * cos(ϕ_disk) - y_2 * sin(ϕ_disk)
            y_3 = x_2 * sin(ϕ_disk) + y_2 * cos(ϕ_disk)
            z_3 = z_2

            # then I translate back to angular coordinates
            θ = acos(z_3)

            if (y_3 > 0)
                ϕ = acos(x_3 / sin(θ))
            elseif (y_3 < 0)
                ϕ = 2 * π - acos(x_3 / sin(θ))
            elseif ((x_3 > 0) && (y_3 == 0))
                ϕ = 0.0
            elseif ((x_3 < 0) && (y_3 == 0))
                ϕ = π
            elseif ((x_3 == 0) && (y_3 == 0))
                ϕ = 0.0 # I prefer to explicitely break the pole issue, but it should never happen anyways
            else
                println("There is something wrong : $(x_3), $(y_3), $(z_3)")
            end

            index = (i - 1) * N + j
            L[index] = (θ, ϕ)
        end
    end

    # create a Cluster_initial object with all these infos, by calling the standard initializer of the structure
    Etot_t = Inf # t is for temporary -- This is a mock energy, the correct one will be computed just afterwards
    Ltot_t = Inf # t is for temporary -- This is a mock momentum, the correct one will be computed just afterwards
    
    clust = Cluster_initial(G, M, 
                        n_disks, N, lmax, 
                        m_min, m_max, γₘ,
                        a_min, a_max, γₐ,
                        e_min, e_max, γₑ,
                        σ, κ,
                        K, L,
                        Etot_t, Ltot_t)
    
    # compute the total energy and total momentum of the cluster
    update_Ltot!(clust)
    update_Etot_v3!(clust) # v3 is the fastest of our three algorithms for computing the energy of a cluster

    return(clust)
end

function draw_from_power_law(y_min::Float64, y_max::Float64, γ::Float64)
    #### Returns a Float64. It is a number drawn at random from a power law of index γ between y_min and y_max ####
    #### It assumes that y_min is greater than y_max                                                           ####
    # first draw at random uniformly between 0 and 1
    r = rand()

    # then transform that random variable R into random variable Y that follows the right power law
    y = 0.0 # This is just to help the compiler with scopes
    if (γ > - 1.0)
        x = y_min^(γ + 1) + r * (y_max^(γ + 1) - y_min^(γ + 1))
        y = x^(1 / (γ + 1))
    elseif (γ == -1.0)
        x = log(y_min) + r * log(y_max / y_min)
        y = exp(x)
    else
        x = y_max^(γ + 1) + r * (y_min^(γ + 1) - y_max^(γ + 1))
        y = x^(1 / (γ + 1))
    end

    return(y)
end

function update_Ltot!(clust::Cluster_initial) 
    #### Updates the total angular momentum (in norm) of a cluster ####
    # First let's compute the vector total angular momentum
    Ltot_vect = [0.0 0.0 0.0] # Will store the vector total angular momentum

    for i = 1:clust.n_disks
        for j = 1:clust.N
            index = (i - 1) * clust.N + j

            m, a, e = clust.K[index]
            θ, ϕ = clust.L[index]

            x = sin(θ) * cos(ϕ)
            y = sin(θ) * sin(ϕ)
            z = cos(θ)
            l = m * sqrt(clust.G * clust.M * a * (1 - e^(2)))

            Ltot_vect[1] += l * x
            Ltot_vect[2] += l * y
            Ltot_vect[3] += l * z
        end
    end

    # Then let's imagine that we rotate the sphere so that the vector total angular momentum is positively aligned with z
    # in these conditions, we only need to remember the norm of the total angular momentum
    Ltot = sqrt(Ltot_vect[1]^(2) + Ltot_vect[2]^(2) + Ltot_vect[3]^(2))
    clust.Ltot = Ltot
end

function update_Etot!(clust::Cluster_initial) 
    #### Updates the total energy of a cluster using the formula Etot = - 0.5 ∑ᵢ ∑ⱼ ∑ₗ (2l + 1) / (4 π) × Jₗ[Kᵢ, Kⱼ] Pₗ(Lᵢ . Lⱼ) ####
    #### and complete coupling coefficients.                                                                                   ####
    # Passing by reference the clsuter's fields we will need
    n_disks, N, lmax = clust.n_disks, clust.N, clust.lmax
    L, K = clust.L, clust.K

    # Getting the constants we will need
    lmax_half = div(lmax, 2)
    n_bins = 10^3

    # Preallocating some arrays we will need
    Pₗ = zeros(lmax_half) # Will hold the Legendre polynomials
    Pl = zeros(lmax_half) # Read the comments in coupling_coefficients_v2 for this one
    Cₗ = zeros(lmax_half) # Read the comments in coupling_coefficients_v2 for this one
    R₁ = zeros(n_bins) # Read the comments in coupling_coefficients_v2 for this one
    R₂ = zeros(n_bins) # Read the comments in coupling_coefficients_v2 for this one
    Iₗ = zeros(lmax_half) # Read the comments in coupling_coefficients_v2 for this one
    Jₗ = zeros(lmax_half) # Will hold the coupling coefficients

    # Let's do the triple sum
    Etot = 0.0

    for star_1 = 1:(n_disks * N)
        θ₁, ϕ₁ = L[star_1]
        x₁ = sin(θ₁) * cos(ϕ₁)
        y₁ = sin(θ₁) * sin(ϕ₁)
        z₁ = cos(θ₁)  
        
        K₁ = K[star_1]

        for star_2 = 1:star_1 # we use the symetry of the dot product and of the coupling coefficients
            θ₂, ϕ₂ = L[star_2]
            x₂ = sin(θ₂) * cos(ϕ₂)
            y₂ = sin(θ₂) * sin(ϕ₂)
            z₂ = cos(θ₂)

            K₂ = K[star_2]

            x = x₁ * x₂ + y₁ * y₂ + z₁ * z₂ # This is L₁ . L₂

            Legendre_polynomials!(x, lmax, Pₗ = Pₗ) # Pₗ now holds the (even) Legendre polynomials in L₁ . L₂
            Coupling_coefficients_v2!(K₁, K₂, lmax, Pₗ = Pl, Cₗ = Cₗ, R₁ = R₁, R₂ = R₂, Iₗ = Iₗ, Jₗ = Jₗ) # Jₗ now holds the coupling coefficents between K₁ and K₂

            for l_half = 1:lmax_half
                l = 2 * l_half
                pref = -0.5 *  (2 * l + 1) / (4 * π)

                if (star_1 == star_2) # because we use the symmetry, we have to account twice for each non-diagonal term and once for the diagonal terms
                    Etot += pref * Jₗ[l_half] * Pₗ[l_half]
                else
                    Etot += 2 * pref * Jₗ[l_half] * Pₗ[l_half]
                end
            end
        end
    end

    clust.Etot = Etot
end

function update_Etot_v2!(clust::Cluster_initial) 
    #### Updates the total energy of a cluster using the formula                                                                               #### 
    #### Etot = - 0.5 ∑ᵢ,ⱼ ∑ₗ,ₘ Yₗₘ(Lᵢ) Yₗₘ(Lⱼ) × (G mᵢ mⱼ 4 π Pₗ(0)² / (2 l + 1)) × ∫dEᵢ / π ∫dEⱼ / π μ(rᵢ) μ(rⱼ) Min(rᵢ, rⱼ)ˡ / Max(rᵢ, rⱼ)ˡ⁺¹  ####
    #### and the multipole algorithm. Should theoretically be the fastest, but with this implementation it isn't                               ####
    #### Please feel free to improve this implementation and use it                                                                            ####

    # Passing by reference the cluster's fields we will need
    n_disks, N, lmax, G = clust.n_disks, clust.N, clust.lmax, clust.G
    K, L = clust.K, clust.L

    # getting the constants we will need
    lmax_half = div(lmax, 2) # number of interesting modes
    n_stars = n_disks * N # number of stars
    n_bins = 10^3 # Number of Riemann nodes in the integrals
    n_nodes = n_stars * n_bins # Number of nodes in the multipole algorithm

    # preallocating the arrays we will need
    Pₗ = zeros(lmax_half) # Will store the Legendre polynomials
    γₗ = zeros(lmax_half) # Will store the prefactors γₗ = 4 π Pₗ(0)² / (2 l + 1)
    ordered_n = Array{Tuple{Float64,Int64}, 1}(undef, n_nodes) # Will store the nodes (r, star) of the multipole algorithm
    Cos_E = zeros(n_bins) # Will store the nodes of the Riemann sums
    P = zeros(2 * lmax + 1, lmax_half) # Will hold the partial sums P of the multipole algorithm
    Q = zeros(2 * lmax + 1, lmax_half) # Will hold the partial sums Q of the multipole algorithm
    δP = zeros(2 * lmax + 1, lmax_half) # Will hold the steps δP of the multipole algorithm
    δQ = zeros(2 * lmax + 1, lmax_half) # Will hold the steps δQ of the multipole algorithm
    S = zeros(2 * lmax + 1, lmax_half) # Will hold the final sums of the multipole algorithm
    Yₗₘ = zeros(2 * lmax + 1, lmax_half) # Will hold the spherical harmonics
    Pₗₘ = zeros(lmax + 1, lmax_half) # Read the comments in spherical_harmonics! for this one



    # Let's compute the prefactors
    prefactor = G / (n_bins^2)

    Legendre_polynomials!(0.0, lmax, Pₗ = Pₗ) # Pₗ now hold the (even) Legendre polynomials in 0
    for l_half = 1:lmax_half
        l = 2 * l_half
        γₗ[l_half] = -0.5 * 4 * π * Pₗ[l_half]^(2) / (2 * l + 1)
    end

    # Let's compute the nodes of the Riemann sums
    for k = 1:n_bins
        Cos_E[k] = cos(π * (2 * k - 1) / (2 * n_bins))  # we use the mid-point rule for faster convergence
    end

    # Let's create the unordered list of nodes (r, star)
    for star = 1:n_stars
        m, a, e = K[star]

        for k = 1:n_bins
            index = (star - 1) * n_bins + k
            cos_E = Cos_E[k]

            r = a * (1 - e * cos_E)

            ordered_n[index] = (r, star)
        end
    end

    # finally let's order the list of nodes -- I think it uses the quicksort or Pivot algorithm
    sort!(ordered_n, by = x -> x[1])



    ### Let's compute the energy

    # let's start with P
    r_min, mass_min, μ_min = 0.0, 0.0, 0.0 # this is just to help with scope
    Yₗₘ_min = zeros(2 * lmax + 1, lmax_half) # this is just to help with scope

    for n = 1:n_nodes
        # acess the parameters of this node
        r, star = ordered_n[n]
        mass, a, e = K[star]
        θ, ϕ = L[star]

        # get μ for this node
        μ = r / a

        # get the spherical harmonics for this node
        Spherical_harmonics!(θ, ϕ, lmax, Yₗₘ = Yₗₘ, Pₗₘ = Pₗₘ)

        for l_half = 1:lmax_half
            l = 2 * l_half

            for m = (-l):l
                Y = Yₗₘ[lmax + 1 + m, l_half]
                Y_min = Yₗₘ_min[lmax + 1 + m, l_half]

                if (n > 1)
                    P[lmax + 1 + m, l_half] *= Y * μ * mass / (Y_min * μ_min * mass_min)
                    P[lmax + 1 + m, l_half] *= (r_min / r)^(l + 1)
                end

                δP[lmax + 1 + m, l_half] = (Y * μ * mass)^(2) / r
                P[lmax + 1 + m, l_half] += δP[lmax + 1 + m, l_half]
                S[lmax + 1 + m, l_half] += P[lmax + 1 + m, l_half]
            end
        end

        r_min, mass_min, μ_min = r, mass, μ  
        Yₗₘ_min[:, :] = Yₗₘ
    end

    # let's do the same for Q
    r_sup, mass_sup, μ_sup = 0.0, 0.0, 0.0 # this is just to help with scope
    Yₗₘ_sup = zeros(2 * lmax + 1, lmax_half) # this is just to help with scope

    for n = n_nodes:(-1):1
        # acess the parameters of this node
        r, star = ordered_n[n]
        mass, a, e = K[star]
        θ, ϕ = L[star]

        # get μ for this node
        μ = r / a

        # get the spherical harmonics for this node -- 5 ms
        Spherical_harmonics!(θ, ϕ, lmax, Yₗₘ = Yₗₘ, Pₗₘ = Pₗₘ)

        for l_half = 1:lmax_half
            l = 2 * l_half

            for m = (-l):l
                Y = Yₗₘ[lmax + 1 + m, l_half]
                Y_sup = Yₗₘ_sup[lmax + 1 + m, l_half]

                if (n < n_nodes)
                    Q[lmax + 1 + m, l_half] *= Y * μ * mass / (Y_sup * μ_sup * mass_sup)
                    Q[lmax + 1 + m, l_half] *= (r / r_sup)^(l)

                    δQ[lmax + 1 + m, l_half] = Y * μ * mass * Y_sup * μ_sup * mass_sup * (r / r_sup)^(l) * (1 / r_sup)
                end

                Q[lmax + 1 + m, l_half] += δQ[lmax + 1 + m, l_half]
                S[lmax + 1 + m, l_half] += Q[lmax + 1 + m, l_half]
            end
        end

        r_sup, mass_sup, μ_sup = r, mass, μ  
        Yₗₘ_sup[:, :] .= Yₗₘ
    end


    # now we have all we need to compute the energy
    Etot = 0.0

    for l_half = 1:lmax_half
        l = 2 * l_half
        γ_l = γₗ[l_half]

        for m = (-l):l
            Etot += prefactor * γ_l * S[lmax + 1 + m, l_half]
        end
    end


    # Let's update the energy attribute of the cluster (a Float cannot be passed by reference)
    clust.Etot = Etot
end

function update_Etot_v3!(clust::Cluster_initial) 
    #### Updates the total energy of a cluster using the formula Etot = - 0.5 ∑ᵢ ∑ⱼ ∑ₗ (2l + 1) / (4 π) × Jₗ[Kᵢ, Kⱼ] Pₗ(Lᵢ . Lⱼ) ####
    #### and interpolated coupling coefficients. This is the fastest of my three methods (for now).                            ####
    # Passing by reference the cluster's fields that we will need
    n_disks, N, lmax = clust.n_disks, clust.N, clust.lmax
    L, K = clust.L, clust.K

    # getting the constants we will need
    lmax_half = div(lmax, 2)
    n_stars = n_disks * N

    # Preallocating an array we will need
    Pₗ = zeros(lmax_half) # will hold the Legendre polynomials

    # Let's do the triple sum
    Etot = 0.0

    for star_1 = 1:n_stars
        K₁ = K[star_1]

        θ₁, ϕ₁ = L[star_1]
        x₁ = sin(θ₁) * cos(ϕ₁)
        y₁ = sin(θ₁) * sin(ϕ₁)
        z₁ = cos(θ₁)

        for star_2 = 1:star_1 # we use the symetry of the dot product and of the coupling coefficients
            K₂ = K[star_2]

            θ₂, ϕ₂ = L[star_2]
            x₂ = sin(θ₂) * cos(ϕ₂)
            y₂ = sin(θ₂) * sin(ϕ₂)
            z₂ = cos(θ₂)

            x = x₁ * x₂ + y₁ * y₂ + z₁ * z₂ # This is L₁ . L₂

            Legendre_polynomials!(x, lmax, Pₗ = Pₗ) # Pₗ now holds the (even) Legendre polynomials at point L₁ . L₂

            for l_half = 1:lmax_half
                l = 2 * l_half
        
                pref = -0.5 * (2 * l + 1) / (4 * π)
                Jl = Coupling_coefficient_interpolated_v2(l, K₁, K₂) # this is 3 times slower than theory, I don't know why
                Pl = Pₗ[l_half]

                if (star_1 == star_2) # because we use the symmetry, we have to account twice for each non-diagonal term and once for the diagonal terms
                    Etot += pref * Jl * Pl
                else
                    Etot += 2 * pref * Jl * Pl
                end
            end
        end
    end

    # let's update Etot
    clust.Etot = Etot
end

function get_ES_from_clust(clust::Cluster_initial)
    #### Returns a tuple of 2 Float64. They are the normalized total energy and the total angular momentum of the cluster ####
    # getting the constants we will need
    n_stars = length(clust.K)

    # Getting Etot_norm
    J_typ = clust.G * clust.m_min^2 / clust.a_min # typical value of J (circular orbits and twice the same population) -- To be improved
    Etot_norm = clust.Etot / (J_typ * n_stars^2) # Our normalization is Etot_norm = Etot / (J_typ N^2)
    
    # Getting s (= Ltot_norm)
    Ltot_max = 0.0
    for k = 1:n_stars
        mₖ, aₖ, eₖ = clust.K[k]
        lₖ = mₖ * sqrt(clust.G * clust.M * aₖ * (1 - eₖ^2))
        Ltot_max += lₖ
    end
    s = clust.Ltot / Ltot_max # our normalization is s = Ltot / Ltot_max

    # Returning the tuple
    return((Etot_norm, s))
end

function plot_that_cluster(clust::Cluster_initial) 
    #### Returns a plot object. It represents the positions of the stars on the unit sphere, the color being their momentum ####
    plotlyjs()
    
    ### plot the sphere
    X(z,θ) = sqrt(1 - z^2) * cos(θ)
    Y(z,θ) = sqrt(1 - z^2) * sin(θ)
    Z(z,θ) = z
    C(z, θ) = 0.0 # to impose a uniform color

    zs = range(-1, 1, length = 50)
    θs = range(0, 2 * pi, length = 50)


    my_plot = surface(X.(zs', θs), Y.(zs', θs), Z.(zs', θs), 
                        color = "grey", alpha = 0.4, 
                        showaxis = false, grid = false)

    ### Plot the stars

    # first we need to compute their cartesian position, and momentums (for the colors)
    x = zeros(clust.n_disks * clust.N)
    y = zeros(clust.n_disks * clust.N)
    z = zeros(clust.n_disks * clust.N)
    Momentums = zeros(clust.n_disks * clust.N)

    for i = 1:clust.n_disks
        for j = 1:clust.N
            index = (i - 1) * clust.N + j

            m, a, e = clust.K[index]
            θ, ϕ = clust.L[index]

            x[index] = sin(θ) * cos(ϕ)
            y[index] = sin(θ) * sin(ϕ)
            z[index] = cos(θ)
            Momentums[index] = log10(m * sqrt(clust.G * clust.M * a * (1 - e^(2))))
        end
    end

    # now we can plot the stars as a scatter plot
    L_min = log10(clust.m_min * sqrt(clust.G * clust.M * clust.a_min * (1 - clust.e_max^(2))))
    L_max = log10(clust.m_max * sqrt(clust.G * clust.M * clust.a_max * (1 - clust.e_min^(2))))

    my_plot = scatter!(x, y, z,
                    xlims = (-1.0, 1.0), ylims = (-1.0, 1.0), zlims = (-1.0, 1.0),
                    zcolor = Momentums, color = :thermal, clims = (L_min, L_max),
                    markersize = 1.0, markerstrokewidth = 0,
                    cbar = false, leg = false)

    # now let's return the figure
    return(my_plot)
end