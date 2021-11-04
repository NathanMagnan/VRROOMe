function create_Bins!(x_min::Float64, x_max::Float64, n_pop::Int64, Bins::Array{Tuple{Float64,Float64}, 1}, type::String)
    #### Fill the arrays Bins with tuples of 2 Floats. Each tuple contains the lower and upper limit of the bin ####
    #### type can be either "log" of "lin", to create logarithmic (m, a) or linear (e) bins                     ####
    for k = 1:n_pop
        if (type == "log")
            x_down = 10^(log10(x_min) + (log10(x_max) - log10(x_min)) * (k - 1) / n_pop)
            x_up = 10^(log10(x_min) + (log10(x_max) - log10(x_min)) * k / n_pop)
        elseif (type == "lin")
            x_down = x_min + (x_max - x_min) * (k - 1) / n_pop
            x_up = x_min + (x_max - x_min) * k / n_pop
        else
            return("type sould be either lin or log")
        end

        Bins[k] = (x_down, x_up)
    end
end

function create_KLN(m_min::Float64, m_max::Float64, γₘ::Float64, n_pop_m::Int64, 
                        a_min::Float64, a_max::Float64, γₐ::Float64, n_pop_a::Int64, 
                        e_min::Float64, e_max::Float64, γₑ::Float64, n_pop_e::Int64;
                        G::Float64 = 1.0, M::Float64 = 1.0)
    #### Returns a tuple of 6 Arrays. They are boundaries of the bins of stellar population used by the distribution, #### 
    #### and the K, L, N Arrays corresponding to these bins                                                           ####
    # getting the constants we will need
    n_pop = n_pop_m * n_pop_a * n_pop_e # Total number of populations

    # Preallocating the arrays we will need
    Bins_m = Array{Tuple{Float64,Float64}, 1}(undef, n_pop_m) # Will hold the bottom and top of the bins in mass
    Bins_a = Array{Tuple{Float64,Float64}, 1}(undef, n_pop_a) # Will hold the bottom and top of the bins in semi-major axis
    Bins_e = Array{Tuple{Float64,Float64}, 1}(undef, n_pop_e) # Will hold the bottom and top of the bins in eccentricity
    K = Array{Tuple{Float64,Float64,Float64}, 1}(undef, n_pop) # Will hold the m, a, e of each stellar populations
    L = zeros(n_pop) # Will hold the angular momentum of each stellar populations
    N = zeros(n_pop) # Will hold the number of particles in each stellar population

    # Let's compute the boundaries of the bins
    create_Bins!(m_min, m_max, n_pop_m, Bins_m, "log")
    create_Bins!(a_min, a_max, n_pop_a, Bins_a, "log")
    create_Bins!(e_min, e_max, n_pop_e, Bins_e, "lin")

    # Let's fill the bins with stars
    for k_m = 1:n_pop_m
        m_down, m_up = Bins_m[k_m]
        mₖ = 10^((log10(m_down) + log10(m_up)) / 2)
        vₘₖ = m_up - m_down # Volume of the bin in mass
        if (vₘₖ == 0.0)
            vₘₖ = 1.0
        end

        for k_a = 1:n_pop_a
            a_down, a_up = Bins_a[k_a]
            aₖ = 10^((log10(a_down) + log10(a_up)) / 2)
            vₐₖ = a_up - a_down # Volume of the bin in sma
            if (vₐₖ == 0.0)
                vₐₖ = 1.0
            end

            for k_e = 1:n_pop_e
                e_down, e_up = Bins_e[k_e]
                eₖ = (e_down + e_up) / 2
                vₑₖ = e_up - e_down # volume of the bin in eccentricity
                if (vₑₖ == 0.0)
                    vₑₖ = 1.0
                end

                index = k_e + n_pop_e * (k_a - 1) + n_pop_e * n_pop_a * (k_m - 1)

                K[index] = (mₖ, aₖ, eₖ)
                L[index] = mₖ * sqrt(G * M * aₖ * (1 - eₖ^2))
                N[index] = vₘₖ * vₐₖ * vₑₖ * mₖ^(γₘ) * aₖ^(γₐ) * eₖ^(γₑ) # population = volume * density
            end
        end
    end
    N /= sum(N) # We normalize all distributions to have 1 star in total

    # Return
    return((Bins_m, Bins_a, Bins_e, K, L, N))
end

function get_EL_for_dist(Etot_norm::Float64, s::Float64, 
                        L::Array{Float64, 1}, N::Array{Float64, 1}, m_min::Float64, a_min::Float64;
                        G::Float64 = 1.0)
    #### Return a tuple of 2 Floats. They are the the unnormalized total energy and momentum that can be used by a distribution ####
    #### NOTE : the cluster provided normalized energy and momentum, but the distribution works with unnormalized ones          ####
    # getting the constants we will need
    n_pop = length(L)

    # E
    J_typ = G * m_min^2 / a_min # typical value of J (circular orbits and twice the same population) -- To be improved
    Etot = Etot_norm * J_typ # Etot_norm = Etot / (J N^2) but a distribution only has 1 star

    # L
    Ltot_max_dist = 0.0
    for k = 1:n_pop
        Ltot_max_dist += L[k] * N[k]
    end
    Ltot = s * Ltot_max_dist # Ltot_norm = s = Ltot / Ltot_max

    # Return
    return((Etot, Ltot))
end