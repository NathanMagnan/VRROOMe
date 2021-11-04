include("./Legendre_polynomials.jl")

function Spherical_harmonics_a!(x::Float64, lmax::Int64;
                                Yₗ = nothing, Pₗ = nothing) 
    #### Returns an array of size lmax_half. Each cell contains one (even axisymmetric) spherical harmonic at point cos(θ) = x ####
    #### The "a" in the name is here to outline that this function only gives the "a"xisymmetric harmonics (i.e. m = 0).         ####
    # getting the constants we will need
    lmax_half = div(lmax, 2)

    # let's preallocate the arrays that haven't been given
    if isnothing(Yₗ)
        Yₗ = zeros(lmax_half) # WIll hold the spherical harmonics
    end
    if isnothing(Pₗ)
        Pₗ = zeros(lmax_half) # WIll hold the Legendre polynomials
    end

    # let's comptute the (even) Legendre polynomials at point x
    Legendre_polynomials!(x, lmax, Pₗ = Pₗ)

    # let's compute the spherical harmonics
    for l_half = 1:lmax_half
        l = 2 * l_half
        Yₗ[l_half] = Pₗ[l_half] * sqrt((2 * l + 1) / (4 * π)) # Yₗ₀ = Pₗ₀ * sqrt((2 l + 1) / 4 π)
    end

    return(Yₗ)
end

function Spherical_harmonics!(θ::Float64, ϕ::Float64, lmax::Int64; 
                                Yₗₘ = nothing, Pₗₘ = nothing) 
    #### Returns an array of size lmax_half * (2 lmax + 1). Each cell contains one (even) spherical harmonic at point (θ, ϕ) ####
    #### This function has no "a" in its title because it returns ALL the spherical harmonics (i.e. m = 0 and m != 0)        ####
    # getting the constants we will need
    lmax_half = div(lmax, 2)

    # let's preallocate the arrays that haven't been given
    if isnothing(Yₗₘ)
        Yₗₘ = zeros(2 * lmax + 1, lmax_half) # Will hold the spherical harmonics
    end
    if isnothing(Pₗₘ)
        Pₗₘ = zeros(lmax + 1, lmax_half) # Will hold the normalized associated Legendre polynomials
    end

    # Let's compute the normalzed associated Legendre polynomials
    Normalized_associated_Legendre_polynomials!(cos(θ), lmax, P = Pₗₘ)

    # Let's compute the spherical harmonics
    for l_half = 1:lmax_half
        l = 2 * l_half

        for m = 0:l
            # The passage from Pₗₘ to Yₗₘ is not the same if m is positive, null or negative. 
            # Cf Wikipedia : https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
            if (m != 0)
                Yₗₘ[lmax + 1 - m, l_half] = (-1)^(m) * sqrt(2) * Pₗₘ[m + 1, l_half] * sin(m * ϕ)
                Yₗₘ[lmax + 1 + m, l_half] = (-1)^(m) * sqrt(2) * Pₗₘ[m + 1, l_half] * cos(m * ϕ)
            else
                Yₗₘ[lmax + 1, l_half] = Pₗₘ[1, l_half]
            end
        end
    end

    return(Yₗₘ)
end