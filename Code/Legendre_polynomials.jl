function Legendre_polynomials!(x::Float64, lmax::Int64; Pₗ = nothing) 
    #### returns an array of size lmax_half. Each cell contains one (even) Legendre polynomials evaluated at point x ####
    # getting all the constant we will need
    lmax_half = div(lmax, 2)

    # Preallocating Pₗ if it hasn't been given
    if isnothing(Pₗ)
        Pₗ = zeros(lmax_half) # will store all the Legendre polynomials
    end

    # Let's start the recurrence relation -- We use Bonnet's recurrence formula
    P_l = nothing # these 3 lines just to help the interpreter with scopes
    P_lmin = nothing
    P_lminmin = nothing

    for l = 0:lmax
        if (l == 0)
            P_l = 1
        elseif (l == 1)
            P_l = x
        else
            P_l = x * P_lmin * (2 * l - 1) / l - P_lminmin * (l - 1) / l # l Pₗ = x * (2 l - 1) Pₗ₋₁ + (l - 1) Pₗ₋₂
        end

        if ((l % 2 == 0) && (l > 0)) # we only store the even harmonics, and l = 0 is useless because it brings a constant to the energy
            l_half = div(l, 2)
            Pₗ[l_half] = P_l
        end

        P_lminmin = P_lmin
        P_lmin = P_l
    end
end

function Normalized_associated_Legendre_polynomials!(x::Float64, lmax::Int64; P = nothing) 
    #### Returns an array of size lmax_half * (lmax + 1). Each cell contains one normalized associated Legendre polynomails Pᴺₗₘ(x)        ####
    #### The definition is Pᴺₗₘ = √((2l + 1) / (4 π)) × √(l - m)! / (l + m)! × Pₗₘ where Pₗₘ is the standard associated Legendre polynomial. ####
    # getting all the constants we will need
    lmax_half = div(lmax, 2)

    # Preallocating P if it hasn't been given
    if isnothing(P)
        P = zeros(lmax + 1, lmax_half) # Will store the normalized associated Legendre polynomials
    end

    # Let's start the recurrence -- We use the method from Numerical Recipes (eq. 6.7.9)
    Pₘₘ = nothing # this is just to help the compilator with scopes

    for m = 0:lmax
        # First let's get Pₘₘ by a first recurrence relation
        if (m == 0)
            Pₘₘ = sqrt(1 / (4 * π)) # P₀₀ = sqrt(1 / 4π)
        else
            Pₘₘ = -1.0 * sqrt((2 * m + 1) / (2 * m)) * sqrt(1 - x^(2)) * Pₘₘ # Pₘ,ₘ = - sqrt((2m + 1) / 2m) * sqrt(1 - x²) * Pₘ₋₁,ₘ₋₁
        end
        
        # The we go from m,m to m,l step by step
        Pₗₘ = nothing # this is just to help the compilator with scopes
        Pₗₘ_min = nothing
        Pₗₘ_minmin = nothing

        for l = m:lmax
            if (l == m)
                Pₗₘ = Pₘₘ
            elseif (l == (m + 1))
                Pₗₘ = x * sqrt(2 * m + 3) * Pₘₘ # Pₘ,ₘ₊₁ = x sqrt(2m + 3) Pₘₘ
            else
                Pₗₘ = sqrt((4 * l^2 - 1) / (l^2 - m^2)) * ( x * Pₗₘ_min - sqrt( ((l - 1)^2 - m^2) / (4 * (l - 1)^2 -1)) * Pₗₘ_minmin )
                # Pₘₗ = sqrt((4l² - 1) / (l² - m²)) * [x Pₘ,ₗ₋₁ - sqrt(((l - 1)² - m²) / (4 (l - 1)² - 1)) Pₘₗ₋₂]
            end

            if (((l % 2) == 0) && (l > 0)) # We store the coefficients only if l is even
                l_half = div(l, 2)
                P[m + 1, l_half] = Pₗₘ
            end

            Pₗₘ_minmin = Pₗₘ_min
            Pₗₘ_min = Pₗₘ
        end
    end

    return(P)
end