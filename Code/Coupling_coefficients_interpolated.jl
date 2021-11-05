using Interpolations
using HDF5
using StaticArrays

try
    global const DIRECTORY_SL = dirname(dirname(@__FILE__)) * "/Data/" # Directory in which the sl are stored
    namefile = DIRECTORY_SL * "data_sl_l_" * string(2) * ".hf5"
    file = h5open(namefile,"r")
    close(file)
catch
    println("File Coupling_coefficients_interpolated - ERROR - could not load the sl coefficients from the file")
end

const ALPHA_MIN = 10^(-3) # Minimum value of α over which we interpolated the sₗ
const ALPHA_MAX = 10^(0) # Maximum value of α over which we interpolated the sₗ
const ECC_MIN = 0.0 # Minimum value of e over which we interpolated the sₗ
const ECC_MAX = 0.99 # Maximum value of e over which we interpolated the sₗ
const LMAX = 50 # largest mode we can compute from the interpolated sₗ
const LMAX_HALF = div(50, 2)


function Interpolating_function_sl(l::Int64)
    #### Returns one interpolating function sₗ ####
    #Getting the name of the file
    namefile = DIRECTORY_SL*"data_sl_l_"*string(l)*".hf5" # Name of the file

    # Reading the file
    file = h5open(namefile,"r") # Opening the file in read-only mode
    tab_logalpha = read(file,"tab_logalpha") # Range of log(alpha)
    tab_ein  = read(file,"tab_ein")  # Range of ein
    tab_eout = read(file,"tab_eout") # Range of eout
    tab_sl = read(file,"tab_sl") # Table of sl
    close(file) # Closing the file

    # getting the bounds and steps
    logalpha_min, delta_logalpha, logalpha_max = tab_logalpha # Values for the range in log(alpha)
    ein_min,  delta_ein,  ein_max  = tab_ein  # Values for the range in ein
    eout_min, delta_eout, eout_max = tab_eout # Values for the range in eout

    # getting the lists of points for which we have exact values
    range_logalpha = logalpha_min:delta_logalpha:logalpha_max # Range in log(alpha)
    range_ein = ein_min:delta_ein:ein_max # Range in ein
    range_eout = eout_min:delta_eout:eout_max # Range in eout

    # getting the interpolation function -- 229 allocations / 100 MiB
    INT_sl = scale(interpolate(tab_sl, BSpline(Cubic(Line(OnGrid())))),
                   range_logalpha,range_ein,range_eout) # Making the 3D gridded interpolation

    # Returning the interpolation function
    return INT_sl
end


const TAB_INT_SL = SVector{LMAX_HALF}([Interpolating_function_sl(l) for l = 2:2:LMAX]) # This is a list of the interpolating function. It is constant that can be accessed by all scripts


function Coupling_coefficients_interpolated(K::Array{Tuple{Float64,Float64,Float64},1}, lmax::Int64; G::Float64 = 1.0) 
    #### returns an array of size n_pop * n_pop * lmax_half. Each cell contains an (even axisymmetric) coupling coefficient ####
    # getting all the constants we will need
    n_pop = length(K)
    lmax_half = div(lmax, 2)

    # preallocating all the arrays we will need
    J = zeros(n_pop, n_pop, lmax_half) # will store all the coupling coefficients

    # Let's filling up J
    for star_1 = 1:n_pop
        m, a, ecc = K[star_1] # (m,a,e) of the unprimed particle

        for star_2 = 1:star_1
            mp, ap, eccp = K[star_2]  # (m,a,e) of the primed particle
            
            if (a <= ap) # Unprimed particle is inner, and primed is outer
                ain, ein = a, ecc # Inner particle
                aout, eout = ap, eccp # Outer particle
            else # Primed particle is inner, and unprimed is outer
                ain, ein = ap, eccp # Inner particle
                aout, eout = a, ecc # Outer particle
            end
    
            alpha = ain / aout # Ratio of the semi-major axes, with alpha <= 1.0

            # Checking that (alpha,ein,eout) are within the interpolated ranges
            @assert ALPHA_MIN <= alpha <= ALPHA_MAX "get_Jl: alpha not in the interpolated range"
            @assert ECC_MIN   <= ein   <= ECC_MAX   "get_Jl: ein not in the interpolated range"
            @assert ECC_MIN   <= eout  <= ECC_MAX   "get_Jl: eout not in the interpolated range"

            for l_half = 1:lmax_half
                l = 2* l_half

                # Checking that l is within the interpolated range
                @assert 2 <= l <= LMAX "get_Jl: l not in the interpolation range"

                # Reading sl
                logalpha = log(alpha) # Value of log(alpha), i.e. the variable w.r.t. which the interpolation was performed
                l_half = div(l, 2) # index of mode l
                sl = TAB_INT_SL[l_half](logalpha, ein, eout) # Getting the value from the interpolated functions
            
                # Computing the coupling coefficient
                Jl = (sl * G * m * mp) / aout # Jₗ[K₁, K₂] = G m₁ m₂ sₗ[K₁, K₂] / a_out by definition

                # Storing the coupling coefficient -- We use the symmetry
                J[star_1, star_2, l_half] = Jl
                J[star_2, star_1, l_half] = Jl
            end
        end
    end

    return(J)
end

function Coupling_coefficient_interpolated_v2(l::Int64,
                                                K::Tuple{Float64,Float64,Float64}, Kp::Tuple{Float64,Float64,Float64};
                                                G::Float64 = 1.0)
    #### Returns a single (even axisymmetric) coupling coefficient ####

    # Giving a name to the annuli parameters
    m,  a,  ecc  = K  # (m,a,e) of the unprimed particle
    mp, ap, eccp = Kp # (m,a,e) of the primed particle
    
    if (a <= ap) # Unprimed particle is inner, and primed is outer
        ain, ein = a, ecc # Inner particle
        aout, eout = ap, eccp # Outer particle
    else # Primed particle is inner, and unprimed is outer
        ain, ein = ap, eccp # Inner particle
        aout, eout = a, ecc # Outer particle
    end
    
    alpha = ain / aout # Ratio of the semi-major axes, with alpha <= 1.0

    ##### Checking that (alpha,ein,eout) are within the interpolated ranges
    @assert ALPHA_MIN <= alpha <= ALPHA_MAX "get_Jl: alpha not in the interpolated range"
    @assert ECC_MIN   <= ein   <= ECC_MAX   "get_Jl: ein not in the interpolated range"
    @assert ECC_MIN   <= eout  <= ECC_MAX   "get_Jl: eout not in the interpolated range"
    @assert 2 <= l <= LMAX "get_Jl: l not in the interpolation range"

    # Reading sl
    logalpha = log(alpha) # Value of log(alpha), i.e. the variable w.r.t. which the interpolation was performed
    l_half = div(l, 2) # index of mode l
    sl = TAB_INT_SL[l_half](logalpha, ein, eout) # Getting the value from the interpolated functions

    # Computing the coupling coefficient
    Jl = (sl * G * m * mp) / aout # Jₗ[K₁, K₂] = G m₁ m₂ sₗ[K₁, K₂] / a_out by definition

    # Returning the coupling coefficient
    return Jl # Output
end
