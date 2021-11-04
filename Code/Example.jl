using Random
Random.seed!(1)


### First let's create an initial cluster
include("./Cluster.jl")

lmax = 10

n_disks = 16
n_stars_per_disk = 64
σ = 0.006

m_min = 1.0
m_max = 100.0
γₘ = -2.0

a_min = 1.0
a_max = 100.0
γₐ = 0.0

e_min = 0.0
e_max = 0.3
γₑ = 1.0

clust = Cluster_initial(n_disks, n_stars_per_disk, lmax, 
                        σ = σ,
                        m_min = m_min, m_max = m_max, γₘ = γₘ,
                        a_min = a_min, a_max = a_max, γₐ = γₐ,
                        e_min = e_min, e_max = e_max, γₑ = γₑ)

# Then we can compute the (normalized) binding energy and spin of this cluster
Etot_norm, s = get_ES_from_clust(clust) # Here, Etot_norm = Etot / (J N^2) = -E_binding
                                        # J = G mₘᵢₙ² / aₘᵢₙ

# We can also plot it, it should reproduce the figure "Initial_cluster.png"
using Plots

my_plot_1 = plot_that_cluster(clust)

display(my_plot_1)



### Now let's create a distribution function that will mimick this cluster

# First we create the m, a, e bins, and fill them with stars i.e. create the N(K) vector
include("./Miscellaneous.jl")

n_pop_m = 10
n_pop_a = 3
n_pop_e = 1

Bins_m, Bins_a, Bins_e, K, L, N = create_KLN(m_min, m_max, γₘ, n_pop_m, 
                                            a_min, a_max, γₐ, n_pop_a, 
                                            e_min, e_max, γₑ, n_pop_e)

# Then we create the matrix of coupling coefficients J
include("./Coupling_coefficients_interpolated.jl")

J = Coupling_coefficients_interpolated(K, lmax)

# Then we translate the normalized binding energy and spin into the pure energy and angular momentum
include("./Distribution.jl")

Etot, Ltot = get_EL_for_dist(Etot_norm, s,
                                L, N, m_min, a_min)



### Now let's optimize this distribution's entropy, which is equivalent to relaxing the cluster
include("./Solve_by_energy_sequence.jl")
using LinearAlgebra

dist = solve_by_energy_sequence(J, L, N, Etot, Ltot)
println( norm(dist.cost) < 10^(-6) ) # This is to check that the optimization method converged



### Finally, we can plot this distribution's segregation in mass, it reproduce the figure "Relaxed_cluster.png"
my_plot_2 = plot_that_dist("m", # "m" means that we want to plot the segregation in mass. 
                                # "a" or "e" would give the segregations in semi-major axis and eccentricity
                        Etot_norm, s, lmax,
                        dist,
                        Bins_m, n_pop_m,
                        Bins_a, n_pop_a,
                        Bins_e, n_pop_e)

display(my_plot_2)