"""
    FoodWebGenerator

Generate and analyse model food webs with the Preferential Preying Model (PPM)
or the Niche Model.

# Conventions
- `adj[i, j] == 1` means species `i` **consumes** species `j` (rows = consumers).
- Basal species have an all-zero row and trophic level `0.0`.
- Connectance is `C = L / S²`, where `L` is the number of links.

# Exports
`FoodWeb`, `generate_food_web`, `trophic_levels`, `trophic_coherence`

# References
- Williams, R. and Martinez, N. (2000). Simple rules yield complex food webs. *Nature*, 404:180–183.
- Johnson, S., Domínguez-García, V., Donetti, L., and Muñoz, M. A. (2014). Trophic coherence
  determines food-web stability. *PNAS*, 111(50):17923–17928.
- Klaise, J. and Johnson, S. (2016). From neurons to epidemics: How trophic coherence affects
  spreading processes. *Chaos*, 26(6):065310.
"""
module FoodWebGenerator

using LinearAlgebra, Random, Distributions, StatsBase
import Base: show

include("niche_model.jl")
include("ppm_model.jl")

export FoodWeb, generate_food_web, trophic_levels, trophic_coherence

# --------------------------
#   Core Data Structures
# --------------------------
"""
    FoodWeb

Food web with adjacency matrix and derived structural properties.

# Fields
- `adj_matrix::AbstractMatrix`   : square consumer–resource matrix; `adj_matrix[i, j] == 1` if `i` eats `j`.
- `trophic_levels::AbstractArray`: trophic level of each species (basal = 0).
- `coherence::Real`              : trophic coherence `q` (0 = maximally coherent).
- `connectance::AbstractFloat`   : `C = L / S²`.
- `basal_species::AbstractArray` : indices of species with no prey.

# Constructors
    FoodWeb(adj, tl, q, conn, basal)   # all fields supplied and validated
    FoodWeb(adj, tl)                   # computes q, connectance and basal species
    FoodWeb(adj)                       # also computes trophic levels

# Invariants (checked by the full constructor)
1. `adj` is square.
2. `length(tl) == size(adj, 1)`.
3. `0 < connectance ≤ 1` (exactly `0.0` when `S == 1`).

# Examples
```julia
fw = FoodWeb([0 0 0; 1 0 0; 0 1 0])
fw.trophic_levels   # [0.0, 1.0, 2.0]
fw.basal_species    #[1]
```
"""
struct FoodWeb
    adj_matrix::AbstractMatrix
    trophic_levels::AbstractArray
    coherence::Real
    connectance::AbstractFloat
    basal_species::AbstractArray

    function FoodWeb(adj, tl, q, conn, basal)
        size(adj, 1) == size(adj,2) || throw(ArgumentError("Adjacency matrix must be square"))
        length(tl) == size(adj, 1) || throw(DimensionMismatch("Trophic levels vs adjacency matrix size"))
        if size(adj, 1) == 1
            (conn == 0.0) || throw(ArgumentError("For S=1, connectance must be 0.0"))
        else
            (0 < conn <= 1.0) || throw(ArgumentError("Connectance must be in (0,1] for S>1"))
        end
        new(adj, tl, q, conn, basal)
    end

    function FoodWeb(adj, tl)
        size(adj, 1) == size(adj,2) || throw(ArgumentError("Adjacency matrix must be square"))
        length(tl) == size(adj, 1) || throw(DimensionMismatch("Trophic levels vs adjacency matrix size"))
        
        q = trophic_coherence(adj, tl)
        c = sum(adj) / (size(adj, 1) ^ 2)
        basal = findall(==(0), sum(adj, dims=2)[:])

        new(adj, tl, q, c, basal)
    end

    function FoodWeb(adj)
        size(adj, 1) == size(adj,2) || throw(ArgumentError("Adjacency matrix must be square"))
        
        tl = trophic_levels(adj)
        q = trophic_coherence(adj, tl)
        c = sum(adj) / (size(adj, 1) ^ 2)
        basal = findall(==(0), sum(adj, dims=2)[:])

        new(adj, tl, q, c, basal)
    end
end

# ----------------------------------
#    Trophic Level Calculations
# ----------------------------------
"""
    trophic_levels(adj_matrix::AbstractMatrix) -> Vector{Float64}

Compute trophic levels by solving the linear system of Levine (1980) as used in
Klaise & Johnson (2016): `(D - A) h = D·1`, where `D = diag(max(k_in, 1))`.
The result is shifted by −1 so that basal species sit at `0.0`; values below
`1e-2` are set to `0.0` to remove numerical noise.

# Arguments
- `adj_matrix`: `adj_matrix[i, j] == 1` means `i` consumes `j`.

# Returns
Vector of trophic levels (basal = 0.0, a consumer of basal species only = 1.0).

# Throws
- `ArgumentError` if the diagonal contains non-zero entries (self-loops).
- `ArgumentError` if there are no basal species (no solution exists).

# Example
```julia
trophic_levels([0 1; 0 0])   # [1.0, 0.0]
```
"""
function trophic_levels(adj_matrix::AbstractMatrix)
    # Validate: no self-loops
    if any(adj_matrix[i,i] != 0 for i in axes(adj_matrix, 1))
        throw(ArgumentError("Adjacency matrix contains self-loops"))
    end

    in_degree = sum(adj_matrix, dims=2)[:]
    basal_species = findall(==(0), in_degree)

    if isempty(basal_species)
        throw(ArgumentError("Food web constains no basal species (cycle detected)"))
    end

    # Ah = b; h: trophic levels
    b = copy(in_degree)
    b[basal_species] .= 1
    A = diagm(b) .- adj_matrix

    tls = A \ b .- 1.0  # Adjust to have basal species at TL=0

    threshold = 1e-2  # TL < 1e-2 -> TL = 0
    tls[tls .< threshold] .= 0.0
    return tls
end

"""
    generate_food_web(S, C, seed, model; T=0.25, basal=nothing)
    generate_food_web(S, C, seed; model=:niche, T=0.25, basal=nothing)
    generate_food_web(S, C, basal, seed; model=:ppm, T=0.25)
    generate_food_web(adj_matrix)
    generate_food_web(adj_matrix, trophic_levels)

Build a [`FoodWeb`](@ref) from a model or from a user-supplied matrix.

# Arguments
- `S::Int`                       : species richness (≥ 1).
- `C::AbstractFloat`             : target connectance `L/S²`, in `(0, 1]` (exactly `0.0` if `S == 1`).
  The niche model additionally needs `C < 0.5`.
- `seed::Int`                    : seed for a fresh `Xoshiro` generator (reproducible output).
- `model::Symbol`                : `:ppm` or `:niche`.
- `basal::Int`                   : number of basal species, `1 ≤ basal ≤ S`. **Required for `:ppm`**, ignored for `:niche`.
- `T::AbstractFloat`             : PPM temperature; low `T` gives more coherent webs (default `0.25`). Ignored for `:niche`.
- `adj_matrix::AbstractMatrix`   : user adjacency matrix (no self-loops, ≥ 1 basal species).
- `trophic_levels::AbstractArray`: user trophic levels (basal = 0).

# Defaults
The `(S, C, seed; …)` method uses `model=:niche`; the `(S, C, basal, seed; …)` method uses `model=:ppm`.

# Returns
A `FoodWeb`.

# Throws
`ArgumentError` for an invalid `model`, `S`, `C` or `basal`, or for a missing `basal` with `:ppm`.

# Examples
```julia
fw = generate_food_web(20, 0.15, 3, 42)                          # PPM, 3 basal species, seed 42
fw = generate_food_web(20, 0.15, 42; model=:ppm, basal=3, T=0.5)  # same model, keyword form
fw = generate_food_web(20, 0.15, 42, :niche)                      # niche model
fw = generate_food_web([0 0 0; 1 0 0; 0 1 0])                     # from a matrix
```
"""
function generate_food_web(S::Int, C::AbstractFloat, seed::Int, model::Symbol; T=0.25, basal=nothing)
    # Parameter validation
    (model in [:niche, :ppm]) || throw(ArgumentError("`model` can take the values `:niche` or `:ppm`."))
    (S >= 1)  || throw(ArgumentError("Species richness must be ≥ 1"))
    if S == 1
        (C == 0.0) || throw(ArgumentError("For S=1, connectance must be 0.0"))
    else
        (0 < C <= 1.0) || throw(ArgumentError("Connectance must be in (0,1] for S>1"))
    end

    rng = Xoshiro(seed)
    if model == :ppm
        (!isnothing(basal)) || throw(ArgumentError("Basal species number (`basal`) must be provided when using the Preferential Preying Model (PPM)"))
        (1 <= basal <= S) || throw(ArgumentError("Invalid basal species number"))
        adj_matrix = ppm(S, C, basal, rng; T=T)
    else
        adj_matrix = niche(S, C, rng)
    end
    
    # Post-generation validation and packaging
    tl_final      = trophic_levels(adj_matrix)
    q             = trophic_coherence(adj_matrix, tl_final)
    basal_species = findall(==(0), sum(adj_matrix, dims=2)[:])

    return FoodWeb(adj_matrix, tl_final, q, sum(adj_matrix) / S^2, basal_species)

end

function generate_food_web(S::Int, C::AbstractFloat, seed::Int; model=:niche, T=0.25, basal=nothing)
    # Parameter validation
    (model in [:niche, :ppm]) || throw(ArgumentError("`model` can take the values `:niche` or `:ppm`."))
    (S >= 1)  || throw(ArgumentError("Species richness must be ≥ 1"))
    if S == 1
        (C == 0.0) || throw(ArgumentError("For S=1, connectance must be 0.0"))
    else
        (0 < C <= 1.0) || throw(ArgumentError("Connectance must be in (0,1] for S>1"))
    end

    rng = Xoshiro(seed)
    if model == :ppm
        (!isnothing(basal)) || throw(ArgumentError("Basal species number (`basal`) must be provided when using the Preferential Preying Model (PPM)"))
        (1 <= basal <= S) || throw(ArgumentError("Invalid basal species number"))
        adj_matrix = ppm(S, C, basal, rng; T=T)
    else
        adj_matrix = niche(S, C, rng)
    end
    
    # Post-generation validation and packaging
    tl_final      = trophic_levels(adj_matrix)
    q             = trophic_coherence(adj_matrix, tl_final)
    basal_species = findall(==(0), sum(adj_matrix, dims=2)[:])

    return FoodWeb(adj_matrix, tl_final, q, sum(adj_matrix) / S^2, basal_species)
end

function generate_food_web(S::Int, C::AbstractFloat, basal::Int, seed::Int; model=:ppm, T=0.25)
    # Parameter validation
    (model in [:niche, :ppm]) || throw(ArgumentError("`model` can take the values `:niche` or `:ppm`."))
    (S >= 1)  || throw(ArgumentError("Species richness must be ≥ 1"))
    if S == 1
        (C == 0.0) || throw(ArgumentError("For S=1, connectance must be 0.0"))
    else
        (0 < C <= 1.0) || throw(ArgumentError("Connectance must be in (0,1] for S>1"))
    end

    rng = Xoshiro(seed)
    if model == :ppm
        (!isnothing(basal)) || throw(ArgumentError("Basal species number (`basal`) must be provided when using the Preferential Preying Model (PPM)"))
        (1 <= basal <= S) || throw(ArgumentError("Invalid basal species number"))
        adj_matrix = ppm(S, C, basal, rng; T=T)
    else
        adj_matrix = niche(S, C, rng)
    end
    
    # Post-generation validation and packaging
    tl_final      = trophic_levels(adj_matrix)
    q             = trophic_coherence(adj_matrix, tl_final)
    basal_species = findall(==(0), sum(adj_matrix, dims=2)[:])

    return FoodWeb(adj_matrix, tl_final, q, sum(adj_matrix) / S^2, basal_species)
end

function generate_food_web(adj_matrix::AbstractMatrix)
    return FoodWeb(adj_matrix)
end

function generate_food_web(adj_matrix::AbstractMatrix, trophic_levels::AbstractArray)
    return FoodWeb(adj_matrix, trophic_levels)
end

# --------------------------------
#     Coherence Calculations
# --------------------------------
"""
    trophic_coherence(adj_matrix, trophic_levels) -> Float64

Trophic coherence `q` (Johnson et al. 2014):

    q = sqrt(mean((h_i - h_j - 1)²))   over all links i → j

i.e. the root-mean-square deviation of trophic distances from 1.
`q = 0` is a perfectly coherent web; larger `q` is less coherent.

If `minimum(trophic_levels) < 1`, levels are assumed to be basal = 0 and are shifted by +1
internally. Returns `0.0` when there are no links.
"""
function trophic_coherence(adj_matrix::AbstractMatrix, trophic_levels::AbstractArray)
    # Handle trophic levels vector with tl[basal_species] == 0.0 (as generated by generate_food_web)
    if minimum(trophic_levels) < 1.0
        trophic_levels = trophic_levels .+ 1.0
    end

    sum_sq_diff = 0.0
    link_count = 0

    for i in axes(adj_matrix, 1)
        for j in axes(adj_matrix, 2)
            if adj_matrix[i, j] > 0
                diff = trophic_levels[i] - trophic_levels[j] - 1
                sum_sq_diff += diff^2
                link_count += 1
            end
        end
    end

    if link_count == 0
        return 0.0  # No links, coherence is zero by convention
    end

    q = sqrt(sum_sq_diff / link_count)

    return q
end

function show(io::IO, fw::FoodWeb)
    println(io, "FoodWeb:")
    println(io, "  Species richness (S): ", size(fw.adj_matrix, 1))
    println(io, "  Connectance (C): ", fw.connectance)
    println(io, "  Trophic coherence (q): ", fw.coherence)
    println(io, "  Basal species indices: ", fw.basal_species)
    println(io, "  Trophic levels: ", fw.trophic_levels)
    println(io, "  Adjacency matrix (top left corner):")
    S = size(fw.adj_matrix, 1)
    max_S = min(S, 10) # Show up to 10x10 matrix
    for i in 1:max_S
        println(io, "    ", fw.adj_matrix[i, 1:max_S])
    end
    if S > 10
        println(io, "    ...")
    end
end
end