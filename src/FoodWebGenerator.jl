"""
    FoodWebGenerator
Module for generating and analyzing food webs using the Preferential Preying Model.
References:
- Klaise, J. and Johnson, S. (2016). From neurons to epidemics: How trophic coherence affects spreading processes. Chaos: An Interdisciplinary Journal of Nonlinear Science, 26(6):065310.
- Johnson, S., Domínguez-García, V., Donetti, L., and Muñoz, M. A. (2014). Trophic coherence determines food-web stability. 
Proceedings of the National Academy of Sciences, 111(50):17923–17928.
"""
module FoodWebGenerator

using LinearAlgebra, Random, Distributions, StatsBase
import Base: show

export FoodWeb, generate_food_web, trophic_levels, trophic_coherence

# --------------------------
#   Core Data Structures
# --------------------------

"""
    FoodWeb
Structure representing a food web with validated adjacency matrix and properties.

 Fields
- `adj_matrix::AbstractMatrix`: Consumer-resource interactions matrix
- `trophic_levels::AbstractArray`: Trophic levels of each species
- `coherence::Real`: Trophic coherence measure (q)
- `connectance:AbstractFloat`: C = L/S²  (L = number of links)
- `basal_species::AbstractArray`: Indices of basal species (TL = 0)

# Invariants
1. adjacency matrix must be square
2. trophic_levels length matches adjacency matrix dimensions
3. 0 < connectance ≤ 1
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
    trophic_levels(adj_matrix::AbstractMatrix{Int})
Calculate species trophic levels using linear system solution.

# Arguments
- `adj_matrix::AbstractMatrix`: Adjacency matrix where adj_matrix[i,j] == 1 indicates i consumes j

# Returns
- Vector of trophic levels where basal species = 0.0

# Notes
Implements the method form Klaise & Johnson (2016)
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

    return A \ b .- 1.0  # Adjust to have basal species at TL=0
end

"""
    generate_food_web(S::Int, C::AbstractFloat, basal::Int, seed::Int; T=0.25)
    generate_food_web(adj_matrix::AbstractMatrix)
    generate_food_web(adj_matrix::AbstractMatrix, trophic_levels::AbstractArray)
Generate a food web with specified richness, connectance, and basal species.

# Arguments
- `S::Int`: Number of species (nodes) in the network.
- `C::AbstractFloat`: Connectance (fraction of realized links).
- `basal`: Number of basal species (1 ≤ basal < S)
- `seed::Int`: Random seed for reproducibility. *(Method 1 only)*
- `adj_matrix::AbstractMatrix`: Predefined adjacency matrix. *(Method 2 and 3 only)*
- `trophic_levels::AbstractArray`: Predefined trophic levels. *(Method 3 only)*
- `T::AbstractFloat`: Temperature parameter controlling trophic coherence (default: 0.25).

# Returns
- `FoodWeb`: A validated food web structure with adjacency matrix, trophic levels, coherence, connectance, and basal species indices.
# Example
```
julia>fw = generate_food_web(10, 0.15, 2, 42) # Uses seed 42
julia>fw = generate_food_web([0 0; 0 1]) # Uses predefined adjacency matrix
julia>fw = generate_food_web([0 0; 0 1], [0.0, 1.0]) # Uses predefined adjacency matrix and trophic levels
```
"""
function generate_food_web(S::Int, C::AbstractFloat, basal::Int, seed::Int; T=0.25)
    # Parameter validation
    (S >= 1)  || throw(ArgumentError("Species richness must be ≥ 1"))
    if S == 1
        (C == 0.0) || throw(ArgumentError("For S=1, connectance must be 0.0"))
    else
        (0 < C <= 1.0) || throw(ArgumentError("Connectance must be in (0,1] for S>1"))
    end
    (1 <= basal <= S) || throw(ArgumentError("Invalid basal species number"))

    # Random.seed!(seed)
    rng = Xoshiro(seed)
    
    adj_matrix = zeros(Int64, S, S)
    tl = zeros(S)
    tl[1:basal] .= 1.0

    L = round(Int, C * S^2)
    β = ((S^2 - basal^2) / (2 * L)) - 1

    for i in (basal + 1) : S
        # Preferential attachment core
        # connect one random node to consumer i
        n = i - 1  # ensures no self-loop

        j = rand(rng, 1:n)  # initial prey
        # j = rand(1:n)  # initial prey

        adj_matrix[i, j] = 1.0

        # Calculate niche probabilities
        # find available prey excluding itself
        available = findall(iszero, @view adj_matrix[i, 1:n])

        # calculate number of additional prey
        x = rand(rng, Beta(1, β))
        # x = rand(Beta(1, β)) 

        κ = min(round(Int, x * n), length(available))

        if κ > 0 && !isempty(available)
            # Calculate probabilities using existing trophic levels
            P = exp.(-abs.(tl[j] .- tl[available]) ./ T)

            # Weighted selection without replacement
            # selected = sample(available, Weights(P), κ; replace=false)
            selected = sample(rng, available, Weights(P), κ; replace=false)

            adj_matrix[i, selected] .= 1.0
        end

        # Update trophic levels incrementally
        k_in = sum(adj_matrix, dims=2)  # in-degree
        tl[i] = 1 + sum(adj_matrix[i, 1:i] .* tl[1:i]) / k_in[i]
    end

    # Post-generation validation and packaging
    tl_final = trophic_levels(adj_matrix)
    q = trophic_coherence(adj_matrix, tl_final)
    basal_species = findall(==(0), sum(adj_matrix, dims=2)[:])

    return FoodWeb(adj_matrix, tl_final, q, C, basal_species)
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
    trophic_coherence(adj::AbstractMatrix, tl::AbstractArray)
Calculate trophic coherence (q) from adjacency matrix and trophic levels.

# Arguments
- `adj_matrix::AbstractMatrix`: Food web adjacency matrix
- `tl::AbstractArray`: Precomputed trophic levels

# Returns
- Trophic coherence measure q (q=0 maximally coherent)

# Notes
Implements the formula from Johnson et al. (2014):
q = √(⟨(h_i - h_j - 1)^2⟩) for all i→j links
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