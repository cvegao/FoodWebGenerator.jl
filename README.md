# FoodWebGenerator.jl
--
[![CI](https://github.com/cvegao/FoodWebGenerator.jl/actions/workflows/test.yml/badge.svg)](https://github.com/cvegao/FoodWebGenerator.jl/actions/workflows/test.yml)
--

A Julia module for generating and analyzing ecological food webs using the **Preferential Preying Model** (default) or the **Niche Model**. If using the Preferential Preying Model, the package produces networks with controllable species richness, connectance, number of basal species and trophic coherence, following the methodology of Klaise & Johnson (2016) and Johnson et al. (2014). The methodology presented in Williams and Martinez (200) is used if the Niche Model is chosen.

## Features

- Generate food webs with specified species richness, connectance, and basal species count (if PPM is chosen)
- Compute **trophic levels** via linear system solution (Klaise & Johnson, 2016)
- Compute **trophic coherence** (q) for any food web (Johnson et al., 2014)
- Validated `FoodWeb` struct with automatic structural consistency checks
- Reproducible generation via seed or custom `AbstractRNG`
- Lightweight dependency footprint: `LinearAlgebra`, `Random`, `Distributions`, `StatsBase`

## Installation

This module is registered in the Julia General Registry.

```julia
using Pkg
Pkg.add("FoodWebGenerator")
```

Or, in the Julia REPL package mode (`]`):

```
pkg> add FoodWebGenerator
```

## Quick Start

```julia
using FoodWebGenerator

# Generate a food web with 20 species, connectance 0.15, 3 basal species
fw = generate_food_web(20, 0.15, 42, basal=3)   # seed = 42

println(fw)
# FoodWeb:
#   Species richness (S): 20
#   Connectance (C):      0.15
#   Trophic coherence (q): 0.18...
#   Basal species indices: [1, 2, 3]
#   Trophic levels: [0.0, 0.0, 0.0, ...]
#   Adjacency matrix (top left corner):
#   ...

# Generate a food web with 20 species and connectance 0.15 using the Niche Model
fw = generate_food_web(20, 0.15, 42, :niche)   # seed = 42

println(fw)
# FoodWeb:
#   Species richness (S): 20
#   Connectance (C):      0.15
#   Trophic coherence (q): 0.9090163521404221...
#   Basal species indices: [1, 3, 9, 10, 12, 18]
#   Trophic levels: [0.0, 1.0, 0.0, 1.0, 2.0, ...]
#   Adjacency matrix (top left corner):
#   ...

# Access fields directly
fw.adj_matrix      # S×S Int adjacency matrix
fw.trophic_levels  # Vector{Float64} of trophic levels
fw.coherence       # Float64, trophic coherence q
fw.connectance     # Float64, C = L / S²
fw.basal_species   # Vector{Int}, indices of basal species
```

## API Reference

### `FoodWeb`

```julia
struct FoodWeb
    adj_matrix     :: AbstractMatrix
    trophic_levels :: AbstractArray
    coherence      :: Real
    connectance    :: AbstractFloat
    basal_species  :: AbstractArray
end
```

A validated structure representing a food web. Three constructors are available:

| Constructor | Description |
|---|---|
| `FoodWeb(adj, tl, q, conn, basal)` | Full constructor with all fields |
| `FoodWeb(adj, tl)` | Computes `q`, `connectance`, and `basal_species` automatically |
| `FoodWeb(adj)` | Computes all derived quantities from the adjacency matrix |

**Invariants enforced at construction:**
- Adjacency matrix must be square
- Length of `trophic_levels` must match matrix dimensions
- Connectance must be in `(0, 1]` (or exactly `0.0` for `S=1`)

---

### `generate_food_web`

```julia
generate_food_web(S, C, seed, method=:ppm; T=0.25, basal=nothing)  -> FoodWeb
generate_food_web(adj_matrix)                                      -> FoodWeb
generate_food_web(adj_matrix, trophic_levels)                      -> FoodWeb
```

Generate a food web using the Preferential Preying Model.

| Parameter       | Type            | Description                                                                      |
|-----------------|-----------------|----------------------------------------------------------------------------------|
| `S`             | `Int`           | Number of species (nodes)                                                        |
| `C`             | `AbstractFloat` | Target connectance — fraction of realized links (`L / S²`)                       |
| `basal`         | `Int`           | Number of basal species (producers), `1 ≤ basal < S` (when `model=:ppm`)         |
| `seed`          | `Int`           | Integer seed for a fresh `Xoshiro` RNG (method 1)                                |
| `model`         | `Symbol`        | Model to be used to generate the adjacency matrix (method 1; `:ppm` or `:niche`) |
| `T`             | `AbstractFloat` | Temperature parameter controlling trophic coherence (default: `0.25`)            |
| `adj_matrix`    | `AbstractMatrix`| Predefined binary adjacency matrix (method 2 & 3)                                |
| `trophic_level` | `AbstractArray` | Predefined trophic levels (TL ≥ 0) (method 3).                                   |

**Returns:** a validated `FoodWeb` object.

```julia
# Method 1
fw = generate_food_web(50, 0.10, 123; basal=10)        # PPM
fw = generate_food_web(50, 0.10, 123, :ppm; basal=10)  # PPM
fw = generate_food_web(50, 0.10, 123, :niche)          # Niche model 

# Adjust trophic coherence via temperature (only PPM)
fw_coherent    = generate_food_web(50, 0.10, 1; T=0.1 basal=5)   # more coherent
fw_incoherent  = generate_food_web(50, 0.10, 1; T=1.0, basal=5)  # less coherent

# Method 2
adj_matrix = [0 0 0; 1 0 0; 0 1 0]
fw         = generate_food_web(adj_matrix)

# Method 3
adj_matrix     = [0 0 0; 1 0 0; 0 1 0]
trophic_levels = [0, 1, 2]
fw             = generate_food_web(adj_matrix, trophic_levels)
```

---

### `trophic_levels`

```julia
trophic_levels(adj_matrix::AbstractMatrix) -> AbstractArray
```

Compute trophic levels for all species from an adjacency matrix, where `adj_matrix[i, j] == 1` indicates species `i` consumes species `j`. Basal species are assigned trophic level `0.0`. Throws `ArgumentError` if no basal species are found (i.e., a cycle is detected).

---

### `trophic_coherence`

```julia
trophic_coherence(adj_matrix::AbstractMatrix, tl::AbstractArray) -> AbstractFloat
```

Compute trophic coherence `q` — the standard deviation of trophic distances across all links:

$$q = \sqrt{\langle (h_i - h_j - 1)^2 \rangle}$$

where the average is taken over all directed links `i → j`. A value of `q = 0` indicates a maximally coherent food web.

## Algorithm

`generate_food_web` implements the **Preferential Preying Model** (PPM):

1. Basal species are seeded at trophic level 1.
2. Each new consumer `i` is connected to a random initial prey `j` (no self-loops).
3. Additional prey are sampled *without replacement* using weights proportional to `exp(-|h_j - h_prey| / T)`, which biases prey choice toward species at similar trophic levels.
4. The number of additional prey follows a `Beta(1, β)` distribution, where β is derived from `S`, `C`, and the number of basal species.
5. After all species are placed, trophic levels are recomputed exactly via linear algebra, and the final `FoodWeb` object is validated and returned.

The temperature parameter `T` directly controls trophic coherence: lower `T` produces more coherent webs (links tend to span exactly one trophic level).

## Dependencies

| Package | Role |
|---|---|
| `LinearAlgebra` | Trophic level linear system solve |
| `Random` | RNG interface and `Xoshiro` |
| `Distributions` | `Beta` distribution for prey count sampling |
| `StatsBase` | Weighted sampling without replacement |

## References
- Johnson, S., Domínguez-García, V., Donetti, L., and Muñoz, M. A. (2014). *Trophic coherence determines food-web stability.* Proceedings of the National Academy of Sciences, 111(50):17923–17928.
- Klaise, J. and Johnson, S. (2016). *From neurons to epidemics: How trophic coherence affects spreading processes.* Chaos: An Interdisciplinary Journal of Nonlinear Science, 26(6):065310.
- Williams, R. and Martinez, N. (2000). *Simple Rules Yield Complex Food Webs*. Nature 404, n.º 6774: 180-83.

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
