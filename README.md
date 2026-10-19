# FoodWebGenerator.jl

[![CI](https://github.com/cvegao/FoodWebGenerator.jl/actions/workflows/test.yml/badge.svg)](https://github.com/cvegao/FoodWebGenerator.jl/actions/workflows/test.yml)

Generate and analyse model food webs in Julia with the **Preferential Preying Model** (PPM; Johnson et al., 2014; Klaise & Johnson, 2016) or the **Niche Model** (Williams & Martinez, 2000).

The PPM gives control over species richness, connectance, number of basal species and trophic coherence. The Niche Model gives control over richness and connectance.

## Features

- Two generative models: `:ppm` and `:niche`
- Trophic levels by linear-system solution (Klaise & Johnson, 2016)
- Trophic coherence `q` for any web (Johnson et al., 2014)
- Validated `FoodWeb` struct (square matrix, consistent lengths, valid connectance)
- Reproducible: every generator takes an integer seed (`Xoshiro`)
- Import your own adjacency matrix (and optionally trophic levels)
- Light dependencies: `LinearAlgebra`, `Random`, `Distributions`, `StatsBase`

## Installation

```julia
using Pkg
Pkg.add("FoodWebGenerator")   # if registered
# otherwise:
Pkg.add(url="[https://github.com/cvegao/FoodWebGenerator.jl](https://github.com/cvegao/FoodWebGenerator.jl)")
```

## Conventions

- `adj[i, j] = 1` means **species `i` eats species `j`** (rows are consumers).
- Basal species have an empty row and trophic level `0.0`.
- Connectance is `C = L / S²`.
- Trophic coherence `q = 0` is a perfectly coherent web.

## Quick start

```julia
using FoodWebGenerator

# PPM: S = 20, C = 0.15, 3 basal species, seed 42
fw = generate_food_web(20, 0.15, 3, 42)

# Niche model: S = 20, C = 0.15, seed 42
fw = generate_food_web(20, 0.15, 42, :niche)

println(fw)
# FoodWeb:
#   Species richness (S): 20
#   Connectance (C): …
#   Trophic coherence (q): …
#   Basal species indices: […]
#   Trophic levels: […]
#   Adjacency matrix (top left corner):
#   …

fw.adj_matrix      # S×S Int matrix
fw.trophic_levels  # trophic levels (basal = 0)
fw.coherence       # trophic coherence q
fw.connectance     # C = L / S²
fw.basal_species   # indices of basal species
```

## API

### `generate_food_web`

```julia
generate_food_web(S, C, seed, model; T=0.25, basal=nothing)   # (1)
generate_food_web(S, C, seed; model=:niche, T=0.25, basal=nothing)  # (2)
generate_food_web(S, C, basal, seed; model=:ppm, T=0.25)      # (3)
generate_food_web(adj_matrix)                                 # (4)
generate_food_web(adj_matrix, trophic_levels)                 # (5)
```

| Argument | Type | Description |
|---|---|---|
| `S` | `Int` | Species richness (≥ 1) |
| `C` | `AbstractFloat` | Target connectance in `(0, 1]` (`0.0` if `S = 1`); `C < 0.5` for the niche model |
| `seed` | `Int` | Seed for a new `Xoshiro` RNG |
| `model` | `Symbol` | `:ppm` or `:niche` |
| `basal` | `Int` | Number of basal species, `1 ≤ basal ≤ S`; **required for `:ppm`**, ignored by `:niche` |
| `T` | `AbstractFloat` | PPM temperature controlling coherence (default `0.25`); ignored by `:niche` |
| `adj_matrix` | `AbstractMatrix` | Own adjacency matrix (no self-loops, ≥ 1 basal species) |
| `trophic_levels` | `AbstractArray` | Own trophic levels (basal = 0) |

> **Defaults.** Method (2) defaults to `model=:niche`. For PPM, either use method (3) or pass `model=:ppm, basal=…` explicitly.

```julia
generate_food_web(50, 0.10, 5, 123)                         # PPM (3)
generate_food_web(50, 0.10, 123; model=:ppm, basal=5)       # PPM (2)
generate_food_web(50, 0.10, 123, :ppm; basal=5)             # PPM (1)
generate_food_web(50, 0.10, 123, :niche)                    # Niche (1)

# Tune coherence with T (PPM only)
fw_coherent   = generate_food_web(50, 0.10, 5, 1; T=0.1)    # lower q
fw_incoherent = generate_food_web(50, 0.10, 5, 1; T=1.0)    # higher q

# Own data
adj = [0 0 0; 1 0 0; 0 1 0]
fw  = generate_food_web(adj)
fw  = generate_food_web(adj, [0.0, 1.0, 2.0])
```

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

| Constructor | Description |
|---|---|
| `FoodWeb(adj, tl, q, conn, basal)` | All fields supplied; validated |
| `FoodWeb(adj, tl)` | Computes `q`, connectance and basal species |
| `FoodWeb(adj)` | Computes trophic levels and everything else |

Invariants: square matrix; `length(tl) == S`; `0 < connectance ≤ 1` (`0.0` if `S = 1`).

### `trophic_levels`

```julia
trophic_levels(adj_matrix::AbstractMatrix) -> Vector{Float64}
```

Solves `(D − A)h = D·1` with `D = diag(max(k_in, 1))` and shifts by −1 so basal species are `0.0` (values < 1e-2 are set to 0). Throws `ArgumentError` for self-loops or if there are no basal species.

### `trophic_coherence`

```julia
trophic_coherence(adj_matrix, trophic_levels) -> Float64
```

$$q = \sqrt{\left\langle (h_i - h_j - 1)^2 \right\rangle}$$

The average runs over all links `i → j`. This is the root-mean-square deviation of trophic distances from 1. If the minimum level is below 1, levels are shifted by +1 internally. It returns `0.0` for a web with no links.

## Algorithms

**Preferential Preying Model**
1. Basal species are seeded at level 1.
2. Each consumer `i` picks a random initial prey `j < i`.
3. It draws `x ~ Beta(1, β)` with `β = (S² − basal²)/(2L) − 1` and `L = round(C·S²)`, and adds `κ = min(round(x·(i−1)), available)` further prey.
4. Further prey are sampled without replacement with weights `exp(−|h_j − h_k| / T)`.
5. The consumer's level becomes `1 + mean(prey levels)`.
6. After all species are placed, trophic levels and `q` are recomputed exactly for the final web.

Low `T` concentrates links on prey at similar trophic level, so webs are more coherent.

**Niche Model**
1. Each species gets a niche value `n_i ~ U(0,1)`.
2. Its feeding range is `r_i = n_i · x`, with `x ~ Beta(1, 1/(2C) − 1)`.
3. The range centre is drawn from `U(r_i/2, n_i)`, and `i` eats all species whose niche value lies in its range.
4. The species with the smallest `n_i` has `r_i = 0`, so there is always a basal species.
5. Isolated species are redrawn until none remain.

## Testing

```julia
using Pkg; Pkg.test("FoodWebGenerator")
```

The tests cover the generators, the structural properties of the niche webs (no self-loops, binary matrix, at least one basal species, no isolated species), and mean connectance against the target.

## Dependencies

| Package | Role |
|---|---|
| `LinearAlgebra` | Trophic-level linear solve |
| `Random` | RNG and `Xoshiro` |
| `Distributions` | `Beta` and `Uniform` sampling |
| `StatsBase` | Weighted sampling without replacement |

## References

- Johnson, S., Domínguez-García, V., Donetti, L., and Muñoz, M. A. (2014). Trophic coherence determines food-web stability. *PNAS*, 111(50):17923–17928.
- Klaise, J. and Johnson, S. (2016). From neurons to epidemics: How trophic coherence affects spreading processes. *Chaos*, 26(6):065310.
- Williams, R. and Martinez, N. (2000). Simple rules yield complex food webs. *Nature*, 404(6774):180–183.

## License

MIT. See [LICENSE](LICENSE).