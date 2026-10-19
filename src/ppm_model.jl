"""
    ppm(S, C, basal, rng=Random.default_rng(); T=0.25) -> Matrix{Int}

Preferential Preying Model (Johnson et al. 2014). Species `1:basal` are basal (internal level 1).
Each consumer `i > basal`, in order:
1. picks a random initial prey `j < i`;
2. draws `x ~ Beta(1, β)`, with `β = (S² − basal²)/(2L) − 1` and `L = round(C·S²)`, and adds
   `κ = min(round(x·(i−1)), available)` more prey;
3. samples these prey without replacement with weights `exp(−|h_j − h_k| / T)`;
4. sets its level to `1 + mean(prey levels)`.

Small `T` favours prey close in trophic level to the first prey, so the web is more coherent.
Internal; use `generate_food_web(…; model=:ppm)`.
"""
function ppm(S::Int, C::AbstractFloat, basal::Int, rng=Random.default_rng(); T=0.25)
    # variables definition
    adj_matrix   = zeros(Int64, S, S)
    tl           = zeros(S)
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

    return adj_matrix
end