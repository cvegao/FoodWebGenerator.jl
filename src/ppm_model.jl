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