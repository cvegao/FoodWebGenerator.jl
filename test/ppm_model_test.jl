include("../src/ppm_model.jl")

@testset "Preferential Preying Model Tests" begin
    @testset "Connectance matches target C" begin

        rng = MersenneTwister(1234)

        # (S, C) combinations to test across a range of richness and connectance targets
        combos = [
            (30, 0.05),
            (30, 0.15),
            (50, 0.10),
            (50, 0.25),
            (100, 0.05),
            (100, 0.20),
        ]

        B = 4  # basal species

        n_reps = 200   # number of replicate webs per (S, C) combo
        atol_C = 0.04  # absolute tolerance allowed between mean empirical C and target C

        for (S, C) in combos
            empirical_C = zeros(n_reps)

            for r in 1:n_reps
                adj = ppm(S, C, B, rng)
                empirical_C[r] = sum(adj) / S^2
            end

            mean_C = mean(empirical_C)

            @testset "S=$S, C=$C" begin
                @test isapprox(mean_C, C; atol=atol_C)
            end
        end
    end

    @testset "Structural properties" begin

        rng = MersenneTwister(42)
        S, C, B = 40, 0.12, 8
        n_reps = 50

        for r in 1:n_reps
            adj = ppm(S, C, B, rng)

            # No species is its own predator/prey (no self-loops)
            @test all(iszero, diag(adj))

            # Matrix has the correct shape
            @test size(adj) == (S, S)

            # Entries are binary (0/1)
            @test all(x -> x == 0 || x == 1, adj)

            # At least one basal species exists (a species with no prey, i.e. an all-zero column)
            n_prey = vec(sum(adj, dims=1))
            @test any(iszero, n_prey)

            # No fully isolated species (every species has at least one prey or one predator)
            n_predators = vec(sum(adj, dims=2))
            @test all(i -> n_prey[i] > 0 || n_predators[i] > 0, 1:S)
        end
    end

    @testset "Edge cases" begin
        rng = MersenneTwister(7)

        # S = 1: no possible links, should return an all-zero 1x1 matrix
        adj1 = ppm(1, 0.1, 1, rng)
        @test size(adj1) == (1, 1)
        @test adj1[1, 1] == 0

        # Very low connectance still returns a valid, fully-populated adjacency matrix
        adj_low = ppm(20, 0.01, 2, rng)
        @test size(adj_low) == (20, 20)
        @test all(iszero, diag(adj_low))
    end
end