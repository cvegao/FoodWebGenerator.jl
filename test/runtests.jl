using FoodWebGenerator

using Test

@testset "FoodWebGenerator Tests" begin
    @testset "generate_food_web" begin
        @testset "Method 1: generate_food_web(S, C, basal, seed)" begin
            # Test 1: Generate a food web with a fixed seed and check properties
            S = 10
            C = 0.15
            basal = 2
            seed = 42
            fw = generate_food_web(S, C, basal, seed)
            
            @test size(fw.adj_matrix) == (S, S)
            @test length(fw.trophic_levels) == S
            @test 0 < fw.connectance <= 1.0
            @test all(fw.trophic_levels .>= 0.0)
        end

        @testset "Method 2: generate_food_web(adj_matrix)" begin
            # Test 2: Generate a food web with a predefined adjacency matrix
            adj_matrix = [0 0; 1 0]
            fw2 = generate_food_web(adj_matrix)
            
            @test size(fw2.adj_matrix) == (2, 2)
            @test length(fw2.trophic_levels) == 2
            @test fw2.connectance == sum(adj_matrix) / (2^2)

            @testset "self-loops" begin
                adj_matrix_with_self_loop = [1 0; 0 0]
                @test_throws ArgumentError generate_food_web(adj_matrix_with_self_loop)
            end
        end

        
        @testset "Method 3: generate_food_web(adj_matrix, trophic_levels)" begin
            # Test 3: Generate a food web with predefined adjacency matrix and trophic levels
            adj_matrix = [0 0; 1 0]
            tl         = [0.0, 1.0]
            fw3        = generate_food_web(adj_matrix, tl)
            
            @test size(fw3.adj_matrix) == (2, 2)
            @test length(fw3.trophic_levels) == 2
            @test fw3.connectance == sum(adj_matrix) / (2^2)
        end
    end

    @testset "trophic_coherence" begin
        # Test 4: Calculate trophic coherence for a simple food web
        adj_matrix = [0 1; 0 0]
        tl = [1.0, 0.0]
        q = trophic_coherence(adj_matrix, tl)
        
        @test q >= 0.0
    end

    @testset "trophic_levels" begin
        # Test 5: Calculate trophic levels for a simple food web
        adj_matrix = [0 1; 0 0]
        tl = trophic_levels(adj_matrix)
        
        @test length(tl) == 2
        @test tl[1] == 2.0
        @test tl[2] == 1.0
    end
end