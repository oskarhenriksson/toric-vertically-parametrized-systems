using Test

function run_tests()
    @testset "Graph theoretic computations" begin
        @testset "Small graph" begin
            N = matrix(QQ, [-1 1 0 0 0; 0 0 0 0 0; 0 0 -1 0 1; 0 0 1 -1 0; 0 0 0 1 -1])
            M = matrix(ZZ, [1 0 0 0 0; 1 2 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1])
            g, _ = reaction_graph(N, M)
            @test length(Graphs.vertices(g)) == 7
            @test length(Graphs.weakly_connected_components(g)) == 3
            @test delta(N, M) == 1
        end

    end
end

run_tests()