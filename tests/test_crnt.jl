using Test

function run_tests()
    @testset "Computing CRNT structure" begin
        @testset "Incomplete stoichiometric matrix" begin
            N = matrix(QQ, [-1 1 0 0 0; 0 0 0 0 0; 0 0 -1 0 1; 0 0 1 -1 0])
            M = matrix(ZZ, [1 0 0 0 0; 1 2 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1])
            @test_throws ArgumentError product_matrix(N, M)
            @test_throws ArgumentError reaction_pairs(N, M)
        end

        @testset "Small network" begin
            N = matrix(QQ, [-1 1 0 0 0; 0 0 0 0 0; 0 0 -1 0 1; 0 0 1 -1 0; 0 0 0 1 -1])
            M = matrix(ZZ, [1 0 0 0 0; 1 2 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1])
            P = product_matrix(N, M)
            reactions = reaction_pairs(N, M)
            @test length(reactions) == 5
            @test reactions[1] == (reactant=M[:, 1], product=P[:, 1])
        end

        @testset "Siphons" begin
            N_siph =matrix(QQ, [
                [0, 0, 1, -1, 1, -1, 2, -2],  # X1
                [-2, 2, 1, -1, -1, 1, 0, 0],  # X2
                [1, -1, -2, 2, -1, 1, -3, 3], # X3
                [1, -1, 0, 0, 1, -1, 1, -1]   # X4
            ])


            M_siph = matrix(ZZ, [
                [0, 0, 0, 1, 0, 1, 0, 2],   # X1
                [2, 0, 0, 1, 1, 0, 0, 0],   # X2
                [0, 1, 2, 0, 1, 0, 3, 0],   # X3
                [0, 1, 0, 0, 0, 1, 0, 1]    # X4
            ])

            @test Set(minimal_siphons(N_siph, M_siph)) == Set([Set([2, 3, 1]), Set([4, 2, 3])])
            @test siphon_test(N_siph, M_siph) == true
        end

    end
end

run_tests();