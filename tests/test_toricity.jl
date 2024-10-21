using Test

function run_tests()
    @testset "Checking toricity" begin

        @testset "Small network" begin
            N = matrix(QQ, [-1 1 0 0 0; 0 0 0 0 0; 0 0 -1 0 1; 0 0 1 -1 0; 0 0 0 1 -1])
            M = matrix(ZZ, [1 0 0 0 0; 1 2 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1])
            A = toric_invariance_space(N, M)
            A_expected = matrix(ZZ, [0 0 1 1 1; 1 1 0 0 0])
            @test rref(A)[2] == rref(A_expected)[2]
        end

        @testset "No invariance" begin
            C_nontoric = matrix(QQ, [-1 1 1 0 0 0; -1 1 0 0 0 1; 0 0 0 1 -1 -1])
            M_nontoric = matrix(ZZ, [1 -1 -2 0 1 0; 2 -2 -1 2 -1 -2; 0 -2 -1 -2 0 -1])
            A = toric_invariance_space(C_nontoric, M_nontoric)
            @test A == zero_matrix(ZZ, 0, 3)
        end

        @testset "Reduced Shinar-Feinberg" begin
            # Reduced Shinar-Feinberg
            N_rsf = matrix(QQ, [[-1, 1, 0, 0, 0, 0, 0, 0], [1, -1, -1, 1, 0, 1, 0, 0], [0, 0, 1, -1, -1, 0, 0, 0], [0, 0, 0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 0, -1, 1, 1], [0, 0, 0, 0, 0, 1, -1, -1]])
            M_rsf = matrix(ZZ, [[1, 0, 0, 0, 0, 0, 0, 1], [0, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1]])
            A_rsf = toric_invariance_space(N_rsf, M_rsf)
            @test A_rsf == matrix(ZZ, [1 1 1 0 1 0; 0 0 0 1 -1 0])
            @test injectivity_test(N_rsf, M_rsf, A_rsf) == true
        end


    end
end

run_tests()