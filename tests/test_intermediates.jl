using Test

function run_tests()
    @testset "Networks with intermediate species" begin

        @testset "Small network with input-1 and input-2 intermediates" begin
            N_small = matrix(QQ, [
                [-1, 1, 0, 0, 0],
                [-1, 1, 0, -1, 1],
                [1, -1, -1, 0, 0],
                [0, 0, 1, 1, -1]
            ])
            M_small = matrix(ZZ, [
                [1, 0, 0, 0, 0],
                [1, 0, 0, 1, 0],
                [0, 1, 1, 0, 0],
                [0, 0, 0, 0, 1]
            ])
            result_small = intermediates(N_small, M_small)
            expected_small = [
                (species=3, input_complexes=[[1,1,0,0]], output_complexes=[[0,1,0,0],[1,1,0,0]]),
                (species=4, input_complexes=[[1,1,0,0],[0,1,0,0]], output_complexes=[[0,1,0,0]])
            ]
            @test result_small == expected_small

            result_small_single_input = single_input_intermediates(N_small, M_small)
            expected_small_single_input = [(species=3, input_complexes=[[1,1,0,0]], output_complexes=[[0,0,0,1],[1,1,0,0]])]
            @test result_small_single_input == expected_small_single_input
        end

        @testset "Cyclic example" begin
            N_cyclic = matrix(QQ, [
                [-1, 0, 1],
                [1, -1, 0],
                [0, 1, -1]
            ])
            M_cyclic = matrix(ZZ, [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]
            ])
            result_cyclic = intermediates(N_cyclic, M_cyclic)
            expected_cyclic = [
                (species=1, input_complexes=[[0,1,0]], output_complexes=[[0,1,0]]),
                (species=3, input_complexes=[[0,1,0]], output_complexes=[[0,1,0]])
            ]
            @test result_cyclic == expected_cyclic

            # Toric invariance space
            A_cyclic = toric_invariance_space(N_cyclic, M_cyclic)
            A_cyclic_expected = matrix(ZZ, [[1, 1, 1]])
            @test rref(A_cyclic)[2] == A_cyclic_expected

            N_cyclic_reduced, M_cyclic_reduced = reduced_network(N_cyclic, M_cyclic, result_cyclic)
            A_cyclic_tilde = toric_invariance_space(N_cyclic_reduced, M_cyclic_reduced)
            A_cyclic_lifted = lift_exponent_matrix(A_cyclic_tilde, M_cyclic, result_cyclic)
            @test row_space(A_cyclic_lifted) == A_cyclic_expected

        end

        @testset "Reversible pairs" begin
            N_reversible_pairs = matrix(QQ, [
                [-1, 1, 0, 0],
                [1, -1, 0, 0],
                [0, 0, -1, 1],
                [0, 0, 1, -1]
            ])
            M_reversible_pairs = matrix(ZZ, [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ])
            result_reversible_pairs = intermediates(N_reversible_pairs, M_reversible_pairs)
            expected_reversible_pairs = [
                (species = 1, input_complexes = [[0, 1, 0, 0]], output_complexes = [[0, 1, 0, 0]]),
                (species = 3, input_complexes = [[0, 0, 0, 1]], output_complexes = [[0, 0, 0, 1]])
            ]
            @test result_reversible_pairs == expected_reversible_pairs

            # Toric invariance space
            A_reversible_pairs = toric_invariance_space(N_reversible_pairs, M_reversible_pairs)
            A_reversible_pairs_expected = matrix(ZZ, [0 0 1 1; 1 1 0 0])
            @test row_space(A_reversible_pairs) == row_space(A_reversible_pairs_expected)

            # Toric invariance space via reduced network
            N_reversible_pairs_reduced, M_reversible_pairs_reduced = reduced_network(N_reversible_pairs, M_reversible_pairs, result_reversible_pairs)
            A_reversible_pairs_tilde = toric_invariance_space(N_reversible_pairs_reduced, M_reversible_pairs_reduced)
            A_reversible_pairs_lifted = lift_exponent_matrix(A_reversible_pairs_tilde, M_reversible_pairs, result_reversible_pairs)
            @test row_space(A_reversible_pairs_lifted) == row_space(A_reversible_pairs_expected)



        end

        @testset "Shinar-Feinberg example" begin
            N_sf = matrix(QQ, [[-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1], [1, -1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, -1, -1, 0, 0, 0, -1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1]])
            M_sf = matrix(ZZ, [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])

            # Toric invariance space for full network
            A_sf = toric_invariance_space(N_sf, M_sf)
            A_sf_expected = matrix(ZZ, [1 1 1 0 1 1 0 1 1; 0 0 0 1 -1 0 0 0 0])
            @test row_space(A_sf) == row_space(A_sf_expected)

            # Find the intermediates
            result_sf = intermediates(N_sf, M_sf)
            expected_sf = [
                (species = 6, input_complexes = [[0, 0, 0, 1, 1, 0, 0, 0, 0]], output_complexes = [[0, 1, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0, 0]]),
                (species = 8, input_complexes = [[0, 0, 1, 0, 0, 0, 1, 0, 0]], output_complexes = [[0, 0, 1, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 1, 0, 0]]),
                (species = 9, input_complexes = [[1, 0, 0, 0, 0, 0, 1, 0, 0]], output_complexes = [[1, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 1, 0, 0]])
            ]
            @test result_sf == expected_sf

            # Find toric invariance space for reduced network and lift to original network
            N_sf_reduced, M_sf_reduced = reduced_network(N_sf, M_sf, result_sf)
            A_sf_tilde = toric_invariance_space(N_sf_reduced, M_sf_reduced)
            A_sf_lifted = lift_exponent_matrix(A_sf_tilde, M_sf, result_sf)
            @test row_space(A_sf_lifted) == row_space(A_sf_expected)
        end

    end
end

run_tests();

# X1 + X2 -> X3
# X3 -> X4
# X4 -> X5
# X5 -> X1 +X2
# X1 + X2 -> X3
# X3 -> X6

# X1 + X2 -> X3 -> X4 -> X5 -> X6 -> X3 -> X7

# c -> Y2 -> Y1 -> Y2 -> c'
# c -> Y2 -> c'
# c -> c'

# c -> Y1 -> c'
# c'' -> Y1 -> c'''

