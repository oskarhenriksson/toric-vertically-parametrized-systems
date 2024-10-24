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
                (species=3, input_reactions=[1], output_reactions=[2, 3]),
                (species=4, input_reactions=[3, 4], output_reactions=[5])
            ]
            @test result_small == expected_small

            result_small_single_input = intermediates(N_small, M_small, only_single_input=true)
            expected_small_single_input = [(species=3, input_reactions=[1], output_reactions=[2, 3])]
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
                (species=1, input_reactions=[3], output_reactions=[1]),
                (species=2, input_reactions=[1], output_reactions=[2]),
                (species=3, input_reactions=[2], output_reactions=[3]),
            ]
            @test result_cyclic == expected_cyclic

            # Toric invariance space
            A_cyclic = toric_invariance_space(N_cyclic, M_cyclic)
            A_cyclic_expected = matrix(ZZ, [[1, 1, 1]])
            @test rref(A_cyclic)[2] == A_cyclic_expected

            intermediate_species_cyclic = [i.species for i in result_cyclic]
            N_cyclic_reduced, M_cyclic_reduced = reduced_network(N_cyclic, M_cyclic, intermediate_species_cyclic)
            A_cyclic_tilde = toric_invariance_space(N_cyclic_reduced, M_cyclic_reduced)
            A_cyclic_lifted = lift_exponent_matrix(A_cyclic_tilde, M_cyclic, result_cyclic)
            @test rref(A_cyclic_lifted)[2] == A_cyclic_expected

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
                (species=1, input_reactions=[2], output_reactions=[1]),
                (species=2, input_reactions=[1], output_reactions=[2]),
                (species=3, input_reactions=[4], output_reactions=[3]),
                (species=4, input_reactions=[3], output_reactions=[4])
            ]
            @test result_reversible_pairs == expected_reversible_pairs

            # Toric invariance space
            A_reversible_pairs = toric_invariance_space(N_reversible_pairs, M_reversible_pairs)
            A_reversible_pairs_expected = matrix(ZZ, [0 0 1 1; 1 1 0 0])
            @test rref(A_reversible_pairs)[2] == rref(A_reversible_pairs_expected)[2]

            # Toric invariance space via reduced network
            intermediate_species_reversible_pairs = [i.species for i in result_reversible_pairs]
            N_reversible_pairs_reduced, M_reversible_pairs_reduced = reduced_network(N_reversible_pairs, M_reversible_pairs, intermediate_species_reversible_pairs)
            A_reversible_pairs_tilde = toric_invariance_space(N_reversible_pairs_reduced, M_reversible_pairs_reduced)
            A_reversible_pairs_lifted = lift_exponent_matrix(A_reversible_pairs_tilde, M_reversible_pairs, result_reversible_pairs)
            @test rref(A_reversible_pairs_lifted)[2] == rref(A_reversible_pairs_expected)[2]



        end

        @testset "Shinar-Feinberg example" begin
            N_sf = matrix(QQ, [[-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1], [1, -1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, -1, -1, 0, 0, 0, -1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1]])
            M_sf = matrix(ZZ, [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])

            # Toric invariance space for full network
            A_sf = toric_invariance_space(N_sf, M_sf)
            A_sf_expected = matrix(ZZ, [1 1 1 0 1 1 0 1 1; 0 0 0 1 -1 0 0 0 0])
            @test rref(A_sf)[2] == A_sf_expected

            # Find the intermediates
            result_sf = intermediates(N_sf, M_sf)
            expected_sf = [
                (species=6, input_reactions=[6], output_reactions=[7, 8]),
                (species=8, input_reactions=[9], output_reactions=[10, 11]),
                (species=9, input_reactions=[12], output_reactions=[13, 14])
            ]
            @test result_sf == expected_sf

            # Find toric invariance space for reduced network and lift to original network
            intermediate_species_sf = [i.species for i in result_sf]
            N_sf_reduced, M_sf_reduced = reduced_network(N_sf, M_sf, intermediate_species_sf)
            A_sf_tilde = toric_invariance_space(N_sf_reduced, M_sf_reduced)
            A_sf_lifted = lift_exponent_matrix(A_sf_tilde, M_sf, result_sf)
            @test rref(A_sf_lifted)[2] == A_sf_expected
        end

    end
end

run_tests();