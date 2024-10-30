using Test

function run_tests()
    @testset "Shinar-Feinberg" begin
        N = matrix(QQ, [[-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1], [1, -1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, -1, -1, 0, 0, 0, -1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1]])
        M = matrix(ZZ, [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])

        # Linearity
        @test is_linear(M) == false

        # Row reduced stoichiometric matrix
        C = rref(N)[2][1:rank(N), :]

        # Conserved quantities
        W = kernel(N, side=:left)

        # Graph structure and deficiency
        g, _ = reaction_graph(N, M)
        m = length(Graphs.vertices(g))
        @test m == 13
        ell = length(Graphs.weakly_connected_components(g))
        @test ell == 4
        @test delta(N, M) == 2

        # Consistency
        @test nonempty_positive_kernel(N) == true

        # Generic nondegeneracy
        @test has_nondegenerate_zero(N, M) == true

        # Reduced network (without intermediate species)
        intermediates_result = intermediates(N, M)
        Ntilde, Mtilde = reduced_network(N, M, intermediates_result)
        @test length(complexes(Ntilde, Mtilde)) == 10


        # Fundamental partition
        FP = fundamental_partition(N)
        @test Set(Set.(FP)) == Set([Set(collect(3:14)), Set([1, 2])])

        # Toric invariance space
        A = toric_invariance_space(C, M)
        @test rref(A)[2] == matrix(ZZ, [1 1 1 0 1 1 0 1 1; 0 0 0 1 -1 0 0 0 0])

        Atilde = toric_invariance_space(Ntilde, Mtilde)
        Atilde_lifted = lift_exponent_matrix(Atilde, M, intermediates_result)
        @test rref(Atilde_lifted)[2] == matrix(ZZ, [1 1 1 0 1 1 0 1 1; 0 0 0 1 -1 0 0 0 0])

        # Generically finitely many cosets?
        @test has_nondegenerate_zero(C, M, QQ.(A)) == true

        # Injectivity test
        @test injectivity_test(C, M, A) == false

        #### Injectivity test for reduced network
        @test injectivity_test(Ntilde, Mtilde, Atilde) == true

        # All steady states nondegenerate?
        @test all_positive_roots_nondegenerate(C, M) == true

        # Mixed volume bound on the number of cosets
        F = coset_counting_system(C, M, A)
        @test HC.mixed_volume(F) == 4

        # Siphon test
        @test siphon_test(N, M) == true

        # Capacity for multistationarity
        @test coset_with_multistationarity(A, W) == false

        # ACR
        @test zero_columns(A) == [7]

        # The deficiency one theorem
        @test convered_by_deficiency_one_theorem(N, M) == false

        # binomiality
        F = vertical_system(C, M)
        binomiality_result = binomiality_check(F, verbose=false)
        @test binomiality_result.generically == true 
        @test binomiality_result.for_all_positive == true

    end
end

run_tests();