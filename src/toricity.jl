using Oscar
using DataStructures

# Fundamental partition

# Transitive closure of a relation
function transitive_closure(relation)
    all_elements = unique(vcat(collect.(relation)...))

    # Initialize DisjointSets structure with elements
    uf = DisjointSets(all_elements)

    # Use the first element of each block as a representative to union with the others
    for block in relation
        representative = first(block)
        for element in block
            union!(uf, representative, element)
        end
    end

    # Create closure by grouping elements by their "root representative"
    closure = Dict{Int,Vector{Int}}()
    for element in all_elements
        representative = find_root(uf, element)
        if !haskey(closure, representative)
            closure[representative] = []
        end
        push!(closure[representative], element)
    end

    return collect(values(closure))
end

function fundamental_partition(C::QQMatrix)

    if ncols(C) == 0
        return Vector{Int}[]
    end

    # Basis of elementary vectors of ker(C)
    G = kernel(C, side=:right)
    E = transpose(rref(transpose(G))[2])

    # Supports of the elementary vectors
    supports = supp.(Vector.(eachcol(E)))

    @req Set(union(supports...)) == Set(1:ncols(C)) "The supports do not cover all columns"

    # Compute the transitive closure
    return transitive_closure(supports)
end


function toric_invariance_space(C::QQMatrix, M::ZZMatrix)
    @req ncols(C) == ncols(M) "C and M need to have the same number of columns"
    if ncols(C) == 0
        return identity_matrix(ZZ, nrows(M))
    end
    FP = fundamental_partition(C)
    directions = Vector{ZZRingElem}[]
    for block in FP
        for i in 1:(length(block)-1)
            push!(directions, M[:, block[i]] - M[:, block[i+1]])
        end
    end
    return kernel(matrix(ZZ, hcat(directions...)), side=:left)
end


# Injectivity wrt ker(A)
function injectivity_matrix(C::QQMatrix, M::ZZMatrix, A::QQMatrix)
    R, mu, alpha = polynomial_ring(QQ, "mu" => 1:ncols(M), "alpha" => 1:nrows(M))
    return vcat(C * diagonal_matrix(R, mu) * transpose(M) * diagonal_matrix(R, alpha), R.(A))
end

function injectivity_test(C::QQMatrix, M::ZZMatrix, A::QQMatrix)
    C = row_space(C)
    @req nrows(C) + nrows(A) == nrows(M) "The injectivity matrix needs to be square"
    inj_matrix = injectivity_matrix(C, M, A)
    inj_polynomial = det(inj_matrix)
    signs = [sign(c) for c in Oscar.coefficients(inj_polynomial)]
    return all(signs .== 1) || all(signs .== -1)
end

injectivity_test(C::QQMatrix, M::ZZMatrix, A::ZZMatrix) = injectivity_test(C, M, QQ.(A))


import HomotopyContinuation as HC

function coset_counting_system(C::QQMatrix, M::ZZMatrix, A::ZZMatrix)
    C = row_space(C)
    C = Matrix(Rational{Int}.(C))
    M = Matrix(Int.(M))
    A = Matrix(Int.(A))
    HC.@var x[1:nrows(M)] k[1:ncols(M)] c[1:nrows(A)]
    steady_state_part = C * (k .* [prod(x .^ m) for m in eachcol(M)])
    linear_part = A * x - c
    return HC.System(vcat(steady_state_part, linear_part), variables=x, parameters=vcat(k, c))
end

function has_nonpositive_entry(C::HC.AbstractSolutionCertificate)
    if !HC.is_certified(C)
        return false
    elseif HC.is_complex(C)
        return true
    else
        return any(i -> HC.Arblib.is_negative(real(HC.Arblib.ref(C.I, i, 1))), 1:length(C.I))
    end
end

function number_of_cosets_for_random_parameters(C::QQMatrix, M::ZZMatrix, A::ZZMatrix;
    certify = true, root_bound::Union{Int,Nothing}=nothing)
    F = coset_counting_system(C, M, A)
    k_star = rand(1:100, ncols(M))
    c_star = Matrix(Int.(A)) * rand(1:100, nrows(M))
    res = HC.solve(F, only_non_zero=true, target_parameters=vcat(k_star, c_star), show_progress=false)
    cert_res = HC.certify(F, res, target_parameters=vcat(k_star, c_star))
    certs = HC.certificates(cert_res)
    if certify
        if isnothing(root_bound)
            root_bound = HC.mixed_volume(F)
        end
        if HC.nsolutions(res) < root_bound
            error("Root bound not attained")
        end
        if HC.ndistinct_certified(cert_res) < root_bound
            error("Failed to certify all solutions")
        end
        number_of_positive = 0
        for c in certs
            if HC.is_positive(c)
                number_of_positive += 1
            elseif !has_nonpositive_entry(c)
                error("Failed to decide on positivity for $(c.I)")
            end
        end
    else
        number_of_positive = count(HC.is_positive, certs)
    end
    return number_of_positive
end


function coset_counting_analysis(N::QQMatrix, M::ZZMatrix, A::ZZMatrix; printing_function=println)
    toricity_flag = false
    finite_flag = false
    @req nrows(N) == nrows(M) "The stoichiometric and kinetic matrix must have the same number of rows"
    C = row_space(N)
    @req rank(C) + rank(A) == nrows(M) "Too small toric invariance space"

    if injectivity_test(C, M, A)
        printing_function("Injectivity test passed for reduced network")
        printing_function("Toricity!")
        toricity_flag = true
        finite_flag = true
        return (toricity=toricity_flag, finite=finite_flag)
    end

    if all_positive_roots_nondegenerate(C, M)
        printing_function("Finitely many cosets")
        finite_flag = true
    end

    if finite_flag
        F = coset_counting_system(C, M, A)
        mv = HC.mixed_volume(F)

        printing_function("Mixed volume bound on the number of cosets: $(mv)")
        grc = HC.nsolutions(HC.solve(F, only_non_zero=true, target_parameters=rand(ComplexF64, length(F.parameters)), show_progress=false))
        printing_function("Estimation of generic root count: $(grc)")

        if mv == 1
            printing_function("Toricity!")
            toricity_flag = true
            return (toricity=toricity_flag, finite=finite_flag)
        end

        # Constant number of cosets?
        N_A = reconstruct_stoichiometric_matrix_with_row_space_and_conserved_quantities(C, QQ.(A))
        if positive_vector_in_rowspace(A) && siphon_test(N_A, M)
            printing_function("Constant number of cosets")
            coset_count = nothing
            try
                coset_count = number_of_cosets_for_random_parameters(C, M, A, certify=true)
                printing_function("Number of cosets for random parameters: $(coset_count)")
                if coset_count == 1
                    printing_function("Toricity!")
                    toricity_flag = true
                end
            catch e
                println(e)
            end
        end
    end
    return (toricity=toricity_flag, finite=finite_flag)
end


function coset_with_multistationarity(A::ZZMatrix, W::QQMatrix)
    n = ncols(A)
    B = kernel(A, side=:right)
    R, alpha = polynomial_ring(QQ, "alpha" => 1:nrows(B))
    multistat_matrix = vcat(R.(W), transpose(B) * diagonal_matrix(alpha))

    if ncols(B) + nrows(W) < n
        return true
    end

    if has_minor_of_constant_sign(multistat_matrix, n)
        return false
    end

    if ncols(B) + nrows(W) == n
        return true
    end

    if ncols(B) + nrows(W) > n
        @info "Inconclusive: check if the minors vanish simultaneously"
        return false
    end

end
