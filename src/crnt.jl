
# General CRNT functions

function product_matrix(N::QQMatrix, M::ZZMatrix)
    @req size(N) == size(M) "Stoichiometric and kinetic matrix must have the same size"
    @req all(is_one, denominator.(N)) "Stoichiometric matrix needs to have integer entries"
    N = ZZ.(N)
    return M + N
end

function reaction_pairs(N::QQMatrix, B::ZZMatrix)
    P = product_matrix(N, B)
    reactions = [(reactant=B[:, j], product=P[:, j]) for j in 1:size(B, 2)]
    return reactions
end

function stoichiometric_and_kinetic_matrix_from_reaction_pairs(reactions::Vector{NamedTuple{(:reactant, :product),Tuple{Vector{T},Vector{T}}}}) where {T}
    @assert all(length(reactions[1].reactant) == length(rxn.reactant) && length(reactions[1].reactant) == length(rxn.product) for rxn in reactions) "All vectors in all pairs must have the same length"
    @assert T == ZZRingElem || T == Int "Entries must be of type ZZRingElem or Int"
    M = matrix(ZZ, hcat([rxn.reactant for rxn in reactions]...))
    N = matrix(QQ, hcat([rxn.product - rxn.reactant for rxn in reactions]...))
    return N, M
end


function minimal_siphons(N::QQMatrix, M::ZZMatrix)
    n = nrows(M)
    reactions = reaction_pairs(N, M)

    # Helper function to check if Z is a siphon
    function is_siphon(Z)
        if isempty(Z)
            return false  # Exclude the empty set
        end
        for rxn in reactions
            for i in Z
                if rxn.product[i] > 0
                    if !any(j -> rxn.reactant[j] > 0, Z)
                        return false
                    end
                end
            end
        end
        return true
    end

    # Helper function for growing siphons
    function extend_to_siphon(candidate, remaining_indices, list_of_siphons)
        # Add current candidate to the list if it is a siphon and does not contain a smaller siphon
        if is_siphon(candidate)
            if !any(s -> issubset(s, candidate) && s != candidate, list_of_siphons)
                push!(list_of_siphons, candidate)
            end
        end
        for i in remaining_indices
            new_candidate = union(candidate, [i])
            new_remaining = filter(x -> x > i, remaining_indices)
            extend_to_siphon(new_candidate, new_remaining, list_of_siphons)
        end
        return list_of_siphons 
    end

    return extend_to_siphon(Set{Int}(), 1:n, Set{Int}[])
end

function siphon_test(N::QQMatrix, M::ZZMatrix)
    siphons = minimal_siphons(N, M)
    W = kernel(N, side=:left)
    for Z in siphons
        if !(positive_vector_in_rowspace(W[:, collect(Z)]))
            return false
        end
    end
    return true
end

function reconstruct_stoichiometric_matrix_with_row_space_and_conserved_quantities(C::QQMatrix, W::QQMatrix)
    @req nrows(W) == rank(W) "Matrix of conserved quantities needs to have full row rank"
    @req nrows(C) == rank(C) "Coefficient matrix needs to have full row rank"
    @req nrows(C) + nrows(W) == ncols(W) "Not enough conserved quantities"
    W = rref(W)[2]
    C = rref(C)[2]
    pivot_colums = [findfirst(row .!= 0) for row in eachrow(W)]
    non_pivot_columns = setdiff(1:ncols(W), pivot_colums)
    N = zero_matrix(QQ,ncols(W),ncols(C))
    N[non_pivot_columns,:] = C 
    N[pivot_colums,:] = -W[:,non_pivot_columns]*C
    # Todo: Remove these checks
    @req row_space(kernel(N, side=:left)) == row_space(W) "Reconstruction failed"
    @req row_space(N) == row_space(C) "Reconstruction failed"
    return N
end
