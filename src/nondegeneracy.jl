
# Generic nondegeneracy checks

# Checks whether (C.diag(k).x^M, L*x-b) has a nondegenerate zero
function has_nondegenerate_zero(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M));
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)
    C = (rref(C)[2])[1:rank(C), :]
    L = rref(L)[2]
    @req ncols(C) == ncols(M) "C and M need to have the same number of columns"
    @req ncols(L) == nrows(M) "L needs to have the same number of columns as M has rows"
    s = rank(C)
    G = kernel(C, side=:right)
    for _ in 1:number_of_attempts
        u = rand(-max_entry_size:max_entry_size, ncols(G))
        nrows(L) == 0 ? h = ones(Int, nrows(M)) : h = rand(-max_entry_size:max_entry_size, nrows(M))
        degeneracy_matrix = vcat(C * diagonal_matrix(G * u) * transpose(M) * diagonal_matrix(h), L)
        if rank(degeneracy_matrix) == s + nrows(L)
            return true
        end
    end
    if certify
        R, u, h = polynomial_ring(QQ, "u" => 1:ncols(G), "h" => 1:nrows(M))
        nrows(L) == 0 ? h = ones(Int, nrows(M)) : nothing
        symbolic_degeneracy_matrix = vcat(C * diagonal_matrix(G * u) * transpose(M) * diagonal_matrix(h), R.(L))
        return !all(is_zero, minors_iterator(symbolic_degeneracy_matrix, s + nrows(L)))
    else
        return false
    end
end

function generic_local_acr(C::QQMatrix, M::ZZMatrix, i::Int; 
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)
    @req has_nondegenerate_zero(C, M) "The system needs to have a nondegenerate zero"
    return !has_nondegenerate_zero(C, M[setdiff(1:nrows(M), i), :], 
        number_of_attempts=number_of_attempts, max_entry_size=max_entry_size, certify=certify)
end


function all_positive_roots_nondegenerate(C::QQMatrix, M::ZZMatrix, L::Union{Nothing,QQMatrix}=nothing; extreme_rays=nothing)
    C = row_space(C)
    if isnothing(L)
        L = zero_matrix(QQ, 0, nrows(M))
    end
    
    # Start with a simplified check
    G = kernel(C, side=:right)
    R, alpha, h = polynomial_ring(QQ, "alpha" => 1:ncols(G), "h" => 1:nrows(M))
    simplified_nondeg_matrix = vcat(C * diagonal_matrix(R, R.(G) * alpha) * transpose(M) * diagonal_matrix(R, h), R.(L))
    if has_minor_of_constant_sign(simplified_nondeg_matrix, rank(C) + nrows(L))
        return true
    end

    # The full check (requires extreme rays)
    if isnothing(extreme_rays)
        extreme_rays = rays_of_nonnegative_kernel(C)
    end
    R, lambda, h = polynomial_ring(QQ, "lambda" => 1:ncols(extreme_rays), "h" => 1:nrows(M))
    nondeg_matrix = vcat(C * diagonal_matrix(extreme_rays * lambda) * transpose(M) * diagonal_matrix(h), R.(L))
    return has_minor_of_constant_sign(nondeg_matrix, rank(C) + nrows(L))
end

all_positive_roots_nondegenerate(C::QQMatrix, M::ZZMatrix, L::ZZMatrix; extreme_rays=nothing) = 
    all_positive_roots_nondegenerate(C, M, QQ.(L), extreme_rays=extreme_rays)
