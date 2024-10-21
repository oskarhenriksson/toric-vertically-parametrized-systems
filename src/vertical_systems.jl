

# General vertical system functions
vertical_system = function (C::QQMatrix, M::ZZMatrix)
    A, a = polynomial_ring(QQ, "a" => 1:ncols(M))
    R, x = polynomial_ring(A, "x" => 1:nrows(M))
    return C * diagonal_matrix(a) * [prod(x .^ m) for m in eachcol(M)]
end

# Checks the matrix of monomials for linearity
function is_linear(B::ZZMatrix)
    for i in 1:ncols(B)
        col = B[:,i]
        if !is_unit_vector(col) && !is_zero(col)
            return false
        end
    end
    return true
end