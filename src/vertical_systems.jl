

# General vertical system functions
vertical_system = function (C::QQMatrix, M::ZZMatrix; a=nothing)
    if isnothing(a)
        A, a = rational_function_field(QQ, "a" => 1:ncols(M))
        R, x = polynomial_ring(A, "x" => 1:nrows(M))
    else
        @req length(a) == ncols(M) "The number of parameters must match the number of columns in M"
        R, x = polynomial_ring(QQ, "x" => 1:nrows(M))
    end
    return C * diagonal_matrix(a) * [prod(x .^ m) for m in eachcol(M)]
end

# Checks the matrix of monomials for linearity
function is_linear(M::ZZMatrix)
    for col in eachcol(M)
        if !is_unit_vector(Vector(col)) && !is_zero(Vector(col))
            return false
        end
    end
    return true
end