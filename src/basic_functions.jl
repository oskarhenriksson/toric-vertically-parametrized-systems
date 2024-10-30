
# Basic functions

# Iterate over rows and columns of a matrix 
# (will be included in future Oscar versions) 
Base.eachrow(a::MatrixElem) = Slices(a, (1, :), (axes(a, 1),))
Base.eachcol(a::MatrixElem) = Slices(a, (:, 1), (axes(a, 2),))


function is_polynomial_of_constant_signs(f::MPolyRingElem)
    coeffs = coefficients(f)
    signs = sign.(coeffs)
    return  all(signs .== 1) || all(signs .== -1)
end

function has_minor_of_constant_sign(M::MatrixElem, r::Int)
    for minor in minors_iterator(M, r)
        if is_polynomial_of_constant_signs(minor)
            return true
        end
    end
    return false
end

# Check if element is nonzero
is_non_zero(x) = !is_zero(x)

# Support of a vector
supp(v::Vector) = findall(is_non_zero, v)

# Row space
function row_space(A::MatrixElem) 
    Arref = rref(A)[2][1:rank(A), :]
    if ncols(Arref) > 0 && nrows(Arref) > 0
        return Arref*sign(Arref[1,1])
    else
        return Arref
    end
end

# Zero columns of matrix
zero_columns(A::MatrixElem) = [i for i in 1:ncols(A) if all(is_zero, A[:, i])]

# Checks whether N has a positive vector in its kernel
# If N is the stoichiometric matrix of a network, this corresponds to checking if the network is consistent
function nonempty_positive_kernel(N::Union{QQMatrix,ZZMatrix})
    inequalities = (-identity_matrix(QQ, ncols(N)), -ones(Int, ncols(N)))
    equalities = (N, zeros(Int, nrows(N)))
    P = polyhedron(inequalities, equalities)
    return is_feasible(P)
end

# Positive vector in row space
function positive_vector_in_rowspace(A::Union{ZZMatrix,QQMatrix})
    A = QQ.(A)
    kerA = kernel(A, side=:right)
    return nonempty_positive_kernel(transpose(kerA))
end


function rays_of_nonnegative_kernel(C::QQMatrix)
    inequalities = (-identity_matrix(QQ, ncols(C)), zeros(Int, ncols(C)))
    equalities = (C, zeros(Int, nrows(C)))
    P = polyhedron(inequalities, equalities)
    return transpose(matrix(QQ, rays(P)))
end

# Function to generate all p-by-p minors of matrix M
# (Will be included in future versions of Oscar!)
function minors_iterator(M::MatrixElem, k::Int)
    row_indices = AbstractAlgebra.combinations(1:nrows(M), k)
    col_indices = AbstractAlgebra.combinations(1:ncols(M), k)
    return (det(M[rows, cols]) for rows in row_indices for cols in col_indices)
end

function is_unit_vector(v::Vector)
    sum(v) == 1 && all(x -> x == 0 || x == 1, v)
end

function unit_vector(n::Int, i::Int)
    @req 1 <= i <= n "Index i must be between 1 and n"
    return [j == i ? 1 : 0 for j in 1:n]
end
