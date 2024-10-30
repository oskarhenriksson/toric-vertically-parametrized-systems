using Oscar

function binomiality_check(F; printing_function=println, verbose::Bool=true)

    # Specialize F at random parameters and check if binomiality can be ruled out
    Qax = parent(first(F))
    Qa = base_ring(Qax)
    @req isa(Qa, AbstractAlgebra.Generic.RationalFunctionField{QQFieldElem, QQMPolyRingElem}) "Coefficient field must be a rational function field"
    Qx, x = polynomial_ring(QQ, symbols(Qax))
    a_star = rand(-100:100, ngens(Qa))
    phi = hom(Qax, Qx, c->evaluate(c,a_star), x)
    F_specialized = phi.(F)
    G_specialized = groebner_basis(ideal(F_specialized), complete_reduction=true)
    if !all(is_binomial, G_specialized)
        verbose && printing_function("Gröbner basis: Not binomial!")
        return (generically = false, for_all_positive = false)
    end

    # Compute a Gröbner basis over the rational function field
    G, T = groebner_basis_with_transformation_matrix(ideal(F), complete_reduction=true)
    if all(is_binomial, G)
        verbose && printing_function("Gröbner basis: Generically binomial!")
    else
        verbose && printing_function("Gröbner basis: Not binomial!")
        return (generically = false, for_all_positive = false)
    end

    # Check that the leading terms do not vanish
    LC = map(g->leading_coefficient(g), G)
    for c in LC   
        signs_of_numerator = sign.(coefficients(numerator(c)))
        if !all(s -> s .== 1 || s .== -1, signs_of_numerator)
            verbose && printing_function("Specialization for all rate constants: Inconclusive")
            verbose && printing_function("Nonvansishing leading coefficient")
            return (generically = true, for_all_positive = false)
        end
    end

    # Check that the all denominors in the transformation matrix are nonzero
    for p in T
        for c in coefficients(p)
            if !is_polynomial_of_constant_signs(denominator(c))
                verbose && printing_function("Specialization for all rate constants: Inconclusive")
                verbose && printing_function("Potential non-vanishing denominator in transformation matrix")
                return (generically = true, for_all_positive = false)
            end
        end
    end

    # Check that each original generator can be expressed in terms of the Gröbner basis
    # with nonvanishing denominators
    for f in F
        Q,r = reduce_with_quotients(f, elements(G))
        for p in Q
            for c in coefficients(p)
                if !is_polynomial_of_constant_signs(denominator(c))
                    verbose && printing_function("Specialization for all rate constants: Inconclusive")
                    verbose && printing_function("Potential non-vanishing denominator: inconclusive!")
                    return (generically = true, for_all_positive = false)
                end
            end
        end
    end

    verbose && printing_function("Specialization for all rate constants: Verified")
    return (generically = true, for_all_positive = true)
end