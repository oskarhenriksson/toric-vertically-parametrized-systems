using Oscar

function binomiality_check(F)
    G, T = groebner_basis_with_transformation_matrix(ideal(F), complete_reduction=true)
    if all(is_binomial, G)
        println("Gröbner basis: Generically binomial!")
    else
        println("Gröbner basis: Not binomial!")
        return (generically = false, for_all_positive = false)
    end

    # Check that the leading terms do not vanish
    LC = map(g->leading_coefficient(g), G)
    for c in LC   
        signs_of_numerator = sign.(coefficients(numerator(c)))
        if !all(s -> s .== 1 || s .== -1, signs_of_numerator)
            println("Specialization for all rate constants: Inconclusive")
            println("Nonvansishing leading coefficient")
            return (generically = true, for_all_positive = false)
        end
    end

    # Check that the all denominors in the transformation matrix are nonzero
    for p in T
        for c in coefficients(p)
            if !is_polynomial_of_constant_signs(denominator(c))
                println("Specialization for all rate constants: Inconclusive")
                println("Potential non-vanishing denominator in transformation matrix")
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
                    println("Specialization for all rate constants: Inconclusive")
                    println("Potential non-vanishing denominator: inconclusive!")
                    return (generically = true, for_all_positive = false)
                end
            end
        end
    end

    println("Specialization for all rate constants: Verified")
    return (generically = true, for_all_positive = true)
end