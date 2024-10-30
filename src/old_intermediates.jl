using Oscar

function intermediates(N::QQMatrix, M::ZZMatrix; only_single_input::Bool=false)

    @req size(N) == size(M) "Stoichiometric and kinetic matrix must have the same size"
    n = size(N, 1)  # Number of species (rows in N)
    r = size(N, 2)  # Number of reactions (columns in N)
    P = M + N       # Product matrix

    intermediates_data = Vector{NamedTuple{(:species, :input_reactions, :output_reactions),Tuple{Int,Vector{Int},Vector{Int}}}}()  # List to store named tuples of intermediates and their reactions

    # Find candidate species for intermediates:
    # Only appears in monomolecular complexes
    # 
    for i in 1:n
        in_rxns = Int[]   # Input reactions where species i is formed
        out_rxns = Int[]  # Output reactions where species i is consumed

        skip_species = false 

        # Loop over all reactions
        for j in 1:r

            # Skip reaction if species i is neither consumed nor produced
            if M[i, j] == 0 && P[i, j] == 0
                continue
            end

            # The species is not an intermediate if it appears with coefficient greater than 1
            if M[i, j] > 1 || P[i, j] > 1
                skip_species = true
                break
            end

            # List out-reactions and check that species i does not react with another species
            if M[i, j] == 1
                for row in 1:n
                    # Not an intermediate if it reacts with another species
                    if row != i && M[row, j] > 0
                        skip_species = true
                        break
                    end
                end
                if skip_species
                    break
                end
                push!(out_rxns, j)
            end

            # List in-reactions and check that species i is not produced with another species
            if P[i, j] == 1
                for row in 1:n
                    # Not an intermediate if it is produced with another species
                    if row != i && P[row, j] > 0
                        skip_species = true
                        break
                    end
                end
                if skip_species
                    break
                end
                push!(in_rxns, j)
            end
        end

        # If species should be skipped, move to the next species
        if skip_species
            continue
        end

        # Species i is an intermediate if it has exactly 1 input reaction and at least 1 output reaction
        if length(in_rxns) > 0 && length(out_rxns) > 0
            if only_single_input && length(in_rxns) > 1
                continue
            end
            push!(intermediates_data, (species=i, input_reactions=in_rxns, output_reactions=out_rxns))
        end
    end

    return intermediates_data
end

function single_input_intermediates(N::QQMatrix, M::ZZMatrix)
    return intermediates(N, M, only_single_input=true)
end

function lift_exponent_matrix(Atilde::ZZMatrix, B::ZZMatrix, intermediates_result::Vector{NamedTuple{(:species, :input_reactions, :output_reactions),Tuple{Int, Vector{Int},Vector{Int}}}})
    n = nrows(B)
    non_intermediates = setdiff(1:n, [I.species for I in intermediates_result])
    A = zero_matrix(ZZ, nrows(Atilde), n)
    A[:, non_intermediates] = Atilde
    filled_in = deepcopy(non_intermediates)

    # Keep working while not all species have been filled into A
    while Set(filled_in) != Set(1:n)
        i = setdiff(1:n, filled_in)[1]
        
        # Find the corresponding index for species i in the intermediates list
        idx = findfirst(I -> I.species == i, intermediates_result)

        # Initialize a chain starting with species i
        chain = [i]
        end_of_chain = false

        # Process the chain until it reaches a species whose dependencies are all filled in
        while !end_of_chain
            # Get the first input reaction for the current intermediate species
            next_complex = B[:, [I.input_reactions[1] for I in intermediates_result][idx]]
            
            # Check if all species involved in this next complex have been filled in
            if Set(supp(next_complex)) ⊆ filled_in
                end_of_chain = true
                
                # Once the chain has reached a fully filled complex, update A for all species in the chain
                for j in chain
                    A[:, j] = A[:, j] + A * next_complex
                end
                
                # Mark all species in the chain as filled in
                filled_in = Set(filled_in) ∪ Set(chain)
            
            # If the next complex involves exactly one species, extend the chain
            elseif length(supp(next_complex)) == 1
                i = supp(next_complex)[1]
                
                # Find the corresponding intermediate index for this species
                idx = findfirst(I -> I.species == i, intermediates_result)
                
                # Check if the species is already in the chain 
                if i in chain
                    # If a cycle is detected, create a new row in A for the chain and stop extending
                    new_row = zero_matrix(ZZ, 1, n)
                    new_row[1, chain] = ones(Int, length(chain)) 
                    A = vcat(A, new_row) 
                    filled_in = Set(filled_in) ∪ Set(chain) 
                    end_of_chain = true
                else
                    # If no cycle is detected, add the species to the chain and continue
                    push!(chain, i)
                end
            else
                # If the complex involves multiple species, raise an error
                error("Error!")
            end
        end
    end
    return A
end


# To do: Decide how to deal with the trivial species (zero-rows in both N and M)
# To do: Check that the intermediates truly are intermediates?
function reduced_network(N::QQMatrix, M::ZZMatrix, intermediate_species::Union{Nothing,Vector{Int}}=nothing)
    reactions = reaction_pairs(N, M)
    reduced_reactions = deepcopy(reactions)
    if isnothing(intermediate_species)
        intermediate_species = [i.species for i in intermediates(N, M)]
    end
    for i in intermediate_species
        in_reactions = [rxn for rxn in reduced_reactions if rxn.product == identity_matrix(ZZ, nrows(M))[:, i]]
        out_reactions = [rxn for rxn in reduced_reactions if rxn.reactant == identity_matrix(ZZ, nrows(M))[:, i]]
        for in_rxn in in_reactions
            for out_rxn in out_reactions
                # Combine only if it gives a nontrivial reaction
                if in_rxn.reactant != out_rxn.product
                    # New reaction bypassing the intermediate
                    push!(reduced_reactions, (reactant=in_rxn.reactant, product=out_rxn.product))
                end
            end
        end
        # Remove old reactions involving the intermediate
        reactions_to_remove = union(in_reactions, out_reactions)
        reduced_reactions = setdiff(reduced_reactions, reactions_to_remove)
    end

    # Deal with the special case of the reduced network being empty
    if length(reduced_reactions) == 0
        N_reduced = zero_matrix(QQ, nrows(M), 0)
        B_reduced = zero_matrix(ZZ, nrows(M), 0)
    else
        N_reduced, B_reduced = stoichiometric_and_kinetic_matrix_from_reaction_pairs(reduced_reactions)
    end

    # Remove the rows corresponding to the intermediate before returning
    new_species = [i for i = 1:nrows(M) if i ∉ intermediate_species]
    return N_reduced[new_species, :], B_reduced[new_species, :]
end


function reduced_network(N::QQMatrix, M::ZZMatrix, 
    intermediates_result::Vector{NamedTuple{(:species, :input_reactions, :output_reactions),Tuple{Int, Vector{Int},Vector{Int}}}})
    return reduced_network(N, M, [i.species for i in intermediates_result])
end




# To do: Decide how to deal with the trivial species (zero-rows in both N and M)
function reduced_network(N::QQMatrix, M::ZZMatrix, intermediate_species::Union{Nothing,Vector{Int}}=nothing)
    reactions = reaction_pairs(N, M)
    reduced_reactions = deepcopy(reactions)
    if isnothing(intermediate_species)
        intermediate_species = [i.species for i in intermediates(N, M)]
    end
    for i in intermediate_species
        in_reactions = [rxn for rxn in reduced_reactions if rxn.product == identity_matrix(ZZ, nrows(M))[:, i]]
        out_reactions = [rxn for rxn in reduced_reactions if rxn.reactant == identity_matrix(ZZ, nrows(M))[:, i]]
        for in_rxn in in_reactions
            for out_rxn in out_reactions
                # Combine only if it gives a nontrivial reaction
                if in_rxn.reactant != out_rxn.product
                    # New reaction bypassing the intermediate
                    push!(reduced_reactions, (reactant=in_rxn.reactant, product=out_rxn.product))
                end
            end
        end
        # Remove old reactions involving the intermediate
        reactions_to_remove = union(in_reactions, out_reactions)
        reduced_reactions = setdiff(reduced_reactions, reactions_to_remove)
    end

    # Deal with the special case of the reduced network being empty
    if length(reduced_reactions) == 0
        N_reduced = zero_matrix(QQ, nrows(M), 0)
        B_reduced = zero_matrix(ZZ, nrows(M), 0)
    else
        N_reduced, B_reduced = stoichiometric_and_kinetic_matrix_from_reaction_pairs(reduced_reactions)
    end

    # Remove the rows corresponding to the intermediate before returning
    new_species = [i for i = 1:nrows(M) if i ∉ intermediate_species]
    return N_reduced[new_species, :], B_reduced[new_species, :]
end