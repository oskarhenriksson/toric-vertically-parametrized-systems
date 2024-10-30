using Oscar
import Graphs

# Given a list of (candidate) intermediate vertices, find their inputs and outputs.
# If an intermediate vertex has no inputs or outputs, it is removed from the list of intermediates,
# and the analysis is re-done
function true_intermediate_vertices_with_inputs_outputs(g::Graphs.DiGraph, intermediate_vertex_candidates::Vector{Int}; only_single_input::Bool=false)

    result = []

    non_intermediate_vertices = setdiff(Graphs.vertices(g), intermediate_vertex_candidates)

    # Find inputs for each intermediate vertex
    for Y in intermediate_vertex_candidates
        inputs = Int[] 

        # For each non-intermediate vertex, run a DFS to find a path to Y
        for X in non_intermediate_vertices
            visited = Set{Int}()
            stack = [X]
            found_path = false

            while !isempty(stack) && !found_path
                current = pop!(stack) # last-in, first-out

                # Skip vertex if already visited
                if current in visited
                    continue
                end
                push!(visited, current)

                # If we have reached Y, record X as an input for Y
                if current == Y
                    push!(inputs, X)
                    found_path = true
                    break
                end

                # Add only intermediate neighbors for further exploration
                for neighbor in Graphs.outneighbors(g, current)
                    if neighbor in intermediate_vertex_candidates
                        push!(stack, neighbor)
                    end
                end
            end
        end

        # If no inputs were found, remove Y from list of intermediates and re-run
        if isempty(inputs)
            new_candidates = setdiff(intermediate_vertex_candidates, [Y])
            return true_intermediate_vertices_with_inputs_outputs(g, new_candidates, only_single_input=only_single_input)
        end

        # Remove if there are multiple inputs and we only wanted single inputs
        if only_single_input && length(inputs) > 1
            new_candidates = setdiff(intermediate_vertex_candidates, [Y])
            return true_intermediate_vertices_with_inputs_outputs(g, new_candidates, only_single_input=only_single_input)
        end

        # Find outputs for each intermediate vertex
        outputs = Int[]  
        visited = Set{Int}()
        stack = [Y]
        while !isempty(stack)
            current = pop!(stack)

            # Skip if already visited
            if current in visited
                continue
            end
            push!(visited, current)

            # If we reach a non-intermediate, record it as an output for Y
            if current in non_intermediate_vertices
                push!(outputs, current)
                continue
            end

            # Continue exploring neighbors
            for neighbor in Graphs.outneighbors(g, current)
                push!(stack, neighbor)
            end
        end

        # If no outputs were found, remove Y from list of intermediates and re-run
        if isempty(outputs)
            new_candidates = setdiff(intermediate_vertex_candidates, [Y])
            return true_intermediate_vertices_with_inputs_outputs(g, new_candidates, only_single_input=only_single_input)
        end

        push!(result, (intermediate_vertex=Y, inputs=inputs, outputs=outputs))
    end

    return result
end


function intermediates(N::QQMatrix, M::ZZMatrix; only_single_input::Bool=false)

    n = nrows(M)
    G, embedding = reaction_graph(N, M)
    vertex_list = keys(embedding)

    # Find the monomoleculuar complexes for which the species only appears in that complex
    intermediate_candidate_vertices = Int[]
    for v in vertex_list
        for i = 1:n
            if embedding[v] == unit_vector(n, i)
                if is_empty([w for w in vertex_list if (w != v && embedding[w][i] > 0)])
                    push!(intermediate_candidate_vertices, v)
                end
            end
        end
    end

    # Check which of these complexes are truly intermediates, and what their input/output vertices are
    intermediate_vertices = true_intermediate_vertices_with_inputs_outputs(G, intermediate_candidate_vertices, only_single_input=only_single_input)

    # Compute the species and the input/output complexes
    result = NamedTuple{(:species, :input_complexes, :output_complexes),Tuple{Int,Vector{Vector{Int}},Vector{Vector{Int}}}}[]
    for I in intermediate_vertices
        species = findfirst(!is_zero, embedding[I.intermediate_vertex])
        input_complexes = [Int.(embedding[v]) for v in I.inputs]
        output_complexes = [Int.(embedding[v]) for v in I.outputs]
        push!(result, (species=species, input_complexes=input_complexes, output_complexes=output_complexes))
    end

    return sort(result, by = x -> x.species)
end

function single_input_intermediates(N::QQMatrix, M::ZZMatrix)
    return intermediates(N, M, only_single_input=true)
end


function reduced_network(N::QQMatrix, M::ZZMatrix,
    intermediates_result::Vector{NamedTuple{(:species, :input_complexes, :output_complexes),Tuple{Int,Vector{Vector{Int}},Vector{Vector{Int}}}}})
    reactions = reaction_pairs(N, M)
    reduced_reactions = deepcopy(reactions)
    for I in intermediates_result
        i = I.species
        inputs = I.input_complexes
        outputs = I.output_complexes
        
        # Remove reactions where i participates
        for rxn in reduced_reactions
            if rxn.reactant[i] > 0 ||   rxn.product[i] > 0 
                reduced_reactions = setdiff(reduced_reactions, [rxn])
            end 
        end

        # Add reactions bypassing i
        for c in inputs
            for cprime in outputs
                if c != cprime
                    new_rxn = (reactant=c, product=cprime)
                    if new_rxn ∉ reduced_reactions
                        push!(reduced_reactions, new_rxn)
                    end 
                end
            end
        end
    end
    
    if length(reduced_reactions) == 0
        N_reduced = zero_matrix(QQ, nrows(M), 0)
        M_reduced = zero_matrix(ZZ, nrows(M), 0)
    else
        N_reduced, M_reduced = stoichiometric_and_kinetic_matrix_from_reaction_pairs(reduced_reactions)
    end

    non_intermediates = setdiff(1:nrows(M), [I.species for I in intermediates_result])

    return N_reduced[non_intermediates, :], M_reduced[non_intermediates, :]
end


function lift_exponent_matrix(Atilde::ZZMatrix, M::ZZMatrix, intermediates_result::Vector{NamedTuple{(:species, :input_complexes, :output_complexes),Tuple{Int,Vector{Vector{Int}},Vector{Vector{Int}}}}})
    n = nrows(M)
    intermediates = [i.species for i in intermediates_result]
    non_intermediates = setdiff(1:n, intermediates)
    A = zero_matrix(ZZ, nrows(Atilde), n)
    A[:, non_intermediates] = Atilde
    for I in intermediates_result
        @req length(I.input_complexes) == 1 "Only single-input intermediates are allowed"
        i = I.species
        A[:, i] = Atilde * first(I.input_complexes)[non_intermediates]
    end
    return A
end