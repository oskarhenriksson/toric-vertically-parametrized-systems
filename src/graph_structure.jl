using Oscar
import Graphs 

# Graph theoretic helper functions
function is_terminal_strongly_connected_component(g, scc)
    for v in scc
        for w in Graphs.outneighbors(g, v)
            if w ∉ scc
                return false
            end
        end
    end
    return true
end

function terminal_strongly_connected_components(g::Graphs.SimpleGraphs.SimpleDiGraph)
    sccs = Graphs.strongly_connected_components(g)
    return [scc for scc in sccs if is_terminal_strongly_connected_component(g, scc)]
end

function is_weakly_reversible(g::Graphs.SimpleGraphs.SimpleDiGraph)
    weak_components = Graphs.weakly_connected_components(g)
    for component in weak_components
        subgraph = Graphs.induced_subgraph(g, component)[1]
        if !Graphs.is_strongly_connected(subgraph)
            return false
        end
    end
    return true
end

# Find the complexes of a network
function complexes(N::QQMatrix, M::ZZMatrix)
    reactions = reaction_pairs(N, M)
    return unique(vcat([rxn.reactant for rxn in reactions], [rxn.product for rxn in reactions]))
end

# Construct an (abstract) reaction graph, together a with a vertext-to-complex map
function reaction_graph(N::QQMatrix, M::ZZMatrix)
    reactions = reaction_pairs(N, M)
    vertices = unique(vcat([rxn.reactant for rxn in reactions], [rxn.product for rxn in reactions]))
    index_to_complex_map = Dict(i => v for (i, v) in enumerate(vertices))
    complex_to_index_map = Dict(v => i for (i, v) in enumerate(vertices))
    edges = [(complex_to_index_map[rxn.reactant], complex_to_index_map[rxn.product]) for rxn in reactions]
    g = Graphs.SimpleDiGraphFromIterator(Graphs.Edge.(edges))
    return g, index_to_complex_map
end

# Reconstruct a network from graph and vertex-to-complex map
function network_from_embedded_graph(g, index_to_complex_map)
    reactions = [(reactant=index_to_complex_map[e.src], product=index_to_complex_map[e.dst]) for e in Graphs.edges(g)]
    return stoichiometric_and_kinetic_matrix_from_reaction_pairs(reactions)
end

# Decompose a network in linkage classes
function linkage_classes(N::QQMatrix, M::ZZMatrix)
    g, complex_to_index_map = reaction_graph(N, M)
    weak_components = Graphs.weakly_connected_components(g)
    component_pairs = [
        [(reactant=complex_to_index_map[i], product=complex_to_index_map[j]) for i in comp for j in Graphs.outneighbors(g, i)]
        for comp in weak_components
    ]
    return stoichiometric_and_kinetic_matrix_from_reaction_pairs.(component_pairs)
end

# Compute the deficiency of a network
function delta(N::QQMatrix, M::ZZMatrix)
    g, _ = reaction_graph(N, M)
    return length(Graphs.vertices(g)) - rank(N) - length(Graphs.weakly_connected_components(g))
end

function covered_by_deficiency_zero_theorem(N::QQMatrix, M::ZZMatrix)
    g, _ = reaction_graph(N, M)
    return (is_weakly_reversible(g) && delta(N, M) == 0)
end

# Check if a network is covered by the deficiency one theorem
function convered_by_deficiency_one_theorem(N::QQMatrix, M::ZZMatrix)
    g, embedding = reaction_graph(N, M)

    # There should be one terminal strong compolnent per weak component    
    t = length(terminal_strongly_connected_components(g))
    ell = length(Graphs.weakly_connected_components(g))
    if t > ell
        return false
    end

    # All linkage classes should have deficiency at most 1
    deficiencies = [delta(Ni, Mi) for (Ni, Mi) in linkage_classes(N, M)]
    if any(d_i > 1 for d_i in deficiencies)
        return false
    end

    # The sum of the deficiencies should be equal to the network deficiency
    if sum(deficiencies) != delta(N, M)
        return false
    end

    # If all the checks above are passes, the network is covered
    return true
end
