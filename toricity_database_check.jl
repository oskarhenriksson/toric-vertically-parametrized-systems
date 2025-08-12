using DataFrames, CSV, Dates
include("src/main.jl");
include("tests/runtests.jl");


# Function for writing to both terminal and file
function write_both(file, msg)
    println(msg)         # Write to terminal
    println(file, msg)   # Write to file
end

# Function for reading files of the ODEbase format
function read_matrix(path)
    file_content = read(path, String)
    cleaned_content = replace(file_content, r"[<>\;]" => "")
    lines = split(cleaned_content, "\n")
    return [parse.(Int, split(line, ",")) for line in lines if !isempty(line)]
end

# Function for saving a csv file
function save_as_txt(vector::Vector, filename::String)
    open(filename, "w") do file
        write(file, join(vector, "\n"))
    end
end

timestamp_str = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")

# Load the models
choice_of_models = "all"
#choice_of_models = "ma"

if choice_of_models == "all"
    const model_directory = "./odebase_2023"
elseif choice_of_models == "ma"
    const model_directory = "./odebase_mass_action_oct2023"
end
list_of_models = readdir(model_directory)
list_of_mass_action_models = readdir("./odebase_mass_action_oct2023")
number_of_models = length(list_of_models)


# Function for collecting the kinetic and stoichiometric matrices given a model id
function reactant_and_stoichiometric_matrices(model_id)
    M_path = joinpath(model_directory, model_id, "kinetic_matrix.txt")
    N_path = joinpath(model_directory, model_id, "reconfigured_stoichiometric_matrix.txt")
    M = matrix(ZZ, read_matrix(M_path))
    N = matrix(QQ, read_matrix(N_path))
    return M, N
end

# Create a data frame for storing the data
df = DataFrame(
    ID = String[],
    OriginallyMassAction = Bool[],
    IrrelevantSpecies = Vector{Int}[],
    NumberOfSpeices = Int[],
    NumberOfReactions = Int[],
    Consistent = Bool[],
    Nondegenerate = Union{Bool,Missing}[],
    ####
    ToricRank = Int[],
    MaximalToricRank = Union{Bool,Missing}[],
    Toric = Union{Bool,Missing}[],
    LocallyToric = Union{Bool,Missing}[],
    TrivialFundamentalPartition = Bool[],
    Quasihomogeneous = Bool[],
    #####
    Multistationarity = Union{Bool,Missing}[],
    ACR = Union{Bool,Missing}[],
    LocalACR = Union{Bool,Missing}[],
    #####
    GenericallyBinomial = Union{Bool,Missing}[],
    Binomial = Union{Bool,Missing}[],
    CoveredByDZT = Bool[],
    CoveredByDOT = Bool[],
)


error_reading_list = []
too_large_list = []
inconsistent_list = []
full_rank_list = []
analyzed_list = []
linear_list = []
dzt_list = []
dot_list = []
nondegenerate_list = []
degenerate_list = []
inf_cosets_list = []
nonconserved_in_invariance_space_list = []
maximal_toric_rank_list = []
skipped_injectivity_list = []
trivial_fundamental_partition_list = []
quasi_homogeneous_list = []
toric_list = []
locally_toric_list = []
multistat_list = []
multistat_precluded_list = []
acr_list = []
local_acr_list = []
non_bin_list = []
generic_bin_list = []
bin_for_all_list = []
has_irrelevant_species_list = []
has_irrelevant_reactions_list = []
gb_skip_list = [
    "BIOMD0000000002",
    "BIOMD0000000028",
    "BIOMD0000000030",
    "BIOMD0000000032",
    "BIOMD0000000048",
    "BIOMD0000000061",
    "BIOMD0000000070",
    "BIOMD0000000085",
    "BIOMD0000000086",
    "BIOMD0000000108",
    "BIOMD0000000123",
    "BIOMD0000000162",
    "BIOMD0000000164",
    "BIOMD0000000165",
    "BIOMD0000000182",
    "BIOMD0000000200",
    "BIOMD0000000237",
    "BIOMD0000000250",
    "BIOMD0000000294",
    "BIOMD0000000409",
    "BIOMD0000000430",
    "BIOMD0000000431",
    "BIOMD0000000500",
    "BIOMD0000000530",
    "BIOMD0000000637",
    "BIOMD0000000638"]

# Run the analysis on each model
mkdir(timestamp_str * "_" * choice_of_models)
open(timestamp_str * "_" * choice_of_models * "/report.txt", "w") do file
    for model_id in list_of_models
        
        # Skip hidden files (file names starting with period)
        if model_id[1] == '.'
            continue
        end
    
        # Print out the model ID
        current_model_id = model_id
        write_both(file, "")
        write_both(file, model_id)

        # Read and parse the kinetic matrix (M) and stoichiometric matrix (N)
        M = nothing
        N = nothing
        try
            M, N = reactant_and_stoichiometric_matrices(model_id)
        catch e
            write_both(file, "Error reading matrices")
            println(e)
            push!(error_reading_list, model_id)
            continue  # Skip to the next model
        end

        # Check for irrelevant species and remove them
        irrelevant_species = [i for i = 1:nrows(M) if (all(iszero, M[i, :]) && all(iszero, N[i, :]))]
        if !isempty(irrelevant_species)
            write_both(file, "Irrelevant species (ignored in subsequent analysis): $irrelevant_species")
            push!(has_irrelevant_species_list, model_id)
        end
        N = N[setdiff(1:nrows(N), irrelevant_species), :]
        M = M[setdiff(1:nrows(M), irrelevant_species), :]

        # TODO: Decide what to do with this!
        irrelevant_reactionns = [i for i = 1:ncols(N) if all(iszero, N[:, i])]
        if !isempty(irrelevant_reactionns)
            write_both(file, "Irrelevant reactions: $irrelevant_reactionns")
            push!(has_irrelevant_reactions_list, model_id)
        end

        # Basic info
        n = nrows(M)
        r = ncols(M)
        s = rank(N)
        d = n - s
        write_both(file, "n=$n,r=$r,d=$d")

        # Row reduced stoichiometric matrix
        C = row_space(N)

        if r > 100
            write_both(file, "Too many reactions")
            push!(too_large_list, model_id)
            continue
        end

        # Consistency
        consistent = nonempty_positive_kernel(C)
        if !consistent
            write_both(file, "Inconsistent")
            push!(inconsistent_list, model_id)
            continue
        end

        if d == 0
            write_both(file, "Full rank")
            push!(full_rank_list, model_id)
            continue
        end

        # Check for linear kinetics
        linear = is_linear(M)
        if linear
            write_both(file, "Linear kinetics")
            push!(linear_list, model_id)
            continue
        end

        push!(analyzed_list, model_id)

        # The deficiency zero theorem
        covered_by_dzt = covered_by_deficiency_zero_theorem(N, M)
        if covered_by_dzt
            write_both(file, "DZT: Covered")
            push!(dzt_list, model_id)
        else
            write_both(file, "DZT: Not covered")
        end

        # The deficiency one theorem
        covered_by_dot = covered_by_deficiency_one_theorem(N, M)
        if covered_by_dot
            write_both(file, "DOT: Covered")
            push!(dot_list, model_id)
        else
            write_both(file, "DOT: Not covered")
        end

        # Conserved quantities
        W = kernel(N, side=:left)

        # Generic nondegeneracy
        nondeg = has_nondegenerate_zero(N, M)
        if nondeg
            write_both(file, "Generically nondegenerate")
            push!(nondegenerate_list, model_id)
        else
            write_both(file, "Degenerate")
            push!(degenerate_list, model_id)
        end

        # Model reduction
        intermediates_result = intermediates(N, M, only_single_input=true)
        Ntilde, Mtilde = reduced_network(N, M, intermediates_result)
        Ctilde = row_space(Ntilde)

        # Fundamental partition
        FP = fundamental_partition(C)
        trivial_fundamental_partition = (length(FP) == 1)
        if trivial_fundamental_partition
            write_both(file, "Trivial fundamental partition")
            push!(trivial_fundamental_partition_list, model_id)
        end

        # Toric invariance space
        Atilde = toric_invariance_group(Ctilde, Mtilde)
        A = toric_invariance_group(C, M)

        # Quasihomogeneity
        quasi_homogeneous_with_positive_toric_rank = is_quasi_homogeneous(C, M, A) && rank(A)>0
        if quasi_homogeneous_with_positive_toric_rank
            write_both(file, "Quasi-homogeneous with respect to A (and positive toric rank)")
            push!(quasi_homogeneous_list, model_id)
        end

        Alifted = lift_exponent_matrix(Atilde, M, intermediates_result)
        if row_space(A) != row_space(Alifted)
            write_both(file, "Error in toric invariance space")
            continue
        end 
        
        
        write_both(file, "Toric rank: $(rank(A))")

        if !all(is_zero, A * N)
            write_both(file, "Non-conserved quantities in toric invariance space")
            push!(nonconserved_in_invariance_space_list, model_id)
        end

        infinitely_many_cosets = (rank(A) < d && nondeg)
        if infinitely_many_cosets
            write_both(file, "Infinitely many cosets in open region")
            push!(inf_cosets_list, model_id)
        end

        maximal_toric_rank = (rank(A) == d)

        if maximal_toric_rank
            write_both(file, "Maximal toric rank")
            push!(maximal_toric_rank_list, model_id)

            toric = false
            all_nondeg = false

            # Analyze the coset counting system for the reduced network
            if rank(Ctilde) + rank(Atilde) == nrows(Mtilde)

                result = coset_counting_analysis(Ntilde, Mtilde, Atilde, printing_function=s -> write_both(file, s))
                toric = result.toricity
                locally_toric = result.finite

            end

            # Analyze the coset counting system for the full network
            if !toric && rank(C) + rank(A) == nrows(M)
                result = coset_counting_analysis(N, M, A, printing_function=s -> write_both(file, s))

                toric = result.toricity
                locally_toric = result.finite

            end

            toric && push!(toric_list, model_id)
            locally_toric && push!(locally_toric_list, model_id)

            # Capacity for multistationarity
            if model_id in ["BIOMD0000000161"]
                write_both(file, "Skipped multistationarity test!")
            else
                if coset_with_multistationarity(A, W)
                    write_both(file, "Coset with multistationarity found")
                    push!(multistat_list, model_id)
                else
                    write_both(file, "No coset with multistationarity found")
                    if toric
                        write_both(file, "Multistationarity precluded")
                        push!(multistat_precluded_list, model_id)
                    end
                end
            end

            # (Local) ACR checks
            local_acr_species = zero_columns(A)
            if !isempty(local_acr_species)
                if toric
                    write_both(file, "Species with ACR: $(local_acr_species)")
                    push!(acr_list, model_id)
                    push!(local_acr_list, model_id)
                elseif locally_toric
                    write_both(file, "Species with local ACR: $(local_acr_species)")
                    push!(local_acr_list, model_id)
                end
            end 
        end

        # Check for binomiality (first trivial cases)
        if all([length(supp(Ctilde[i,:]))==2 for i=1:nrows(Ctilde)]) || ncols(Mtilde) == 0
            write_both(file, "Gröbner basis: Generically binomial!")
            write_both(file, "Specialization for all rate constants: Verified")
            push!(generic_bin_list, current_model_id)
            push!(bin_for_all_list, current_model_id)
        elseif current_model_id in gb_skip_list
            write_both(file, "Gröbner basis: Skipped")
        else
            binomiality_result = binomiality_check(vertical_system(Ctilde, Mtilde), printing_function=s -> write_both(file, s), verbose=true)
            binomiality_result.generically ? push!(generic_bin_list, current_model_id) : push!(non_bin_list, current_model_id)
            binomiality_result.for_all_positive && push!(bin_for_all_list, model_id)
        end

    end


    # Summerize the results
    write_both(file, "\nNumber of analyzed networks: $(length(error_reading_list))")

    write_both(file, "\nError reading matrices:\n$(error_reading_list)")
    write_both(file, length(error_reading_list))

    write_both(file, "\nToo large models:\n$(too_large_list)")
    write_both(file, length(too_large_list))

    write_both(file, "\nInconsistent models:\n$(inconsistent_list)")
    write_both(file, length(inconsistent_list))

    write_both(file, "\nFull rank models:\n$(full_rank_list)")
    write_both(file, length(full_rank_list))

    write_both(file, "\nLinear models:\n$(linear_list)")
    write_both(file, length(linear_list))

    write_both(file, "\nDZT covered models:\n$(dzt_list)")
    write_both(file, length(dzt_list))

    write_both(file, "\nDOT covered models:\n$(dot_list)")
    write_both(file, length(dot_list))

    write_both(file, "\nNondegenerate models:\n$(nondegenerate_list)")
    write_both(file, length(nondegenerate_list))

    write_both(file, "\nDegenerate models:\n$(degenerate_list)")
    write_both(file, length(degenerate_list))

    write_both(file, "\nInfinitely many cosets:\n$(inf_cosets_list)")
    write_both(file, length(inf_cosets_list))

    write_both(file, "\nNonconserved in invariance space:\n$(nonconserved_in_invariance_space_list)")
    write_both(file, length(nonconserved_in_invariance_space_list))

    write_both(file, "\nFull-dimensional invariance space:\n$(maximal_toric_rank_list)")
    write_both(file, length(maximal_toric_rank_list))

    write_both(file, "\nSkipped injectivity test:\n$(skipped_injectivity_list)")
    write_both(file, length(skipped_injectivity_list))

    write_both(file, "\nToricity\n:$(toric_list)")
    write_both(file, length(toric_list))

    write_both(file, "\nLocal toricity:\n$(locally_toric_list)")
    write_both(file, length(locally_toric_list))

    write_both(file, "\nLocally toric but not toric:\n$(setdiff(locally_toric_list, toric_list))")
    write_both(file, length(setdiff(locally_toric_list, toric_list)))

    write_both(file, "\nToric but not covered by dzt or dot:\n$(setdiff(toric_list, union(dzt_list, dot_list)))")
    write_both(file, length(setdiff(toric_list, union(dzt_list, dot_list))))

    write_both(file, "\nMultistationarity:\n$(multistat_list)")
    write_both(file, length(multistat_list))

    write_both(file, "\nMultistationarity precluded:\n$(multistat_precluded_list)")
    write_both(file, length(multistat_precluded_list))

    write_both(file, "\nACR:\n$(acr_list)")
    write_both(file, length(acr_list))

    write_both(file, "\nLocal ACR:\n$(local_acr_list)")
    write_both(file, length(local_acr_list))

    write_both(file, "\nNon-binomial:\n$(non_bin_list)")
    write_both(file, length(non_bin_list))

    write_both(file, "\nGeneric binomiality:\n$(generic_bin_list)")
    write_both(file, length(generic_bin_list))

    write_both(file, "\nBinomiality for all positive:\n$(bin_for_all_list)")
    write_both(file, length(bin_for_all_list))

end

for model_id in analyzed_list

    M, N = reactant_and_stoichiometric_matrices(model_id)
    irrelevant_species = [i for i = 1:nrows(M) if (all(iszero, M[i, :]) && all(iszero, N[i, :]))]
    N = N[setdiff(1:nrows(N), irrelevant_species), :]
    M = M[setdiff(1:nrows(M), irrelevant_species), :]
    C = row_space(N)
    A = toric_invariance_group(C, M)
    rank(A)
    new_entry = (ID=model_id,
        OriginallyMassAction = (model_id in list_of_mass_action_models),
        IrrelevantSpecies = irrelevant_species,
        NumberOfSpeices = nrows(M),
        NumberOfReactions = ncols(M),
        Consistent = !(model_id in inconsistent_list),
        Nondegenerate = (model_id in nondegenerate_list),
        ####
        ToricRank = rank(A),
        MaximalToricRank = (model_id in maximal_toric_rank_list),
        Toric = (model_id in toric_list),
        LocallyToric = (model_id in locally_toric_list),
        TrivialFundamentalPartition = (model_id in trivial_fundamental_partition_list),
        Quasihomogeneous = (model_id in quasi_homogeneous_list),
        ####
        Multistationarity = (model_id in multistat_list),
        ACR = (model_id in acr_list),
        LocalACR = (model_id in local_acr_list),
        ####
        GenericallyBinomial = (model_id in generic_bin_list),
        Binomial = (model_id in bin_for_all_list),
        CoveredByDZT = (model_id in dzt_list),
        CoveredByDOT = (model_id in dot_list)
    )

    push!(df, new_entry)

    if model_id in gb_skip_list
        last(df).GenericallyBinomial = missing
        last(df).Binomial = missing
    end

end

CSV.write(timestamp_str * "_" * choice_of_models * "/results.csv", df)
