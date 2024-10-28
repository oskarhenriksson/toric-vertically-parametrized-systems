include("src/main.jl");
include("tests/runtests.jl");

using Dates

# Function for writing to both terminal and file
function write_both(file, msg)
    println(msg)         # Write to terminal
    println(file, msg)   # Write to file
end

# Function for reading files of the ODE base format
function read_matrix(path)
    file_content = read(path, String)
    cleaned_content = replace(file_content, r"[<>\;]" => "")
    lines = split(cleaned_content, "\n")
    return [parse.(Int, split(line, ",")) for line in lines if !isempty(line)]
end

# Function for saving a csv file
function save_as_csv(vector::Vector, filename::String)
    open(filename, "w") do file
        write(file, join(vector, "\n"))
    end
end

timestamp_str = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")

# Load the models
choice_of_models = "all"
#choice_of_models = "ma"

if choice_of_models == "all"
    const model_directory = "../odebase_2023"
elseif choice_of_models == "ma"
    const model_directory = "../odebase_mass_action_oct2023"
end
list_of_models = readdir(model_directory)
number_of_models = length(list_of_models)


# Function for collecting the kinetic and stoichiometric matrices given a model id
function collect_matrices(model_id)
    B_path = joinpath(model_directory, model_id, "kinetic_matrix.txt")
    N_path = joinpath(model_directory, model_id, "reconfigured_stoichiometric_matrix.txt")
    B = matrix(ZZ, read_matrix(B_path))
    N = matrix(QQ, read_matrix(N_path))
    return B, N
end

error_reading = []
too_large = []
inconsistent = []
full_rank = []
analyzed = []
linear = []
dzt = []
dot = []
nondegenerate = []
degenerate = []
nonconserved_in_invariance_space = []
fulldimensional_invariance_space = []
skipped_injectivity = []
toric = []
finite = []
multistat = []
multistat_precluded = []
acr = []
local_acr = []
generic_bin = []
bin_for_all = []
has_irrelevant_species = []

# Run the analysis on each model
mkdir(timestamp_str * "_" * choice_of_models)
open(timestamp_str * "_" * choice_of_models * "/report.txt", "w") do file
    for model_id in list_of_models
        write_both(file, "")
        write_both(file, model_id)

        # Read and parse the kinetic matrix (B) and stoichiometric matrix (N)
        M = nothing
        N = nothing
        try
            M, N = collect_matrices(model_id)
        catch e
            write_both(file, "Error reading matrices")
            println(e)
            push!(error_reading, model_id)
            continue  # Skip to the next model
        end

        # Check for irrelevant species and remove them
        irrelevant_species = [i for i = 1:nrows(M) if (all(iszero, M[i, :]) && all(iszero, N[i, :]))]
        if !isempty(irrelevant_species)
            write_both(file, "Irrelevant species (ignored in subsequent analysis): $irrelevant_species")
            push!(has_irrelevant_species, model_id)
        end
        N = N[setdiff(1:nrows(N), irrelevant_species), :]
        M = M[setdiff(1:nrows(M), irrelevant_species), :]

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
            push!(too_large, model_id)
            continue
        end

        # Consistency
        if !nonempty_positive_kernel(C)
            write_both(file, "Inconsistent")
            push!(inconsistent, model_id)
            continue
        end

        if d == 0
            write_both(file, "Full rank")
            push!(full_rank, model_id)
            continue
        end

        push!(analyzed, model_id)

        # Check for linear kinetics
        if is_linear(M)
            write_both(file, "Linear kinetics")
            push!(linear, model_id)
            #continue
        end

        # The deficiency zero theorem
        if covered_by_deficiency_zero_theorem(N, M)
            write_both(file, "DZT: Covered")
            push!(dzt, model_id)
        else
            write_both(file, "DZT: Not covered")
        end

        # The deficiency one theorem
        if convered_by_deficiency_one_theorem(N, M)
            write_both(file, "DOT: Covered")
            push!(dot, model_id)
        else
            write_both(file, "DOT: Not covered")
        end

        # Conserved quantities
        W = kernel(N, side=:left)

        # Generic nondegeneracy
        if has_nondegenerate_zero(N, M)
            write_both(file, "Generically nondegenerate")
            push!(nondegenerate, model_id)
        else
            write_both(file, "Degenerate")
            push!(degenerate, model_id)
        end

        # Model reduction
        intermediates_result = intermediates(N, M, only_single_input=true)
        Ntilde, Mtilde = reduced_network(N, M, intermediates_result)
        Ctilde = row_space(Ntilde)

        # Toric invariance space
        Atilde = toric_invariance_space(Ctilde, Mtilde)
        A = lift_exponent_matrix(Atilde, M, intermediates_result)
        write_both(file, "Dimension of toric invariance space: $(rank(A))")

        if !all(is_zero, A * N)
            write_both(file, "Non-conserved quantities in toric invariance space")
            push!(nonconserved_in_invariance_space, model_id)
        end


        if rank(A) == d
            write_both(file, "Full-dimensional toric invariance space")
            push!(fulldimensional_invariance_space, model_id)

            toricity_flag = false
            all_nondeg = false

            if model_id in ["BIOMD0000000161", "BIOMD0000000250", "BIOMD0000000409", "BIOMD0000000428", "BIOMD0000000505", "BIOMD0000000594", "BIOMD0000000835", "BIOMD0000000640"]
                write_both(file, "Too slow for injectivity test! Skipped")
                push!(skipped_injectivity, model_id)
            else

                # Analyze the coset counting system for the reduced network
                if rank(Ctilde) + rank(Atilde) == nrows(Mtilde)

                    result = coset_counting_analysis(Ntilde, Mtilde, Atilde, printing_function=s -> write_both(file, s))
                    toricity_flag = result.toricity
                    finite_flag = result.finite

                end

                # Todo: Conclusion about infinitely many cosets?

                # Analyze the coset counting system for the full network
                if !toricity_flag && rank(C) + rank(A) == nrows(M)
                    result = coset_counting_analysis(N, M, A, printing_function=s -> write_both(file, s))

                    toricity_flag = result.toricity
                    finite_flag = result.finite

                end

                toricity_flag && push!(toric, model_id)
                finite_flag && push!(finite, model_id)

                # Capacity for multistationarity
                if coset_with_multistationarity(A, W)
                    write_both(file, "Coset with multistationarity found")
                    push!(multistat, model_id)
                else
                    write_both(file, "No coset with multistationarity found")
                    if toricity_flag
                        write_both(file, "Multistationarity precluded")
                        push!(multistat_precluded, model_id)
                    end
                end

                # (Local) ACR checks
                local_acr_species = zero_columns(A)
                if !isempty(local_acr_species)
                    if toricity_flag
                        write_both(file, "Species with ACR: $(local_acr_species)")
                        push!(acr, model_id)
                        push!(local_acr, model_id)
                    elseif finite_flag
                        write_both(file, "Species with local ACR: $(local_acr_species)")
                        push!(local_acr, model_id)
                    end
                end
            end
        end

        if model_id ∉ [
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
            "BIOMD0000000161", 
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
            "BIOMD0000000638",
            "BIOMD0000000835"]
            binomiality_result = binomiality_check(vertical_system(C, M))
            binomiality_result.generically && push!(generic_bin, model_id)
            binomiality_result.for_all_positive && push!(bin_for_all, model_id)
        else
            write_both(file, "Gröbner basis: Skipped")
        end

    end

    # Summerize the results
    write_both(file, "Number of analyzed networks: $(length(error_reading))")
    
    write_both(file, "\nError reading matrices:\n$(error_reading)")
    write_both(file, length(error_reading))

    write_both(file, "\nToo large models:\n$(too_large)")
    write_both(file, length(too_large))

    write_both(file, "\nInconsistent models:\n$(inconsistent)")
    write_both(file, length(inconsistent))

    write_both(file, "\nFull rank models:\n$(full_rank)")
    write_both(file, length(full_rank))

    write_both(file, "\nLinear models:\n$(linear)")
    write_both(file, length(linear))

    write_both(file, "\nDZT covered models:\n$(dzt)")
    write_both(file, length(dzt))

    write_both(file, "\nDOT covered models:\n$(dot)")
    write_both(file, length(dot))

    write_both(file, "\nNondegenerate models:\n$(nondegenerate)")
    write_both(file, length(nondegenerate))

    write_both(file, "\nDegenerate models:\n$(degenerate)")
    write_both(file, length(degenerate))

    write_both(file, "\nNonconserved in invariance space:\n$(nonconserved_in_invariance_space)")
    write_both(file, length(nonconserved_in_invariance_space))

    write_both(file, "\nFull-dimensional invariance space:\n$(fulldimensional_invariance_space)")
    write_both(file, length(fulldimensional_invariance_space))

    write_both(file, "\nSkipped injectivity test:\n$(skipped_injectivity)")
    write_both(file, length(skipped_injectivity))

    write_both(file, "\nToricity\n:$(toric)")
    write_both(file, length(toric))

    write_both(file, "\nLocal toricity:\n$(finite)")
    write_both(file, length(finite))

    write_both(file, "\nMultistationarity:\n$(multistat)")
    write_both(file, length(multistat))

    write_both(file, "\nMultistationarity precluded:\n$(multistat_precluded)")
    write_both(file, length(multistat_precluded))

    write_both(file, "\nACR:\n$(acr)")
    write_both(file, length(acr))

    write_both(file, "\nLocal ACR:\n$(local_acr)")
    write_both(file, length(local_acr))

    write_both(file, "\nGeneric binomiality:\n$(generic_bin)")
    write_both(file, length(generic_bin))

    write_both(file, "\nBinomiality for all positive:\n$(bin_for_all)")
    write_both(file, length(bin_for_all))

end

save_as_csv(analyzed, timestamp_str * "_" * choice_of_models * "/analyzed.csv")
save_as_csv(toric, timestamp_str * "_" * choice_of_models * "/toric_networks.csv")
save_as_csv(finite, timestamp_str * "_" * choice_of_models * "/locally_toric_networks.csv")
save_as_csv(multistat, timestamp_str * "_" * choice_of_models * "/multistationary_networks.csv")
save_as_csv(multistat_precluded, timestamp_str * "_" * choice_of_models * "/non_multistationarity_networks.csv")
save_as_csv(acr, timestamp_str * "_" * choice_of_models * "/acr_networks.csv")
save_as_csv(local_acr, timestamp_str * "_" * choice_of_models * "/local_acr_networks.csv")
save_as_csv(generic_bin, timestamp_str * "_" * choice_of_models * "/generic_binomiality_networks.csv")
save_as_csv(bin_for_all, timestamp_str * "_" * choice_of_models * "/binomiality_for_all_positive_networks.csv")
save_as_csv(has_irrelevant_species, timestamp_str * "_" * choice_of_models * "/irrelevant_species.csv")
save_as_csv(error_reading, timestamp_str * "_" * choice_of_models * "/error_reading.csv")
save_as_csv(too_large, timestamp_str * "_" * choice_of_models * "/too_large.csv")
save_as_csv(inconsistent, timestamp_str * "_" * choice_of_models * "/inconsistent.csv")
save_as_csv(full_rank, timestamp_str * "_" * choice_of_models * "/full_rank.csv")
save_as_csv(linear, timestamp_str * "_" * choice_of_models * "/linear.csv")
save_as_csv(dzt, timestamp_str * "_" * choice_of_models * "/dzt.csv")
save_as_csv(dot, timestamp_str * "_" * choice_of_models * "/dot.csv")
save_as_csv(nondegenerate, timestamp_str * "_" * choice_of_models * "/nondegenerate.csv")
save_as_csv(degenerate, timestamp_str * "_" * choice_of_models * "/degenerate.csv")
save_as_csv(nonconserved_in_invariance_space, timestamp_str * "_" * choice_of_models * "/nonconserved_in_invariance_space.csv")
save_as_csv(fulldimensional_invariance_space, timestamp_str * "_" * choice_of_models * "/fulldimensional_invariance_space.csv")
save_as_csv(skipped_injectivity, timestamp_str * "_" * choice_of_models * "/skipped_injectivity.csv")
