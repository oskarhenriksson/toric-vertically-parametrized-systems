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
        write(file, join(vector, ","))
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

        if ncols(M) < 10
            binomiality_result = binomiality_check(vertical_system(C, M))
            binomiality_result.generically && push!(generic_bin, model_id)
            binomiality_result.for_all_positive && push!(bin_for_all, model_id)
        else
            write_both(file, "Gröbner basis: Skipped")
        end

    end
end
