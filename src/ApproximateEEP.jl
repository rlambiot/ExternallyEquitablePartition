"""
    ApproximateEEP.jl

Find approximately externally equitable partitions using Metropolis-Hastings
Monte Carlo optimization.

The "energy" of a partition measures how far it is from being externally 
equitable. The algorithm explores the space of partitions, accepting moves
that decrease energy (or occasionally increase it, to escape local minima).
"""
module ApproximateEEP

using Random
using LinearAlgebra
using Statistics
using Dates
using JSON

export find_approximate_eep,
       compute_energy,
       compute_defect,
       compute_imbalance,
       run_metropolis,
       run_analysis_approximate

#=============================================================================
  Energy/Defect Functions
=============================================================================#

"""
    compute_defect(A::AbstractMatrix, partition::Vector{Vector{Int}})

Compute the external equitability defect of a partition.

The defect measures how far a partition is from being externally equitable.
For each pair of cells (i,j) with i≠j, we compute the variance of the number
of neighbors that vertices in cell i have in cell j.

Returns:
- total_defect: Sum of variances across all cell pairs
- defect_matrix: Matrix where entry (i,j) is the variance for that pair
"""
function compute_defect(A::AbstractMatrix, partition::Vector{Vector{Int}})
    k = length(partition)
    defect_matrix = zeros(Float64, k, k)
    
    for i in 1:k
        cell_i = partition[i]
        if length(cell_i) <= 1
            continue  # Single vertex has zero variance
        end
        
        for j in 1:k
            if i == j
                continue  # Only external cells
            end
            
            cell_j = partition[j]
            
            # Compute neighbor counts for each vertex in cell_i to cell_j
            counts = [sum(A[v, u] for u in cell_j; init=0) for v in cell_i]
            
            # Variance of these counts
            if length(counts) > 1
                defect_matrix[i, j] = var(counts; corrected=false)
            end
        end
    end
    
    total_defect = sum(defect_matrix)
    return total_defect, defect_matrix
end

"""
    compute_imbalance(partition::Vector{Vector{Int}}, n::Int)

Compute the imbalance of a partition.

Imbalance = 1 - EffectiveSize / NumCells

- Imbalance = 0: All cells have equal size (perfectly balanced)
- Imbalance → 1: One cell has almost all vertices (maximally imbalanced)
"""
function compute_imbalance(partition::Vector{Vector{Int}}, n::Int)
    k = length(partition)
    if k <= 1
        return 0.0
    end
    
    # Compute effective size
    H = 0.0
    for cell in partition
        p = length(cell) / n
        if p > 0
            H -= p * log(p)
        end
    end
    eff_size = exp(H)
    
    # Imbalance: how far from perfectly balanced
    return 1.0 - eff_size / k
end

"""
    compute_energy(A::AbstractMatrix, partition::Vector{Vector{Int}}; 
                   α::Float64=1.0, β::Float64=0.0, γ::Float64=0.0)

Compute the total energy of a partition.

Energy = α * defect + β * num_cells + γ * imbalance

Parameters:
- α: Weight for equitability defect (higher = more equitable)
- β: Weight for number of cells (higher = prefer fewer cells)
- γ: Weight for imbalance (higher = prefer balanced cell sizes)

The default (α=1, β=0, γ=0) optimizes purely for equitability.

Note: Imbalance = 1 - EffectiveSize/NumCells, so:
- Imbalance = 0 when all cells are equal size
- Imbalance → 1 when one cell dominates
"""
function compute_energy(A::AbstractMatrix, partition::Vector{Vector{Int}}; 
                        α::Float64=1.0, β::Float64=0.0, γ::Float64=0.0)
    n = size(A, 1)
    
    # Defect term (how far from externally equitable)
    defect, _ = compute_defect(A, partition)
    
    # Number of cells term
    k = length(partition)
    
    # Imbalance term (how far from equal-sized cells)
    imbalance = compute_imbalance(partition, n)
    
    return α * defect + β * k + γ * imbalance
end

#=============================================================================
  Metropolis Moves
=============================================================================#

"""
    propose_move_vertex(partition::Vector{Vector{Int}}, n::Int)

Propose moving a random vertex to a different cell.
Returns (new_partition, move_info) or (nothing, nothing) if move is invalid.
"""
function propose_move_vertex(partition::Vector{Vector{Int}}, n::Int)
    k = length(partition)
    if k < 2
        return nothing, nothing
    end
    
    # Pick a random non-singleton cell to move from
    valid_sources = findall(c -> length(c) > 1, partition)
    if isempty(valid_sources)
        return nothing, nothing
    end
    
    source_idx = rand(valid_sources)
    source_cell = partition[source_idx]
    
    # Pick a random vertex to move
    vertex = rand(source_cell)
    
    # Pick a random different cell to move to
    target_idx = rand(setdiff(1:k, source_idx))
    
    # Create new partition
    new_partition = [copy(c) for c in partition]
    filter!(v -> v != vertex, new_partition[source_idx])
    push!(new_partition[target_idx], vertex)
    
    # Remove empty cells
    filter!(!isempty, new_partition)
    
    return new_partition, (type=:move, vertex=vertex, from=source_idx, to=target_idx)
end

"""
    propose_swap_vertices(partition::Vector{Vector{Int}}, n::Int)

Propose swapping two vertices between different cells.
"""
function propose_swap_vertices(partition::Vector{Vector{Int}}, n::Int)
    k = length(partition)
    if k < 2
        return nothing, nothing
    end
    
    # Pick two different cells
    idx1, idx2 = rand(1:k, 2)
    while idx1 == idx2
        idx2 = rand(1:k)
    end
    
    cell1, cell2 = partition[idx1], partition[idx2]
    if isempty(cell1) || isempty(cell2)
        return nothing, nothing
    end
    
    # Pick random vertices
    v1 = rand(cell1)
    v2 = rand(cell2)
    
    # Create new partition with swapped vertices
    new_partition = [copy(c) for c in partition]
    filter!(v -> v != v1, new_partition[idx1])
    filter!(v -> v != v2, new_partition[idx2])
    push!(new_partition[idx1], v2)
    push!(new_partition[idx2], v1)
    
    return new_partition, (type=:swap, v1=v1, v2=v2, cell1=idx1, cell2=idx2)
end

"""
    propose_merge_cells(partition::Vector{Vector{Int}}, n::Int)

Propose merging two random cells into one.
"""
function propose_merge_cells(partition::Vector{Vector{Int}}, n::Int)
    k = length(partition)
    if k < 3  # Need at least 3 cells to merge and stay non-trivial
        return nothing, nothing
    end
    
    # Pick two different cells to merge
    idx1, idx2 = rand(1:k, 2)
    while idx1 == idx2
        idx2 = rand(1:k)
    end
    
    # Create new partition with merged cells
    new_partition = Vector{Int}[]
    for (i, cell) in enumerate(partition)
        if i == idx1
            push!(new_partition, vcat(partition[idx1], partition[idx2]))
        elseif i != idx2
            push!(new_partition, copy(cell))
        end
    end
    
    return new_partition, (type=:merge, cell1=idx1, cell2=idx2)
end

"""
    propose_split_cell(partition::Vector{Vector{Int}}, n::Int)

Propose splitting a random cell into two.
"""
function propose_split_cell(partition::Vector{Vector{Int}}, n::Int)
    # Find cells with at least 2 vertices
    splittable = findall(c -> length(c) >= 2, partition)
    if isempty(splittable)
        return nothing, nothing
    end
    
    idx = rand(splittable)
    cell = partition[idx]
    
    # Random split
    shuffled = shuffle(cell)
    split_point = rand(1:length(cell)-1)
    cell1 = shuffled[1:split_point]
    cell2 = shuffled[split_point+1:end]
    
    # Create new partition
    new_partition = Vector{Int}[]
    for (i, c) in enumerate(partition)
        if i == idx
            push!(new_partition, cell1)
            push!(new_partition, cell2)
        else
            push!(new_partition, copy(c))
        end
    end
    
    return new_partition, (type=:split, cell=idx)
end

"""
    propose_move(partition::Vector{Vector{Int}}, n::Int; 
                 move_probs::NTuple{4,Float64}=(0.5, 0.2, 0.15, 0.15))

Propose a random move with given probabilities for each move type.
"""
function propose_move(partition::Vector{Vector{Int}}, n::Int; 
                      move_probs::NTuple{4,Float64}=(0.5, 0.2, 0.15, 0.15))
    r = rand()
    cumprob = cumsum(collect(move_probs))
    
    if r < cumprob[1]
        return propose_move_vertex(partition, n)
    elseif r < cumprob[2]
        return propose_swap_vertices(partition, n)
    elseif r < cumprob[3]
        return propose_merge_cells(partition, n)
    else
        return propose_split_cell(partition, n)
    end
end

#=============================================================================
  Metropolis-Hastings Algorithm
=============================================================================#

"""
    run_metropolis(A::AbstractMatrix;
                   n_steps::Int=10000,
                   T_init::Float64=1.0,
                   T_final::Float64=0.01,
                   cooling::Symbol=:exponential,
                   α::Float64=1.0,
                   β::Float64=0.0,
                   γ::Float64=0.0,
                   initial_partition::Union{Nothing,Vector{Vector{Int}}}=nothing,
                   initial_k::Int=5,
                   seed::Union{Nothing,Int}=nothing,
                   verbose::Bool=false,
                   record_history::Bool=false)

Run Metropolis-Hastings Monte Carlo to find approximately equitable partitions.

# Arguments
- `A`: Adjacency matrix
- `n_steps`: Number of Monte Carlo steps
- `T_init`: Initial temperature
- `T_final`: Final temperature (for annealing)
- `cooling`: Cooling schedule (:exponential, :linear, :constant)
- `α, β, γ`: Energy function weights (defect, eff_size, n_cells)
- `initial_partition`: Starting partition (random if nothing)
- `initial_k`: Number of initial cells if starting random
- `seed`: Random seed
- `verbose`: Print progress
- `record_history`: Record energy history for plotting

# Returns
- `best_partition`: Partition with lowest energy found
- `best_energy`: Energy of best partition
- `final_partition`: Partition at end of run
- `stats`: Dictionary with run statistics
"""
function run_metropolis(A::AbstractMatrix;
                        n_steps::Int=10000,
                        T_init::Float64=1.0,
                        T_final::Float64=0.01,
                        cooling::Symbol=:exponential,
                        α::Float64=1.0,
                        β::Float64=0.0,
                        γ::Float64=0.0,
                        initial_partition::Union{Nothing,Vector{Vector{Int}}}=nothing,
                        initial_k::Int=5,
                        seed::Union{Nothing,Int}=nothing,
                        verbose::Bool=false,
                        record_history::Bool=false)
    
    if seed !== nothing
        Random.seed!(seed)
    end
    
    n = size(A, 1)
    
    # Initialize partition
    if initial_partition === nothing
        # Random partition into initial_k cells
        k = min(initial_k, n)
        assignments = rand(1:k, n)
        # Ensure all cells non-empty
        for i in 1:k
            if !(i in assignments)
                assignments[rand(1:n)] = i
            end
        end
        partition = [findall(==(i), assignments) for i in 1:k]
        filter!(!isempty, partition)
    else
        partition = [copy(c) for c in initial_partition]
    end
    
    # Compute initial energy
    current_energy = compute_energy(A, partition; α=α, β=β, γ=γ)
    best_energy = current_energy
    best_partition = [copy(c) for c in partition]
    
    # Statistics
    accepted = 0
    rejected = 0
    history = record_history ? Float64[] : nothing
    
    # Temperature schedule
    if cooling == :exponential
        cooling_rate = (T_final / T_init)^(1 / n_steps)
    elseif cooling == :linear
        cooling_rate = (T_init - T_final) / n_steps
    end
    
    T = T_init
    
    for step in 1:n_steps
        # Propose move
        new_partition, move_info = propose_move(partition, n)
        
        if new_partition === nothing
            rejected += 1
            continue
        end
        
        # Compute new energy
        new_energy = compute_energy(A, new_partition; α=α, β=β, γ=γ)
        
        # Metropolis criterion
        ΔE = new_energy - current_energy
        
        if ΔE < 0 || rand() < exp(-ΔE / T)
            # Accept move
            partition = new_partition
            current_energy = new_energy
            accepted += 1
            
            # Update best
            if current_energy < best_energy
                best_energy = current_energy
                best_partition = [copy(c) for c in partition]
            end
        else
            rejected += 1
        end
        
        # Record history
        if record_history
            push!(history, current_energy)
        end
        
        # Update temperature
        if cooling == :exponential
            T *= cooling_rate
        elseif cooling == :linear
            T = max(T_final, T - cooling_rate)
        end
        # :constant keeps T = T_init
        
        # Progress
        if verbose && step % (n_steps ÷ 10) == 0
            defect, _ = compute_defect(A, partition)
            println("Step $step/$n_steps: E=$(round(current_energy, digits=4)), " *
                    "defect=$(round(defect, digits=4)), T=$(round(T, digits=4)), " *
                    "k=$(length(partition)), accept_rate=$(round(accepted/(accepted+rejected), digits=3))")
        end
    end
    
    # Canonicalize best partition
    best_partition = sort([sort(c) for c in best_partition], by=first)
    
    stats = Dict(
        "n_steps" => n_steps,
        "T_init" => T_init,
        "T_final" => T_final,
        "cooling" => string(cooling),
        "accepted" => accepted,
        "rejected" => rejected,
        "acceptance_rate" => accepted / (accepted + rejected),
        "history" => history
    )
    
    return best_partition, best_energy, partition, stats
end

#=============================================================================
  High-Level Interface
=============================================================================#

"""
    find_approximate_eep(A::AbstractMatrix;
                         n_runs::Int=5,
                         n_steps::Int=10000,
                         T_init::Float64=1.0,
                         T_final::Float64=0.01,
                         target_k::Union{Nothing,Int}=nothing,
                         α::Float64=1.0,
                         β::Float64=0.1,
                         γ::Float64=0.0,
                         seed::Union{Nothing,Int}=nothing,
                         verbose::Bool=false)

Find an approximately externally equitable partition using multiple 
Metropolis-Hastings runs.

# Arguments
- `A`: Adjacency matrix
- `n_runs`: Number of independent runs (best result returned)
- `n_steps`: Steps per run
- `T_init`, `T_final`: Temperature range for simulated annealing
- `target_k`: If specified, initializes with this many cells
- `α`: Weight for equitability defect (higher = more equitable)
- `β`: Weight for number of cells (higher = prefer fewer cells)
- `γ`: Weight for imbalance (higher = prefer equal-sized cells)
- `seed`: Random seed
- `verbose`: Print progress

# Energy function
    E(π) = α * Defect + β * NumCells + γ * Imbalance
    
where Imbalance = 1 - EffectiveSize/NumCells ∈ [0, 1)

# Returns
- `best_partition`: Best partition found
- `defect`: Equitability defect of best partition
- `energy`: Total energy
- `all_results`: Results from all runs
"""
function find_approximate_eep(A::AbstractMatrix;
                              n_runs::Int=5,
                              n_steps::Int=10000,
                              T_init::Float64=1.0,
                              T_final::Float64=0.01,
                              target_k::Union{Nothing,Int}=nothing,
                              α::Float64=1.0,
                              β::Float64=0.1,
                              γ::Float64=0.0,
                              seed::Union{Nothing,Int}=nothing,
                              verbose::Bool=false)
    
    if seed !== nothing
        Random.seed!(seed)
    end
    
    n = size(A, 1)
    initial_k = target_k === nothing ? max(2, n ÷ 5) : target_k
    
    all_results = []
    best_partition = nothing
    best_energy = Inf
    
    for run in 1:n_runs
        if verbose
            println("\n" * "="^50)
            println("Run $run/$n_runs")
            println("="^50)
        end
        
        partition, energy, _, stats = run_metropolis(A;
            n_steps=n_steps,
            T_init=T_init,
            T_final=T_final,
            cooling=:exponential,
            α=α, β=β, γ=γ,
            initial_k=initial_k,
            verbose=verbose
        )
        
        defect, _ = compute_defect(A, partition)
        
        push!(all_results, Dict(
            "partition" => partition,
            "energy" => energy,
            "defect" => defect,
            "n_cells" => length(partition),
            "stats" => stats
        ))
        
        if energy < best_energy
            best_energy = energy
            best_partition = partition
        end
        
        if verbose
            println("Run $run result: E=$(round(energy, digits=4)), " *
                    "defect=$(round(defect, digits=4)), k=$(length(partition))")
        end
    end
    
    best_defect, defect_matrix = compute_defect(A, best_partition)
    
    return best_partition, best_defect, best_energy, all_results
end

"""
    is_approximately_equitable(A::AbstractMatrix, partition::Vector{Vector{Int}}; 
                                tolerance::Float64=0.1)

Check if a partition is approximately externally equitable within a tolerance.
"""
function is_approximately_equitable(A::AbstractMatrix, partition::Vector{Vector{Int}}; 
                                     tolerance::Float64=0.1)
    defect, defect_matrix = compute_defect(A, partition)
    n = size(A, 1)
    
    # Normalize defect by graph size
    normalized_defect = defect / (n * length(partition))
    
    return normalized_defect < tolerance, normalized_defect
end

#=============================================================================
  Analysis and Output
=============================================================================#

"""
    effective_size(partition::Vector{Vector{Int}}, n::Int)

Compute effective number of groups.
"""
function effective_size(partition::Vector{Vector{Int}}, n::Int)
    H = 0.0
    for cell in partition
        p = length(cell) / n
        if p > 0
            H -= p * log(p)
        end
    end
    return exp(H)
end

"""
    run_analysis_approximate(A::AbstractMatrix, graph_name::String;
                             output_dir::String="output",
                             n_runs::Int=5,
                             n_steps::Int=10000,
                             α::Float64=1.0,
                             β::Float64=0.1,
                             seed::Union{Nothing,Int}=nothing,
                             verbose::Bool=true)

Run approximate EEP analysis and save results.
"""
function run_analysis_approximate(A::AbstractMatrix, graph_name::String;
                                  output_dir::String="output",
                                  n_runs::Int=5,
                                  n_steps::Int=10000,
                                  α::Float64=1.0,
                                  β::Float64=0.1,
                                  γ::Float64=0.0,
                                  seed::Union{Nothing,Int}=nothing,
                                  verbose::Bool=true)
    
    mkpath(output_dir)
    
    n = size(A, 1)
    m = sum(A) ÷ 2
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    if verbose
        println("="^60)
        println("Approximate EEP Analysis (Metropolis-Hastings)")
        println("="^60)
        println("Graph: $graph_name")
        println("Vertices: $n, Edges: $m")
        println("Parameters: α=$α, β=$β, γ=$γ")
        println("Runs: $n_runs, Steps/run: $n_steps")
        println("-"^60)
    end
    
    start_time = time()
    partition, defect, energy, all_results = find_approximate_eep(A;
        n_runs=n_runs,
        n_steps=n_steps,
        α=α, β=β, γ=γ,
        seed=seed,
        verbose=verbose
    )
    elapsed = time() - start_time
    
    eff_size = effective_size(partition, n)
    n_cells = length(partition)
    cell_sizes = [length(c) for c in partition]
    
    # Check if it's exactly equitable
    is_exact = defect < 1e-10
    
    if verbose
        println("\n" * "-"^60)
        println("BEST RESULT")
        println("-"^60)
        println("Partition: $n_cells cells")
        println("Cell sizes: $cell_sizes")
        println("Defect: $(round(defect, digits=6))")
        println("Effective size: $(round(eff_size, digits=4))")
        println("Is exactly EEP: $is_exact")
        println("Computation time: $(round(elapsed, digits=2)) seconds")
        println("\nCells:")
        for (i, cell) in enumerate(partition)
            println("  Cell $i: $(sort(cell))")
        end
    end
    
    # Save results
    results = Dict(
        "graph_name" => graph_name,
        "timestamp" => timestamp,
        "parameters" => Dict(
            "n_runs" => n_runs,
            "n_steps" => n_steps,
            "alpha" => α,
            "beta" => β,
            "gamma" => γ,
            "seed" => seed === nothing ? "none" : seed
        ),
        "graph_properties" => Dict(
            "n_vertices" => n,
            "n_edges" => m
        ),
        "best_partition" => Dict(
            "cells" => partition,
            "n_cells" => n_cells,
            "cell_sizes" => cell_sizes,
            "defect" => defect,
            "energy" => energy,
            "effective_size" => eff_size,
            "is_exactly_eep" => is_exact
        ),
        "computation_time_seconds" => elapsed
    )
    
    # Save JSON
    results_file = joinpath(output_dir, "$(graph_name)_approx_$(timestamp)_results.json")
    open(results_file, "w") do f
        JSON.print(f, results, 2)
    end
    
    if verbose
        println("\nResults saved to: $results_file")
    end
    
    return results
end

end # module

