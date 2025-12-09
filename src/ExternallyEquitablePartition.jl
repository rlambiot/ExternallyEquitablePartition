"""
    ExternallyEquitablePartition.jl

A Julia module for finding externally equitable partitions (EEPs) of graphs,
with focus on finding the partition with minimum effective size (entropy).

An externally equitable partition groups vertices such that each vertex in cell Cᵢ
has the same number of neighbors in cell Cⱼ for all i ≠ j.
"""
module ExternallyEquitablePartition

using Random
using LinearAlgebra
using Dates
using JSON

export find_min_entropy_eep,
       find_externally_equitable_partition,
       refine_to_eep,
       effective_size,
       entropy_partition,
       is_externally_equitable,
       is_equitable,
       quotient_matrix,
       run_analysis,
       save_results

#=============================================================================
  Core Partition Functions
=============================================================================#

"""
    entropy_partition(partition::Vector{Vector{Int}}, n::Int)

Compute Shannon entropy of a partition: H = -Σ pᵢ log(pᵢ)

# Arguments
- `partition`: Vector of cells (each cell is a vector of vertex indices)
- `n`: Total number of vertices

# Returns
- Shannon entropy (in nats)
"""
function entropy_partition(partition::Vector{Vector{Int}}, n::Int)
    H = 0.0
    for cell in partition
        p = length(cell) / n
        if p > 0
            H -= p * log(p)
        end
    end
    return H
end

"""
    effective_size(partition::Vector{Vector{Int}}, n::Int)

Compute effective number of groups (perplexity): exp(H)

This is the Hill number of order 1, representing the "equivalent number
of equal-sized groups" that would give the same entropy.

# Arguments
- `partition`: Vector of cells
- `n`: Total number of vertices

# Returns
- Effective size ∈ [1, k] where k is the number of cells
"""
function effective_size(partition::Vector{Vector{Int}}, n::Int)
    return exp(entropy_partition(partition, n))
end

"""
    is_externally_equitable(A::AbstractMatrix, partition::Vector{Vector{Int}})

Check if a given partition is externally equitable.

A partition is externally equitable if every vertex in cell Cᵢ has the same
number of neighbors in cell Cⱼ, for all pairs i ≠ j.

# Arguments
- `A`: Adjacency matrix (n × n)
- `partition`: Partition to check

# Returns
- `true` if externally equitable, `false` otherwise
"""
function is_externally_equitable(A::AbstractMatrix, partition::Vector{Vector{Int}})
    for (i, cell_i) in enumerate(partition)
        for (j, cell_j) in enumerate(partition)
            i == j && continue
            isempty(cell_i) && continue
            counts = [sum(A[v, u] for u in cell_j; init=0) for v in cell_i]
            if !allequal(counts)
                return false
            end
        end
    end
    return true
end

"""
    is_equitable(A::AbstractMatrix, partition::Vector{Vector{Int}})

Check if a given partition is (fully) equitable.

A partition is equitable if every vertex in cell Cᵢ has the same number
of neighbors in cell Cⱼ, for ALL pairs (i, j) including i = j.

This is more restrictive than externally equitable.
"""
function is_equitable(A::AbstractMatrix, partition::Vector{Vector{Int}})
    for (i, cell_i) in enumerate(partition)
        for (j, cell_j) in enumerate(partition)
            isempty(cell_i) && continue
            counts = [sum(A[v, u] for u in cell_j; init=0) for v in cell_i]
            if !allequal(counts)
                return false
            end
        end
    end
    return true
end

# Helper function
allequal(x) = isempty(x) || all(==(first(x)), x)

"""
    quotient_matrix(A::AbstractMatrix, partition::Vector{Vector{Int}})

Compute the quotient matrix for a given partition.

Entry Q[i,j] represents the number of neighbors a vertex in cell i has in cell j.
For externally equitable partitions, this is well-defined for i ≠ j.

# Arguments
- `A`: Adjacency matrix
- `partition`: Partition of vertices

# Returns
- Quotient matrix Q (k × k where k is number of cells)
"""
function quotient_matrix(A::AbstractMatrix, partition::Vector{Vector{Int}})
    k = length(partition)
    Q = zeros(Int, k, k)
    
    for i in 1:k
        isempty(partition[i]) && continue
        v = partition[i][1]  # Representative vertex
        for j in 1:k
            Q[i, j] = sum(A[v, u] for u in partition[j]; init=0)
        end
    end
    
    return Q
end

#=============================================================================
  Partition Refinement Algorithm
=============================================================================#

"""
    refine_to_eep(A::AbstractMatrix, partition::Vector{Vector{Int}})

Refine a partition to the coarsest externally equitable partition that refines it.

Uses iterative partition refinement based on external neighbor signatures.

# Arguments
- `A`: Adjacency matrix
- `partition`: Initial partition to refine

# Returns
- Coarsest EEP that refines the input partition
"""
function refine_to_eep(A::AbstractMatrix, partition::Vector{Vector{Int}})
    partition = [copy(cell) for cell in partition if !isempty(cell)]
    
    changed = true
    while changed
        changed = false
        new_partition = Vector{Int}[]
        
        for cell in partition
            if length(cell) <= 1
                push!(new_partition, cell)
                continue
            end
            
            # Group vertices by their external neighbor signature
            signatures = Dict{Vector{Int}, Vector{Int}}()
            for v in cell
                sig = Int[]
                for other_cell in partition
                    if other_cell !== cell  # External cells only
                        push!(sig, sum(A[v, u] for u in other_cell; init=0))
                    end
                end
                if haskey(signatures, sig)
                    push!(signatures[sig], v)
                else
                    signatures[sig] = [v]
                end
            end
            
            # Split if multiple signatures exist
            if length(signatures) > 1
                changed = true
                for vertices in values(signatures)
                    push!(new_partition, vertices)
                end
            else
                push!(new_partition, cell)
            end
        end
        
        partition = new_partition
    end
    
    return [sort(cell) for cell in partition]
end

"""
    find_externally_equitable_partition(A::AbstractMatrix; 
                                         initial_partition=nothing)

Find an externally equitable partition by refining from an initial partition.

# Arguments
- `A`: Adjacency matrix
- `initial_partition`: Starting partition (default: all vertices in one cell)

# Returns
- `partition`: The resulting EEP
- `Q`: Quotient matrix
"""
function find_externally_equitable_partition(A::AbstractMatrix; 
                                              initial_partition::Union{Nothing, Vector{Vector{Int}}}=nothing)
    n = size(A, 1)
    if initial_partition === nothing
        initial_partition = [collect(1:n)]
    end
    
    partition = refine_to_eep(A, initial_partition)
    
    # Canonicalize
    partition = sort([sort(cell) for cell in partition], by=first)
    
    Q = quotient_matrix(A, partition)
    
    return partition, Q
end

#=============================================================================
  Candidate Generation for Minimum Entropy Search
=============================================================================#

"""
    degree_partition(A::AbstractMatrix)

Create initial partition based on vertex degrees.
"""
function degree_partition(A::AbstractMatrix)
    n = size(A, 1)
    degrees = vec(sum(A, dims=2))
    
    deg_groups = Dict{Int, Vector{Int}}()
    for v in 1:n
        d = Int(degrees[v])
        push!(get!(deg_groups, d, Int[]), v)
    end
    
    return collect(values(deg_groups))
end

"""
    generate_random_partition(n::Int, k::Int)

Generate a random partition of 1:n into k non-empty cells.
"""
function generate_random_partition(n::Int, k::Int)
    k = min(k, n)
    assignments = rand(1:k, n)
    
    # Ensure all groups are non-empty
    for i in 1:k
        if !(i in assignments)
            assignments[rand(1:n)] = i
        end
    end
    
    partition = [Int[] for _ in 1:k]
    for (v, g) in enumerate(assignments)
        push!(partition[g], v)
    end
    
    return filter(!isempty, partition)
end

"""
    generate_binary_partitions(n::Int; max_samples::Int=20)

Generate 2-cell partitions (most likely to give coarse EEPs).
"""
function generate_binary_partitions(n::Int; max_samples::Int=20)
    partitions = Vector{Vector{Int}}[]
    
    for size1 in 1:(n÷2)
        n_possible = binomial(n, size1)
        n_samples = min(max_samples, n_possible)
        
        for _ in 1:n_samples
            cell1 = sort(randperm(n)[1:size1])
            cell2 = sort(setdiff(1:n, cell1))
            push!(partitions, [cell1, cell2])
        end
    end
    
    return partitions
end

"""
    generate_candidate_partitions(A::AbstractMatrix; n_random::Int=50)

Generate diverse candidate initial partitions for minimum entropy search.
"""
function generate_candidate_partitions(A::AbstractMatrix; n_random::Int=50)
    n = size(A, 1)
    candidates = Vector{Vector{Int}}[]
    
    # 1. Degree-based partition
    push!(candidates, degree_partition(A))
    
    # 2. Binary partitions (k=2, most likely to stay coarse)
    append!(candidates, generate_binary_partitions(n; max_samples=n_random÷2))
    
    # 3. Small k partitions (k = 3, 4)
    for k in 3:min(4, n-1)
        for _ in 1:(n_random ÷ 4)
            push!(candidates, generate_random_partition(n, k))
        end
    end
    
    # 4. Neighborhood-based: {v ∪ N(v)} vs rest
    for v in 1:n
        neighbors = findall(x -> x > 0, A[v, :])
        cell1 = sort(unique([v; neighbors]))
        cell2 = sort(setdiff(1:n, cell1))
        if !isempty(cell2) && length(cell1) < n
            push!(candidates, [cell1, cell2])
        end
    end
    
    # 5. Distance-based: vertices within distance 2
    A2 = A * A
    for v in 1:n
        close = findall(i -> A[v,i] > 0 || A2[v,i] > 0 || i == v, 1:n)
        far = setdiff(1:n, close)
        if !isempty(far) && length(close) > 1
            push!(candidates, [sort(close), sort(far)])
        end
    end
    
    # 6. Triangle-based partition
    triangles = zeros(Int, n)
    for v in 1:n
        triangles[v] = sum(A[v, u] * A[u, w] * A[w, v] for u in 1:n for w in u+1:n)
    end
    tri_groups = Dict{Int, Vector{Int}}()
    for v in 1:n
        push!(get!(tri_groups, triangles[v], Int[]), v)
    end
    if length(tri_groups) > 1
        push!(candidates, collect(values(tri_groups)))
    end
    
    return candidates
end

#=============================================================================
  Main Search Function
=============================================================================#

"""
    find_min_entropy_eep(A::AbstractMatrix; 
                         n_random::Int=50, 
                         verbose::Bool=false,
                         seed::Union{Nothing,Int}=nothing)

Find the externally equitable partition with MINIMUM effective size
(excluding the trivial universal partition).

# Arguments
- `A`: Adjacency matrix
- `n_random`: Number of random candidates per strategy
- `verbose`: Print progress information
- `seed`: Random seed for reproducibility

# Returns
- `best_partition`: The EEP with minimum effective size (> 1 cell)
- `best_eff_size`: The effective size (exp of entropy)
- `quotient`: The quotient matrix
- `all_eeps`: Dict of all unique EEPs found with their effective sizes
"""
function find_min_entropy_eep(A::AbstractMatrix; 
                              n_random::Int=50, 
                              verbose::Bool=false,
                              seed::Union{Nothing,Int}=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end
    
    n = size(A, 1)
    candidates = generate_candidate_partitions(A; n_random=n_random)
    
    # Track unique EEPs
    seen = Set{Vector{Vector{Int}}}()
    all_eeps = Dict{Vector{Vector{Int}}, Float64}()
    
    best_partition = nothing
    best_eff_size = Inf
    
    for initial in candidates
        # Refine to EEP
        eep = refine_to_eep(A, initial)
        
        # Skip trivial universal partition
        if length(eep) <= 1
            continue
        end
        
        # Canonicalize
        canonical = sort([sort(cell) for cell in eep], by=first)
        
        if canonical in seen
            continue
        end
        push!(seen, canonical)
        
        # Compute effective size
        eff = effective_size(canonical, n)
        all_eeps[canonical] = eff
        
        if verbose
            println("Found EEP: $(length(canonical)) cells, eff_size = $(round(eff, digits=4))")
        end
        
        if eff < best_eff_size
            best_eff_size = eff
            best_partition = canonical
        end
    end
    
    if best_partition === nothing
        @warn "No non-trivial EEP found. Returning universal partition."
        return [collect(1:n)], 1.0, reshape([sum(A[1,:])], 1, 1), all_eeps
    end
    
    Q = quotient_matrix(A, best_partition)
    
    return best_partition, best_eff_size, Q, all_eeps
end

#=============================================================================
  Analysis and Output Functions
=============================================================================#

"""
    run_analysis(A::AbstractMatrix, graph_name::String;
                 output_dir::String="output",
                 n_random::Int=50,
                 seed::Union{Nothing,Int}=nothing,
                 verbose::Bool=true)

Run complete analysis on a graph and save results to output directory.

# Arguments
- `A`: Adjacency matrix
- `graph_name`: Name identifier for the graph
- `output_dir`: Directory to save results
- `n_random`: Number of random candidates
- `seed`: Random seed for reproducibility
- `verbose`: Print progress

# Returns
- Results dictionary with all computed information
"""
function run_analysis(A::AbstractMatrix, graph_name::String;
                      output_dir::String="output",
                      n_random::Int=50,
                      seed::Union{Nothing,Int}=nothing,
                      verbose::Bool=true)
    
    # Ensure output directory exists
    mkpath(output_dir)
    
    n = size(A, 1)
    m = sum(A) ÷ 2  # Number of edges (assuming undirected)
    
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    if verbose
        println("="^60)
        println("Externally Equitable Partition Analysis")
        println("="^60)
        println("Graph: $graph_name")
        println("Vertices: $n, Edges: $m")
        println("Timestamp: $timestamp")
        println("-"^60)
    end
    
    # Run the search
    start_time = time()
    partition, eff_size, Q, all_eeps = find_min_entropy_eep(A; 
                                                             n_random=n_random, 
                                                             verbose=verbose,
                                                             seed=seed)
    elapsed = time() - start_time
    
    # Compute additional metrics
    entropy = entropy_partition(partition, n)
    n_cells = length(partition)
    cell_sizes = [length(cell) for cell in partition]
    is_also_equitable = is_equitable(A, partition)
    
    if verbose
        println("-"^60)
        println("Results:")
        println("  Best partition has $n_cells cells")
        println("  Cell sizes: $cell_sizes")
        println("  Effective size: $(round(eff_size, digits=4))")
        println("  Entropy: $(round(entropy, digits=4))")
        println("  Is also fully equitable: $is_also_equitable")
        println("  Total unique EEPs found: $(length(all_eeps))")
        println("  Computation time: $(round(elapsed, digits=2)) seconds")
        println("-"^60)
        println("Partition:")
        for (i, cell) in enumerate(partition)
            println("  Cell $i: $cell")
        end
        println("Quotient matrix:")
        display(Q)
        println()
    end
    
    # Prepare results dictionary
    results = Dict(
        "graph_name" => graph_name,
        "timestamp" => timestamp,
        "parameters" => Dict(
            "n_random" => n_random,
            "seed" => seed === nothing ? "none" : seed
        ),
        "graph_properties" => Dict(
            "n_vertices" => n,
            "n_edges" => m,
            "density" => 2m / (n * (n-1))
        ),
        "best_partition" => Dict(
            "cells" => partition,
            "n_cells" => n_cells,
            "cell_sizes" => cell_sizes,
            "effective_size" => eff_size,
            "entropy" => entropy,
            "is_also_equitable" => is_also_equitable
        ),
        "quotient_matrix" => Q,
        "all_eeps_found" => length(all_eeps),
        "all_effective_sizes" => sort(collect(values(all_eeps))),
        "computation_time_seconds" => elapsed
    )
    
    # Save results
    save_results(results, output_dir, graph_name, timestamp)
    
    return results
end

"""
    save_results(results::Dict, output_dir::String, graph_name::String, timestamp::String)

Save analysis results to files in the output directory.
"""
function save_results(results::Dict, output_dir::String, graph_name::String, timestamp::String)
    base_name = "$(graph_name)_$(timestamp)"
    
    # Save parameters as JSON
    params_file = joinpath(output_dir, "$(base_name)_params.json")
    open(params_file, "w") do f
        JSON.print(f, Dict(
            "graph_name" => results["graph_name"],
            "timestamp" => results["timestamp"],
            "parameters" => results["parameters"],
            "graph_properties" => results["graph_properties"]
        ), 2)
    end
    
    # Save full results as JSON
    results_file = joinpath(output_dir, "$(base_name)_results.json")
    
    # Convert matrices and arrays for JSON serialization
    json_results = deepcopy(results)
    json_results["quotient_matrix"] = [collect(row) for row in eachrow(results["quotient_matrix"])]
    
    open(results_file, "w") do f
        JSON.print(f, json_results, 2)
    end
    
    # Save human-readable summary
    summary_file = joinpath(output_dir, "$(base_name)_summary.txt")
    open(summary_file, "w") do f
        println(f, "="^60)
        println(f, "Externally Equitable Partition Analysis Results")
        println(f, "="^60)
        println(f, "")
        println(f, "Graph: $(results["graph_name"])")
        println(f, "Timestamp: $(results["timestamp"])")
        println(f, "")
        println(f, "PARAMETERS")
        println(f, "-"^40)
        println(f, "  n_random: $(results["parameters"]["n_random"])")
        println(f, "  seed: $(results["parameters"]["seed"])")
        println(f, "")
        println(f, "GRAPH PROPERTIES")
        println(f, "-"^40)
        println(f, "  Vertices: $(results["graph_properties"]["n_vertices"])")
        println(f, "  Edges: $(results["graph_properties"]["n_edges"])")
        println(f, "  Density: $(round(results["graph_properties"]["density"], digits=4))")
        println(f, "")
        println(f, "BEST PARTITION (Minimum Effective Size)")
        println(f, "-"^40)
        println(f, "  Number of cells: $(results["best_partition"]["n_cells"])")
        println(f, "  Cell sizes: $(results["best_partition"]["cell_sizes"])")
        println(f, "  Effective size: $(round(results["best_partition"]["effective_size"], digits=6))")
        println(f, "  Entropy: $(round(results["best_partition"]["entropy"], digits=6))")
        println(f, "  Is also fully equitable: $(results["best_partition"]["is_also_equitable"])")
        println(f, "")
        println(f, "  Cells:")
        for (i, cell) in enumerate(results["best_partition"]["cells"])
            println(f, "    Cell $i: $cell")
        end
        println(f, "")
        println(f, "QUOTIENT MATRIX")
        println(f, "-"^40)
        Q = results["quotient_matrix"]
        for row in eachrow(Q)
            println(f, "  $(collect(row))")
        end
        println(f, "")
        println(f, "SEARCH STATISTICS")
        println(f, "-"^40)
        println(f, "  Total unique EEPs found: $(results["all_eeps_found"])")
        println(f, "  All effective sizes: $(round.(results["all_effective_sizes"], digits=4))")
        println(f, "  Computation time: $(round(results["computation_time_seconds"], digits=2)) seconds")
        println(f, "")
        println(f, "="^60)
    end
    
    println("Results saved to:")
    println("  - $params_file")
    println("  - $results_file")
    println("  - $summary_file")
end

end # module

