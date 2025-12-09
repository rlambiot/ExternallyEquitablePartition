#!/usr/bin/env julia
"""
    run_analysis.jl

Main script to run externally equitable partition analysis on example graphs.
Results are saved to the 'output' directory.

Usage:
    julia run_analysis.jl
"""

# Include the module
include("src/ExternallyEquitablePartition.jl")
using .ExternallyEquitablePartition
using LinearAlgebra

#=============================================================================
  Example Graphs
=============================================================================#

"""
    path_graph(n)

Create adjacency matrix for path graph P_n: 1 - 2 - 3 - ... - n
"""
function path_graph(n::Int)
    A = zeros(Int, n, n)
    for i in 1:n-1
        A[i, i+1] = 1
        A[i+1, i] = 1
    end
    return A
end

"""
    cycle_graph(n)

Create adjacency matrix for cycle graph C_n
"""
function cycle_graph(n::Int)
    A = path_graph(n)
    A[1, n] = 1
    A[n, 1] = 1
    return A
end

"""
    star_graph(n)

Create adjacency matrix for star graph S_n (1 center + n leaves)
"""
function star_graph(n::Int)
    A = zeros(Int, n+1, n+1)
    for i in 2:n+1
        A[1, i] = 1
        A[i, 1] = 1
    end
    return A
end

"""
    complete_graph(n)

Create adjacency matrix for complete graph K_n
"""
function complete_graph(n::Int)
    A = ones(Int, n, n) - I
    return A
end

"""
    complete_bipartite_graph(m, n)

Create adjacency matrix for complete bipartite graph K_{m,n}
"""
function complete_bipartite_graph(m::Int, n::Int)
    total = m + n
    A = zeros(Int, total, total)
    for i in 1:m
        for j in m+1:total
            A[i, j] = 1
            A[j, i] = 1
        end
    end
    return A
end

"""
    petersen_graph()

Create adjacency matrix for the Petersen graph (10 vertices, 15 edges)
"""
function petersen_graph()
    A = [
        0 1 0 0 1 1 0 0 0 0;
        1 0 1 0 0 0 1 0 0 0;
        0 1 0 1 0 0 0 1 0 0;
        0 0 1 0 1 0 0 0 1 0;
        1 0 0 1 0 0 0 0 0 1;
        1 0 0 0 0 0 0 1 1 0;
        0 1 0 0 0 0 0 0 1 1;
        0 0 1 0 0 1 0 0 0 1;
        0 0 0 1 0 1 1 0 0 0;
        0 0 0 0 1 0 1 1 0 0
    ]
    return A
end

"""
    barbell_graph(n)

Create adjacency matrix for barbell graph: two complete graphs K_n connected by a bridge
"""
function barbell_graph(n::Int)
    total = 2n
    A = zeros(Int, total, total)
    
    # First clique
    for i in 1:n
        for j in i+1:n
            A[i, j] = 1
            A[j, i] = 1
        end
    end
    
    # Second clique
    for i in n+1:total
        for j in i+1:total
            A[i, j] = 1
            A[j, i] = 1
        end
    end
    
    # Bridge
    A[n, n+1] = 1
    A[n+1, n] = 1
    
    return A
end

"""
    wheel_graph(n)

Create adjacency matrix for wheel graph W_n (center + cycle of n vertices)
"""
function wheel_graph(n::Int)
    total = n + 1
    A = zeros(Int, total, total)
    
    # Center (vertex 1) connected to all others
    for i in 2:total
        A[1, i] = 1
        A[i, 1] = 1
    end
    
    # Outer cycle
    for i in 2:total-1
        A[i, i+1] = 1
        A[i+1, i] = 1
    end
    A[2, total] = 1
    A[total, 2] = 1
    
    return A
end

#=============================================================================
  Main Analysis
=============================================================================#

function main()
    println("╔══════════════════════════════════════════════════════════════╗")
    println("║     Externally Equitable Partition Analysis                  ║")
    println("║     Finding Minimum Effective Size Partitions                ║")
    println("╚══════════════════════════════════════════════════════════════╝")
    println()
    
    # Analysis parameters
    n_random = 100
    seed = 42
    output_dir = "output"
    
    # List of graphs to analyze
    graphs = [
        ("path_P6", path_graph(6)),
        ("cycle_C6", cycle_graph(6)),
        ("star_S5", star_graph(5)),
        ("complete_K5", complete_graph(5)),
        ("bipartite_K3_3", complete_bipartite_graph(3, 3)),
        ("petersen", petersen_graph()),
        ("barbell_3", barbell_graph(3)),
        ("wheel_W6", wheel_graph(6)),
    ]
    
    # Store all results for comparison
    all_results = Dict{String, Any}()
    
    for (name, A) in graphs
        println("\n" * "─"^60)
        results = run_analysis(A, name;
                              output_dir=output_dir,
                              n_random=n_random,
                              seed=seed,
                              verbose=true)
        all_results[name] = results
        println()
    end
    
    # Print summary comparison
    println("\n" * "═"^60)
    println("SUMMARY COMPARISON")
    println("═"^60)
    println()
    println("Graph               │ Vertices │ Cells │ Eff. Size │ Equitable?")
    println("────────────────────┼──────────┼───────┼───────────┼───────────")
    
    for (name, _) in graphs
        r = all_results[name]
        n = r["graph_properties"]["n_vertices"]
        k = r["best_partition"]["n_cells"]
        eff = round(r["best_partition"]["effective_size"], digits=3)
        eq = r["best_partition"]["is_also_equitable"] ? "Yes" : "No"
        
        # Format with padding
        name_pad = rpad(name, 19)
        n_pad = lpad(string(n), 8)
        k_pad = lpad(string(k), 5)
        eff_pad = lpad(string(eff), 9)
        eq_pad = lpad(eq, 10)
        
        println("$name_pad │$n_pad │$k_pad │$eff_pad │$eq_pad")
    end
    
    println()
    println("All results saved to: $output_dir/")
    println()
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

