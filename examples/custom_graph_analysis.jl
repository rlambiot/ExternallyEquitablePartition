#!/usr/bin/env julia
"""
    custom_graph_analysis.jl

Template for analyzing your own graph.
Modify the adjacency matrix A below to analyze any graph.
"""

# Navigate to project root
cd(@__DIR__)
cd("..")

include("src/ExternallyEquitablePartition.jl")
using .ExternallyEquitablePartition

#=============================================================================
  CONFIGURE YOUR ANALYSIS HERE
=============================================================================#

# Graph name (used for output files)
GRAPH_NAME = "my_custom_graph"

# Define your adjacency matrix here
# Example: a simple 5-vertex graph
#
#   1 --- 2
#   |     |
#   3 --- 4 --- 5
#
A = [
    0  1  1  0  0;   # vertex 1: connected to 2, 3
    1  0  0  1  0;   # vertex 2: connected to 1, 4
    1  0  0  1  0;   # vertex 3: connected to 1, 4
    0  1  1  0  1;   # vertex 4: connected to 2, 3, 5
    0  0  0  1  0    # vertex 5: connected to 4
]

# Analysis parameters
N_RANDOM = 100      # Number of random initial partitions to try
SEED = nothing      # Random seed (nothing for random, or an integer for reproducibility)
OUTPUT_DIR = "output"

#=============================================================================
  RUN ANALYSIS (no need to modify below)
=============================================================================#

function main()
    println("╔════════════════════════════════════════════════════════════╗")
    println("║  Custom Graph Analysis                                      ║")
    println("╚════════════════════════════════════════════════════════════╝")
    println()
    
    # Validate adjacency matrix
    n = size(A, 1)
    @assert size(A, 2) == n "Adjacency matrix must be square"
    @assert A == A' "Adjacency matrix must be symmetric (undirected graph)"
    @assert all(diag(A) .== 0) "Diagonal must be zero (no self-loops)"
    
    println("Graph: $GRAPH_NAME")
    println("Vertices: $n")
    println("Edges: $(sum(A) ÷ 2)")
    println()
    
    # Run analysis
    results = run_analysis(A, GRAPH_NAME;
                          output_dir=OUTPUT_DIR,
                          n_random=N_RANDOM,
                          seed=SEED,
                          verbose=true)
    
    println("\n" * "═"^60)
    println("Analysis complete!")
    println("Results saved to: $OUTPUT_DIR/")
    println("═"^60)
    
    return results
end

# Run
results = main()

