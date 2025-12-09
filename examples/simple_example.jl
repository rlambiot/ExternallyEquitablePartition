#!/usr/bin/env julia
"""
    simple_example.jl

A simple example demonstrating how to use the ExternallyEquitablePartition module.
"""

# Navigate to project root if needed
cd(@__DIR__)
cd("..")

# Include the module
include("src/ExternallyEquitablePartition.jl")
using .ExternallyEquitablePartition

println("="^50)
println("Simple Example: Finding Externally Equitable Partitions")
println("="^50)

#=============================================================================
  Example 1: Define your own graph
=============================================================================#

println("\n--- Example 1: Custom Graph ---\n")

# Define a simple graph as an adjacency matrix
# This is a "bowtie" graph: two triangles sharing a vertex
#
#     2       4
#    / \     / \
#   1---3---5---6
#
A = [0 1 1 0 0 0;   # vertex 1
     1 0 1 0 0 0;   # vertex 2
     1 1 0 1 0 0;   # vertex 3 (central vertex)
     0 0 1 0 1 1;   # vertex 4
     0 0 0 1 0 1;   # vertex 5
     0 0 0 1 1 0]   # vertex 6

println("Adjacency matrix:")
display(A)
println()

# Find the minimum entropy EEP
partition, eff_size, Q, all_eeps = find_min_entropy_eep(A; verbose=false)

println("Minimum effective size EEP:")
for (i, cell) in enumerate(partition)
    println("  Cell $i: vertices $cell")
end
println()
println("Effective size: $(round(eff_size, digits=4))")
println("Number of cells: $(length(partition))")
println()
println("Quotient matrix Q (entry Q[i,j] = neighbors of v∈Cᵢ in Cⱼ):")
display(Q)
println()

# Check properties
println("Is externally equitable? ", is_externally_equitable(A, partition))
println("Is also fully equitable? ", is_equitable(A, partition))

#=============================================================================
  Example 2: Using a specific initial partition
=============================================================================#

println("\n" * "-"^50)
println("\n--- Example 2: Starting from a Specific Partition ---\n")

# Start with a custom initial partition
initial = [[1, 2, 3], [4, 5, 6]]  # Left triangle vs right triangle
println("Initial partition: $initial")

partition2, Q2 = find_externally_equitable_partition(A; initial_partition=initial)

println("Refined EEP:")
for (i, cell) in enumerate(partition2)
    println("  Cell $i: vertices $cell")
end
println()
println("Effective size: $(round(effective_size(partition2, 6), digits=4))")

#=============================================================================
  Example 3: Run full analysis with output files
=============================================================================#

println("\n" * "-"^50)
println("\n--- Example 3: Full Analysis with File Output ---\n")

# Run complete analysis and save to output folder
results = run_analysis(A, "bowtie_graph";
                       output_dir="output",
                       n_random=50,
                       seed=123,
                       verbose=true)

println("\nResults dictionary keys: ", keys(results))

#=============================================================================
  Example 4: Comparing all found EEPs
=============================================================================#

println("\n" * "-"^50)
println("\n--- Example 4: All Unique EEPs Found ---\n")

# The search finds multiple EEPs - let's examine them
partition, eff_size, Q, all_eeps = find_min_entropy_eep(A; n_random=100, seed=42)

println("Found $(length(all_eeps)) unique non-trivial EEPs:")
println()

# Sort by effective size
sorted_eeps = sort(collect(all_eeps), by=x->x[2])

for (i, (eep, eff)) in enumerate(sorted_eeps)
    println("EEP $i (eff_size = $(round(eff, digits=4))):")
    for (j, cell) in enumerate(eep)
        println("    Cell $j: $cell")
    end
    println()
end

println("="^50)
println("Examples complete! Check the 'output' folder for saved results.")
println("="^50)

