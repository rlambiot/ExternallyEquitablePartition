#!/usr/bin/env julia
"""
    karate_club_approximate.jl

Find approximately externally equitable partitions of the Karate Club network
using Metropolis-Hastings Monte Carlo optimization.

This approach allows finding partitions that are "almost" EEP, which is useful
for real networks that may not have exact non-trivial EEPs.
"""

# Get the project root directory
const PROJECT_ROOT = dirname(@__DIR__)

# Load packages
using Random
using Statistics

# Load our modules
include(joinpath(PROJECT_ROOT, "src", "ApproximateEEP.jl"))
include(joinpath(PROJECT_ROOT, "src", "ExternallyEquitablePartition.jl"))
using .ApproximateEEP
using .ExternallyEquitablePartition

#=============================================================================
  Zachary's Karate Club Network
=============================================================================#

function karate_club_graph()
    edges = [
        (1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8), (1,9), (1,11), (1,12), (1,13), (1,14), (1,18), (1,20), (1,22), (1,32),
        (2,3), (2,4), (2,8), (2,14), (2,18), (2,20), (2,22), (2,31),
        (3,4), (3,8), (3,9), (3,10), (3,14), (3,28), (3,29), (3,33),
        (4,8), (4,13), (4,14),
        (5,7), (5,11),
        (6,7), (6,11), (6,17),
        (7,17),
        (9,31), (9,33), (9,34),
        (10,34),
        (14,34),
        (15,33), (15,34),
        (16,33), (16,34),
        (19,33), (19,34),
        (20,34),
        (21,33), (21,34),
        (23,33), (23,34),
        (24,26), (24,28), (24,30), (24,33), (24,34),
        (25,26), (25,28), (25,32),
        (26,32),
        (27,30), (27,34),
        (28,34),
        (29,32), (29,34),
        (30,33), (30,34),
        (31,33), (31,34),
        (32,33), (32,34),
        (33,34)
    ]
    
    n = 34
    A = zeros(Int, n, n)
    for (i, j) in edges
        A[i, j] = 1
        A[j, i] = 1
    end
    
    return A
end

#=============================================================================
  Main Analysis
=============================================================================#

function main()
    println("╔══════════════════════════════════════════════════════════════════╗")
    println("║  Karate Club - Approximate EEP via Metropolis-Hastings           ║")
    println("╚══════════════════════════════════════════════════════════════════╝")
    println()
    
    A = karate_club_graph()
    n = size(A, 1)
    
    # Create output directory
    output_dir = joinpath(PROJECT_ROOT, "output")
    mkpath(output_dir)
    
    println("Graph: Zachary's Karate Club (n=$n)")
    println()
    
    # =========================================================================
    # Experiment 1: Pure defect minimization (find most equitable partition)
    # =========================================================================
    println("="^70)
    println("EXPERIMENT 1: Minimize defect only (α=1, β=0)")
    println("="^70)
    println()
    
    partition1, defect1, energy1, _ = find_approximate_eep(A;
        n_runs=5,
        n_steps=20000,
        α=1.0,      # Defect weight
        β=0.0,      # No effective size penalty
        γ=0.0,      # No cell count penalty
        T_init=1.0,
        T_final=0.001,
        seed=42,
        verbose=true
    )
    
    println("\nResult:")
    println("  Cells: $(length(partition1))")
    println("  Defect: $(round(defect1, digits=6))")
    println("  Effective size: $(round(ApproximateEEP.effective_size(partition1, n), digits=4))")
    println("  Partition:")
    for (i, cell) in enumerate(partition1)
        println("    Cell $i: $(sort(cell))")
    end
    
    # =========================================================================
    # Experiment 2: Balance defect and effective size
    # =========================================================================
    println("\n" * "="^70)
    println("EXPERIMENT 2: Balance defect and coarseness (α=1, β=0.5)")
    println("="^70)
    println()
    
    partition2, defect2, energy2, _ = find_approximate_eep(A;
        n_runs=5,
        n_steps=20000,
        α=1.0,      # Defect weight
        β=0.5,      # Effective size penalty (prefer coarser partitions)
        γ=0.0,
        T_init=2.0,
        T_final=0.001,
        seed=42,
        verbose=true
    )
    
    println("\nResult:")
    println("  Cells: $(length(partition2))")
    println("  Defect: $(round(defect2, digits=6))")
    println("  Effective size: $(round(ApproximateEEP.effective_size(partition2, n), digits=4))")
    println("  Partition:")
    for (i, cell) in enumerate(partition2)
        println("    Cell $i: $(sort(cell))")
    end
    
    # =========================================================================
    # Experiment 3: Target specific number of cells
    # =========================================================================
    println("\n" * "="^70)
    println("EXPERIMENT 3: Target ~4 cells (α=1, β=0, γ=0.5)")
    println("="^70)
    println()
    
    partition3, defect3, energy3, _ = find_approximate_eep(A;
        n_runs=5,
        n_steps=20000,
        α=1.0,
        β=0.0,
        γ=0.5,      # Penalize many cells
        target_k=4,
        T_init=2.0,
        T_final=0.001,
        seed=42,
        verbose=true
    )
    
    println("\nResult:")
    println("  Cells: $(length(partition3))")
    println("  Defect: $(round(defect3, digits=6))")
    println("  Effective size: $(round(ApproximateEEP.effective_size(partition3, n), digits=4))")
    println("  Partition:")
    for (i, cell) in enumerate(partition3)
        println("    Cell $i: $(sort(cell))")
    end
    
    # =========================================================================
    # Compare with exact EEP search
    # =========================================================================
    println("\n" * "="^70)
    println("COMPARISON: Exact EEP search")
    println("="^70)
    println()
    
    exact_partition, exact_eff, _, _ = find_min_entropy_eep(A; n_random=200, seed=42)
    exact_defect, _ = compute_defect(A, exact_partition)
    
    println("Exact EEP (minimum effective size):")
    println("  Cells: $(length(exact_partition))")
    println("  Defect: $(round(exact_defect, digits=6))")
    println("  Effective size: $(round(exact_eff, digits=4))")
    
    # =========================================================================
    # Summary comparison
    # =========================================================================
    println("\n" * "="^70)
    println("SUMMARY COMPARISON")
    println("="^70)
    println()
    println("Method                    | Cells | Defect    | Eff. Size")
    println("-"^70)
    println("Approx (defect only)      | $(lpad(length(partition1), 5)) | $(lpad(round(defect1, digits=4), 9)) | $(lpad(round(ApproximateEEP.effective_size(partition1, n), digits=3), 9))")
    println("Approx (balanced)         | $(lpad(length(partition2), 5)) | $(lpad(round(defect2, digits=4), 9)) | $(lpad(round(ApproximateEEP.effective_size(partition2, n), digits=3), 9))")
    println("Approx (target k≈4)       | $(lpad(length(partition3), 5)) | $(lpad(round(defect3, digits=4), 9)) | $(lpad(round(ApproximateEEP.effective_size(partition3, n), digits=3), 9))")
    println("Exact EEP                 | $(lpad(length(exact_partition), 5)) | $(lpad(round(exact_defect, digits=4), 9)) | $(lpad(round(exact_eff, digits=3), 9))")
    println()
    
    println("="^70)
    println("Analysis complete!")
    println("="^70)
end

# Run
main()

