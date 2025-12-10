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
using Graphs
using GraphMakie
using CairoMakie
using Colors
using NetworkLayout

# Load our modules
include(joinpath(PROJECT_ROOT, "src", "ApproximateEEP.jl"))
include(joinpath(PROJECT_ROOT, "src", "ExternallyEquitablePartition.jl"))
using .ApproximateEEP
using .ExternallyEquitablePartition

# =============================================================================
#   Zachary's Karate Club Network
# =============================================================================

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

# =============================================================================
#   Visualization Functions
# =============================================================================

"""
    get_color_palette(n)

Generate a visually distinct color palette for n groups.
"""
function get_color_palette(n::Int)
    if n <= 10
        # Use a carefully chosen palette for small n
        base_colors = [
            colorant"#E63946",  # Red
            colorant"#457B9D",  # Steel blue
            colorant"#2A9D8F",  # Teal
            colorant"#E9C46A",  # Yellow
            colorant"#F4A261",  # Orange
            colorant"#9B5DE5",  # Purple
            colorant"#00F5D4",  # Cyan
            colorant"#F15BB5",  # Pink
            colorant"#6A994E",  # Green
            colorant"#BC6C25",  # Brown
        ]
        return base_colors[1:n]
    else
        # Use HSL color wheel for larger n
        return [HSL(h, 0.7, 0.5) for h in range(0, 360, length=n+1)[1:end-1]]
    end
end

"""
    plot_partition(A, partition, title, filename; layout_seed=42)

Create a visualization of a graph with nodes colored by partition.
"""
function plot_partition(A::Matrix, partition::Vector{Vector{Int}}, 
                        title::String, filename::String;
                        layout_seed::Int=42,
                        subtitle::String="")
    
    n = size(A, 1)
    
    # Create Graphs.jl graph
    g = SimpleGraph(n)
    for i in 1:n
        for j in i+1:n
            if A[i,j] > 0
                add_edge!(g, i, j)
            end
        end
    end
    
    # Assign colors to nodes based on partition
    colors = get_color_palette(length(partition))
    node_colors = Vector{RGB{Float64}}(undef, n)
    
    for (cell_idx, cell) in enumerate(partition)
        for v in cell
            node_colors[v] = convert(RGB{Float64}, colors[cell_idx])
        end
    end
    
    # Create figure
    fig = Figure(size=(900, 800), backgroundcolor=:white)
    
    # Title
    Label(fig[1, 1], title, fontsize=18, font=:bold, tellwidth=false)
    
    # Subtitle
    if !isempty(subtitle)
        Label(fig[2, 1], subtitle, fontsize=14, color=:gray40, tellwidth=false)
    end
    
    # Graph plot
    ax = Axis(fig[3, 1], 
              aspect=DataAspect(),
              backgroundcolor=:white)
    hidedecorations!(ax)
    hidespines!(ax)
    
    # Use spring layout with fixed seed for reproducibility
    Random.seed!(layout_seed)
    layout = Spring(; iterations=100, C=2.0)
    
    # Plot the graph
    graphplot!(ax, g,
               node_color=node_colors,
               node_size=25,
               edge_color=(:gray, 0.4),
               edge_width=1.5,
               nlabels=string.(1:n),
               nlabels_fontsize=9,
               nlabels_color=:white,
               layout=layout)
    
    # Legend - show cells with their sizes
    if length(partition) <= 15
        legend_elements = [MarkerElement(color=colors[i], marker=:circle, markersize=15) 
                          for i in 1:length(partition)]
        legend_labels = ["Cell $i (n=$(length(partition[i]))): $(length(partition[i]) <= 6 ? sort(partition[i]) : "$(length(partition[i])) vertices")" 
                        for i in 1:length(partition)]
        
        Legend(fig[4, 1], legend_elements, legend_labels, 
               "Partition Cells",
               orientation=:horizontal,
               tellwidth=false,
               nbanks=min(4, length(partition)))
    else
        Label(fig[4, 1], "$(length(partition)) cells (too many to display individually)", 
              fontsize=12, color=:gray50, tellwidth=false)
    end
    
    # Save
    save(filename, fig, px_per_unit=2)
    println("Saved: $filename")
    
    return fig
end

"""
    plot_comparison(A, partitions, labels, filename; layout_seed=42)

Create a side-by-side comparison of multiple partitions.
"""
function plot_comparison(A::Matrix, partitions::Vector, labels::Vector{String}, 
                         defects::Vector{Float64}, filename::String;
                         layout_seed::Int=42)
    
    n = size(A, 1)
    n_partitions = length(partitions)
    
    # Create Graphs.jl graph
    g = SimpleGraph(n)
    for i in 1:n
        for j in i+1:n
            if A[i,j] > 0
                add_edge!(g, i, j)
            end
        end
    end
    
    # Compute layout once (shared across all plots)
    Random.seed!(layout_seed)
    layout = Spring(; iterations=100, C=2.0)
    pos = layout(g)
    
    # Create figure
    fig = Figure(size=(400 * n_partitions, 500), backgroundcolor=:white)
    
    # Main title
    Label(fig[1, 1:n_partitions], "Karate Club: Comparison of Partitions", 
          fontsize=20, font=:bold, tellwidth=false)
    
    for (idx, (partition, label, defect)) in enumerate(zip(partitions, labels, defects))
        # Assign colors
        colors = get_color_palette(length(partition))
        node_colors = Vector{RGB{Float64}}(undef, n)
        
        for (cell_idx, cell) in enumerate(partition)
            for v in cell
                node_colors[v] = convert(RGB{Float64}, colors[cell_idx])
            end
        end
        
        # Axis
        ax = Axis(fig[2, idx],
                  title="$label\n($(length(partition)) cells, defect=$(round(defect, digits=4)))",
                  aspect=DataAspect(),
                  backgroundcolor=:white)
        hidedecorations!(ax)
        hidespines!(ax)
        
        # Plot
        graphplot!(ax, g,
                   node_color=node_colors,
                   node_size=20,
                   edge_color=(:gray, 0.3),
                   edge_width=1.0,
                   nlabels=string.(1:n),
                   nlabels_fontsize=8,
                   nlabels_color=:white,
                   layout=_->pos)  # Use precomputed positions
    end
    
    save(filename, fig, px_per_unit=2)
    println("Saved: $filename")
    
    return fig
end

# =============================================================================
#   Main Analysis
# =============================================================================

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
    
    # Store results for visualization
    all_partitions = []
    all_labels = String[]
    all_defects = Float64[]
    
    # =========================================================================
    # Experiment 1: Pure defect minimization (find most equitable partition)
    # =========================================================================
    α1, β1, γ1 = 1.0, 0.0, 0.0
    
    println("="^70)
    println("EXPERIMENT 1: Minimize defect only (α=$α1, β=$β1, γ=$γ1)")
    println("="^70)
    println()
    
    partition1, defect1, energy1, _ = find_approximate_eep(A;
        n_runs=5,
        n_steps=20000,
        α=α1,
        β=β1,
        γ=γ1,
        T_init=1.0,
        T_final=0.001,
        seed=42,
        verbose=true
    )
    
    println("\nResult:")
    println("  Cells: $(length(partition1))")
    println("  Defect: $(round(defect1, digits=6))")
    println("  Effective size: $(round(ApproximateEEP.effective_size(partition1, n), digits=4))")
    
    # =========================================================================
    # Experiment 2: Balance defect and effective size
    # =========================================================================
    α2, β2, γ2 = 1.0, 0.5, 0.0
    
    println("\n" * "="^70)
    println("EXPERIMENT 2: Balance defect and coarseness (α=$α2, β=$β2, γ=$γ2)")
    println("="^70)
    println()
    
    partition2, defect2, energy2, _ = find_approximate_eep(A;
        n_runs=5,
        n_steps=20000,
        α=α2,
        β=β2,
        γ=γ2,
        T_init=2.0,
        T_final=0.001,
        seed=42,
        verbose=true
    )
    
    println("\nResult:")
    println("  Cells: $(length(partition2))")
    println("  Defect: $(round(defect2, digits=6))")
    println("  Effective size: $(round(ApproximateEEP.effective_size(partition2, n), digits=4))")
    
    label2 = "α=$α2, β=$β2, γ=$γ2"
    push!(all_partitions, partition2)
    push!(all_labels, label2)
    push!(all_defects, defect2)
    
    # =========================================================================
    # Experiment 3: Target ~3-4 cells with moderate defect tolerance
    # =========================================================================
    α3, β3, γ3 = 1.0, 0.2, 0.1
    
    println("\n" * "="^70)
    println("EXPERIMENT 3: Target moderate coarseness (α=$α3, β=$β3, γ=$γ3)")
    println("="^70)
    println()
    
    partition3, defect3, energy3, _ = find_approximate_eep(A;
        n_runs=10,
        n_steps=30000,
        α=α3,
        β=β3,
        γ=γ3,
        target_k=4,
        T_init=3.0,
        T_final=0.001,
        seed=123,
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
    
    label3 = "α=$α3, β=$β3, γ=$γ3"
    push!(all_partitions, partition3)
    push!(all_labels, label3)
    push!(all_defects, defect3)
    
    # =========================================================================
    # Experiment 4: Find partition with ~5 cells
    # =========================================================================
    α4, β4, γ4 = 1.0, 0.1, 0.05
    
    println("\n" * "="^70)
    println("EXPERIMENT 4: Target ~5 cells (α=$α4, β=$β4, γ=$γ4)")
    println("="^70)
    println()
    
    partition4, defect4, energy4, _ = find_approximate_eep(A;
        n_runs=10,
        n_steps=30000,
        α=α4,
        β=β4,
        γ=γ4,
        target_k=5,
        T_init=3.0,
        T_final=0.001,
        seed=456,
        verbose=true
    )
    
    println("\nResult:")
    println("  Cells: $(length(partition4))")
    println("  Defect: $(round(defect4, digits=6))")
    println("  Effective size: $(round(ApproximateEEP.effective_size(partition4, n), digits=4))")
    println("  Partition:")
    for (i, cell) in enumerate(partition4)
        println("    Cell $i: $(sort(cell))")
    end
    
    label4 = "α=$α4, β=$β4, γ=$γ4"
    push!(all_partitions, partition4)
    push!(all_labels, label4)
    push!(all_defects, defect4)
    
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
    
    push!(all_partitions, exact_partition)
    push!(all_labels, "Exact EEP")
    push!(all_defects, exact_defect)
    
    # =========================================================================
    # Summary comparison
    # =========================================================================
    println("\n" * "="^70)
    println("SUMMARY COMPARISON")
    println("="^70)
    println()
    println("Method                         | Cells | Defect    | Eff. Size")
    println("-"^75)
    println("Exp1 (α=$α1, β=$β1, γ=$γ1)       | $(lpad(length(partition1), 5)) | $(lpad(round(defect1, digits=4), 9)) | $(lpad(round(ApproximateEEP.effective_size(partition1, n), digits=3), 9))")
    println("Exp2 (α=$α2, β=$β2, γ=$γ2)       | $(lpad(length(partition2), 5)) | $(lpad(round(defect2, digits=4), 9)) | $(lpad(round(ApproximateEEP.effective_size(partition2, n), digits=3), 9))")
    println("Exp3 (α=$α3, β=$β3, γ=$γ3)     | $(lpad(length(partition3), 5)) | $(lpad(round(defect3, digits=4), 9)) | $(lpad(round(ApproximateEEP.effective_size(partition3, n), digits=3), 9))")
    println("Exp4 (α=$α4, β=$β4, γ=$γ4)    | $(lpad(length(partition4), 5)) | $(lpad(round(defect4, digits=4), 9)) | $(lpad(round(ApproximateEEP.effective_size(partition4, n), digits=3), 9))")
    println("Exact EEP                      | $(lpad(length(exact_partition), 5)) | $(lpad(round(exact_defect, digits=4), 9)) | $(lpad(round(exact_eff, digits=3), 9))")
    println()
    
    # =========================================================================
    # Generate visualizations
    # =========================================================================
    println("="^70)
    println("GENERATING VISUALIZATIONS")
    println("="^70)
    println()
    
    # Individual partition plots
    plot_partition(A, partition2, 
                   "Approximate EEP ($label2)",
                   joinpath(output_dir, "karate_approx_balanced.png");
                   subtitle="Cells: $(length(partition2)), Defect: $(round(defect2, digits=4)), Eff.Size: $(round(ApproximateEEP.effective_size(partition2, n), digits=2))")
    
    plot_partition(A, partition3, 
                   "Approximate EEP ($label3)",
                   joinpath(output_dir, "karate_approx_moderate.png");
                   subtitle="Cells: $(length(partition3)), Defect: $(round(defect3, digits=4)), Eff.Size: $(round(ApproximateEEP.effective_size(partition3, n), digits=2))")
    
    plot_partition(A, partition4, 
                   "Approximate EEP ($label4)",
                   joinpath(output_dir, "karate_approx_5cell.png");
                   subtitle="Cells: $(length(partition4)), Defect: $(round(defect4, digits=4)), Eff.Size: $(round(ApproximateEEP.effective_size(partition4, n), digits=2))")
    
    plot_partition(A, exact_partition, 
                   "Exact EEP (Minimum Effective Size)",
                   joinpath(output_dir, "karate_exact_eep.png");
                   subtitle="Cells: $(length(exact_partition)), Defect: 0.0, Eff.Size: $(round(exact_eff, digits=2))")
    
    # Comparison plot
    plot_comparison(A, all_partitions, all_labels, all_defects,
                    joinpath(output_dir, "karate_partition_comparison.png"))
    
    println()
    println("="^70)
    println("Analysis complete!")
    println("="^70)
    println()
    println("Output files:")
    println("  - karate_approx_balanced.png")
    println("  - karate_approx_moderate.png")
    println("  - karate_approx_5cell.png")
    println("  - karate_exact_eep.png")
    println("  - karate_partition_comparison.png")
end

# Run
main()
