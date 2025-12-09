#!/usr/bin/env julia
"""
    karate_club_analysis.jl

Analyze the Zachary Karate Club network for externally equitable partitions
and visualize the three best (minimum effective size) partitions.
"""

# Get the project root directory
const PROJECT_ROOT = dirname(@__DIR__)

# Load packages
using Graphs
using GraphMakie
using CairoMakie
using Colors
using Random
using NetworkLayout

# Load our module
include(joinpath(PROJECT_ROOT, "src", "ExternallyEquitablePartition.jl"))
using .ExternallyEquitablePartition

#=============================================================================
  Zachary's Karate Club Network
=============================================================================#

"""
    karate_club_graph()

Create the famous Zachary Karate Club network (34 nodes, 78 edges).
This network represents friendships between members of a karate club
that eventually split into two groups.
"""
function karate_club_graph()
    # Edge list for Zachary's Karate Club
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
  Visualization Functions
=============================================================================#

"""
    get_color_palette(n)

Generate a visually distinct color palette for n groups.
"""
function get_color_palette(n::Int)
    if n <= 8
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
                        layout_seed::Int=42)
    
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
    fig = Figure(size=(800, 700), backgroundcolor=:white)
    
    # Title
    Label(fig[1, 1], title, fontsize=20, font=:bold, tellwidth=false)
    
    # Graph plot
    ax = Axis(fig[2, 1], 
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
               edge_color=(:gray, 0.5),
               edge_width=1.5,
               nlabels=string.(1:n),
               nlabels_fontsize=10,
               nlabels_color=:white,
               layout=layout)
    
    # Legend
    legend_elements = [MarkerElement(color=colors[i], marker=:circle, markersize=15) 
                       for i in 1:length(partition)]
    legend_labels = ["Cell $i: $(sort(partition[i]))" for i in 1:length(partition)]
    
    Legend(fig[3, 1], legend_elements, legend_labels, 
           "Partition Cells",
           orientation=:horizontal,
           tellwidth=false,
           nbanks=min(3, length(partition)))
    
    # Save
    save(filename, fig, px_per_unit=2)
    println("Saved: $filename")
    
    return fig
end

#=============================================================================
  Main Analysis
=============================================================================#

function main()
    println("╔══════════════════════════════════════════════════════════════╗")
    println("║  Karate Club - Externally Equitable Partition Analysis       ║")
    println("╚══════════════════════════════════════════════════════════════╝")
    println()
    
    # Create output directory
    output_dir = joinpath(PROJECT_ROOT, "output")
    mkpath(output_dir)
    
    # Get Karate Club adjacency matrix
    A = karate_club_graph()
    n = size(A, 1)
    m = sum(A) ÷ 2
    
    println("Graph: Zachary's Karate Club")
    println("Vertices: $n")
    println("Edges: $m")
    println()
    
    # Run analysis with more random samples to find diverse partitions
    println("Searching for externally equitable partitions...")
    println("-"^60)
    
    partition, eff_size, Q, all_eeps = find_min_entropy_eep(A; 
                                                             n_random=200, 
                                                             verbose=false,
                                                             seed=42)
    
    println("Found $(length(all_eeps)) unique non-trivial EEPs")
    println()
    
    # Sort all EEPs by effective size
    sorted_eeps = sort(collect(all_eeps), by=x->x[2])
    
    # Get the three best (minimum effective size)
    n_best = min(3, length(sorted_eeps))
    best_eeps = sorted_eeps[1:n_best]
    
    println("="^60)
    println("TOP $n_best PARTITIONS (Minimum Effective Size)")
    println("="^60)
    
    for (rank, (eep, eff)) in enumerate(best_eeps)
        println()
        println("--- Rank $rank ---")
        println("Effective size: $(round(eff, digits=4))")
        println("Number of cells: $(length(eep))")
        println("Cell sizes: $([length(c) for c in eep])")
        
        # Check if also fully equitable
        is_eq = is_equitable(A, eep)
        println("Is also fully equitable: $is_eq")
        
        println("Cells:")
        for (i, cell) in enumerate(eep)
            println("  Cell $i: $(sort(cell))")
        end
        
        # Compute quotient matrix
        Q_eep = quotient_matrix(A, eep)
        println("Quotient matrix:")
        for row in eachrow(Q_eep)
            println("  $(collect(row))")
        end
    end
    
    # Generate visualizations
    println()
    println("="^60)
    println("GENERATING VISUALIZATIONS")
    println("="^60)
    println()
    
    for (rank, (eep, eff)) in enumerate(best_eeps)
        title = "Karate Club - Rank $rank EEP\n" *
                "Effective size: $(round(eff, digits=3)), " *
                "Cells: $(length(eep)), " *
                "Sizes: $([length(c) for c in eep])"
        
        filename = joinpath(output_dir, "karate_club_partition_rank$(rank).png")
        
        plot_partition(A, eep, title, filename; layout_seed=123)
    end
    
    # Also save full analysis results
    println()
    results = run_analysis(A, "karate_club";
                          output_dir=output_dir,
                          n_random=200,
                          seed=42,
                          verbose=false)
    
    # Save summary of all effective sizes found
    eff_sizes_file = joinpath(output_dir, "karate_club_all_effective_sizes.txt")
    open(eff_sizes_file, "w") do f
        println(f, "All Effective Sizes Found for Karate Club Network")
        println(f, "="^50)
        println(f, "Total unique EEPs: $(length(sorted_eeps))")
        println(f, "")
        println(f, "Rank | Eff.Size | #Cells | Cell Sizes")
        println(f, "-"^50)
        for (i, (eep, eff)) in enumerate(sorted_eeps)
            sizes = [length(c) for c in eep]
            println(f, "$i | $(round(eff, digits=4)) | $(length(eep)) | $sizes")
        end
    end
    println("Saved: $eff_sizes_file")
    
    println()
    println("="^60)
    println("ANALYSIS COMPLETE")
    println("="^60)
    println()
    println("Output files in '$output_dir/':")
    println("  - karate_club_partition_rank1.png")
    println("  - karate_club_partition_rank2.png")
    println("  - karate_club_partition_rank3.png")
    println("  - karate_club_*_params.json")
    println("  - karate_club_*_results.json")
    println("  - karate_club_*_summary.txt")
    println("  - karate_club_all_effective_sizes.txt")
    println()
end

# Run
main()

