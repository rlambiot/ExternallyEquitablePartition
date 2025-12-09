# Externally Equitable Partitions

A Julia package for finding **externally equitable partitions** (EEPs) of graphs, with a focus on identifying the partition with **minimum effective size** (entropy-based measure).

## What is an Externally Equitable Partition?

### Definition

An **externally equitable partition** of a graph G = (V, E) is a partition of the vertex set V into cells C₁, C₂, ..., Cₖ such that:

> For all i ≠ j, every vertex in cell Cᵢ has the **same number of neighbors** in cell Cⱼ.

This means the number of edges from any vertex v ∈ Cᵢ to cell Cⱼ depends only on i and j, not on the specific choice of v.

### Comparison with Equitable Partitions

| Property | Equitable Partition | Externally Equitable Partition |
|----------|--------------------|---------------------------------|
| Constraint | ∀ i, j (including i = j) | Only for i ≠ j |
| Internal structure | Constrained | Unconstrained |
| Strictness | More restrictive | Less restrictive |

**Key difference**: Externally equitable partitions relax the constraint on *internal* neighbors (within the same cell), allowing vertices in the same cell to have different degrees within that cell.

### The Quotient Matrix

For an externally equitable partition, we can define a **quotient matrix** Q where:

```
Q[i,j] = number of neighbors that any vertex in Cᵢ has in Cⱼ
```

For i ≠ j, this is well-defined by the EEP property. The quotient matrix captures the "coarse-grained" structure of the graph.

## Effective Size (Entropy-Based Measure)

### Definition

Given a partition with cells of sizes n₁, n₂, ..., nₖ (total n vertices), define:

- **Cell probabilities**: pᵢ = nᵢ / n
- **Shannon entropy**: H = -Σᵢ pᵢ ln(pᵢ)
- **Effective size**: N_eff = exp(H)

The effective size represents the "equivalent number of equal-sized groups" and satisfies:

- **Minimum** (N_eff = 1): Universal partition (all vertices in one cell)
- **Maximum** (N_eff = n): Discrete partition (each vertex alone)

### Why Minimize Effective Size?

Finding the EEP with **minimum effective size** (excluding the trivial universal partition) identifies the most parsimonious non-trivial grouping that respects the external connectivity structure. This is useful for:

1. **Coarse-graining** dynamical systems on networks
2. Finding **structural equivalence classes**
3. **Dimensionality reduction** while preserving external interactions
4. Network **controllability analysis**

## The Algorithm

### Overview

The algorithm uses **iterative partition refinement** to find externally equitable partitions:

1. Start with an initial partition
2. Compute the "external signature" of each vertex (neighbor counts in other cells)
3. Split cells where vertices have different external signatures
4. Repeat until no more splits occur

The result is the **coarsest EEP that refines the initial partition**.

### Multi-Start Search for Minimum Effective Size

Since different initial partitions can lead to different EEPs, we use a multi-start heuristic:

1. **Generate diverse candidates**: degree-based, neighborhood-based, random k-partitions
2. **Refine each** to an EEP
3. **Track all unique EEPs** found
4. **Return the one with minimum effective size** (excluding trivial)

### Pseudocode

```
function find_min_entropy_eep(A):
    candidates = generate_candidate_partitions(A)
    best = nil
    best_eff_size = ∞
    
    for initial in candidates:
        eep = refine_to_eep(A, initial)
        
        if |eep| > 1:  # exclude trivial
            eff = effective_size(eep)
            if eff < best_eff_size:
                best = eep
                best_eff_size = eff
    
    return best

function refine_to_eep(A, partition):
    repeat:
        for each cell C in partition:
            signatures = {}
            for v in C:
                sig = [count neighbors of v in C' for C' ≠ C]
                signatures[sig].add(v)
            
            if |signatures| > 1:
                split C into sub-cells by signature
    until no changes
    
    return partition
```

## Installation

### Requirements

- Julia 1.6 or later
- Packages: `JSON` (for saving results)

### Setup

```julia
# From the project directory
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Or install dependencies manually:

```julia
using Pkg
Pkg.add("JSON")
```

## Usage

### Basic Usage

```julia
include("src/ExternallyEquitablePartition.jl")
using .ExternallyEquitablePartition

# Define adjacency matrix
A = [0 1 0 0 0 0;
     1 0 1 0 0 0;
     0 1 0 1 0 0;
     0 0 1 0 1 0;
     0 0 0 1 0 1;
     0 0 0 0 1 0]

# Find minimum entropy EEP
partition, eff_size, Q, all_eeps = find_min_entropy_eep(A; verbose=true)

println("Best partition: ", partition)
println("Effective size: ", eff_size)
println("Quotient matrix:")
display(Q)
```

### Run Full Analysis with Output

```julia
include("src/ExternallyEquitablePartition.jl")
using .ExternallyEquitablePartition

A = [0 1 1 1 1;
     1 0 0 0 0;
     1 0 0 0 0;
     1 0 0 0 0;
     1 0 0 0 0]

results = run_analysis(A, "star_graph_S4";
                       output_dir="output",
                       n_random=100,
                       seed=42,
                       verbose=true)
```

This saves:
- `output/<name>_<timestamp>_params.json` - Analysis parameters
- `output/<name>_<timestamp>_results.json` - Full results in JSON
- `output/<name>_<timestamp>_summary.txt` - Human-readable summary

### Running the Demo

```bash
julia run_analysis.jl
```

## Output Files

All results are saved to the `output/` directory:

| File | Contents |
|------|----------|
| `*_params.json` | Input parameters and graph properties |
| `*_results.json` | Complete results including partition, quotient matrix, all EEPs found |
| `*_summary.txt` | Human-readable summary report |

### Example Output Structure

```
output/
├── petersen_2024-01-15_14-30-22_params.json
├── petersen_2024-01-15_14-30-22_results.json
└── petersen_2024-01-15_14-30-22_summary.txt
```

## API Reference

### Main Functions

| Function | Description |
|----------|-------------|
| `find_min_entropy_eep(A)` | Find EEP with minimum effective size |
| `find_externally_equitable_partition(A)` | Find EEP from default/custom initial partition |
| `refine_to_eep(A, partition)` | Refine partition to coarsest EEP |
| `run_analysis(A, name)` | Full analysis with file output |

### Utility Functions

| Function | Description |
|----------|-------------|
| `is_externally_equitable(A, partition)` | Check if partition is EEP |
| `is_equitable(A, partition)` | Check if partition is fully equitable |
| `effective_size(partition, n)` | Compute effective size (exp of entropy) |
| `entropy_partition(partition, n)` | Compute Shannon entropy |
| `quotient_matrix(A, partition)` | Compute quotient matrix |

## Examples

### Path Graph P₆

```
1 — 2 — 3 — 4 — 5 — 6
```

Minimum entropy EEP: `[{1, 6}, {2, 5}, {3, 4}]`
- Effective size ≈ 3.0 (three equal-sized cells)

### Star Graph S₄

```
    2
    |
5 — 1 — 3
    |
    4
```

Minimum entropy EEP: `[{1}, {2, 3, 4, 5}]`
- Effective size ≈ 1.72 (center vs leaves)

### Cycle Graph C₆

Multiple EEPs exist, including bipartitions and antipodal groupings.

## Mathematical Background

### Relation to Graph Automorphisms

The **orbit partition** under graph automorphisms is always equitable (and hence externally equitable). However, externally equitable partitions can be coarser than the orbit partition.

### Relation to Weisfeiler-Lehman

The partition refinement algorithm is related to the **1-dimensional Weisfeiler-Lehman** algorithm for graph isomorphism testing, but focuses on external neighbor counts only.

### Applications

1. **Network dynamics**: EEPs identify vertex groups with identical external influence patterns
2. **Model reduction**: Quotient graphs preserve certain spectral and dynamical properties
3. **Symmetry detection**: Coarse EEPs reveal approximate symmetries
4. **Controllability**: Related to structural controllability of networked systems

## References

- Godsil, C., & Royle, G. (2001). *Algebraic Graph Theory*. Springer.
- Weisfeiler, B., & Leman, A. (1968). *The reduction of a graph to canonical form*.
- O'Clery, N., Yuan, Y., Stan, G. B., & Barahona, M. (2013). Observability and coarse-graining of consensus dynamics through the external equitable partition. *Physical Review E*.

## License

MIT License

