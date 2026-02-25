# cell_list.jl
# Cell-linked list for alchemical GCMC in Julia.
# Adaptations from Allen & Tildesley link_list_module.f90:
#   - Julia 1-based indexing throughout
#   - Fractional (alchemical) particle lives OUTSIDE the list at all times
#   - Variable N ∈ [0, N_max] is safe; N=0 is a no-op in all functions
#   - Small boxes (sc < 3, e.g. 4-atom test cases): sc is clamped to 1,
#     which puts every atom in a single cell — correct but O(N²), same as
#     naive. No special code paths needed anywhere else.
#   - Pre-allocated neighbour scratch buffer avoids per-call heap allocation

using StaticArrays

# ════════════════════════════════════════════════════════════════════════════
# Struct
# ════════════════════════════════════════════════════════════════════════════

mutable struct CellList
    sc       ::Int           # cells per dim (≥1); sc=1 → single-cell / naive fallback
    N_max    ::Int           # maximum normal-particle count
    head     ::Array{Int,3}  # head[cx,cy,cz], 1-based, shape (sc,sc,sc)
    list     ::Vector{Int}   # list[i] → next atom in same cell, 0 = end-of-chain
    cell     ::Matrix{Int}   # cell[dim,i] = 1-based cell index of atom i, shape (3,N_max)
    nbr_buf  ::Vector{Int}   # pre-allocated scratch for neighbours!()
    nbr_n    ::Int           # number of valid entries currently in nbr_buf
end

# ════════════════════════════════════════════════════════════════════════════
# Constructor
# ════════════════════════════════════════════════════════════════════════════

"""
    CellList(N_max, r_cut_box)

Allocate a CellList for up to `N_max` normal particles.
`r_cut_box = r_cut / L` (box units, so r_cut_box = r_cut_σ / L_σ).

If the box is too small for ≥ 3 cells per dimension (e.g. a tiny 4-atom test
system), `sc` is clamped to 1 so that everything lands in one cell and the
code degenerates to O(N²) without any change to the calling code.

Initialise with:
    cl = CellList(sim.N_max + 1, sim.r_cut_box)
    make_list!(cl, μstate.N, μstate.r_box)
"""
function CellList(N_max::Int, r_cut_box::Float64)
    sc_raw = floor(Int, 1.0 / r_cut_box)
    if sc_raw < 3
        @warn "Box too small for ≥3 cells/dim (sc=$sc_raw). Using sc=1 (naive O(N²) fallback). Fine for small test systems."
        sc = 1
    else
        sc = sc_raw
        println("CellList: sc=$sc cells per dimension  (N_max=$N_max)")
    end

    head    = zeros(Int, sc, sc, sc)
    list    = zeros(Int, N_max)
    cell    = zeros(Int, 3, N_max)
    nbr_buf = zeros(Int, N_max)   # worst-case: all N_max atoms are neighbours

    return CellList(sc, N_max, head, list, cell, nbr_buf, 0)
end

# ════════════════════════════════════════════════════════════════════════════
# Cell-index mapping
# ════════════════════════════════════════════════════════════════════════════

"""
    c_index(cl, ri_box) → SVector{3,Int}

Map a position in box units (components in −0.5 … 0.5) to a 1-based cell triple.
"""
@inline function c_index(cl::CellList, ri::AbstractVector)
    @assert all(abs.(ri) .<= 0.5 + 1e-10) "Atom outside main box: $ri"
    c1 = clamp(floor(Int, (ri[1] + 0.5) * cl.sc) + 1, 1, cl.sc)
    c2 = clamp(floor(Int, (ri[2] + 0.5) * cl.sc) + 1, 1, cl.sc)
    c3 = clamp(floor(Int, (ri[3] + 0.5) * cl.sc) + 1, 1, cl.sc)
    return SVector{3,Int}(c1, c2, c3)
end

# ════════════════════════════════════════════════════════════════════════════
# Low-level list operations
# ════════════════════════════════════════════════════════════════════════════

"""
    create_in_list!(cl, i, ci)

Insert atom `i` as the new head of cell `ci`.
"""
@inline function create_in_list!(cl::CellList, i::Int, ci::SVector{3,Int})
    cl.list[i]                    = cl.head[ci[1], ci[2], ci[3]]
    cl.head[ci[1], ci[2], ci[3]] = i
    cl.cell[1, i] = ci[1]
    cl.cell[2, i] = ci[2]
    cl.cell[3, i] = ci[3]
end

"""
    destroy_in_list!(cl, i, ci)

Remove atom `i` from cell `ci`.
"""
function destroy_in_list!(cl::CellList, i::Int, ci::SVector{3,Int})
    this = cl.head[ci[1], ci[2], ci[3]]

    if this == i
        cl.head[ci[1], ci[2], ci[3]] = cl.list[i]
        return
    end

    while true
        nxt = cl.list[this]
        if nxt == i
            cl.list[this] = cl.list[i]   # link over i
            return
        elseif nxt == 0
            error("destroy_in_list!: atom $i not found in cell $ci")
        else
            this = nxt
        end
    end
end

"""
    move_in_list!(cl, i, ci_new)

Move atom `i` from its current cell to `ci_new`. No-op if cell unchanged.
"""
@inline function move_in_list!(cl::CellList, i::Int, ci_new::SVector{3,Int})
    ci_old = SVector{3,Int}(cl.cell[1,i], cl.cell[2,i], cl.cell[3,i])
    ci_old == ci_new && return
    destroy_in_list!(cl, i, ci_old)
    create_in_list!(cl, i, ci_new)
end

# ════════════════════════════════════════════════════════════════════════════
# Build / rebuild
# ════════════════════════════════════════════════════════════════════════════

"""
    make_list!(cl, N, r_box)

Rebuild the entire cell list from scratch for `N` normal particles.
Safe for N = 0 (just clears the head array).
"""
function make_list!(cl::CellList, N::Int, r_box::Vector{<:AbstractVector})
    fill!(cl.head, 0)
    @inbounds for i in 1:N
        ci = c_index(cl, r_box[i])
        create_in_list!(cl, i, ci)
    end
end

# ════════════════════════════════════════════════════════════════════════════
# Neighbour shell — writes into pre-allocated buffer, zero heap allocation
# ════════════════════════════════════════════════════════════════════════════

# 27 displacement vectors for the 3×3×3 neighbourhood, listed with inversion
# symmetry: columns 1:13 = lower half-shell, column 14 = (0,0,0) self cell,
# columns 15:27 = upper half-shell.
const _NK = 13
const _D  = let
    rows = [
        (-1,-1,-1), ( 0,-1,-1), ( 1,-1,-1),
        (-1, 1,-1), ( 0, 1,-1), ( 1, 1,-1),
        (-1, 0,-1), ( 1, 0,-1), ( 0, 0,-1),
        ( 0,-1, 0), ( 1,-1, 0), (-1,-1, 0),
        (-1, 0, 0), ( 0, 0, 0), ( 1, 0, 0),   # col 14 = self cell (0,0,0)
        ( 1, 1, 0), (-1, 1, 0), ( 0, 1, 0),
        ( 0, 0, 1), (-1, 0, 1), ( 1, 0, 1),
        (-1,-1, 1), ( 0,-1, 1), ( 1,-1, 1),
        (-1, 1, 1), ( 0, 1, 1), ( 1, 1, 1),
    ]
    m = Matrix{Int}(undef, 3, 2*_NK+1)
    for (col, (dx,dy,dz)) in enumerate(rows)
        m[1,col] = dx; m[2,col] = dy; m[3,col] = dz
    end
    m
end

"""
    neighbours!(cl, i, ci; half=false)

Fill `cl.nbr_buf[1:cl.nbr_n]` with neighbour indices of a probe at cell `ci`,
skipping index `i` (pass `i=0` for the fractional particle).

- `half=false` (default): search all 27 cells. Use for trial moves and the
  fractional particle — the probe position need not be registered in the list.
- `half=true`: Newton half-shell (13 cells + down-list in own cell). Only
  valid when atom `i` is currently registered at `ci`. Use for total-energy sums.

Zero heap allocation — results are in `cl.nbr_buf[1:cl.nbr_n]`.
"""
function neighbours!(cl::CellList, i::Int, ci::SVector{3,Int}; half::Bool=false)
    sc = cl.sc

    if half
        (cl.cell[1,i]==ci[1] && cl.cell[2,i]==ci[2] && cl.cell[3,i]==ci[3]) ||
            error("neighbours!: cell mismatch for atom $i stored=$(cl.cell[:,i]) given=$ci")
        col_start = _NK + 2   # upper half-shell starts at column 15
        col_end   = 2*_NK + 1
    else
        col_start = 1
        col_end   = 2*_NK + 1
    end

    n = 0

    @inbounds for col in col_start:col_end
        cj1 = mod1(ci[1] + _D[1,col], sc)
        cj2 = mod1(ci[2] + _D[2,col], sc)
        cj3 = mod1(ci[3] + _D[3,col], sc)

        # Half-shell self-cell: start from list[i] to skip already-counted pairs
        j = (col == _NK + 1 && half) ? cl.list[i] : cl.head[cj1, cj2, cj3]

        while j != 0
            if j != i
                n += 1
                cl.nbr_buf[n] = j
            end
            j = cl.list[j]
        end
    end

    cl.nbr_n = n
    return nothing
end

# ════════════════════════════════════════════════════════════════════════════
# Consistency check (debugging only)
# ════════════════════════════════════════════════════════════════════════════

function check_list(cl::CellList, N::Int, r_box::Vector{<:AbstractVector})
    for i in 1:N
        ci = c_index(cl, r_box[i])
        (cl.cell[1,i]==ci[1] && cl.cell[2,i]==ci[2] && cl.cell[3,i]==ci[3]) ||
            error("check_list: atom $i wrong cell — expected $ci got $(cl.cell[:,i])")
    end
    for cx in 1:cl.sc, cy in 1:cl.sc, cz in 1:cl.sc
        j = cl.head[cx, cy, cz]
        while j != 0
            (cl.cell[1,j]==cx && cl.cell[2,j]==cy && cl.cell[3,j]==cz) ||
                error("check_list: atom $j stored in wrong cell")
            j = cl.list[j]
        end
    end
    println("check_list: OK (N=$N, sc=$(cl.sc))")
end

############## ORIGINAL

using StaticArrays

# ────────────────────────────────────────────────────────────────────────────
# Mutable container that owns the cell-list state.
# One instance lives for the whole simulation; pass it everywhere.
# ────────────────────────────────────────────────────────────────────────────
mutable struct CellList
    sc    ::Int                        # number of cells along each dimension
    head  ::Array{Int,3}               # head[cx,cy,cz]  (1-based, sc×sc×sc)
    list  ::Vector{Int}                # list[i] = next atom after i in same cell; 0 = end
    cell  ::Matrix{Int}                # cell[dim, i] = 1-based cell index of atom i  (3 × N_max)
end

# ────────────────────────────────────────────────────────────────────────────
# Construction / initialisation
# ────────────────────────────────────────────────────────────────────────────

"""
    CellList(N_max, r_cut_box)

Allocate a CellList for up to `N_max` normal particles.
`r_cut_box = r_cut / L` (box units).
The fractional particle is NOT stored in the list.
"""
function CellList(N_max::Int, r_cut_box::Float64)
    sc = floor(Int, 1.0 / r_cut_box)
    sc >= 3 || error("Box too small for cell list: sc = $sc (need ≥ 3)")
    println("CellList: sc = $sc subcells per dimension")

    head = zeros(Int, sc, sc, sc)
    list = zeros(Int, N_max)
    cell = zeros(Int, 3, N_max)
    return CellList(sc, head, list, cell)
end

# ────────────────────────────────────────────────────────────────────────────
# Low-level helpers
# ────────────────────────────────────────────────────────────────────────────

"""
    c_index(cl, ri_box) -> SVector{3,Int}

Map a position in box units (−0.5 … 0.5) to a 1-based cell triple.
"""
@inline function c_index(cl::CellList, ri::AbstractVector)
    any(abs.(ri) .> 0.5) && error("Atom outside main box: $ri")
    ci = floor.(Int, (ri .+ 0.5) .* cl.sc) .+ 1   # shift to 1-based
    ci = clamp.(ci, 1, cl.sc)
    return SVector{3,Int}(ci)
end

"""
    create_in_list!(cl, i, ci)

Insert atom `i` as the new head of cell `ci` (1-based triple).
"""
@inline function create_in_list!(cl::CellList, i::Int, ci::SVector{3,Int})
    cl.list[i]                      = cl.head[ci[1], ci[2], ci[3]]
    cl.head[ci[1], ci[2], ci[3]]   = i
    cl.cell[:, i]                   = ci
end

"""
    destroy_in_list!(cl, i, ci)

Remove atom `i` from cell `ci`.
"""
function destroy_in_list!(cl::CellList, i::Int, ci::SVector{3,Int})
    this = cl.head[ci[1], ci[2], ci[3]]

    if this == i                               # atom i is the head
        cl.head[ci[1], ci[2], ci[3]] = cl.list[i]
        return
    end

    while true                                 # walk the chain
        nxt = cl.list[this]
        if nxt == i
            cl.list[this] = cl.list[i]        # link over i
            return
        elseif nxt == 0
            error("destroy_in_list!: atom $i not found in cell $ci")
        else
            this = nxt
        end
    end
end

"""
    move_in_list!(cl, i, ci_new)

Move atom `i` from its current cell to `ci_new` (no-op if unchanged).
"""
@inline function move_in_list!(cl::CellList, i::Int, ci_new::SVector{3,Int})
    ci_old = SVector{3,Int}(cl.cell[:, i])
    ci_old == ci_new && return
    destroy_in_list!(cl, i, ci_old)
    create_in_list!(cl, i, ci_new)
end

# ────────────────────────────────────────────────────────────────────────────
# Build / rebuild the full list
# ────────────────────────────────────────────────────────────────────────────

"""
    make_list!(cl, N, r_box)

Rebuild the entire cell list from scratch for `N` normal particles.
Call this after loading a configuration or after a box-size change.
"""
function make_list!(cl::CellList, N::Int, r_box::Vector{<:AbstractVector})
    fill!(cl.head, 0)
    @inbounds for i in 1:N
        ci = c_index(cl, r_box[i])
        create_in_list!(cl, i, ci)
    end
end

# ────────────────────────────────────────────────────────────────────────────
# Neighbour iteration
# ────────────────────────────────────────────────────────────────────────────

# The 13 half-shell displacement vectors (plus the 13 mirrors = full 26 neighbours + self cell).
# Listed with inversion symmetry so that [1:13] is the "upper" half-shell and
# [14] is (0,0,0) (self), matching Allen & Tildesley convention.
const _NK = 13
const _D = let rows = [
    (-1,-1,-1), ( 0,-1,-1), ( 1,-1,-1),
    (-1, 1,-1), ( 0, 1,-1), ( 1, 1,-1),
    (-1, 0,-1), ( 1, 0,-1), ( 0, 0,-1),
    ( 0,-1, 0), ( 1,-1, 0), (-1,-1, 0),
    (-1, 0, 0), ( 0, 0, 0), ( 1, 0, 0),   # index 14 = self cell (k=0 in A&T)
    ( 1, 1, 0), (-1, 1, 0), ( 0, 1, 0),
    ( 0, 0, 1), (-1, 0, 1), ( 1, 0, 1),
    (-1,-1, 1), ( 0,-1, 1), ( 1,-1, 1),
    (-1, 1, 1), ( 0, 1, 1), ( 1, 1, 1),
]
    # pack into 3 × 27 matrix, columns -13:13 mapped to indices 1:27
    m = Matrix{Int}(undef, 3, 2*_NK+1)
    for (col, (dx,dy,dz)) in enumerate(rows)
        m[:,col] = [dx, dy, dz]
    end
    m
end

"""
    neighbours(cl, i, ci; half=true) -> Vector{Int}

Return list of neighbour atom indices for atom `i` in cell `ci`.

- `half=true`  : Newton's-third-law half-shell (13 neighbour cells + down-list in own cell).
                 Requires that `ci == cl.cell[:,i]` (atom has not moved since list was built).
- `half=false` : Full 26 surrounding cells + own cell (all atoms except `i` itself).
                 Use this when `i` is a *probe* position not stored in the list (e.g. a
                 proposed insertion site or the fractional particle).
"""
function neighbours(cl::CellList, i::Int, ci::SVector{3,Int}; half::Bool=true)
    sc = cl.sc

    if half
        SVector{3,Int}(cl.cell[:, i]) == ci ||
            error("neighbours: cell mismatch for atom $i")
        k1, k2 = _NK+1+1, 2*_NK+1   # columns for k=1..13 in _D  (upper half-shell)
    else
        k1, k2 = 1, 2*_NK+1         # all 27 columns
    end

    buf = Int[]
    sizehint!(buf, 64)

    @inbounds for col in k1:k2
        # neighbour cell index (1-based, periodic)
        cj1 = mod1(ci[1] + _D[1,col], sc)
        cj2 = mod1(ci[2] + _D[2,col], sc)
        cj3 = mod1(ci[3] + _D[3,col], sc)

        # which atom to start from in this cell?
        if col == _NK+1 && half           # self cell, half-shell: start down-list from i
            j = cl.list[i]
        else
            j = cl.head[cj1, cj2, cj3]
        end

        while j != 0
            j != i && push!(buf, j)
            j = cl.list[j]
        end
    end

    return buf
end

# ────────────────────────────────────────────────────────────────────────────
# Consistency check (for debugging)
# ────────────────────────────────────────────────────────────────────────────

function check_list(cl::CellList, N::Int, r_box::Vector{<:AbstractVector})
    for i in 1:N
        ci = c_index(cl, r_box[i])
        SVector{3,Int}(cl.cell[:,i]) == ci ||
            error("check_list: atom $i has wrong cell stored")
    end

    for cx in 1:cl.sc, cy in 1:cl.sc, cz in 1:cl.sc
        j = cl.head[cx, cy, cz]
        while j != 0
            all(cl.cell[:,j] .== [cx,cy,cz]) ||
                error("check_list: atom $j stored in wrong cell")
            j = cl.list[j]
        end
    end
    println("check_list: OK")
end
