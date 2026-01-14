# PT Ester-Core Placement Strategy Implementation

## Overview
Implemented PT (pentaerythritol ester) molecule placement using an ester-core-to-core strategy with minimal disruption to existing code structure.

## Key Implementation Details

### 1. PT Ester-Core Point Computation
**Function**: `compute_pt_ester_core(elems, X, k=8)`
- Identifies the k oxygen atoms closest to the molecular COM
- Returns the centroid of these k O atoms as the "ester-core point"
- Falls back to molecular COM if insufficient oxygen atoms present
- Configurable via `cfg.pt_k` (default: 8)

### 2. PT-to-PT Core Placement
**Function**: `place_pt_to_pt_core(mono_elems, mono_X, anchor_X, anchor_core, ...)`
- Places a new PT molecule using the ester-core strategy
- Randomly rotates the molecule
- Positions new PT's core at shell distance (dmin–dmax) from anchor core
- Adds small lateral offset for diversity
- Returns placed coordinates and metadata

### 3. Global PT Core Centroid
**Function**: `compute_global_pt_core_centroid(cluster_list, molecule_indices, molecules_data, pt_idx, pt_k=8)`
- Computes average of all placed PT core points
- Used as reference anchor for non-PT molecule placement
- Ensures non-PT molecules orient around the central ester-region cluster

### 4. Three-Mode Placement Strategy

#### Mode 1: Pure PT Only
- Detects via `_check_only_pt_molecules(molecules_data)`
- Uses core-to-core placement for all PT molecules
- Optimizes ester-region clustering without constraint

#### Mode 2: Mixed (PT + Non-PT)
**Phase 1: Place all PT molecules**
- Uses core-to-core strategy for PT↔PT anchoring
- Falls back to head-insert if PT anchors non-PT initially

**Phase 2: Place non-PT molecules**
- Computes global PT core centroid
- Anchors non-PT molecules to PT molecules (preferred)
- Uses tighter distances (75-85% of standard) to keep close to PT
- Non-PT placement uses existing head-insert with adaptive scaling

#### Mode 3: No PT (Standard Behavior)
- Falls back to original head-insert algorithm
- Preserves backward compatibility
- Cycles through molecule types in order

## Minimal Disruption Design

### Preserved Elements
- All existing utilities: `head_vector()`, `place_head_insert()`, clash checking
- Config parameters: `dmin`, `dmax`, `lateral`, `jitter_deg`, `max_trials_add`
- CLI interface, output format, and result structure unchanged
- Deterministic under same random seed

### Added Code
- 3 new helper functions (~150 lines total)
- Modified `_generate_single_seed()` with branching logic (~90 lines net)
- No changes to I/O, assembly, or restraint generation

### Integration Points
- `pt_idx` tracking in placement loop
- `molecule_indices` list distinguishes PT vs non-PT
- `placed_by_type` counter manages sequential expansion
- Reuses existing RNG seeding and clash detection

## Constraint Handling

### For PT-to-PT Placement
- Shell distance (dmin–dmax): Controls core-to-core separation
- Lateral offset: Adds 3D diversity around shell
- Random rotation: Full spherical orientation
- Jitter: Small additional rotation for variation

### For Non-PT Placement
- Distance scaling: `tighter_dmin = dmin × 0.75`, `tighter_dmax = dmax × 0.85`
- Lateral scaling: `tighter_lateral = lateral × 0.75`
- Maintains clash checking unchanged
- Head vector orientation preserved

## Configuration

### Default Parameters
- `pt_k`: 8 (oxygen atoms for core definition)
- Distance scaling: 0.75–0.85 ratio for non-PT
- All other sampling parameters unchanged

### Backward Compatibility
- Systems without PT use standard algorithm
- Legacy configs work without modification
- `cfg.pt_k` optional (defaults to 8 if missing)

## Testing Recommendations

1. **Pure PT Cluster**: Verify ester regions cluster tightly
2. **Mixed Cluster (PT + H2SO4)**: Confirm H2SO4 centers near PT core
3. **Non-PT Only**: Validate standard placement unchanged
4. **Determinism**: Same seed produces identical results
5. **Clash Detection**: Verify no overlaps with tighter spacing

## Code Structure Summary

```
compute_pt_ester_core()           → Identifies ester-core from k closest O atoms
  ↓
place_pt_to_pt_core()             → PT-to-PT core-based rigid placement
  ↓
compute_global_pt_core_centroid() → Average PT core centroid for non-PT anchor
  ↓
_generate_single_seed() [modified]
  ├─→ Mode 1 (Pure PT):     Core-to-core placement
  ├─→ Mode 2 (Mixed):       PT cores first, then non-PT near global core
  └─→ Mode 3 (No PT):       Standard head-insert
```

## Performance Impact
- Minimal overhead: O(k) per PT placement (k ≈ 8 O atoms)
- No impact on non-PT or PT-free systems
- Clash checking remains O(N²) as before
- Overall complexity unchanged from original
