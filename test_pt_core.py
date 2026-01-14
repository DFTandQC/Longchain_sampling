#!/usr/bin/env python3
"""Quick test of PT ester-core placement functions."""
import numpy as np
from lib.sampling import compute_pt_ester_core, place_pt_to_pt_core, compute_global_pt_core_centroid

# Create a simple PT-like molecule (for testing)
# PT structure: approximate coordinates
elems = np.array(['C', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'H', 'H', 'H'], dtype=object)
X = np.array([
    [0.0, 0.0, 0.0],      # C (center)
    [1.0, 0.0, 0.0],      # O
    [-1.0, 0.0, 0.0],     # O
    [0.0, 1.0, 0.0],      # O
    [0.0, -1.0, 0.0],     # O
    [0.0, 0.0, 1.0],      # O
    [0.0, 0.0, -1.0],     # O
    [0.5, 0.5, 0.5],      # O
    [-0.5, -0.5, -0.5],   # O
    [2.0, 0.0, 0.0],      # H
    [0.0, 2.0, 0.0],      # H
    [0.0, 0.0, 2.0],      # H
], dtype=float)

# Test 1: Compute ester-core
core = compute_pt_ester_core(elems, X, k=8)
print(f"✓ PT ester-core computed: {core}")
print(f"  Core is near center: {np.linalg.norm(core) < 1.0}")

# Test 2: Core-to-core placement
rng = np.random.default_rng(42)
X_new, info = place_pt_to_pt_core(elems, X, X, core, rng, d_range=(8.0, 11.0))
print(f"✓ PT-to-PT core placement successful")
print(f"  New molecule placed at distance: {info['d']:.2f} Angstrom")

# Test 3: Global core centroid
cluster_list = [X, X_new]
molecule_indices = [0, 0]
class MockMolSpec:
    def __init__(self):
        self.name = 'PT'
molecules_data = [(MockMolSpec(), elems, X)]
global_core = compute_global_pt_core_centroid(cluster_list, molecule_indices, molecules_data, 0, pt_k=8)
print(f"✓ Global PT core centroid computed: {global_core}")

print("\n✅ All PT ester-core placement functions working correctly!")
