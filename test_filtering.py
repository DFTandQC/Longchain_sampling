#!/usr/bin/env python3
"""Quick test of filtering functions."""

import numpy as np
from lib.sampling import (
    radius_of_gyration,
    count_inter_molecular_contacts,
    compute_rmsd_fingerprint,
    filter_accept
)
from lib.config import ClusterConfig, MoleculeSpec

# Test 1: Radius of gyration
X1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
rg1 = radius_of_gyration(X1)
print(f"Test 1 - Radius of gyration: {rg1:.3f} (expect ~0.707)")
assert 0.6 < rg1 < 0.8, f"Rg out of range: {rg1}"

# Test 2: Inter-molecular contacts
X2 = np.array([[0, 0, 0], [0.1, 0, 0], [10, 0, 0], [10.1, 0, 0]])
mol_indices2 = [0, 0, 1, 1]
contacts2 = count_inter_molecular_contacts(X2, ['C', 'C', 'C', 'C'], mol_indices2, contact_cut=3.8)
print(f"Test 2 - Inter-molecular contacts (dist=10): {contacts2} (expect 0)")
assert contacts2 == 0, f"Expected 0 contacts, got {contacts2}"

X3 = np.array([[0, 0, 0], [0.1, 0, 0], [2, 0, 0], [2.1, 0, 0]])
mol_indices3 = [0, 0, 1, 1]
contacts3 = count_inter_molecular_contacts(X3, ['C', 'C', 'C', 'C'], mol_indices3, contact_cut=3.8)
print(f"Test 3 - Inter-molecular contacts (dist=2): {contacts3} (expect 4)")
assert contacts3 == 4, f"Expected 4 contacts, got {contacts3}"

# Test 4: Fingerprinting
fp1 = compute_rmsd_fingerprint(X1)
fp2 = compute_rmsd_fingerprint(X1 + 0.001)  # Tiny shift
print(f"Test 4 - RMSD fingerprint: same=True if equal, {fp1 == fp2} (expect True)")

# Test 5: Filter accept (basic config)
cfg = ClusterConfig(enable_filter=False)  # Disabled should always pass
result, reason = filter_accept([], [], [], cfg)
print(f"Test 5 - Filter disabled: {result} (expect True)")
assert result is True, f"Expected True when filter disabled, got {result}"

# Test 6: Filter with max_rg
cfg_with_filter = ClusterConfig(
    enable_filter=True,
    max_rg=1.0,
    min_contacts=1,
    contact_cut=3.8
)
# Create mock cluster
X_cluster = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], 
                      [10, 0, 0], [10.1, 0, 0], [10, 1, 0], [10, 0, 1]])
molecule_indices = [0, 0, 0, 0, 1, 1, 1, 1]
molecules_data = [
    (MoleculeSpec("A", "test.xyz", 1), np.array(['C']*4), np.array([[0,0,0]]*4)),
    (MoleculeSpec("B", "test.xyz", 1), np.array(['C']*4), np.array([[0,0,0]]*4))
]
cluster_list = [X_cluster[:4], X_cluster[4:]]
result, reason = filter_accept(cluster_list, molecule_indices, molecules_data, cfg_with_filter)
print(f"Test 6 - Large Rg rejection: passes={result}, reason={reason}")
# This should fail due to large Rg (>1.0)
assert result is False and "rg_too_large" in str(reason), f"Expected Rg rejection, got {reason}"

print("\nAll filtering tests passed!")
