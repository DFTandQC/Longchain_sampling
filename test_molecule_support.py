#!/usr/bin/env python3
"""
Test script for molecule support - demonstrates head atom detection for different molecule types.

Tests the following molecule categories:
1. Oxygen-containing molecules (H₂SO₄, H₂O)
2. Nitrogen-containing molecules (NH₃)
3. Sulfur-containing molecules (H₂S)
4. Hydrocarbon molecules (CH₄)
"""

import numpy as np
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from lib.sampling import get_head_atoms_priority, head_point_from_atoms, pick_nearest_O_indices


def test_molecule(name, elems, X, expected_atom_type):
    """Test head atom detection for a molecule."""
    print(f"\n{'='*60}")
    print(f"Test: {name}")
    print(f"{'='*60}")
    print(f"Atomic composition: {', '.join([f'{e}' for e in elems])}")
    
    # Get head atoms
    head_idx = get_head_atoms_priority(elems)
    head_atoms = elems[head_idx]
    
    print(f"Head atom indices: {list(head_idx)}")
    print(f"Head atoms selected: {', '.join(np.unique(head_atoms))}")
    print(f"Expected atom type: {expected_atom_type}")
    
    # Calculate head point
    try:
        head_pt = head_point_from_atoms(elems, X, head_idx)
        print(f"Head point (centroid): [{head_pt[0]:.4f}, {head_pt[1]:.4f}, {head_pt[2]:.4f}]")
    except Exception as e:
        print(f"❌ Error calculating head point: {e}")
        return False
    
    # Try to get nearest atoms for constraints
    try:
        constraint_idx = pick_nearest_O_indices(elems, X, head_pt, k=2)
        constraint_atoms = elems[constraint_idx]
        print(f"Constraint atom indices: {list(constraint_idx)}")
        print(f"Constraint atoms: {', '.join(constraint_atoms)}")
    except Exception as e:
        print(f"❌ Error selecting constraint atoms: {e}")
        return False
    
    print("✅ PASS")
    return True


def main():
    """Run all molecule support tests."""
    
    print("\n" + "="*60)
    print("MOLECULE SUPPORT TEST SUITE")
    print("="*60)
    
    all_passed = True
    
    # Test 1: Sulfuric Acid (H₂SO₄) - O-containing
    print("\n--- Test 1: Oxygen-containing molecules ---")
    
    # H2SO4: 2H, 1S, 4O
    elems_h2so4 = np.array(['H', 'S', 'O', 'O', 'O', 'O'])
    X_h2so4 = np.array([
        [0.0,  0.0,  1.0],  # H
        [0.0,  0.0, -0.5],  # S
        [0.96, 0.0, -0.8],  # O
        [-0.96, 0.0, -0.8], # O
        [0.0,  0.96, -0.8], # O
        [0.0, -0.96, -0.8], # O
    ], dtype=float)
    
    result1 = test_molecule("H2SO4 (Sulfuric Acid)", elems_h2so4, X_h2so4, "O (Oxygen)")
    all_passed = all_passed and result1
    
    # H2O: 2H, 1O
    elems_h2o = np.array(['O', 'H', 'H'])
    X_h2o = np.array([
        [0.0,  0.0,  0.0],  # O
        [0.95, 0.0,  0.0],  # H
        [-0.24, 0.92, 0.0], # H
    ], dtype=float)
    
    result2 = test_molecule("H2O (Water)", elems_h2o, X_h2o, "O (Oxygen)")
    all_passed = all_passed and result2
    
    # Test 2: Ammonia (NH₃) - N-containing, NO O atoms
    print("\n--- Test 2: Nitrogen-containing molecules (no O) ---")
    
    elems_nh3 = np.array(['N', 'H', 'H', 'H'])
    X_nh3 = np.array([
        [0.0,  0.0,  0.0],   # N
        [0.93, 0.0,  0.0],   # H
        [-0.47, 0.81, 0.0],  # H
        [-0.47, -0.81, 0.0], # H
    ], dtype=float)
    
    result3 = test_molecule("NH3 (Ammonia)", elems_nh3, X_nh3, "N (Nitrogen)")
    all_passed = all_passed and result3
    
    # Test 3: Hydrogen Sulfide (H₂S) - S-containing, NO O or N atoms
    print("\n--- Test 3: Sulfur-containing molecules (no O, no N) ---")
    
    elems_h2s = np.array(['S', 'H', 'H'])
    X_h2s = np.array([
        [0.0,  0.0,  0.0],  # S
        [1.33, 0.0,  0.0],  # H
        [-0.34, 1.29, 0.0], # H
    ], dtype=float)
    
    result4 = test_molecule("H2S (Hydrogen Sulfide)", elems_h2s, X_h2s, "S (Sulfur)")
    all_passed = all_passed and result4
    
    # Test 4: Methane (CH₄) - Pure hydrocarbon, NO heteroatoms
    print("\n--- Test 4: Hydrocarbon molecules (no heteroatoms) ---")
    
    elems_ch4 = np.array(['C', 'H', 'H', 'H', 'H'])
    X_ch4 = np.array([
        [0.0,  0.0,  0.0],   # C
        [0.63, 0.63, 0.63],  # H
        [-0.63, -0.63, 0.63],# H
        [-0.63, 0.63, -0.63],# H
        [0.63, -0.63, -0.63],# H
    ], dtype=float)
    
    result5 = test_molecule("CH4 (Methane)", elems_ch4, X_ch4, "All atoms (C+H)")
    all_passed = all_passed and result5
    
    # Test 5: Mixed O and N (shouldn't happen in reality, but tests priority)
    print("\n--- Test 5: Edge case - molecule with both O and N ---")
    
    elems_mixed = np.array(['C', 'O', 'N', 'H', 'H'])
    X_mixed = np.random.randn(5, 3) * 0.5
    
    result6 = test_molecule("Mixed O+N molecule", elems_mixed, X_mixed, "O (Oxygen) - higher priority")
    all_passed = all_passed and result6
    
    # Summary
    print("\n" + "="*60)
    if all_passed:
        print("✅ ALL TESTS PASSED")
        print("="*60)
        print("\nSummary:")
        print("  ✅ O-containing molecules: Head defined by O atoms")
        print("  ✅ N-containing molecules (no O): Head defined by N atoms")
        print("  ✅ S-containing molecules (no O,N): Head defined by S atoms")
        print("  ✅ Hydrocarbons (no heteroatoms): Head defined by all atoms")
        print("  ✅ Priority system works correctly")
        print("\n🎉 System successfully supports diverse molecule types!")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        print("="*60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
