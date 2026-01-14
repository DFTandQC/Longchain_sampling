#!/usr/bin/env python3
"""
Configuration Comparison Tool - Test multiple parameter sets quickly

用法: python test_configurations.py
"""

import subprocess
import json
import sys
from pathlib import Path

# Define parameter sets to test
TEST_CONFIGS = {
    "compact": {
        "description": "Tight molecules (紧凑)",
        "params": {
            "dmin": 4.0,
            "dmax": 8.0,
            "lateral": 1.0,
            "jitter_deg": 20.0
        }
    },
    "balanced": {
        "description": "Balanced spacing (平衡)",
        "params": {
            "dmin": 6.0,
            "dmax": 10.0,
            "lateral": 2.0,
            "jitter_deg": 25.0
        }
    },
    "dispersed": {
        "description": "Loose molecules (分散)",
        "params": {
            "dmin": 10.0,
            "dmax": 15.0,
            "lateral": 3.0,
            "jitter_deg": 25.0
        }
    }
}

def run_test(config_name, params, nseeds=10):
    """Run sampling with given parameters."""
    print(f"\n{'='*70}")
    print(f"Testing: {config_name.upper()} - {TEST_CONFIGS[config_name]['description']}")
    print(f"{'='*70}")
    print(f"Parameters:")
    for key, val in params.items():
        print(f"  {key:15} = {val}")
    
    cmd = [
        "python", "main.py",
        "--molecules", '[{"name":"H2SO4","file":"monomer/opt-cisSA-B97-3c.xyz","count":1},'
                       '{"name":"PT","file":"monomer/opt-PT-B97-3c.xyz","count":2}]',
        "--dmin", str(params["dmin"]),
        "--dmax", str(params["dmax"]),
        "--lateral", str(params["lateral"]),
        "--jitter_deg", str(params["jitter_deg"]),
        "--nseeds", str(nseeds),
        "--out", f"out_{config_name}"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout + result.stderr
    
    # Parse results
    if "SUCCESS" in output or "Generated" in output:
        print("\n✅ SUCCESS - All seeds generated")
        # Count how many seeds were generated
        if "Generated" in output:
            for line in output.split('\n'):
                if "Generated" in line:
                    print(f"   {line.strip()}")
    else:
        print("\n❌ FAILED - See errors below:")
        print(output[-500:] if len(output) > 500 else output)
    
    return result.returncode == 0

def main():
    print("\n" + "="*70)
    print("CLUSTER GEOMETRY CONFIGURATION COMPARISON TOOL")
    print("="*70)
    print("\nThis tool tests different parameter sets for molecule spacing.")
    print("Generate different styles of clusters to find your preferred configuration.\n")
    
    # Run all tests
    results = {}
    for config_name, config in TEST_CONFIGS.items():
        success = run_test(config_name, config["params"], nseeds=15)
        results[config_name] = success
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    for config_name, success in results.items():
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{config_name:15} {status}")
    
    print(f"\n{'='*70}")
    print("OUTPUT LOCATIONS:")
    print(f"{'='*70}")
    print("Compact:   out_compact/seeds_H2SO4_combined.xyz")
    print("Balanced:  out_balanced/seeds_H2SO4_combined.xyz")
    print("Dispersed: out_dispersed/seeds_H2SO4_combined.xyz")
    print("\n💡 Open these files in a molecular viewer (Jmol, Avogadro, etc.)")
    print("   to visually compare the different configurations.\n")
    
    return 0 if all(results.values()) else 1

if __name__ == "__main__":
    sys.exit(main())
