import sys
from pathlib import Path

# Ensure project root (sampling/) is on sys.path so 'lib' package imports work
proj_root = str(Path(__file__).resolve().parents[1])
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)

from lib.config import ClusterConfig
from lib.sampling import read_xyz, _generate_single_seed

# Minimal smoke test: build one seed with tightened mixed-cluster distances
cfg = ClusterConfig()
cfg.nseeds = 1
cfg.random_seed = 42
cfg.dmin = 7.0
cfg.dmax = 11.0
cfg.lateral = 1.5

# Load molecules defined in cfg.molecules
molecules_data = []
for spec in cfg.molecules:
    elems, X, _ = read_xyz(spec.file)
    molecules_data.append((spec, elems, X))

res = _generate_single_seed(1, cfg, molecules_data, cfg.random_seed * 1000)
if res is None:
    print('SMOKE_TEST: FAIL - _generate_single_seed returned None')
    sys.exit(2)
else:
    print('SMOKE_TEST: OK')
    print('composition:', res['composition'])
    print('n_atoms total:', res['X_all'].shape[0])
    print('first 3 coords:\n', res['X_all'][:3])
    print('pairs count:', len(res['pairs']))
    print('pairs sample:', res['pairs'][:10])
