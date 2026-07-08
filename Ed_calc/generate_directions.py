#!/usr/bin/env python3
"""Generate N random unit vectors uniformly distributed on the unit sphere.

Usage: python3 generate_directions.py [N] [seed]
Output: one line per direction: index dx dy dz
"""
import sys
import numpy as np

n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 42
np.random.seed(seed)

for i in range(n):
    v = np.random.randn(3)
    v = v / np.linalg.norm(v)
    print(f"{i+1:4d}  {v[0]: .8f}  {v[1]: .8f}  {v[2]: .8f}")
