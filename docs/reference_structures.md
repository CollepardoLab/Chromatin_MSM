# Reference Structures

Each subdirectory contains a reference structure for a tetra-nucleosome system at a given linker DNA length (in bp, e.g. `25`, `30`). Suffixes `_A` and `_AT` denote variant DNA sequences. 

## How reference structures were generated

### 1. Feature encoding

Each tetra-nucleosome conformation was encoded using **distance-symmetrized features** computed from the six pairwise inter-nucleosome distances (d12, d13, d14, d23, d24, d34) between the four nucleosome core particles.

The symmetrization respects the two-fold pseudo-symmetry of the tetra-nucleosome, collapsing the six raw distances into six symmetric features. The method follows the approach of [Ding and Zhang](https://www.nature.com/articles/s41467-021-21377-z):

```python
def calculate_symmetrized_features(dists):
    s1 = dists[:,0] + dists[:,5]                              # d12 + d34
    s2 = dists[:,1] + dists[:,4]                              # d13 + d24
    s3 = dists[:,2]                                           # d14
    s4 = dists[:,3]                                           # d23
    s5 = dists[:,0]*dists[:,1] + dists[:,4]*dists[:,5]        # d12*d13 + d24*d34
    s6 = (dists[:,0]**2)*dists[:,1] + (dists[:,4]**2)*dists[:,5]  # d12²*d13 + d24²*d34

    return np.column_stack((s1, s2, s3, s4, s5, s6))
```

`dists` is an `(N, 6)` array of distances in column order: d12, d13, d14, d23, d24, d34.

### 2. Clustering

K-means clustering (k = 3) was applied to the symmetrized feature vectors using the [deeptime](https://deeptime-ml.github.io/) implementation:

```python
from deeptime.clustering import KMeans

estimator = KMeans(n_clusters=3)
clustering = estimator.fit(features).fetch_model()
```

### 3. Selecting the reference structure

The most populated cluster (largest number of assigned frames) was identified. The frame closest to that cluster's centroid in feature space was selected as the reference structure and saved as `<linker-length>rep.pdb`.
