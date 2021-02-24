# Assignment 2 Tutorial

RSA – 20 conditions in behaving animal data

## Question 1

Map ISC of representational dissimilarity matrices (RDMs) over searchlights.

Here we need calculate the RDM matrix for each subject at each searchlight, and
then calculated the average pairwise correlation between those RDM across
subejects. ISC statnds for inter-subject correlation. The RDM can be calculated
using Pearson correlation distance, which is 1 - Pearson r, where the Pearson r
is calculated between every pair of condition patterns. Because we have 20
conditions, the full correlation matrix can be calculated as a 20 by 20 square
symetrical matrix with ones along the diagonal. Subtracting this matrix from
1 gives you a correlation distance matrix, which is also symetrical, but has
zeros along the diagonal. this can be done in Python using Numpy like this:

```python
import numpy as np

# make up some fake data with 20 rows and 200 columns:
samples = np.random.random((20,200))

r_matrix = np.coorcoef(samples)
r_matrix.shape # this should yield (20,20)

dist = 1 - r_matrix
```
However, we only need half of this matrix without the diagonal, either the upper
or lower triangle, which could be derived from this square symmetrical matrix.
But there is an easier way provided by the SciPy Python Library. This library
provides lots of useful functions for doing RSA analysis, and cluster analyses.
For the RDM, we can use pdist (note that the functionality for pdist in SciPy is
almost identical to the pdist function in Matlab). Usign Pdist in this context
we can do the following:

```python
from scipy.spatial.distance import pdist

dist = pdist(samples, metric='correlation')
dist.shape # this should yield (190,)
```
**pdist** gives us a vector of length (n * (n-1)) / 2, for n equaling the number
of conditions or rows in the input. Note that as a default pdist uses Euclidean
distance as the distance metric, so it necessary to set the metric value to
"correlation" if you want correlation distance. 


 

## Question 2

Calculate RDMs and make MDS plots for searchlights centered on EBA, FFA, pSTS,
OP (LIP?), IPS, PM coordinates (get coordinates from from Neurosynth).

## Question 3

Calculate correlations of RDMs between these 5 loci then plot as 2nd order DSM
and MDS.

