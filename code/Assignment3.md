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
symettrical matrix with ones along the diagonal. Subtracting 1 this matrix from
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


First lets define 

## Question 2

Calculate RDMs and make MDS plots for searchlights centered on EBA, FFA, pSTS,
OP (LIP?), IPS, PM coordinates (get coordinates from from Neurosynth).

## Question 3

Calculate correlations of RDMs between these 5 loci then plot as 2nd order DSM
and MDS.

