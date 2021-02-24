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
For the RDM, we can use **pdist** (note that the functionality for **pdist** in
SciPy is almost identical to the **pdist** function in Matlab). Using **pdist**
in this context we can do the following:

```python
import scipy.spatial.distance as dist

rdm = dist.pdist(samples, metric='correlation')
rdm.shape # this should yield (190,)
```
**pdist** gives us a vector of length (n * (n-1)) / 2, for n equaling the number
of conditions or rows in the input. Note that as a default **pdist** uses Euclidean
distance as the distance metric, so it necessary to set the metric value to
"correlation" if you want correlation distance. 

In order to run an **RDM** searchlight in the context of our MVPA toolset, we
can define an **RDM Measure**, which is a member of the abstract set of
**dataset measures**. The concept of a **dataset measure** in PyMVPA,
CoSMoMVPA, and in our own aMVPA, is an abstraction that defines the behavior of a
function that takes a **Dataset** as input, and returns a single value, or a set
of values, that represent the outcome of some measure that desribes the dataset
as a whole (as opposed to say, a single feature, sample, or whatever). For
Lab Assignment 1, we defined the function **cross_validated_classification**.
This function is an example of a **dataset measure** because it takes a dataset
as input, and returns the classification accuracy across data-folds for some
classification procedure, e.g. LinearSVC from **scikit-learn**. In order to make
the concept of a **dataset measure** useful in different contexts, we also need
to know what kind of data structure will be returned be a measure. For this
let's stipulate that any **dataset measure** returns values in the form of a
column vector. So within the **amvpa.py**, we can define an **rdm measure** as
follows:

```python
def rdm_measure(ds, metric="correlation"):
    rdm = dist.pdist(ds.samples, metric=metric)
    return rdm.reshape((rdm.shape[0],-1)) # enforce column vector output
```

Now that we have defined the RDM measure, we can use it directly to run a
searchlight on a given dataset. I have written a helper function to load a
dataset for the 20-condition Behaving Animals dataset that we have been using.
This is just so that I don't have to keep rewriting the code for each new
analysis. See [load_ds.py](https://github.com/andycon/mvpc/blob/master/code/load_ds.py).

We can use
[load_ds.py](https://github.com/andycon/mvpc/blob/master/code/load_ds.py) to
load a dataset for subject '01':

```python
from load_ds import load

s = '01'
ds = load(s)
```

This will load a dataset that has contains 200 rows, with ten chunks
corresponding to the ten runs from the experiment, thus we have ten sets of 20
patterns. But for RDM analysis, we need just one pattern for each of the 20
conditions. To get this we can take the samples for each chunk and average those
together. I have written a helper function for this and added it to
[amvpa.py](https://github.com/andycon/mvpc/blob/master/code/amvpa.py#"class%20Dataset").



 

## Question 2

Calculate RDMs and make MDS plots for searchlights centered on EBA, FFA, pSTS,
OP (LIP?), IPS, PM coordinates (get coordinates from from Neurosynth).

## Question 3

Calculate correlations of RDMs between these 5 loci then plot as 2nd order DSM
and MDS.

