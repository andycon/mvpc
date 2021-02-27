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
analysis. See [load_ds.py][1].

We can use [load_ds.py][1] to load a dataset for subject '01':

```python
from load_ds import load

s = '01'
ds = load(s)
ds.shape # This should yield (200, 34928)
```

This will load a dataset that has contains 200 rows, with ten chunks
corresponding to the ten runs from the experiment, thus we have ten sets of 20
patterns. But for RDM analysis, we need just one pattern for each of the 20
conditions. To get this we can take the samples for each chunk and average those
together. I have written a helper function for this and added it to
[amvpa.py][2] as a class function within the Dataset class definition:

```python
ds = ds.mean_samples_across_chunks()
ds.shape # should now yield (20, 34928)
```
Alternatively, if you have updated your preprocessed data to include the
"allruns" GLM, you can skip the last step by passing keyword argument
**allruns** set to True to the load function:

```python
ds = load(s, allruns=True)
```

Because we have already loaded a precomputed searchlight map as part of our
**load** routine, running an RDM searchlight can be done like this:

```python
rdm_sl_ds = searchlight(ds, rdm_measure)
rdm_sl_ds.shape # This should yield (190, 34928)
```

This will take a while a while to compute. Note that the shape of the dataset
returned by the searchlight call is (190, 34298). This is because each column of
the samples is one vectorized RDM (lower/upper triangle). This could be useful,
but it is not yet what we need to find the solution for Q1, which requires the
ISC between RDMs at each searchlight. So, we need to do this for every subject,
then iterate over all features (columns) to combine searchlight-RDMs across
subjects to compute the average Pearson correlation between all pairs of
subject-specific RDMs.

One way to do this is to define a new function that can do the work for us. This
new function will be another example of a **dataset measure**. Because all of
our subjects are in the same space, and the features correspond across subjects,
we can load each subject and combine them into a single dataset by stacking the
samples vertically, and setting the chunks sample attribute to reflect the
subject IDs. Our new function will take this type of dataset as input. First
load the combined dataset:

```python
subs = [1, 12, 17, 24, 27, 31, 32, 33, 34, 36, 37, 41] 

ds_all = None
for s in subs:
    ds = load("{:02}".format(s), allruns=True)
    if ds_all is None:
        ds_all = ds
    else:
        ds_all.append(ds)

ds_all.set_sa("chunks",np.repeat(subs,20))
```
Let's define the **dataset measure** function called
**inter_chunk_rdm_correlation**:

```python
def inter_chunk_rdm_correlation(ds, metric="correlation"):
    """The average correlation between vectorized distance matrices computed for
    each chunk in the dataset.

    This **dataset measure** function breaks the input dataset up by chunks,
    computes the RDM for the samples in each chunk. The vectorized RDMs are
    vertically stacked as a set of row vectors and then the average correlation
    between RDMs is returned. 

    Note that if the chunks attribute designates individual subjects, then this
    function computes the inter-subject correlation (ISC) for RDMs.

    Parameters
    ----------
    ds : amvpa.Dataset object
        Chunks attribute must be set, and each chunk must have the same number
        of rows.

    metric: str, default "correlation"
        This sets the pdist metric parameter for the RDM distance metric. It
        does not affect the correlation between RDMs computed between chunks,
        which is set to return the average Pearson correlation.

    Returns
    -------
    mu_corr : np.array
        This will be a single value, but to conform to **dataset measure**
        output, this single value is wrapped in a 1 x 1 column vector.
    """
    rsa_by_chunks = None
    for ch in np.unique(ds.chunks):
        ds_ch = ds.select_chunk(ch)
        rsa_ch = rdm(ds_ch, metric=metric).reshape((1,-1)) # make row vector
        if rsa_by_chunks is None:
            rsa_by_chunks = rsa_ch
        else:
            rsa_by_chunks = np.vstack((rsa_by_chunks, rsa_ch))
    # Use 1 minus pdist, which will convert correlation distance to Pearson r,
    # and average the result:
    mu_corr = np.mean(1-dist.pdist(rsa_by_chunks,metric="correlation"))

    # return value as a numpy column vector
    return np.array([mu_corr]).reshape((1,1))
```

In order to run the ISC RDM searchlight we can put it all together in a
[script][3]:

```python
from load_ds import load
from amvpa import *

nproc = 12 # number of processors to use.
save_filename = "../Q1_isc_rdm_sl_allruns.nii.gz"

subs = [1, 12, 17, 24, 27, 31, 32, 33, 34, 36, 37, 41]

ds_all = None
for s in subs:
    ds = load("{:02}".format(s), allruns=True)
    if ds_all is None:
        ds_all = ds
    else:
        ds_all.append(ds)

ds_all.set_sa("chunks",np.repeat(subs,20))

res = searchlight(ds_all, inter_subject_rdm_correlation, nproc=nproc)

res.save_to_nifti(save_filename)
```
[1]: https://github.com/andycon/mvpc/blob/master/code/load_ds.py
[2]: https://github.com/andycon/mvpc/blob/master/code/amvpa.py
[3]: https://github.com/andycon/mvpc/blob/master/code/lab2_Q1.py

## Question 2

Calculate RDMs and make MDS plots for searchlights centered on EBA, FFA, pSTS,
OP (LIP?), IPS, PM coordinates (get coordinates from from Neurosynth).

First we need to look up some MNI coordinates for centers of ROIs.  Note that
the set below is neither a complete nor an authoritative set of MNI coordinates
for answering the question. These are just for illustrating how to solve the
problem.

```python
i# Define a dictionary for the coordinates
coords = {
    'EBA':  {   'R': (48, -72, 2),
                'L': (-48, -74, 4)
                },
    'FFA':  {   'R': (42, -50, -20),
                'L': (-40, -52, -20)
                },
    'STS':  {   'R': (56, -34, 4),
                'L': (-60, -34, 4)
                },
    'IPS':  {   'L': (-34, -50, 50),
                'R': (34, -54, 50)
                },
    'IPL':  {   'R': (58, -50, 41),
                'L': (-59, -50, 37)
                }
    }
```

The affine transformation matrix for our data, which is stored in the header of
of Nifti files, and in the 'a' attribute of our MVPA datasets is used to map the
coodinates of voxels in our data grid into MNI space. In order to map the MNI
coordinates for our ROIs (from NeuroSynth), we need to use the inverse of our
affine matrix. [Click here for a nice tutorial][4] on affine spatial
transformations.

I have added the following function to [amvpa.py][1] for this purpose:

```python
def get_f_indx_for_mni_coord(mni_coord, grid, aff): 
    x,y,z = grid 
    vol_idx = np.arange(x*y*z).reshape((x,y,z)) 
    indcs = np.indices((x,y,z)) 
    mni = np.ones((4,1)) 
    mni[:3,0] = mni_coord 
    c = np.round(np.dot(np.linalg.inv(aff),mni)) 
    idx = vol_idx[int(c[0]), int(c[1]), int(c[2])] 
    return idx 
```
Using this function and our **masked_neighborhood** function, we can define ROIs
for based on our coordinates:

```python
ds = load('01', allruns=True)

eba_l_idx = get_f_indx_for_mni_coord(   coord['EBA']['L'], 
                                        ds.fa['f_indx'],
                                        ds.a['aff'])
print(eba_l_idx) # should yield f_index: 86351

# Use this feature index to get masked_neighborhood with 200 voxels:
eba_l_roi = masked_neighborhood(eba_l_idx, ds.fa['f_indx'], grid, nvox=200)

print(len(eba_l_roi)) # should yield: 200
print(eba_l_roi[:5])  # should yield: [86351, 81038, 86282, 86352, 86420]

# Now define a new dataset by selecting just the features in eba_l_roi:
eba_l_ds = ds.select_features(eba_l_roi)

eba_l_ds.shape # should yield: (20, 200)
```

As a sanity check we can use our ROI dataset to make a mask and see what it
looks like using AFNI (or whatever). To do this we make a new dataset with the
same dataset attributes 'a' and feature attributes 'fa' as the ROI dataset, and
set the samples to a single row of ones, and save the mask dataset to a nifti
file:
```python
roi_mask = Dataset(np.ones((1,eba_l_ds.shape[1])))
roi_mask.a = eba_l_ds.a
roi_mask.fa = eba_l_ds.fa

# save to subject directory
roi_mask.save_to_nifti("../sub-rid000001/eba_l_mask.nii")
```
Here is a this mask on a SUMA surface:

![Left EBA on MNI Surface][5]


Finally to get the RDM for the left EBA:

```python
rdm = rdm_measure(eba_l_ds)
```
Putting this all together we can calculate the average RDM for a single ROI as
follows. (This code is in [lab2_Q2.py][6]):

```python
# Two additional imports are needed:
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

roi = "EBA"
hemi = "L"
nvox = 200

subs = [1, 12, 17, 24, 27, 31, 32, 33, 34, 36, 37, 41]

rdm_all = None
for s in subs:
    ds = load("{:02}".format(s), allruns=True)
    idx = get_f_indx_for_mni_coord(coords[roi][hemi],ds.a['grid'],ds.a['aff'])
    roi_nbh = masked_neighborhood(idx, ds.fa['f_indx'], ds.a['grid'],nvox=nvox)
    roi_ds = ds.select_features(roi_nbh)
    rdm = rdm_measure(ds)

    if rdm_all is None:
        rdm_all = rdm
    else:
        rdm_all = rdm_all + rdm

rdm_all = rdm_all/len(subs)

Ax = plt.matshow(squareform(rdm_all[:,0]))
plt.colorbar()
h = "Left"
if hemi == 'R':
    h = "Right"
plt.title("Average RDM for {} {}".format(h,roi))
plt.savefig("avg_{}_{}.png".format(h,roi))
```

The resulting RDM looks like this:

![Average RDM for Left EBA][7]

Note: MDS examples are yet to come.... stay tuned.



[4]: https://www.cs.utexas.edu/users/fussell/courses/cs384g-fall2011/lectures/lecture07-Affine.pdf
[5]: https://github.com/andycon/mvpc/blob/master/code/eba_l.png
[6]: https://github.com/andycon/mvpc/blob/master/code/lab2_Q2.py
[7]: https://github.com/andycon/mvpc/blob/master/code/avg_Left_EBA.png

## Question 3

Calculate correlations of RDMs between these 5 loci then plot as 2nd order DSM
and MDS.

