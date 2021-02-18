from load_ds import load
from amvpa import *
import sys

subs = [1, 12, 17, 24, 27, 31, 32, 33, 34, 36, 37, 41] 

ds_all = None
for s in subs:
    ds = load("{:02}".format(s))
    ds = ds.mean_samples_across_chunks()
    if ds_all is None:
        ds_all = ds
    else:
        ds_all.append(ds)

ds_all.set_sa("chunks",np.repeat(subs,20))

ds_all.samples = ds_all.samples[:,:100]
ds_all.fa['f_indx'] = ds_all.fa['f_indx'][:100]
ds_all.fa['sl_map'] = ds_all.fa['sl_map'][:100]

res = searchlight(ds_all, inter_subject_rdm_correlation)


res.save_to_nifti("../group_isc_rsa_searchlight.nii.gz")


