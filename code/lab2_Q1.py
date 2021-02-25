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



