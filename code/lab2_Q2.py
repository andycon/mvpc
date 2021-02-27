from load_ds import load
from amvpa import *
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform

# Some MNI coordinates for centers of ROIs
# Note that this is neither a complete nor an authoritative set of MNI
# coordinates for answering the question. These are just for illustrating how to
# solve the problem.
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





