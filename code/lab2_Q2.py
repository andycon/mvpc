from amvpa import *

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


# The affine transformation matrix for our data, which is stored in the header
# of of Nifti files, and in the 'a' attribute of our MVPA datasets is used to
# map the coodinates of voxels in our data grid into MNI space. In order to map
# the MNI coordinates for our ROIs (from NeuroSynth), we need to use the inverse
# of our affine matrix.

