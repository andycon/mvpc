from amvpa import *

def load(sub, allruns=False):
    print("Loading subject ====> {}".format(sub))
    dpath = "../../preproc/sub-rid0000{}/".format(sub)
    opath = "../sub-rid0000{}/".format(sub)

    tasks = ["beh","tax"]
    data_fn = dpath+"Qtstats_{}_run-{}.nii.gz"
    allruns_fn = dpath+"Qtstats_allruns.nii.gz"
    mask_fn = dpath+"glasser_masks.nii.gz"
    animals = ['bird','insect','primate','reptile','ungulate']
    behaviors = ['eating','fighting','running','swimming']
    twenty_conds = ['{}_{}'.format(a,b) for a in animals for b in behaviors]
    animals = np.repeat(animals,4) # all animals for one run
    behaviors = np.tile(behaviors,5) # all behaviors for each run
        
    ## Pre-calculated searchlight map
    sl_map_fn = "glassGM_Searchlight.txt"

    ds = None
    if allruns:
        ds = Dataset(allruns_fn, mask=mask_fn)
        ds.set_sa('targets',twenty_conds)
        ds.set_sa('animals', animals)
        ds.set_sa('behaviors', behaviors)
    else:

        for task in tasks:
            for r in range(1,6):
                if ds is None:
                    ds = Dataset(data_fn.format(task, r), mask=mask_fn)
                else:
                    ds.append(Dataset(data_fn.format(task,r), mask=mask_fn))

        ds.set_sa('targets', np.tile(twenty_conds,10))
        ds.set_sa('chunks', np.repeat(range(10),20))
        ds.set_sa('animals', np.tile(animals,10))
        ds.set_sa('behaviors', np.tile(behaviors, 10))

    # Load a pre-computed searchlight space

    sl_map = eval(open(sl_map_fn,'r').read())
    ds.set_searchlight_map(sl_map)

    return ds

