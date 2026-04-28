# pixi run -e main python -u src/notebooks/spatial/kennedi_Xenium/tools/run_cnmf.py &> outs/spatial/tools/cnmf.out

from pathlib import Path
import numpy as np
import pickle
import multiprocessing
from cnmf import cNMF


if __name__ == "__main__":
    N_ITERS = 20 # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
    N_GENES = 2000 ## Number of over-dispersed genes to use for running the actual factorizations
    SEED = 123
    CORES = 10
    RUN_NAME = 'kennedi_xenium'
    HOMEDIR = Path(".")
    DATADIR = HOMEDIR / "data" / "processed" / "spatial" / "Xenium" / "kennedi_flu"
    FILE_NAME = "integrated.h5ad"
    SAVEDIR = DATADIR / "tools" / "cnmf"
    SAVEDIR.mkdir(parents=True, exist_ok=True)
    Ks = range(5,20())

    ## Specify the Ks to use as a space separated list in this case "5 6 7 8 9 10"
    K = ' '.join([str(i) for i in Ks])

    countfn = str(DATADIR / FILE_NAME)

    ## Initialize the cnmf object that will be used to run analyses
    cnmf_obj = cNMF(output_dir=SAVEDIR, name=RUN_NAME)
    cnmf_obj.prepare(counts_fn=countfn, components=np.array(Ks), n_iter=N_ITERS, seed=SEED, num_highvar_genes=N_GENES)

    # mediate partial fail
    cnmf_obj.update_nmf_iter_params()

    # Function to be executed by each process
    def factorize_worker(worker_i):
        cnmf_obj.factorize(worker_i=worker_i, total_workers=CORES, skip_completed_runs=True)

    # Create a pool of workers and distribute the tasks
    if __name__ == "__main__":
        with multiprocessing.Pool(processes=CORES) as pool:
            pool.map(factorize_worker, range(CORES))
    
    cnmf_obj.combine()

    cnmf_obj.k_selection_plot(close_fig=False)

    # Save the cnmf_obj as a pickle file
    with open(SAVEDIR / "cnmf_obj.pkl", "wb") as f:
        pickle.dump(cnmf_obj, f)