import os
import magnetopost as mp

# Prior to running this script, execute the following from the directory of
# this script:
#   mkdir runs && cd runs
#   wget -r -np -nH -R "index.html*" --cut-dirs=2 http://mag.gmu.edu/git-data/ccmc/SWPC_SWMF_052811_2/
#   cd SWPC_SWMF_052811_2; find . -name "*.cdf.gz" | xargs -J{} gunzip {}

# Get directory of this script
data_dir = os.path.split(os.path.abspath(__file__))[0]

# TODO: Get rCurrents from file
info = {
        "model": "SWMF",
        "run_name": "SWPC_SWMF_052811_2",
        "rCurrents": 4.0,
        "file_type": "cdf",
        "dir_run": os.path.join(data_dir, "SWPC_SWMF_052811_2"),
        "dir_plots": os.path.join(data_dir, "SWPC_SWMF_052811_2.plots"),
        "dir_derived": os.path.join(data_dir, "SWPC_SWMF_052811_2.derived"),
        "deltaB_files": {
            "YKC": os.path.join(data_dir, "SWPC_SWMF_052811_2", "2006_YKC_pointdata.txt")
        }
}

# Locations to compute B. See config.py for list of known points.
points  = ["YKC", "GMpoint1"]

compute = True
plot    = True
n_steps = None  # If None, process and plot all files. Use 3 or so for debugging.

# Create output dirs, if needed, and list of files to process
mp.util.setup(info)

if compute:
    mp.postproc.job_ie(info, points, n_steps=n_steps, stitch_only=compute==False)
    mp.postproc.job_ms(info, points, n_steps=n_steps, stitch_only=compute==False)

if plot:
    mp.plot.surf_point(info, points, n_steps=n_steps)
    mp.plot.msph_point(info, points, n_steps=n_steps)
