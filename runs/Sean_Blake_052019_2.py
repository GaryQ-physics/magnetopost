import os
import magnetopost as mp

# Prior to running this script, execute the following from the directory of
# this script:
#
#   wget -r -np -nH --cut-dirs=2 -R "index.html*" http://mag.gmu.edu/git-data/sblake/Sean_Blake_052019_2/
#
#   Generated files and plots are available by replacing trailing slash in above URL
#   with .plots or .derived.

# Get directory of this script
data_dir = os.path.split(os.path.abspath(__file__))[0]

# TODO: Get rCurrents from file
info = {
        "model": "SWMF",
        "run_name":  "Sean_Blake_052019_2",
        "rCurrents": 1.8,
        "file_type": "out",
        "dir_run":     os.path.join(data_dir, "Sean_Blake_052019_2", "RAW"),
        "dir_plots":   os.path.join(data_dir, "Sean_Blake_052019_2.plots"),
        "dir_derived": os.path.join(data_dir, "Sean_Blake_052019_2.derived")
}

# Locations to compute B. See config.py for list of known points.
points  = ["Colaba", "GMpoint1"]

compute = True
plot    = True
n_steps = None  # If None, process and plot all files. Use 2 or so for debugging.

# Create output dirs, if needed, and list of files to process
mp.util.setup(info)

if compute:
    mp.postproc.job_ie(info, points, n_steps=n_steps, stitch_only=compute==False)
    mp.postproc.job_ms(info, points, n_steps=n_steps, stitch_only=compute==False)

if plot:
    mp.plot.surf_point(info, points, n_steps=n_steps)
    mp.plot.msph_point(info, points, n_steps=n_steps)
