import magnetopost as mp

# Prior to running this script, execute the following from the directory of
# this script
#
# mkdir /Big/Disk/runs/DIPTSUR2; cd runs; ln -s /Big/Disk/runs/DIPTSUR2
#
# cd runs
# wget -r -np -nH --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/DIPTSUR2
#
# Then modify dir_run in the info dictionary below
#
# Optional:
#   OS-X
#       find . -name "index.*" | xargs =i rm -f {}
#   Linux
#      find . -name "index.*" | xargs rm -f {}
#

# TODO: Get rCurrents from file
info = {
        "model": "SWMF",
        "run_name": "DIPTSUR2",
        "rCurrents": 1.8,
        "file_type": "out",
        "dir_run": "/Users/weigel/git/magnetopost/runs/DIPTSUR2",
        "dir_plots": "/Users/weigel/git/magnetopost/runs/DIPTSUR2.plots"
}

# Locations to compute B. See config.py for list of known points.
points_surf  = ["colaba"]
points_msph  = ["GMpoint1"]

compute = False
plot    = True
n_steps = None  # If None, process all files

# Create output dirs if needed and list of files to process
mp.util.setup(info)

if compute:
    mp.postproc.job_ie(info, points_surf, n_steps=n_steps)
    mp.postproc.job_ms(info, points_msph, n_steps=n_steps)

if plot:
    for point in points_surf:
        mp.plot.surf_point(info, point, n_steps=n_steps)
    for point in points_msph:
        mp.plot.msph_point(info, point, n_steps=n_steps)
