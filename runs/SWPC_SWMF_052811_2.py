import magnetopost as mp

# Prior to running this script, execute the following from the directory of
# this script
#
# mkdir runs
#  or
# mkdir /Big/Disk/runs; ln -s /Big/Disk/runs
#
# cd runs
# wget -r -np -nH --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/SWPC_SWMF_052811_2
#
# Then modify dir_run in the following dictionary

# TODO: Get rCurrents from file
info = {
        "model": "SWMF",
        "run_name": "SWPC_SWMF_052811_2",
        "rCurrents": 4.0,
        "file_type": "cdf",
        "dir_run": "/Users/weigel/git/magnetopost/runs/SWPC_SWMF_052811_2",
        "dir_plots": "/Users/weigel/git/magnetopost/runs/SWPC_SWMF_052811_2.plots"
}

# Locations to compute B. See config.py for list of known points.
points_surf  = ["YKC"]
points_msph  = ["GMpoint1"]

compute = True
plot    = True
n_steps = None  # If None, process all files

# Create output dirs if needed and list of files to process
mp.util.setup(info)

if compute:
    #mp.postproc.job_ie(info, points_surf, n_steps=n_steps)
    mp.postproc.job_ms(info, points_msph, n_steps=n_steps)

if plot:
    if points_surf is not None:
        for point in points_surf:
            mp.plot.surf_point(info, point, n_steps=n_steps)
    if points_msph is not None:
        for point in points_msph:
            mp.plot.msph_point(info, point, n_steps=n_steps)
