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
        "dir_run": "/Users/weigel/git/magnetopost/runs/SWPC_SWMF_052811_2"
}

points = ["YKC"] # Locations to compute B. See config.py for list of known points.
plot_only = False   # If True, re-create plots only
n_steps = 10        # If None, process all files

if plot_only == False:
    # Create output dirs and list of files to process
    mp.util.setup(info)

    # Do calculations
    mp.postproc.job_ie(info, points, n_steps=n_steps)
    mp.postproc.job_ms(info, points, n_steps=n_steps)

for point in points:
    mp.plot.surf_point(info, point, n_steps=n_steps)
    mp.plot.msph_point(info, point, n_steps=n_steps)

