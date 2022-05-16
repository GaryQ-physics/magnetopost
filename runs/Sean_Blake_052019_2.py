import magnetopost as mp

# Prior to running this script, execute the following from the directory of
# this script
#
# mkdir runs
#  or
# mkdir /Big/Disk/runs; ln -s /Big/Disk/runs
#
# cd runs
# wget -r -np -nH --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/Sean_Blake_052019_2/
#
# Then modify dir_run in the following dictionary

# TODO: Get rCurrents from file
info = {
        "model": "SWMF",
        "run_name": "Sean_Blake_052019_2",
        "rCurrents": 1.8,
        "file_type": "out",
        "dir_run": "/Users/weigel/git/magnetopost/runs/Sean_Blake_052019_2/RAW"
}

points = ["colaba"] # Locations to compute B. See config.py for list of known points.
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
