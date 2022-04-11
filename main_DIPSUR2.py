import util
import magnetopost as mp

# Prior to running this script, execute the following from the directory of
# this script
#
# mkdir runs
#  or
# mkdir /Big/Disk/runs; ln -s /Big/Disk/runs
#
# cd runs
# wget -r -np -nH --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/DIPTSUR2/GM/
# wget -r -np -nH --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/DIPTSUR2/IE/
#
# Then modify the following line
dir_run = "/Users/weigel/git/magnetopost/runs/DIPTSUR2"

# TODO: Get rCurrents from file
info = {
        "model": "SWMF",
        "run_name": "DIPTSUR2",
        "rCurrents": 1.8
}

util.setup(dir_run, info)

n_steps = None

points = ["colaba"]
mp.postproc.job_ie(points, dir_run, n_steps=n_steps)
mp.postproc.job_ms(points, dir_run, n_steps=n_steps)

from magnetopost.plotting import extract_magnetometer_data as emd
emd.surface_point(dir_run, "colaba", n_steps=n_steps)

#mp.plot.surf_point(points, dir_run, n_steps=n_steps)
#mp.plot.msph_point(points, dir_run, n_steps=n_steps)
