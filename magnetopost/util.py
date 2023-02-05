import os
import numpy as np

def Tstr(time, length=6):
    return '%.4d%.2d%.2dT%.2d%.2d%.2d'%(time[:6])


def setup(info):
    import magnetopost as mp

    assert os.path.exists(info["dir_run"]), "dir_run = " + info["dir_run"] + " not found"

    dir_steps = os.path.join(info["dir_derived"], "timeseries", "timesteps")

    if not os.path.exists(info["dir_plots"]):
        os.mkdir(info["dir_plots"])
        mp.logger.info("Created " + info["dir_plots"])

    if not os.path.exists(info["dir_derived"]):
        os.mkdir(info["dir_derived"])
        mp.logger.info("Created " + info["dir_derived"])
    
    if not os.path.exists(dir_steps):
        os.makedirs(os.path.join(dir_steps))
        mp.logger.info("Created " + dir_steps)

    info['files'] = {}

    if info['file_type'] == 'cdf':

        for subdir in ["GM_CDF", "IONO-2D_CDF"]:
            if subdir == "GM_CDF":
                # Note lower case 'cdf' in _GM_cdf_list
                file = open(os.path.join(info['dir_run'], subdir, info['run_name'] + '_GM_cdf_list'), 'r')
                key = "magnetosphere"
                prefix = ""
            if subdir == "IONO-2D_CDF":
                # Note upper case 'CDF' in _GM_CDF_list
                file = open(os.path.join(info['dir_run'], subdir, info['run_name'] + '_IE_CDF_list'), 'r')
                key = "ionosphere"
                prefix = info['run_name'] + '.swmf.' 

            info['files'][key] = {}

            lines = file.readlines()

            for line in lines:
                line = line.strip()
                linea = line.split(" ")
                if linea[0].endswith('.cdf') == False:
                    continue

                datea = linea[3].split("/")
                timea = linea[5].split(":")
                time = (int(datea[0]), int(datea[1]), int(datea[2]), int(timea[0]), int(timea[1]), int(timea[2]))

                info['files'][key][time] = os.path.join(info['dir_run'], subdir, prefix + linea[0])

    if info['file_type'] == 'out':

        generate_filelist_txts(info)

        info['files']['magnetosphere'] = {}
        with open(os.path.join(info["dir_derived"], 'magnetosphere_files.txt'), 'r') as f:
            for line in f.readlines():
                items = line.split(' ')
                time = tuple([int(ti) for ti in items[:6]])
                info['files']['magnetosphere'][time] = os.path.join(info['dir_run'], items[-1][:-1])

        info['files']['ionosphere'] = {}
        with open(os.path.join(info["dir_derived"], 'ionosphere_files.txt'), 'r') as f:
            for line in f.readlines():
                items = line.split(' ')
                time = tuple([int(ti) for ti in items[:6]])
                info['files']['ionosphere'][time] = os.path.join(info['dir_run'], items[-1][:-1])


def generate_filelist_txts(info):

    import os
    import re
    import json

    import magnetopost as mp

    dir_run = info["dir_run"]

    fn = os.path.join(info["dir_derived"], 'run.info.py')
    with open(fn, 'w') as outfile:
        outfile.write(json.dumps(info))

    mp.logger.info("Wrote {}".format(fn))

    if 'dir_magnetosphere' in info:
        dir_data = os.path.join(dir_run, info['dir_magnetosphere'])
    else:
        dir_data = os.path.join(dir_run, 'GM/IO2')

    magnetosphere_outs = sorted(os.listdir(dir_data))

    fn = os.path.join(info["dir_derived"], 'magnetosphere_files.txt')
    k = 0
    with open(fn,'w') as fl:
        regex = r"3d__.*\.out$"
        for fname in magnetosphere_outs:
            if re.search(regex, fname):
                k = k + 1
                assert(fname[:4] == '3d__')
                assert(fname[9] == '_')
                if fname[10] == 'e':
                    Y = int(fname[11:15])
                    M = int(fname[15:17])
                    D = int(fname[17:19])
                    assert(fname[19] == '-')
                    h = int(fname[20:22])
                    m = int(fname[22:24])
                    s = int(fname[24:26])
                    assert(fname[26] == '-')
                    mil = int(fname[27:30])
                    assert(fname[30:] == '.out')
                    fl.write(f'{Y} {M} {D} {h} {m} {s} {mil} {dir_data}/{fname}\n')

    mp.logger.info("Wrote {} file names to {}".format(k, fn))

    if 'dir_ionosphere' in info:
        dir_data = os.path.join(dir_run, info['dir_ionosphere'])
    else:
        dir_data = os.path.join(dir_run, 'IE/ionosphere')

    ionosphere_outs = sorted(os.listdir(dir_data))

    fn = os.path.join(info["dir_derived"], 'ionosphere_files.txt')
    k = 0
    with open(fn,'w') as fl:
        regex = r"i_.*\.tec$"

        for fname in ionosphere_outs:
            if re.search(regex, fname):
                k = k + 1
                assert(fname[:2] == 'i_')
                if fname[2] == 'e':
                    Y = int(fname[3:7])
                    M = int(fname[7:9])
                    D = int(fname[9:11])
                    assert(fname[11] == '-')
                    h = int(fname[12:14])
                    m = int(fname[14:16])
                    s = int(fname[16:18])
                    assert(fname[18] == '-')
                    mil = int(fname[19:22])
                    assert(fname[22:] == '.tec')
                    fl.write(f'{Y} {M} {D} {h} {m} {s} {mil} IE/ionosphere/{fname}\n')

    mp.logger.info("Wrote {} file names to {}".format(k, fn))


def GetMagnetometerCoordinates(magnetometer, time, csys, ctype, hxway=True):

    from magnetopost.config import defined_magnetometers

    if isinstance(magnetometer, str):
        magnetometer = defined_magnetometers[magnetometer]

    if hxway:
        from hxform import hxform as hx
        return hx.transform(np.array(magnetometer.coords), time, magnetometer.csys, csys, ctype_in=magnetometer.ctype, ctype_out=ctype)
    else: 
        import spacepy.coordinates as sc
        from spacepy.time import Ticktock
        cvals = sc.Coords(magnetometer.coords, magnetometer.csys, magnetometer.ctype)
        t_str = '%04d-%02d-%02dT%02d:%02d:%02d'%(time[:6])
        cvals.ticks = Ticktock(t_str, 'ISO')
        _ = cvals.convert('MAG', 'car')
        newcoord = cvals.convert(csys, ctype)
        return newcoord.data[0, :]
