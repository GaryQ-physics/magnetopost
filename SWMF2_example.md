one thing to note is that SWPC run already has a file listing all cdfs with their times. its called `SWPC_SWMF_052811_2_GM_cdf_list` and `SWPC_SWMF_052811_2_IE_CDF_list`. All we do is reformat it. In the DIPTSUR2 case, there was no such files, but we can generate the list by parsing the filenames

note also, for SWPC this reformatting is done with awk, so that needs to be installed on system. (it should be, I believe its part of POSIX)

# Instructions

## 0. install repos

cd to installation directory, then install

```

```

## 1. Get files

cd to empty directory, then

```
wget -r -np -nH --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/SWPC_SWMF_052811_2/GM_CDF/
wget -r -np -nH --cut-dirs=2 http://mag.gmu.edu/git-data/sblake/SWPC_SWMF_052811_2/IONO-2D_CDF/
```

(you dont need the raw\_output/ directory)

## 2. Generate lists in textfile

```
cd SWPC_SWMF_052811_2;
mkdir derived; 
python -c 'from magnetopost.model_patches import SWMF2; SWMF2.generate_filelist_txts()'
```

## 3. write info file

manually go into text editor and write `derived/run.info.py`.
so you get the following:

```
$ cat derived/run.info.py
{
  "model"                 : "SWMF2",
  "run_name"              : "SWPC",
  "rCurrents"             : 4.0
}

$ ls derived/
 magnetosphere_files.txt
 ionosphere_files.txt
 run.info.py
```

## 4. run the postprocessing

```
mkdir -p derived/timeseries/slices
python -c 'from magnetopost import postproc; postproc.job(("YKC", "GMpoint6"))'
```

## 5. download the relevant magnetometer file
```
wget -P derived/ http://mag.gmu.edu/git-data/GaryQ-Physics/2006_YKC_pointdata.txt
```

## 6. run plotting script

now create the following scipt `do` so `$ cat do` yield:

```
#!/usr/bin/env python
from magnetopost.plotting import extract_magnetometer_data as emd
emd.surface_point('SWPC_SWMF_052811_2','YKC')
emd.msph_point('SWPC_SWMF_052811_2','GMpoint6')
```

then run as

```
./do --rootdir=$(pwd)/../
```

may need to chmod +x
