# Magnetosphere

```
$ cd /media/.../sblake/DIPTSUR2/
$ ls
 derived/
 GM/
 Param.in
 PostProcessConfig.txt
$ cat PostProcess.conf.py
{
  "model"                 : "SWMF"
  "run_name"              : "DIPTSUR2"
  "magnetosphere_files"   : "./GM/IO2/"
  "ionosphere_files"      : "./IE/ionosphere/"
  "magnetometer_files"    : "./GM/IO2/"

  "rCurrents"             : 1.8
 }

$ cat run.info.py
{
  "model"                 : "SWMF"
  "run_name"              : "DIPTSUR2"
  "rCurrents"             : 1.8
 }

$ cat magnetosphere_files.txt
2019 09 02 04 10 00 000 ./GM/IO2/3d__var_2_e20190902-041000-000.out
2019 09 02 04 11 00 000 ./GM/IO2/3d__var_2_e20190902-041100-000.out
2019 09 02 04 12 00 016 ./GM/IO2/3d__var_2_e20190902-041200-016.out
2019 09 02 04 13 00 000 ./GM/IO2/3d__var_2_e20190902-041300-000.out
2019 09 02 04 14 00 000 ./GM/IO2/3d__var_2_e20190902-041400-000.out
2019 09 02 04 15 00 031 ./GM/IO2/3d__var_2_e20190902-041500-031.out
2019 09 02 04 16 00 000 ./GM/IO2/3d__var_2_e20190902-041600-000.out
2019 09 02 04 17 00 000 ./GM/IO2/3d__var_2_e20190902-041700-000.out
...

$ cat ionosphere_files.txt
2019 09 02 04 10 00 000 ./IE/ionosphere/i_e20190902-041000-000.tec
2019 09 02 04 11 00 000 ./IE/ionosphere/i_e20190902-041100-000.tec
2019 09 02 04 12 00 016 ./IE/ionosphere/i_e20190902-041200-016.tec
2019 09 02 04 13 00 000 ./IE/ionosphere/i_e20190902-041300-000.tec
2019 09 02 04 14 00 000 ./IE/ionosphere/i_e20190902-041400-000.tec
2019 09 02 04 15 00 031 ./IE/ionosphere/i_e20190902-041500-031.tec
2019 09 02 04 16 00 000 ./IE/ionosphere/i_e20190902-041600-000.tec
2019 09 02 04 17 00 000 ./IE/ionosphere/i_e20190902-041700-000.tec
...


$ python -c 'from magnetopost import postproc; postproc.job(("colaba", "GMpoint1"), do_summary=True)'
```
