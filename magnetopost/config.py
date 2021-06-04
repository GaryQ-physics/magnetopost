from collections import namedtuple
Magnetometer = namedtuple('Magnetometer', ['name','csys','ctype','coords'])


defined_magnetometers = {
    'YKC'      : Magnetometer(name='YKC'     ,csys='GEO',ctype='sph',coords=(1.,  62.480 ,  245.518 ) ),
    'FRN'      : Magnetometer(name='FRN'     ,csys='GEO',ctype='sph',coords=(1.,  37.0913, -119.7193) ),
    'FUR'      : Magnetometer(name='FUR'     ,csys='GEO',ctype='sph',coords=(1.,  48.17  ,  11.28   ) ),
    'colaba'   : Magnetometer(name='colaba'  ,csys='GEO',ctype='sph',coords=(1.,  18.907 ,  72.815  ) ),
    'GMpoint1' : Magnetometer(name='GMpoint1',csys='GSM',ctype='car',coords=(0.71875,0.09375,-3.71875) ),
    }



########################################################################
import os
import resource

if os.path.exists('/home/gary/'):
    #https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
    soft, hard = int(13*2**30), int(13*2**30)
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))
elif os.path.exists('/home/gquaresi/'):
    soft, hard = 90*2**30, 90*2**30
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))
