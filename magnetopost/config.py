from collections import namedtuple
Magnetometer = namedtuple('Magnetometer', ['name','csys','ctype','coords'])

defined_magnetometers = {
    'center'    : Magnetometer(name='center'  ,csys='GSM',ctype='car',coords=(0.,0.,0.) ),
    'Colaba'    : Magnetometer(name='Colaba'  ,csys='GEO',ctype='sph',coords=(1.,  18.907 ,  72.815  ) ), # r, lat, lon (in degrees)

    'YKC'       : Magnetometer(name='YKC'     ,csys='GEO',ctype='sph',coords=(1.,  62.480 ,  245.518 ) ),
    'YKC_N'     : Magnetometer(name='YKC_N'   ,csys='GEO',ctype='sph',coords=(1.,  63.480 ,  245.518 ) ),
    'YKC_S'     : Magnetometer(name='YKC_S'   ,csys='GEO',ctype='sph',coords=(1.,  61.480 ,  245.518 ) ),
    'FRN'       : Magnetometer(name='FRN'     ,csys='GEO',ctype='sph',coords=(1.,  37.0913, -119.7193) ),
    'FUR'       : Magnetometer(name='FUR'     ,csys='GEO',ctype='sph',coords=(1.,  48.17  ,  11.28   ) ),
    'OTT'       : Magnetometer(name='OTT'     ,csys='GEO',ctype='sph',coords=(1.,  45.40  ,  284.45  ) ),#!! from paper, todo: check acuracy
    'FRD'       : Magnetometer(name='FRD'     ,csys='GEO',ctype='sph',coords=(1.,  38.20  ,  282.63  ) ),#!! from paper, todo: check acuracy

    'gridpnt1'  : Magnetometer(name='gridpnt1',csys='GEO',ctype='sph',coords=(1., 18*(175./174.), 72.) ),
    'gridpnt2'  : Magnetometer(name='gridpnt2',csys='GEO',ctype='sph',coords=(1., 18*(175./174.), 72.+180.) ),
    'GMpoint1'  : Magnetometer(name='GMpoint1',csys='GSM',ctype='car',coords=(0.71875,0.09375,-3.71875) ),
    'GMpoint2'  : Magnetometer(name='GMpoint2',csys='GSM',ctype='car',coords=(-146.,-14.,-14.) ),
    'GMpoint3'  : Magnetometer(name='GMpoint3',csys='GSM',ctype='car',coords=(-156.,-124.,-124.) ),
    'GMpoint4'  : Magnetometer(name='GMpoint4',csys='GSM',ctype='car',coords=(0.71875,0.09375,-2.21875) ),
    'GMpoint5'  : Magnetometer(name='GMpoint5',csys='GSM',ctype='car',coords=(0.71875,0.09375,-1.71875) ),
    'GMpoint6'  : Magnetometer(name='GMpoint6',csys='GSM',ctype='car',coords=(0.,0.,-10.) ),
    'GMpoint7'  : Magnetometer(name='GMpoint7',csys='GSM',ctype='car',coords=(0.,0.,-6.) ),
    'GMpoint8'  : Magnetometer(name='GMpoint8',csys='GSM',ctype='car',coords=(0.,0.,-5.) ),
    'GMpoint9'  : Magnetometer(name='GMpoint9',csys='GSM',ctype='car',coords=(6.,0.,0.) ),
    'GMpoint10' : Magnetometer(name='GMpoint10',csys='GSM',ctype='car',coords=(4.,0.,-4.) ),

    'zerolat'   : Magnetometer(name='zerolat',csys='SM',ctype='sph',coords=(2., 0., 33.) ),
    'zerolon'   : Magnetometer(name='zerolon',csys='SM',ctype='sph',coords=(2., 33., 0.) ),
    'zeroboth'  : Magnetometer(name='zeroboth',csys='SM',ctype='sph',coords=(2., 0., 0.) ),
    'otherside' : Magnetometer(name='otherside',csys='SM',ctype='sph',coords=(2., 0., 180.) ),
}
