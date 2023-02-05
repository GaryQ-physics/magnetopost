import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
import cutplane_plot as cp
import cxtransform as cx

time = [2003, 11, 20, 7, 0, 0, 0]
cp.plot(time, 'p', 'xz', logz=True, dx=0.5, dy=0.5,
        xlims=[-30, 30], ylims=[-20, 20])

# Manually colobar limits 
cp.plot(time, 'p', 'xz', logz=True, dx=0.5, dy=0.5,
        xlims=[-30, 30], ylims=[-20, 20], zlims=[0, 100])

if False:    
    cp.plot(time, 'p', 'xy', dx=0.5, dy=0.5,
            xlims=[-30, 30], ylims=[-20, 20], zlims=[0, 10])
    cp.plot(time, 'p', 'yz', dx=0.5, dy=0.5,
            xlims=[-30, 30], ylims=[-20, 20], zlims=[0, 10])

if False:
    mag = [1.01, 57.50, 176.] # r, lat, long in MAG coordinate system
    cp.plot(time, mag, None, 'p', dx=0.5, dy=0.5,
            xlims=[-30, 30], ylims=[-20, 20], png=False)

    # Convert GSM position to MAG position
    mag = cx.GSMtoMAG([-1, 0, 0.7], time, 'car', 'sph')
    cp.plot(time, mag, None, 'p', dx=0.5, dy=0.5,
            xlims=[-30, 30], ylims=[-20, 20], png=False)

if False:
    cp.plot(time, 'p', mag, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])
    cp.plot(time, 'p', [[1, 0, 0], [0, 0, 1]], field_lines=mag, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])
    cp.plot(time, 'p', 'xz', field_lines=mag, 
            dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20], debug=True)

if False:
    dtheta = 0.1
    import numpy as np
    angles = (np.pi/2.)*np.arange(-1, 1+dtheta, dtheta)
    
    field_lines = []
    r = 1.
    for angle in angles:
        field_lines.append([r*np.sin(angle), 0., r*np.cos(angle)])
        field_lines[-1] = cx.GSMtoMAG(field_lines[-1], time, 'car', 'sph')

    cp.plot(time, 'p', 'xz', field_lines=field_lines, dx=0.5, dy=0.5, xlims=[-30, 30], ylims=[-20, 20])    
