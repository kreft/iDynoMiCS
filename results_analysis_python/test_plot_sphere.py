#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import numpy
import toolbox_plotting
import toolbox_schematic_new as toolbox_schematic


fig = toolbox_plotting.SlideFigure()
axis = fig.add_subplot('', 111, frameon=True, projection='3d')

#X, Y = numpy.meshgrid([0, 1], [0, 1])
#axis.plot_surface(X, Y, 0.3, color='k', zorder=1)
Y, Z = numpy.meshgrid([0, 1], [0, 0.5])
axis.plot_surface(1.0, Y, Z, color='k', zorder=1)

sphere = toolbox_schematic.Sphere()
sphere.set_defaults(edgecolor='none', facecolor='r', zorder=0)
sphere.set_points((0.5, 0.5, 0.5), 0.5)
sphere.draw(axis)

axis.set_xlabel('x')
axis.set_ylabel('y')
axis.set_zlabel('z')
