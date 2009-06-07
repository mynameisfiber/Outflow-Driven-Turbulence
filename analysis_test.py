#!/usr/bin/env python

"""This script demonstrates how one can script MayaVi and use its
contour related modules.  Notice the magic line at the top.
"""
# Author: Prabhu Ramachandran <prabhu_r@users.sf.net>
# Copyright (c) 2005-2008, Enthought, Inc.
# License: BSD Style.

# Standard library imports
from os.path import join, abspath
import numpy

# Enthought library imports
from enthought.mayavi.scripts import mayavi2
from enthought.mayavi.sources.array_source import ArraySource
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
from enthought.mayavi.modules.outline import Outline
from enthought.mayavi.modules.grid_plane import GridPlane
from enthought.mayavi.modules.contour_grid_plane import ContourGridPlane
from enthought.mayavi.modules.iso_surface import IsoSurface
from enthought.mayavi.modules.scalar_cut_plane import ScalarCutPlane

@mayavi2.standalone                        
def contour(i,x):
    """The script itself.  We needn't have defined a function but
    having a function makes this more reusable.
    """
    # 'mayavi' is always defined on the interpreter.
    # Create a new scene.
    mayavi.new_scene()

    # Read in datacube
    fd = open("output/output-cube-%.8d-%.3d"%(i,0))
    rho, rhovx, rhovz, rhovy = numpy.fromfile(file=fd,dtype=numpy.float32).reshape((4,x,x,x))
    src = ArraySource(transpose_input_array=False)
    src.scalar_data = numpy.log10(rho).T.copy()
    mayavi.add_source(src)

    # Create an outline for the data.
    o = Outline()
    mayavi.add_module(o)

    # Create one ContourGridPlane normal to the 'x' axis.
    cgp = ContourGridPlane()
    mayavi.add_module(cgp)
    # Set the position to the middle of the data.
    cgp.grid_plane.position = 0

    # Another normal to 'y' axis.
    cgp = ContourGridPlane()
    mayavi.add_module(cgp)
    # Set the axis and position to the middle of the data.
    cgp.grid_plane.axis = 'y'
    cgp.grid_plane.position = 0
    
    # Another normal to 'z' axis.
    cgp = ContourGridPlane()
    mayavi.add_module(cgp)
    # Set the axis and position to the middle of the data.
    cgp.grid_plane.axis = 'z'
    cgp.grid_plane.position = 0

    # An isosurface module.
    iso = IsoSurface(compute_normals=False)
    mayavi.add_module(iso)
    iso.contour.contours = [220.0]

    # An interactive scalar cut plane.
    cp = ScalarCutPlane()
    mayavi.add_module(cp)
    cp.implicit_plane.normal = 0,0,1


if __name__ == '__main__':
    contour(1086,200)
