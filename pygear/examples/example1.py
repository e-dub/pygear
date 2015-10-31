#!/usr/bin/python
# coding: utf-8

"""
Example problem 1 for pyGear

changelog:
    2011/07/20  version 0.10    jdamerau@arcor.de
    2011/12/01  version 0.11    jdamerau@arcor.de

license:
This code is published under the terms of the GNU General Public License v3
http://www.gnu.org/licenses/gpl-3.0.html

project home:
https://sourceforge.net/projects/pygear/

The code has been tested with Python 2.6, http://www.python.org
required modules (tested with version):
1. pyGear (0.21), https://sourceforge.net/projects/pygear/
"""

from __future__ import print_function

from pygear import *
from example_data import *

# EXAMPLE 1: CREATE GEAR GEOMETRY AND STEP-EXPORT FROM GEAR DATA
geardata = extgear_1  # select gear data here
mygear = CylindricalGearWheel(geardata)  # create cylindrical gear wheel instance
print(mygear)  # print gear data
mygear_solid = mygear.makeOCCSolid()                  # create OCC-3d-solid of gear wheel
# writeOCCShape(mygear_solid, 'extgear_1.stp', 'step')      # write 3d-solid of gear to STEP-file (in working directory)
# writeOCCShape(mygear_solid, 'extgear_1.igs', 'iges')      # write 3d-solid of gear to IGES-file (in working directory), IGES interface has bug!
# writeOCCShape(mygear_solid, 'extgear_1.wrl', 'vrml')      # write 3d-solid of gear to VRML-file (in working directory)
displayOCCShape(mygear_solid)
