#!/usr/bin/python
# coding: utf-8

"""Example problem 1 for pyGear

Summary
-------
Create gear geometry and a export STEP file

"""

from __future__ import print_function

import pygear.core
import pygear.examples.example_data

geardata = pygear.examples.example_data.extgear_1  # select gear data here
mygear = pygear.core.CylindricalGearWheel(geardata)  # create cylindrical gear wheel instance
print(mygear)  # print gear data
mygear_solid = mygear.makeOCCSolid()  # create OCC-3d-solid of gear wheel

# write 3d-solid of gear to STEP-file (in working directory)
pygear.core.writeOCCShape(mygear_solid, 'extgear_1.stp', 'step')

# write 3d-solid of gear to IGES-file (in working directory), IGES interface has bug!
pygear.core.writeOCCShape(mygear_solid, 'extgear_1.igs', 'iges')

# write 3d-solid of gear to VRML-file (in working directory)
pygear.core.writeOCCShape(mygear_solid, 'extgear_1.wrl', 'vrml')
pygear.core.displayOCCShape(mygear_solid)
