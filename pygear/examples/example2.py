#!/usr/bin/python
# coding: utf-8

"""Example problem 2 for pyGear

Summary
-------
CREATE GEAR PAIR FROM GEAR DATA
EXPORT STEP-FILES,
CALCULATE INERTIA-TENSORS
DISPLAY ONE OF THE GEARS
"""

from __future__ import print_function

import matplotlib.pyplot as plt

import pygear.core
import pygear.examples.example_data


geardata_1 = pygear.examples.example_data.extgear_1  # select gear data for 1st gear here
geardata_2 = pygear.examples.example_data.extgear_2  # select gear data for 2nd gear here

# create cylindrical gear wheel instance for 1st gear
mygear_1 = pygear.core.CylindricalGearWheel(pygear.examples.example_data.extgear_1)

# create cylindrical gear wheel instance for 2nd gear
mygear_2 = pygear.core.CylindricalGearWheel(pygear.examples.example_data.extgear_2)

# create gear pair instance (two gears that can mate are required)
mypair = pygear.core.CylindricalGearPair({'a': 243}, mygear_1, mygear_2)

print(mypair)  # print gear pair data
mygear_1_solid = mygear_1.makeOCCSolid()  # create OCC-3d-solid of gear wheel 1
mygear_2_solid = mygear_2.makeOCCSolid()  # create OCC-3d-solid of gear wheel 2

# write 3d-solid of 1st gear to STEP-file (in working directory)
# write_occ_shape(mygear_1_solid, 'mygear_1.stp', 'step')

# write 3d-solid of 2nd gear to STEP-file
# write_occ_shape(mygear_2_solid, 'mygear_2.stp', 'step')

density = 7.85e-6  # define mass density of gears (steel in this example, take care of unit consistency)
[m_1, J_1, r_g_1] = mygear_1.calcInertiaProperties(density)  # calculate inertia properties of 1st gear
[m_2, J_2, r_g_2] = mygear_2.calcInertiaProperties(density)  # calculate inertia properties of 2nd gear

# print inertia properties of both gears in nice format
print('\nGear 1:')
print('mass:\t', m_1)
print('inertia tensor:')
print(J_1[0])
print(J_1[1])
print(J_1[2])
print('\nGear 2:')
print('mass:\t', m_2)
print('inertia tensor:')
print(J_2[0])
print(J_2[1])
print(J_2[2])

# convert tooth form coordinates to NumPy-array (for plotting)
np_formcoords = pygear.core.pythonocc_array_to_numpy_array(mygear_1.formcoords)

# plot tooth form coordinates
plt.plot(np_formcoords[:, 0], np_formcoords[:, 1], 'x-')

# display plot (close window to continue!)
plt.show()

# display solid
pygear.core.display_occ_shape([mygear_1_solid])
