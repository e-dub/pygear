#!/usr/bin/python
# coding: utf-8

"""Example problem 3 for pyGear

Summary
-------
Create gear through manufacturing simulation (hobbing)

"""

from __future__ import print_function, division

import scipy

import pygear.core
import matplotlib.pyplot as plt
import pygear.examples.example_data


tooldata = pygear.examples.example_data.prottool  # select tool data here
machinedata = pygear.examples.example_data.hobber  # select machine data here
mytool = pygear.core.ToothedRackTool(pygear.examples.example_data.stantool)  # create tool instance
myblank = pygear.core.Blank()  # create blank instance

# set blank to be arc segment (half a tooth is sufficient, must begin at pi/2 for external gears) with 100 points.
# This is necessary in most cases as the tool does usually not cut the tip of the tooth.
myblank.set_circular(125.0, scipy.pi / 2, scipy.pi / 2 * (1 + 1 / pygear.examples.example_data.hobber.get('z')), 100)

print(myblank)  # print blank data
mymachine = pygear.core.GearHobber(machinedata, mytool, myblank)  # create (manufacturing) machine instance
print(mymachine)  # print machine data

# The following section demonstrates the visualization of the hobbing process
# (this is just for demonstration and not the original purpose of pygear)

# get limits of curve parameters of tool (that is internally represented
# by a parametrized and arc-length normalized curve)
tool_param_limit = mytool._get_parameter_limit()

tool_points = 201  # set the resolution of the tool profiles in the plot
tool = scipy.zeros([tool_points, 2])  # allocate array for 2D-points of tool profile
param_step = tool_param_limit / (tool_points - 1) * 2  # calculate step size for curve parameter
angle_points = 61  # set the number of angular positions of the cutter to be displayed
angle_step = 0.8 / (angle_points / 2 - 1)  # calculate angular distance of positions
plt.figure(1)  # create new figure
plt.subplot(121, aspect='equal')  # left subplot is active

# loop over all angular positions
for angle_index in range(int(-(angle_points - 1) / 2), int((angle_points - 1) / 2 + 1)):
    phi = angle_index * angle_step  # calculate actual angle of cutter
    for point_index in range(0, tool_points):  # loop over all points on tool's profile
        point_toolCS = mymachine._getToolPoint((point_index - (tool_points - 1) / 2) * param_step)

        # use internal method of machine for evaluation of tool profile point
        tool[point_index, :] = mymachine._transformToMachineCoords(phi, point_toolCS, elemtype='point')
    plt.plot(tool[:, 0], tool[:, 1], color=str(0.5 * phi + 0.5))  # plot actual cutter position
    plt.hold(True)  # allow multiple curves in same plot

# Execute manufacturing simulation
plt.subplot(122, aspect='equal')  # right subplot is active

# simulate gear manufacturing to create tooth shape (blank is accounted for)
[formcoords, geardata] = mymachine.create_tooth_shape()

# convert tooth form coordinates to NumPy-array (for plotting)
np_formcoords = pygear.core.pythonocc_array_to_numpy_array(formcoords)

plt.plot(np_formcoords[:, 0], np_formcoords[:, 1], 'bx-')  # plot tooth form coordinates
print(geardata)  # print gear data
plt.show()  # display plot (close window to continue!)
geardata.update({'b': 30.0, 'd_s': 70.0})  # add tooth width (mandatory) and shaft diameter to gear data
mygear = pygear.core.CylindricalGearWheel(geardata, None, formcoords)  # create cylindrical gear wheel instance
print(mygear)  # print data of gear wheel
mygear_solid = mygear.makeOCCSolid()  # create OCC-3d-solid of gear wheel
# write_occ_shape(mygear_solid, 'hobbed_gear.stp', 'step')      # write 3d-solid of gear to STEP-file (in working directory)
# write_occ_shape(mygear_solid, 'hobbed_gear.igs', 'iges')      # write 3d-solid of gear to IGES-file (in working directory), IGES interface has bug!
# write_occ_shape(mygear_solid, 'hobbed_gear.wrl', 'vrml')      # write 3d-solid of gear to VRML-file (in working directory)
pygear.core.display_occ_shape([mygear_solid])  # display solid
