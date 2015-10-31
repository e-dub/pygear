"""
Example problem 3 for pyGear

changelog:
    2011/07/20  version 0.10    jdamerau@arcor.de
    2012/03/24  version 0.11    jdamerau@arcor.de
    2012/04/30  version 0.12    jdamerau@arcor.de

license:
This code is published under the terms of the GNU General Public License v3
http://www.gnu.org/licenses/gpl-3.0.html

project home:
https://sourceforge.net/projects/pygear/

The code has been tested with Python 2.6, http://www.python.org
required modules (tested with version):
1. pyGear (0.22), https://sourceforge.net/projects/pygear/
2. matplotlib (1.01), http://matplotlib.sourceforge.net/
3. NumPy (1.4.0-rc1), http://numpy.scipy.org
4. SciPy (0.7.1), http://numpy.scipy.org
"""

from __future__ import print_function

from __future__ import division
from scipy import *
# from numpy.linalg import norm
from pygear import *
from pygear import PythonOCCArrayToNumPyArray
import matplotlib.pyplot as plt
from example_data import *

# EXAMPLE 3: CREATE GEAR THROUGH MANUFACTURING SIMULATION (HOBBING)

tooldata = prottool  # select tool data here
machinedata = hobber  # select machine data here
mytool = ToothedRackTool(stantool)  # create tool instance
myblank = Blank()  # create blank instance
myblank.setCircular(125.0, pi / 2, pi / 2 * (1 + 1 / hobber.get('z')), 100)  # set blank to be arc segment (half a tooth is sufficient, must begin at pi/2 for external gears) with 100 points. This is necessary in most cases as the tool does usually not cut the tip of the tooth.
print(myblank)  # print blank data
mymachine = GearHobber(machinedata, mytool, myblank)  # create (manufacturing) machine instance
print(mymachine)  # print machine data

# The following section demonstrates the visualization of the hobbing process
# (this is just for demonstration and not the original purpose of pygear)
tool_param_limit = mytool._getParameterLimit()  # get limits of curve parameters of tool (that is internally represented by a parametrized and arc-length normalized curve)
tool_points = 201  # set the resolution of the tool profiles in the plot
tool = zeros([tool_points, 2])  # allocate array for 2D-points of tool profile
param_step = tool_param_limit / (tool_points - 1) * 2  # calculate step size for curve parameter
angle_points = 61  # set the number of angular positions of the cutter to be displayed
angle_step = 0.8 / (angle_points / 2 - 1)  # calculate angular distance of positions
plt.figure(1)  # create new figure
plt.subplot(121, aspect='equal')  # left subplot is active
for angle_index in range(int(-(angle_points - 1) / 2), int((angle_points - 1) / 2 + 1)):  # loop over all angular positions
    phi = angle_index * angle_step  # calculate actual angle of cutter
    for point_index in range(0, tool_points):  # loop over all points on tool's profile
        point_toolCS = mymachine._getToolPoint((point_index - (tool_points - 1) / 2) * param_step)
        tool[point_index,:] = mymachine._transformToMachineCoords(phi, point_toolCS, elemtype='point')  # use internal method of machine for evaluation of tool profile point
    plt.plot(tool[:, 0], tool[:, 1], color=str(0.5 * phi + 0.5))  # plot actual cutter position
    plt.hold(True)  # allow multiple curves in same plot

# Execute manufacturing simulation
plt.subplot(122, aspect='equal')  # right subplot is active
[formcoords, geardata] = mymachine.createToothShape()  # simulate gear manufactoring to create tooth shape (blank is accounted for)
np_formcoords = PythonOCCArrayToNumPyArray(formcoords)  # convert tooth form coordinates to NumPy-array (for plotting)
plt.plot(np_formcoords[:, 0], np_formcoords[:, 1], 'bx-')  # plot tooth form coordinates
print(geardata)  # print gear data
plt.show()  # display plot (close window to continue!)
geardata.update({'b': 30.0, 'd_s': 70.0})  # add tooth width (mandatory) and shaft diameter to gear data
mygear=CylindricalGearWheel(geardata, None, formcoords)  # create cylindrical gear wheel instance
print(mygear)  # print data of gear wheel
mygear_solid = mygear.makeOCCSolid()  # create OCC-3d-solid of gear wheel
# writeOCCShape(mygear_solid, 'hobbed_gear.stp', 'step')      # write 3d-solid of gear to STEP-file (in working directory)
# writeOCCShape(mygear_solid, 'hobbed_gear.igs', 'iges')      # write 3d-solid of gear to IGES-file (in working directory), IGES interface has bug!
# writeOCCShape(mygear_solid, 'hobbed_gear.wrl', 'vrml')      # write 3d-solid of gear to VRML-file (in working directory)
displayOCCShape([mygear_solid])  # display solid
