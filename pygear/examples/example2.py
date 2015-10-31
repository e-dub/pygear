"""
Example problem 2 for pyGear

changelog:
    2011/07/30  version 0.10    jdamerau@arcor.de
    2011/12/01  version 0.11    jdamerau@arcor.de

license:
This code is published under the terms of the GNU General Public License v3
http://www.gnu.org/licenses/gpl-3.0.html

project home:
https://sourceforge.net/projects/pygear/

The code has been tested with Python 2.6, http://www.python.org
required modules (tested with version):
1. pythonOCC (0.5), http://www.pythonocc.org
2. wxPython (2.8), http://www.wxpython.org
3. NumPy (1.6.6), http://numpy.scipy.org
4. SciPy (0.10.0rc1), http://numpy.scipy.org
5. pyGear (0.21), https://sourceforge.net/projects/pygear/
"""

from __future__ import print_function

# from scipy import *
from pygear import *
from pygear import PythonOCCArrayToNumPyArray
import matplotlib.pyplot as plt
from example_data import *

# EXAMPLE 2: CREATE GEAR PAIR FROM GEAR DATA, EXPORT STEP-FILES,
#            CALCULATE INERTIA-TENSORS AND DISPLAY ONE OF THE GEARS

geardata_1 = extgear_1  # select gear data for 1st gear here
geardata_2 = extgear_2  # select gear data for 2nd gear here
mygear_1 = CylindricalGearWheel(extgear_1)  # create cylindrical gear wheel instance for 1st gear
mygear_2 = CylindricalGearWheel(extgear_2)  # create cylindrical gear wheel instance for 2nd gear
mypair = CylindricalGearPair({'a': 243}, mygear_1, mygear_2)  # create gear pair instance (two gears that can mate are required)
print(mypair)  # print gear pair data
mygear_1_solid = mygear_1.makeOCCSolid()  # create OCC-3d-solid of gear wheel 1
mygear_2_solid = mygear_2.makeOCCSolid()  # create OCC-3d-solid of gear wheel 2
# writeOCCShape(mygear_1_solid, 'mygear_1.stp', 'step')  # write 3d-solid of 1st gear to STEP-file (in working directory)
# writeOCCShape(mygear_2_solid, 'mygear_2.stp', 'step')  # write 3d-solid of 2nd gear to STEP-file
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
np_formcoords = PythonOCCArrayToNumPyArray(mygear_1.formcoords)

# plot tooth form coordinates
plt.plot(np_formcoords[:, 0], np_formcoords[:, 1], 'x-')

# display plot (close window to continue!)
plt.show()

# display solid
displayOCCShape([mygear_1_solid])
