=================================================
pyGear
=================================================

project home:
https://sourceforge.net/projects/pygear/


What is pyGear?
--------------
PyGear is open-source software for generating gear geometry and calculating
gear properties. It is intended to be used as a pre-processor for Computer-
Aided Design and Computer-Aided Engineering.
It relies on industry standards in order to create the most
exact shapes possible. At the same time it is easy to use, as it let's you
generate geometries in many different ways, even with a minimum of input
data provided. The gear-geometry can be exported as 2d-coordinates or in
3d as STEP, IGES or VRML-file.

PyGear relies on Python and additionally the modules pythonOCC, wxPython,
NumPy and SciPy are required.


pyGear structure
----------------
Currently pyGear consists of the following files:

  README.txt
  setup.py
  pygear.py
  license.txt
  

License information
-------------------
The pyGear-code is published under the terms of the GNU General Public License v3
http://www.gnu.org/licenses/gpl-3.0.html
See the file "license.txt" for terms & conditions for usage, and a DISCLAIMER
OF ALL WARRANTIES.


ChangeLog
-------------------
0.12:   - gear-generator for cylindrical gear wheels (spur, helical, internal, external)
        - 2d-export (x-y-coordinates)
        - 3d-export (STEP, IGES, VRML) and display
        - basic gear pair calculations
0.14:   - tooth thickness allowance --> exact flank
        - addendum allowance
        - acceptance backlash calculation
        - tip chamfer enhancements
        - bug fixes
        - renamed classes:
            "SpurGearWheel" --> "CylindricalGearWheel"
            "SpurGearPair" --> "CylindricalGearPair"
0.20:   - gear hobbing simulation --> exact geometry (external gears only)
        - multi-step manufacturing by defnition of blank
        - new classes:
            "Tool"
            "ToothedRackTool"
            "Machine"
            "GearHobber"
            "Blank"
0.21:   - removed redundancies from code
        - modified to work with newer versions of pythonOCC, NumPy and SciPy
        - removed trochoidal root shape generation as it was not correct for helical gears
          --> use hobbing simulation for exact root fillet shape instead
0.22:	- novel algorithm for envelope computation in manufacturing simulation
        - manufacturing simulation enhanced: replaced linear interpolations for
          edge-points with exact algorithm (within numerical precision)
        - manufacturing simulation more precise and stable
        - bug fixes
0.23    - improved algorithm for envelope computation in manufacturing simulation
        - manufacturing simulation more accurate and robust
        - manufacturing simulation much faster now
        - manufacturing of gears with pointed tip works now
        - some convergence problems fixed
            

Installation
---------------------
make sure you have Python (http://www.python.org) installed.

make sure that the following packages are installed:
1. pythonOCC, http://www.pythonocc.org
2. wxPython, http://www.wxpython.org
3. NumPy, http://numpy.scipy.org
4. SciPy, http://numpy.scipy.org
5. matplotlib, http://matplotlib.sourceforge.net (required for example problems only)

then copy the file pygear.py to 
<python path>\Lib\site-packages

The code has been tested with Python 2.6, http://www.python.org
required modules (tested with version):
1. pythonOCC (0.5), http://www.pythonocc.org
2. wxPython (2.8), http://www.wxpython.org
3. NumPy (1.6.6), http://numpy.scipy.org
4. SciPy (0.10.0rc1), http://numpy.scipy.org
5. matplotlib (1.01), http://matplotlib.sourceforge.net (for example problems)
it is strongly recommended to use the same versions of the modules as pyGear likely won't work
with older versions.


Documentation
------------------
sorry, not yet. See the comments in the source code of file pygear.py for some documentation.

            
Examples
------------------
There are some example problems:

  example_data.py
  example1.py
  example2.py
  example3.py