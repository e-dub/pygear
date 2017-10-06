"""
Classes representing involute gear wheels

changelog:
    2009/12/27  version 0.10    jdamerau@arcor.de
    2010/10/10  version 0.11    jdamerau@arcor.de
    2010/10/28  version 0.12    jdamerau@arcor.de
    2011/01/19  version 0.13    jdamerau@arcor.de
    2011/05/30	version 0.14	jdamerau@arcor.de
    2011/06/01	version 0.20	jdamerau@arcor.de
    2011/12/01  version 0.21    jdamerau@arcor.de
    2012/03/24  version 0.22    jdamerau@arcor.de
    2012/04/30  version 0.23    jochen.damerau@gmail.com
    2015/11/18  version 0.24    jochen.damerau@gmail.com

license:
This code is published under the terms of the GNU General Public License v3
http://www.gnu.org/licenses/gpl-3.0.html

project home:
https://sourceforge.net/projects/pygear/

The code has been tested with Python 2.7, http://www.python.org
required modules (tested with version):
1. pythonOCC (0.16), http://www.pythonocc.org
2. wxPython (3.0.2), http://www.wxpython.org
3. NumPy (1.9.2), http://numpy.scipy.org
4. SciPy (0.16.0), http://numpy.scipy.org
"""

# Imports
from __future__ import division
from math import sin, asin, cos, acos, tan, atan, pi, degrees, radians, sqrt
from copy import *
from sys import maxint
from warnings import *
from OCC.TColgp import TColgp_Array1OfPnt2d, TColgp_Array1OfPnt
from OCC.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCC.Geom2d import Geom2d_Circle, Geom2d_Line, Geom2d_CartesianPoint, Geom2d_OffsetCurve, \
                       Geom2d_BSplineCurve, Geom2d_Curve, Handle_Geom2d_BSplineCurve
from OCC.Geom import Geom_Plane, Geom_Circle, Geom_Line, Geom_BSplineCurve, Geom_OffsetCurve
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_InterCurveCurve, Geom2dAPI_ProjectPointOnCurve
from OCC.GeomAPI import GeomAPI_PointsToBSpline, GeomAPI_ProjectPointOnCurve
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Dir, gp_Dir2d, gp_Ax2, gp_Ax22d, gp_Circ,   \
                   gp_Pln, gp_Vec, gp_Vec2d, gp_Ax1, gp_Trsf, gp_XYZ
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeEdge,  \
                               BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakePolygon,  \
                               BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing,        \
                               BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Transform
from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakePipeShell
from OCC.TopoDS import TopoDS_Wire
from OCC.TopoDS import topods as TopoDS_Cast
from OCC.Geom2dGcc import Geom2dGcc_QualifiedCurve
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GccEnt import *
from OCC.Geom2dGcc import Geom2dGcc_Circ2d2TanRad
from OCC.GCPnts import GCPnts_UniformAbscissa
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.TopAbs import *
from OCC.ShapeFix import *
from OCC.ShapeExtend import ShapeExtend_WireData
from OCC.BRepCheck import *
from OCC.TopExp import TopExp_Explorer
from OCC.BRep import BRep_Tool
from OCC.GProp import GProp_GProps
from OCC.BRepGProp import brepgprop
from scipy import *
from numpy.linalg import svd, det
from scipy.optimize import newton, fsolve, fmin_cg
from scipy.linalg import norm
from scipy.interpolate import InterpolatedUnivariateSpline

__version__ = "0.24"
__all__ = ['GearWheel','CylindricalGearWheel','GearPair','CylindricalGearPair','writeCoords',  \
           'displayOCCShape','writeOCCShape','Tool','ToothedRackTool','Machine','GearHobber',  \
           'Blank']


# Functions needed by classes
def inv(angle):
    """
    Involute of a circle

    INPUT parameter:
    angle : angle of tangent to base circle (numeric)[degrees]

    OUTPUT:
    value of involute-function at angle (numeric)
    """
    
    return tan(radians(angle))-radians(angle)


def sign(number):
    """
    Sign of a number

    INPUT parameter:
    number : number that's sign is to be calculated (numeric)

    OUTPUT:
    value of sign(number) (-1, 0, 1)
    """

    if number>0.0:
         return 1
    elif number<0.0:
         return -1
    else:
         return 0


def NumPyArrayToPythonOCCArray(numpy_array):
    """
    transform array from nparray (NumPy) to TColgp_Array1OfPnt2d (pythonOCC) format

    INPUT parameter:
    numpy_array : array containing points in rows and coordinates in columns (Nx2-nparray, NumPy)

    OUTPUT:
    pythonOCC_array: array containing points in rows and coordinates in columns (TColgp_Array1OfPnt2d, pythonOCC)
    """
    
    # create arrays of points holding coordinates
    pythonOCC_array = TColgp_Array1OfPnt2d(1, size(numpy_array, 0))
    # set entries
    for index in range(0, size(numpy_array, 0)):
         pythonOCC_array.SetValue(index+1, gp_Pnt2d(numpy_array[index,0], numpy_array[index,1]))

    return pythonOCC_array


def PythonOCCArrayToNumPyArray(pythonOCC_array):
    """
    transform array from TColgp_Array1OfPnt2d (pythonOCC) to nparray (NumPy) format

    INPUT parameter:
    pythonOCC_array: array containing points in rows and coordinates in columns (TColgp_Array1OfPnt2d, pythonOCC)

    OUTPUT:
    numpy_array : array containing points in rows and coordinates in columns (Nx2-nparray, NumPy)
    """
    
    # create arrays of points holding coordinates
    numpy_array = zeros([pythonOCC_array.Length(), 2])
    # set entries
    for index in range(1, pythonOCC_array.Length()+1):
         numpy_array[index-1,:] = array([pythonOCC_array.Value(index).X(), pythonOCC_array.Value(index).Y()])

    return numpy_array


def CartesianCoordinatesToPolarCoordinates(x, y):
    """
    convert a tuple that represents a vector in cartesian coordinates to polar coordinates.
    The zero angle of the polar representation corresponds to x-direction of cartesian representation.

    INPUT parameters:
    x: x-value of vector in cartesian coordinates (numeric)
    y: y-value of vector in cartesian coordinates (numeric)
    
    OUTPUT:
    r:      radial coordinate of vector in polar coordinates (numeric)
    phi:    anglular coordinate of vector in polar coordinates (numeric)[radians] 
    """

    r = sqrt(x**2+y**2)
    if x>0.0:   # arctangent is not unique
         phi = atan(y/x)
    elif x<0.0:
         phi = atan(y/x)+pi*sign(y)
    else:   # absolute value of x/y is infinity
         if y>0.0:
              phi = pi/2
         elif y<0.0:
              phi = -pi/2
         else:
              phi = 0.0     # this is arbitrary

    return r, phi


def PolarCoordinatesToCartesianCoordinates(r, phi):
    """
    convert a tuple that represents a vector in polar coordinates to cartesian coordinates.
    The zero angle of the polar representation corresponds to x-direction of cartesian representation.

    INPUT parameters:
    r:      radial coordinate of vector in polar coordinates (numeric)
    phi:    anglular coordinate of vector in polar coordinates (numeric)[radians] 
    
    OUTPUT:
    x: x-value of vector in cartesian coordinates (numeric)
    y: y-value of vector in cartesian coordinates (numeric)
    """

    x = r*cos(phi)
    y = r*sin(phi)

    return x, y
    

def writeCoords(coords, outfile):
    """
    Output python list formatted as array (rows and columns)

    INPUT parameter:
    coords  : array containing points in rows and coordinates in columns (TColgp_Array1OfPnt2d, pythonOCC)
    outfile : filename to output array (string)
    """

    fd = file(outfile, 'w')
    # upper and lower index of point-array
    upper_index = coords.Upper()
    lower_index = coords.Lower()
    for index in range(lower_index, upper_index+1):
         fd.write(str(coords.Value(index).X()).ljust(20)+'\t'+  \
                  str(coords.Value(index).Y()).ljust(20)+'\n')
    fd.close()


def displayOCCShape(shape):
    """
    Display OpenCascade-Shape-object.

    INPUT parameter:
    shape : to be displayed (OpenCascade-Shape-object)
    """

    from OCC.Display.SimpleGui import init_display
    display, start_display, add_menu, add_function_to_menu = init_display()
    
    display.DisplayShape(shape, update=True)
    start_display()
    

def writeOCCShape(shape, outfile, outformat='step'):
    """
    Convert OpenCascade-Shape-object and save to STEP-file.

    INPUT parameter:
    shape     : to be saved to file (OpenCascade-Shape-object)
    outfile   : full path and filename (string)
    outformat : output format (string)
                optional - one of the following: 'iges', 'step', 'stl' or 'vrml'
    """
    
    from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCC.IGESControl import IGESControl_Writer
    from OCC.StlAPI import StlAPI_Writer
    #from OCC.VrmlAPI import VrmlAPI_Writer

    if outformat=='step':
        export = STEPControl_Writer()
        export.Transfer(shape,STEPControl_AsIs)
        export.Write(outfile)
        
    elif outformat=='iges':  # still buggy
        export = IGESControl_Writer('mm',1)
        export.AddShape(shape)
        export.Write(outfile)
        
    elif outformat=='stl':  # This works fine now
        stl_ascii_format = True
        export = StlAPI_Writer()
        export.Write(shape, outfile, stl_ascii_format)
        #export.write_file()

    elif outformat=='vrml':  # No implementation
        pass
        #export = VrmlAPI_Writer()
        #export.Write(shape, outfile)
        
    else:
        raise ValueError, 'output format not known (allowed values: \'iges\', \'step\', \'stl\' or \'vrml\')'  
    

# Classes
class GearWheel:
    """
    Parent Class for all gear wheels.
    """

    # Attributes: settings for geometry construction
    points_flank   = 20  # points along the involute from root form to tip form circle
    points_fillet  = 15  # points on fillet from root to root form circle
    points_tip     = 5   # points along tip circle (half tooth)
    points_root    = 5   # points along root circle from center of gap to
                         # beginning of fillet (half tooth)
    points_shaft   = 20  # points on inner radius of gearwheel (shaft diameter)
    points_chamfer = 4   # points on tip chamfer
    points_ext     = 8   # points on extension of involute beyond base circle (if applicable)(should be at least 8)
    points_width   = 10  # resolution along width

    # Attributes: default settings for parameters
    _x_default       =  0.0    # default value for addendum modification
    _alpha_n_default = 20.0    # default value for pressure angle (DIN 867)
    _beta_default    =  0.0    # default value for helix angle (spur gear)
    _rho_f_default   =  0.38   # default value for fillet radius (divided bei module)(DIN 867)
    _c_default       =  0.25   # default value for tip clearance (divided bei module)(DIN 867)
    _k_default       =  0.0    # default value for tip height modification
    _d_s_default     =  0.0    # default value for shaft diameter (inner gear wheel diameter)
    _h_k_default     =  0.0    # default value for radial value of tip chamfer
    _tol_default     =  1e-6   # default tolerance for comparisons
    _A_s_default     =  0.0    # default value for tooth thickness allowance (DIN 3967)

    # Attributes: gear data
    data          = None   # dictionary containing all gear parameters for macro-geometry
    modifications = None   # dictionary containing all gear parameters for micro-geometry (flank profile modifications)
    formcoords    = None   # list of 2D-coordinates of points of half a tooth profile (TColgp_Array1OfPnt2d, pythonOCC)
    _formwire     = None   # wire of half a tooth profile (TopoDS_Wire, pythonOCC)


    def setResolution(self, curvename, value):
        """
        Set resolution for tooth form representation

        INPUT parameters:
        curvename : segment of tooth flank (string)
                    one of the following: flank, fillet, tip, root, shaft, width
        value     : new value for number of points to represent segment
        """
        
        if curvename=='flank':
             self.points_flank = value
        elif curvename=='fillet':
             self.points_fillet = value
        elif curvename=='tip':
             self.points_tip = value
        elif curvename=='root':
             self.points_root = value
        elif curvename=='shaft':
             self.points_shaft = value
        elif curvename=='width':
             self.points_width = value


    def getResolution(self, curvename):
        """
        Get resolution for tooth form representation

        INPUT parameters:
        curvename : segment of tooth flank (string)
                    one of the following: flank, fillet, tip, root, shaft, width

        OUTPUT:
        number of points used to represent requested segment
        """
        
        if curvename=='flank':
             return self.points_flank
        elif curvename=='fillet':
             return self.points_fillet
        elif curvename=='tip':
             return self.points_tip
        elif curvename=='root':
             return self.points_root
        elif curvename=='shaft':
             return self.points_shaft
        elif curvename=='width':
             return self.points_width


    def setFormCoords(self, formcoords, formwire):
        """
        Set tooth form coordinates

        INPUT parameters:
        coords : list of 2d-coordinate points (TColgp_Array1OfPnt2d, pythonOCC)
        wire   : wire describing the tooth form (TopoDS_Wire, pythonOCC)
        """
        
        # form coordinates: type and value check (at least two points for defining a
        # tooth form (straight flanks) and two coordinates per point)
        if not isinstance(formcoords, TColgp_Array1OfPnt2d):
             raise TypeError, 'instance of TColgp_Array1OfPnt2d expected'
        if formcoords.Length()<2:
             raise TypeError, 'too few points for tooth form'

        if formwire and not isinstance(formwire, TopoDS_Wire):
             raise TypeError, 'instance of TopoDS_Wire expected'

        self.formcoords = formcoords
        self._formwire  = formwire


    def getFormCoords(self):
        """
        Get tooth form coordinates

        INPUT parameters:
        -

        OUTPUT:
        tooth form coordinates (TColgp_Array1OfPnt2d, pythonOCC)
        tooth form (TopoDS_Wire, pythonOCC)
        """

        return self.formcoords, self._formwire
             

    def _makeUnique(self, coords):
        """
        Remove redundant entries from coordinate array

        INPUT parameter:
        coords : list of 2d-coordinate points (TColgp_Array1OfPnt2d, pythonOCC)

        OUTPUT:
        unique_coords : list of unique coordinates (TColgp_Array1OfPnt2d, pythonOCC)
        """
        
        # tolerance for comparisons
        tol = self._tol_default*self.data.get('m_n')

        # upper and lower index of point-array
        upper_index = coords.Upper()
        lower_index = coords.Lower()
    
        # remove redundant entries
        uniques = list()
        for index in range(lower_index, upper_index+1):
             unique = True
             for unique_point in uniques:
                  if abs(coords.Value(index).X()-unique_point[0])<tol and  \
                     abs(coords.Value(index).Y()-unique_point[1])<tol:
                       unique = False
             if unique:
                  uniques.append([coords.Value(index).X(), coords.Value(index).Y()])

        # copy list entries into coordinate array
        length_uniques = len(uniques)
        unique_coords  = TColgp_Array1OfPnt2d(lower_index, lower_index+length_uniques-1)
        for index in range(lower_index, lower_index+length_uniques):
             if abs(uniques[index-1][0])>tol:
                  unique_x = uniques[index-1][0]
             else:
                  unique_x = 0.0
             if abs(uniques[index-1][1])>tol:
                  unique_y = uniques[index-1][1]
             else:
                  unique_y = 0.0
             unique_coords.SetValue(index, gp_Pnt2d(unique_x, unique_y))
        return unique_coords


    def _make3dGeomFrom2dGeom(self, geom2d_object, z_coord_of_plane=0.0):
        """
        Make an OpenCascade 3d-geometry-object from a 2d-geometry-object by projecting
        it on a plane parallel to the x-y-plane of the global cartesian coordinate system
        of the 3d-space. Returns a native 3d-object built from the scratch (the corresponding
        OCC operations seem to be buggy). Only for unbounded objects.

        INPUT parameters:
        geom2d_object    : the OpenCascade 2d-geometry-object to project. Can be one of
                           the following: gp_Pnt2d, Geom2d_Circle, Geom2d_Line, Geom2d_Curve
        z_coord_of_plane : the coordinate where the projection plane cuts the global
                           z-axis of 3d-space (numeric).

        OUTPUT:
        geom3d_object    : the result of the projection as 3d-geometry object. Can be one of
                           the following: gp_Pnt, Geom_Circle, Geom_Line, Geom_BSplineCurve
        """

        if not (type(z_coord_of_plane)==type(0.0) or type(z_coord_of_plane)==type(0)):
            raise TypeError, 'z-coordinate must be numeric'

        # convert point
        if isinstance(geom2d_object, gp_Pnt2d):
            geom3d_object = gp_Pnt(geom2d_object.X(), geom2d_object.Y(), z_coord_of_plane)

        # convert line    
        elif isinstance(geom2d_object, Geom2d_Line):
            location2d    = geom2d_object.Location()
            location3d    = gp_Pnt(location2d.X(), location2d.Y(), z_coord_of_plane)
            direction2d   = geom2d_object.Direction()
            direction3d   = gp_Dir(direction2d.X(), direction2d.Y(), z_coord_of_plane)
            geom3d_object = Geom_Line(location3d, direction3d)

        # convert circle    
        elif isinstance(geom2d_object, Geom2d_Circle):
            circle2d      = geom2d_object.Circ2d()
            location2d    = circle2d.Location()
            location3d    = gp_Pnt(location2d.X(), location2d.Y(), z_coord_of_plane)
            x_direction3d = gp_Dir(circle2d.XAxis().Direction().X(), circle2d.XAxis().Direction().Y(), 0.0)
            y_direction3d = gp_Dir(circle2d.YAxis().Direction().X(), circle2d.YAxis().Direction().Y(), 0.0)
            z_direction3d = gp_Dir(0.0, 0.0, x_direction3d.X()*y_direction3d.Y()-  \
                                   x_direction3d.Y()*y_direction3d.X())   # apply cross-product to conserve orientation
            axis3d        = gp_Ax2(location3d, z_direction3d, x_direction3d)
            radius        = circle2d.Radius()
            geom3d_object = Geom_Circle(axis3d, radius)

        # convert bspline          
        elif isinstance(geom2d_object, Geom2d_BSplineCurve):
            pole_number = geom2d_object.NbPoles()
            poles3d     = TColgp_Array1OfPnt(1, pole_number)
            for index in range(1, pole_number+1):
                pole2d = geom2d_object.Pole(index)
                poles3d.SetValue(index, self._make3dGeomFrom2dGeom(pole2d))
            knot_number    = geom2d_object.NbKnots()
            knots          = TColStd_Array1OfReal(1, knot_number)
            geom2d_object.Knots(knots)
            multiplicities = TColStd_Array1OfInteger(1, knot_number)
            geom2d_object.Multiplicities(multiplicities)
            degree         = geom2d_object.Degree()
            periodic       = geom2d_object.IsPeriodic()
            geom3d_object  = Geom_BSplineCurve(poles3d, knots, multiplicities, degree, periodic)

        # convert offset curve          
        elif isinstance(geom2d_object, Geom2d_OffsetCurve):
            basis_curve2d   = geom2d_object.BasisCurve()
            basis_bspline2d = Handle_Geom2d_BSplineCurve().DownCast(basis_curve2d).GetObject()
            basis_bspline3d = self._make3dGeomFrom2dGeom(basis_bspline2d).GetHandle()
            offset          = geom2d_object.Offset()
            reference_dir   = gp_Dir(0.0, 0.0, 1.0)
            geom3d_object   = Geom_OffsetCurve(basis_bspline3d, offset, reference_dir)
            
            
        else:
            raise TypeError, 'method cannot be applied to given object type'

        return geom3d_object


    def __str__(self):
        """
        Define string conversion of GearWheel objects

        INPUT parameter:
        -

        OUTPUT:
        string representation of class
        """
        
        outstr = 'gear wheel data:\n'
        # output gear data
        for date in self.data:
             outstr += date.ljust(10)+':\t'+str(self.data.get(date))+'\n'

        # output modification data
        if self.modifications:
             outstr += '\nflank modifications:\n'
             for date in self.modifications:
                  outstr += date.ljust(10)+':\t'+str(self.modifications.get(date))+'\n'

        # output tooth form coordinates
        if self.formcoords:
             # upper and lower index of point-array
             upper_index = self.formcoords.Upper()
             lower_index = self.formcoords.Lower()
             outstr += '\ntooth form coordinates:\n'
             for index in range(lower_index, upper_index+1):
                  outstr += str(self.formcoords.Value(index).X()).ljust(20)+'\t'+  \
                            str(self.formcoords.Value(index).Y()).ljust(20)+'\n'

        return outstr
    

    def __init__(self, geardata, flankmods=None, formcoords=None):
        """
        Initialization of GearWheel-object
        Should be overwritten in derived classes

        INPUT parameter:
        geardata   : data of gear wheel (dictionary)
        flankmods  : data of flank modifications (dictionary)
        formcoords : list of 2d-coordinate points (list, list(len=2), numeric)
        """
        
        self.data           = deepcopy(geardata)
        self.modifications  = deepcopy(flankmods)
        self.formcoords     = self.setFormCoords(formcoords, None)

    
    def getGearData(self):
        """
        Return data-attribute of class

        OUTPUT:
        data attribute of class (dictionary)
        """
        
        return self.data


    def setGearData(self, geardata):
        """
        Set data-attribute of class, overwrite current value

        INPUT parameter:
        geardata : dictionary, containing geometric data of gear
                   for content, see method __init__
        """
        
        self.__init__(geardata, self.modifications, self.formcoords)


    def updateGearData(self, geardata):
        """
        Set data-attribute of class, update current value

        INPUT parameter:
        geardata : dictionary, containing geometric data of gear
                   for content, see method __init__
        """
        
        tempdata = self.data.copy()
        tempdata.update(geardata)
        self.__init__(geardata, self.modifications, self.formcoords)

    def getFlankModifications(self):
        """
        Return modifications-attribute of class

        OUTPUT:
        data attribute of class (dictionary)
        """
        
        return self.modifications


    def setFlankModifications(self, flankmods):
        """
        Set modifications-attribute of class, overwrite current value

        INPUT parameter:
        flankmods : dictionary, containing flank modification data of gear
                    for content, see method __init__
        """
        
        self.__init__(self.data, flankmods, self.formcoords)


    def updateFlankModifications(self, flankmods):
        """
        Set modifications-attribute of class, update current value

        INPUT parameter:
        flankmods : dictionary, containing flank modification data of gear
                    for content, see method __init__
        """
        
        tempmods = self.modifications.copy()
        tempmods.update(flankmods)
        self.__init__(self.data, tempmods, self.formcoords)
    

class CylindricalGearWheel(GearWheel):
    """
    Class representing a spur wheel or a helical gear wheel. Applicable for external and internal gears.
    Derived from GearWheel-class
    """
    
    def _toothThickness(self, d_y):
        """
        Tooth thickness in transverse cross-section (chord-length)

        INPUT parameter:
        d_y : two times coordinate of tooth flank in radial direction
              (diameter of y-cylinder)

        OUTPUT:
        s_y  : chord length of tooth thickness at d_y (numeric)
        d_yc : cutting point of diameter through tooth center and chord (numeric)
        """

        # necessary due to numerical rounding errors
        if self.data.get('d')/d_y*cos(radians(self.data.get('alpha_t')))>1.0:
             alpha_yt = 0.0
        else:
             alpha_yt = degrees(acos(self.data.get('d')/d_y*  \
                        cos(radians(self.data.get('alpha_t')))))
        s_yt = d_y*((pi+4*self.data.get('x_E')*tan(radians(self.data.get('alpha_n'))))/  \
               2/self.data.get('z')+inv(self.data.get('alpha_t'))-inv(alpha_yt))	
        s_y  = d_y*(sin(s_yt/d_y))   # tooth thickness (chord-length)
        d_yc = d_y*(cos(s_yt/d_y))   # diameter at center of tooth (cut with chord)

        return s_y, d_yc
    

    def _circle(self, x_C, y_C, r, phi):
        """
        circle in cartesian coordinates

        INPUT parameter:
        x_C : x-value of center point (numeric)
        y_C : y-value of center point (numeric)
        r   : radius (numeric, positive)
        phi : angle parameter of circle point (numeric)[radians]
              zero value refers to radius parallel to x-axes pointing in positive x-direction
        
        OUTPUT:
        x   : x-value of circle point (numeric)
        y   : y-value of circle point (numeric)
        """

        # calculate point on circle
        x = x_C + r*cos(phi)
        y = y_C + r*sin(phi)

        return x, y


    def _analyzeFormcoords(self):
        # ONLY FOR EXTERNAL GEARS SO FAR !!!
        """
        analyze tooth form coordinates in order to get necessary information for
        geometry generator.

        INPUT parameters:
        formcoords : 2D cartesian coordinates of points on the
                     toothflank, describing a half tooth (TColgp_Array1OfPnt2d, pythonOCC)

        OUTPUT:
        suppdata :  supplement data for geardata dictionary (dictionary)
                    the dictionary contains at least the following keys:
                    d_f      : root circle diameter (numeric)
                    d_a      : tip diameter (numeric)
                    d_Ff     : root form diameter (numeric)
                    d_Fa     : tip form diameter (numeric)
                    z        : number of teeth (numeric, integer)
        """

        # transform formcoords to NumPy-array
        half_tooth = PythonOCCArrayToNumPyArray(self.formcoords)

        # convert to polar coordinates
        half_tooth_polar = zeros([size(half_tooth, 0)-1, 2])  
        for index in range(0, size(half_tooth, 0)-1):
             [r, phi] = CartesianCoordinatesToPolarCoordinates(half_tooth[index+1, 0], half_tooth[index+1, 1])
             half_tooth_polar[index, 0] = r
             half_tooth_polar[index, 1] = phi
        d_f = 2*min(half_tooth_polar[:,0])    # minimum radius --> root circle
        d_a = 2*max(half_tooth_polar[:,0])    # maximum radius --> tip circle
        tau = 2*(max(half_tooth_polar[:,1])-min(half_tooth_polar[:,1]))  # pitch angle [radians]
        z   = int(round(2*pi/tau))            # number of teeth from pitch angle

        # for finding form diameters, it is checked if the points are part of the flank involute
        # the limiting points of the flank involute define the form diameters
        if self.data.has_key('alpha_n') and self.data.has_key('alpha_t') and self.data.has_key('x_E'):
             tol = self._tol_default*self.data.get('m_n')    # tolerance for comparisons
             point_on_flank = False
             first_limit_diameter = None
             second_limit_diameter = None
             for point in range(0, size(half_tooth_polar, 0)):
                  [x, y] = self._toothThickness(2*half_tooth_polar[point, 0])
                  if abs(x+2*half_tooth[point+1, 0])<tol and abs(y-2*half_tooth[point+1, 1])<tol:       # numerical round off error may prevent this condition from becoming true!
                       if not point_on_flank:
                            first_limit_diameter = 2*half_tooth_polar[point, 0]
                       point_on_flank = True
                  else:
                       if point_on_flank:
                            second_limit_diameter = 2*half_tooth_polar[point, 0]
                       point_on_flank = False
             if second_limit_diameter==None:
                  second_limit_diameter = 2*half_tooth_polar[point, 0]
             if first_limit_diameter==second_limit_diameter:
                  raise ValueError, 'tooth form coordinate analysis failed'
             if first_limit_diameter>second_limit_diameter:
                  d_Fa = first_limit_diameter
                  d_Ff = second_limit_diameter
             else:
                  d_Fa = second_limit_diameter
                  d_Ff = first_limit_diameter
             if self.data.has_key('d_Ff'):          # use user-parameter if supplied
                  d_Ff = self.data.get('d_Ff')
             if self.data.has_key('d_Fa'):
                  d_Ff = self.data.get('d_Fa')
             return {'d_f':d_f, 'd_a':d_a, 'd_Ff':d_Ff, 'd_Fa':d_Fa, 'z':z}

        else:
             return {'d_f':d_f, 'd_a':d_a, 'z':z}
        

    def __init__(self, geardata, flankmods=None, formcoords=None):
        """
        Initialization of GearWheel-object.
        All parameters in accordance to DIN 3960 and DIN 3967.

        INPUT parameters:
        z        : number of teeth (numeric, integer)
        m_n      : normal module (numeric, positive)
        d        : pitch diameter (numeric)
                   two of the three parameters z, m_n, d, must be supplied
        b        : tooth width (numeric, positive)
        d_f      : root circle diameter (numeric)
                   optional - calculated if not supplied
        d_a      : tip diameter (numeric)
                   optional - calculated if not supplied
        d_Ff     : root form diameter (numeric)
                   optional - will be estimated if not supplied
        d_Fa     : tip form diameter (numeric)
                   optional - set equal da if not supplied (no chamfer)
        rho_f    : fillet radius (numeric)
                   optional - set equal 0.38*mn if not supplied
        x        : addendum modification factor (numeric)
                   optional - set equal 0.0 if not supplied
        alpha_n  : pressure angle (numeric, positive)[degrees]
                   optional - set equal 20.0 if not supplied
        beta     : helix angle (numeric)[degrees]
                   optional - set equal 0.0 if not supplied
        a        : addendum (numeric)
                   optional - no estimation
        c        : tip clearance (numeric, positive, 0.1...0.3*mn)
                   optional - set equal 0.25*mn if not supplied
        alpha_wt : service pressure angle (numeric, positive)[degrees]
                   optional - calculated from z_2 or d_w
        d_w      : service pitch diameter (numeric)
                   optional - calculated from alpha_wt or z_2
        h_k      : radial size of tip chamfer (numeric)
                   optional - set equal d_a-d_Fa or 0.0 if not supplied
        s_aK     : remaining tooth thickness at tip, chord-length (numeric)
                   optional - set equal s_a-2*h_k if not supplied
        z_2      : number of teeth of counter gear (numeric, integer)
                   optional - calculated from alpha_wt or d_w
        d_s      : shaft diameter, inner gear wheel diameter (numeric)
                   optional - set equal 0.0 if not supplied
        A_s      : tooth thickness allowance in normal cross-section (numeric, negative)
                   optional - set equal 0.0 if not supplied 

        All input parameters above are arranged in a dictionary. The keys are
        the names of the parameters as listed above.
                 
        formcoords : 2D cartesian coordinates of points on the
                     toothflank, describing a half tooth (TColgp_Array1OfPnt2d, pythonOCC)

        There are several possibilities for defining a complete gearwheel:
        1) z, m_n, b, (beta), formcoords
        2) z, m_n, b, (beta), d_f, d_a, d_Ff, d_Fa, rho_f
        3) z, m_n, b, (beta), alpha_n, alpha_wt, x, a, rho_f
        4) z, m_n, b, (beta), alpha_n, z_2, x, a, rho_f
        Some parameters can be left out, but the result might differ
        from your real gear. Missing parameters are estimated if
        possible. The helix angle beta doesn't have to be supplied
        for a spur gear.
        The constructor does not check for unit consistency. The user is
        responsible for supplying all values with consistent units.
        """

        self.data = deepcopy(geardata)
        self.modifications = deepcopy(flankmods)

        # form coordinates: value check (at least two points for defining a
        # tooth form (straight flanks) and two coordinates per point)      
        if formcoords:
             self.setFormCoords(formcoords, None)
             self.data.update(self._analyzeFormcoords())

        # number of teeth: value check
        if self.data.has_key('z') and not type(self.data.get('z'))==type(1):
             raise TypeError, 'number of teeth not integer'

        # module: value check
        if self.data.has_key('m_n') and not self.data.get('m_n')>=0:   
             raise ValueError, 'module non-positive'

        # helix angle: set to default if not supplied   
        if not self.data.has_key('beta'):
             self.data.update({'beta':self._beta_default})

        # calculate remaining parameter from number of teeth, module and pitch
        # diameter if only two of three are supplied (error if less are supplied)
        if self.data.has_key('tau') and not self.data.has_key('z'):
             self.data.update({'z':int(2*pi/self.data.get('tau'))})
        if self.data.has_key('z') and self.data.has_key('m_n'):
             self.data.update({'d':self.data.get('m_n')*self.data.get('z')/  \
                               cos(radians(self.data.get('beta')))})
        elif self.data.has_key('z') and self.data.has_key('d'):
             self.data.update({'m_n':self.data.get('d')*  \
                               cos(radians(self.data.get('beta')))/  \
                               self.data.get('z')})
        elif self.data.has_key('m_n') and self.data.has_key('d'):
             self.data.update({'z':int(self.data.get('d')*  \
                               cos(radians(self.data.get('beta')))/  \
                               self.data.get('m_n'))})
        else:
             raise AttributeError, 'insufficient data supplied'

        # calculate transverse pitch angle
        if not self.data.has_key('tau'):
             self.data.update({'tau':degrees(2*pi/self.data.get('z'))})

        # indicator if gear is external (number of teeth positive) or internal
        isexternal = sign(self.data.get('z'))
        # pitch diameter: value check (same sign as number of teeth)
        if not sign(self.data.get('d'))==isexternal:
             raise ValueError, 'sign of pitch diameter'

        # calculate module in transverse cross-section                
        self.data.update({'m_t':self.data.get('m_n')/cos(radians(  \
                          self.data.get('beta')))})

        # tooth width: check for existence and value check        
        if not self.data.has_key('b'):
             raise AttributeError, 'tooth width not supplied'
        if not self.data.get('b')>=0:
             raise ValueError, 'tooth width non-positive'

        # pressure angle: value check, set to default if not supplied
        if self.data.has_key('alpha_n'):
             if self.data.get('alpha_n')<0:
                  raise ValueError, 'pitch angle non-positive'
        else:
             self.data.update({'alpha_n':self._alpha_n_default})

        # addendum modification factor: set to default if not supplied     
        if not self.data.has_key('x'):
             self.data.update({'x':self._x_default})

        # tooth thickness allowance: set to default if not supplied
        if not self.data.has_key('A_s'):
             self.data.update({'A_s':self._A_s_default})
        # tooth thickness allowance: value check
        else:
             if not self.data.get('A_s')<=0:
                  raise ValueError, 'tooth thickness allowance positive'
        # calculate generating addendum modification coefficient:    
        self.data.update({'x_E':self.data.get('x')+self.data.get('A_s')/  \
                          2/tan(radians(self.data.get('alpha_n')))/self.data.get('m_n')})

        # calculate pressure angle in transverse cross-section
        self.data.update({'alpha_t':degrees(atan(tan(radians(self.data.get(  \
                          'alpha_n')))/cos(radians(self.data.get('beta')))))})

        # service pitch diameter: value check
        # calculate service pressure angle from service pitch diameter if not supplied
        if self.data.has_key('d_w') and not self.data.has_key('alpha_wt'):
             if not sign(self.data.get('d_w'))==isexternal:
                  raise ValueError, 'sign of service pitch diameter'
             self.data.update({'alpha_wt':degrees(acos(self.data.get('d')/  \
                               self.data.get('d_w')*cos(radians(  \
                               self.data.get('alpha_t')))))})
             
        # service pitch diameter: calculate from service pressure angle if possible
        if self.data.has_key('alpha_wt') and not self.data.has_key('d_w'):
             self.data.update({'d_w':self.data.get('d')*  \
                               cos(radians(self.data.get('alpha_t')))/  \
                               cos(radians(self.data.get('alpha_wt')))})

        # calculate base circle diameter
        self.data.update({'d_b':self.data.get('d')*cos(radians(self.data.get(  \
                          'alpha_t')))})

        # get further data from form coordinate analysis     
        if formcoords:
             self.data.update(self._analyzeFormcoords())
                
        if not formcoords:
             # tip clearance: value check, set to default if not supplied
             if self.data.has_key('c'):
                  if self.data.get('c')<0.1*self.data.get('m_n') or  \
                     self.data.get('c')>0.3*self.data.get('m_n'):
                       raise ValueError, 'tip clearance out of bounds'
             else:
                  self.data.update({'c':self._c_default*self.data.get('m_n')})

             # fillet radius: value check, set to default if not supplied                       
             if not self.data.has_key('rho_f'):
                  self.data.update({'rho_f':self._rho_f_default*self.data.get('m_n')})
             else:
                  if self.data.get('rho_f')<0:
                       raise ValueError, 'fillet radius negative'

             # CAUTION: THE FOLLOWING SECTION OF CODE WILL BE REMOVED IN FUTURE RELEASES!
             # tool fillet radius: value check                       
             if self.data.has_key('rho_fP'):
                  if self.data.get('rho_fP')<0:
                       raise ValueError, 'tool fillet radius negative'
                  if not self.data.get('beta')==0:
                       raise ValueError, 'fillet trochoid cannot be generated for helical gears'
             # END OF CODE SECTION TO BE REMOVED

             # calculate tip height modification factor if possible (else set to default)
             # various attempts are made
             if self.data.has_key('a'):
                  if self.data.has_key('alpha_wt') and not self.data.has_key('z_2'):
                       if self.data.get('alpha_wt')<0:
                            raise ValueError, 'service pressure angle non-positive'
                       self.data.update({'a_d':self.data.get('a')*cos(radians(   \
                                         self.data.get('alpha_wt')))/cos(radians(  \
                                         self.data.get('alpha_t')))})
                       self.data.update({'z_2':int(self.data.get('a_d')*2/self.data.get(  \
                                         'm_t')-self.data.get('z'))})
                  
                  elif self.data.has_key('z_2') and not self.data.has_key('alpha_wt'):
                       if not type(self.data.get('z_2'))==type(1):
                            raise TypeError, 'number of teeth of counter gear not integer'
                       if self.data.get('z')<0 and self.data.get('z_2')<0:
                            raise ValueError, 'two internal wheels cannot be paired'
                       self.data.update({'a_d':self.data.get('m_t')*(self.data.get('z')+  \
                                         self.data.get('z_2'))/2})
                       self.data.update({'alpha_wt':degrees(acos(self.data.get('a_d')/  \
                                        self.data.get('a')*cos(radians(self.data.get(  \
                                        'alpha_t')))))})
                  if self.data.has_key('alpha_wt') and self.data.has_key('z_2'):
                       x_2 = (inv(self.data.get('alpha_wt'))-inv(self.data.get('alpha_t')))*  \
                             (self.data.get('z')+self.data.get('z_2'))/2/tan(radians(  \
                             self.data.get('alpha_n')))-self.data.get('x')
                       self.data.update({'k':(self.data.get('a')-self.data.get('a_d'))/  \
                                         self.data.get('m_n')-(self.data.get('x')+x_2)})
                  else:
                       self.data.update({'k':self._k_default})
             else:     
                  self.data.update({'k':self._k_default})

             # root circle diameter: value check, calculate if not supplied
             if self.data.has_key('d_f'):
                  if self.data.get('d_f')>self.data.get('d'):
                       raise ValueError, 'root circle diameter greater than pitch diameter'
                  if not sign(self.data.get('d_f'))==isexternal:
                       raise ValueError, 'sign of root circle diameter'
             else:
                  self.data.update({'d_f':self.data.get('d')+2*self.data.get('x_E')*  \
                                    self.data.get('m_n')-2*(self.data.get('m_n')+  \
                                    self.data.get('c'))})

             # tip diameter: value check, calculate if not supplied
             if self.data.has_key('d_a'):
                 #if self.data.get('d_a')<self.data.get('d'):
                     #raise ValueError, 'tip diameter less than pitch diameter'
                 if not sign(self.data.get('d_a'))==isexternal:
                     raise ValueError, 'sign of tip diameter'
             else:
                  self.data.update({'d_a':self.data.get('d')+2*self.data.get('x')*  \
                                    self.data.get('m_n')+2*self.data.get('m_n')+2*  \
                                    self.data.get('k')*self.data.get('m_n')})

             # radial value of tip chamfer: value check, calculate or set to default
             # if not supplied
             if self.data.has_key('h_k'):
                  if self.data.get('h_k')<0:
                       raise ValueError, 'value of tip chamfer negative'
             elif self.data.has_key('d_Fa'):
                  self.data.update({'h_k':abs(self.data.get('d_a')-self.data.get('d_Fa'))/2})
             else:
                  self.data.update({'h_k':self._h_k_default})

             # remaining tooth thickness: value check, set to default if not supplied
             s_a, d_ac = self._toothThickness(self.data.get('d_a'))
             if not self.data.has_key('s_aK'):
                  self.data.update({'s_aK':s_a-2*self.data.get('h_k')})
             if self.data.get('s_aK')<0:
                  raise ValueError, 'remaining tooth thickness at tip negative'
             if self.data.get('s_aK')>s_a:
                  raise ValueError, 'remaining tip tooth thickness greater than tooth thickness'
			 
                  
             # root form diameter: value check
             if self.data.has_key('d_Ff'):
                  if self.data.get('d_Ff')>self.data.get('d'):
                       raise ValueError, 'root form diameter greater than pitch diameter'
                  if self.data.get('d_Ff')<self.data.get('d_f'):
                       raise ValueError, 'root form diameter less than root circle diameter'
                  if not sign(self.data.get('d_Ff'))==isexternal:
                       raise ValueError, 'sign of root form diameter'

             # tip form diameter: value check       
             if self.data.has_key('d_Fa'):
                  if self.data.get('d_Fa')<self.data.get('d'):
                       raise ValueError, 'tip form diameter less than pitch diameter'
                  if self.data.get('d_Fa')>self.data.get('d_a'):
                       raise ValueError, 'tip form diameter greater than tip diameter'
                  if not sign(self.data.get('d_Fa'))==isexternal:
                       raise ValueError, 'sign of tip form diameter'
             else:
                  self.data.update({'d_Fa':self.data.get('d_a')-2*self.data.get('h_k')})

        # shaft diameter: set to default if not supplied
        if not self.data.has_key('d_s'):
             self.data.update({'d_s':self._d_s_default})
        if abs(self.data.get('d_s'))>self._tol_default:
             if not sign(self.data.get('d_s'))==isexternal:
                  raise ValueError, 'sign of shaft diameter'
             if not self.data.get('d_s')<self.data.get('d_f'):
                  raise ValueError, 'shaft diameter greater than root circle diameter'
        
        # calculate tooth form coordinates if not supplied
        if not self.formcoords:
             self._makeFormCoords()
        else:
             self.formcoords = self._makeUnique(self.formcoords)

        # calculate tooth form wire if form coordinates user-supplied
        if not self._formwire:
             self._makeFormWire()

        # cleanup temporary data (should not be attribute of single gearwheel)       
        if self.data.has_key('z_2'):
             self.data.pop('z_2')
        if self.data.has_key('a'):
             self.data.pop('a')
        if self.data.has_key('a_d'):
             self.data.pop('a_d')


    def _makeFormWire(self):
        """
        Make PythonOCC wire of tooth form from user-supplied form coordinates.
        old form wire (if existend) will be replaced!
        This method is used only if user-supplied form coordinates are
        present.

        INPUT parameter:
        -
        """

        # module imports
        from numpy.linalg import norm
        
        # tolerance for comparisons
        tol = self._tol_default*self.data.get('m_n')

        # indicator whether gear is external (number of teeth positive) or internal
        isexternal    = sign(self.data.get('z'))
        
        # delete old form wire if existend
        if self._formwire:
             del(self._formwire)
        # make empty wire for tooth form (edges are added below sequentially)
        toothform_wire = ShapeExtend_WireData()

        # check if form coordinates are present
        if not self.formcoords:     
             _self.makeFormCoords()

        # transform formcoords to NumPy-array
        half_tooth = PythonOCCArrayToNumPyArray(self.formcoords)

        # convert to polar coordinates
        half_tooth_polar = zeros([size(half_tooth, 0), 2])  
        for index in range(0, size(half_tooth, 0)):
             [r, phi] = CartesianCoordinatesToPolarCoordinates(half_tooth[index, 0], half_tooth[index, 1])
             half_tooth_polar[index, 0] = r
             half_tooth_polar[index, 1] = phi

        # find points on root circle
        root_circle_condition = abs(abs(half_tooth_polar[:,0])-abs(self.data.get('d_f')/2))<tol

        # find points on root rounding
        root_rounding_condition = (~root_circle_condition) & \
                                  (isexternal*abs(half_tooth_polar[:,0])<=isexternal*(abs(self.data.get('d_Ff')/2.0+tol)))
        root_rounding_condition[nonzero(root_circle_condition)[0][-1]] = True  # edge points shall belong to both segments
        root_rounding_condition[0] = False    # the center of the gear wheel is not part of the root shape
        
        # find points on involute
        involute_condition = (zeros([size(half_tooth, 0)])==1)
        for index in range(1, size(half_tooth, 0)):
            s_yt, d_yc = self._toothThickness(2*isexternal*abs(half_tooth_polar[index, 0]))     
            if norm(array([-s_yt/2, d_yc/2])-half_tooth[index,:])<tol and isexternal*abs(half_tooth_polar[index,0])>=isexternal*abs(self.data.get('d_Ff')/2-tol):
                involute_condition[index] = True

        # find points on tip chamfer (or tip rounding)
        tip_chamfer_condition = (isexternal*abs(half_tooth_polar[:,0])>=isexternal*abs(self.data.get('d_Fa')/2-tol)) & \
                                (abs(abs(half_tooth_polar[:,0])-abs(self.data.get('d_a')/2))>tol)
        tip_chamfer_condition[0] = False    # the center of the gear wheel is not part of the tip chamfer
        tip_chamfer_condition[nonzero(involute_condition)[0][-1]] = True    # edge points shall belong to both segments

        # find points on tip circle
        tip_circle_condition = abs(abs(half_tooth_polar[:,0])-abs(self.data.get('d_a')/2))<tol
        tip_chamfer_condition[nonzero(tip_circle_condition)[0][0]] = True   # edge points shall belong to both segments

        for condition in [root_circle_condition, root_rounding_condition, involute_condition,  \
                          tip_chamfer_condition, tip_circle_condition]:
            # interpolate b-spline for segment
            # check if there are more than two points --> if not discard segment
            if count_nonzero(condition)>1:
                # determinate degree of spline interpolation depending on segment (necessary to avoid overshooting)
                if (condition==involute_condition).all():
                    mindeg = 3
                    maxdeg = 8
                else:
                    mindeg = 2
                    maxdeg = 5
     
                curve = Geom2dAPI_PointsToBSpline(NumPyArrayToPythonOCCArray(half_tooth[nonzero(condition)][:]), mindeg, maxdeg)
                curve_qualified = Geom2dGcc_QualifiedCurve(Geom2dAdaptor_Curve(curve.Curve()), GccEnt_unqualified)
                # make limiting vertices of segment curve
                start_point      = gp_Pnt2d(half_tooth[nonzero(condition)[0][0]][0],   \
                                            half_tooth[nonzero(condition)[0][0]][1])
                start_vertex     = BRepBuilderAPI_MakeVertex(self._make3dGeomFrom2dGeom(start_point)).Vertex()
                end_point        = gp_Pnt2d(half_tooth[nonzero(condition)[0][-1]][0],   \
                                            half_tooth[nonzero(condition)[0][-1]][1])
                end_vertex       = BRepBuilderAPI_MakeVertex(self._make3dGeomFrom2dGeom(end_point)).Vertex()
                
                # make edge from segment curve
                edge = BRepBuilderAPI_MakeEdge(self._make3dGeomFrom2dGeom(curve.Curve().GetObject()).GetHandle(), \
                                               start_vertex, end_vertex).Edge()
                toothform_wire.AddOriented(edge, TopAbs_FORWARD)
            
        self._formwire  = toothform_wire.Wire()


        
    def _makeFormCoords(self):
        """
        Tooth form coordinates in transverse cross-section (half tooth and half gap)
        points returned in 2D-cartesian coordinates, origin on wheel axis
        old form coordinates (if existend) will be replaced!
        This method should be used only if no user-supplied form coordinates are
        present.

        INPUT parameter:
        -
        """

        # module imports
        from scipy.optimize import fsolve
        from numpy.linalg import norm
        from numpy import insert

        # tolerance for comparisons
        tol = self._tol_default*self.data.get('m_n')
        
        # delete old form coordinates if existend
        if self.formcoords:     
             del(self.formcoords)
        if self._formwire:
             del(self._formwire)

        
        # indicator whether gear is external (number of teeth positive) or internal
        isexternal    = sign(self.data.get('z'))
        inv_extension = False
             
        # indices for adressing parts of tooth form
        lower_index          = 0
        start_rootcirc_index = lower_index+1        # one entry reserved for origin
        end_rootcirc_index   = start_rootcirc_index+self.points_root-1
        start_fillet_index   = end_rootcirc_index
        end_fillet_index     = start_fillet_index+self.points_fillet-1
        start_involute_index = end_fillet_index+1   # can differ from end of fillet
        end_involute_index   = start_involute_index+self.points_flank-1
        start_chamfer_index  = end_involute_index+1
        end_chamfer_index    = start_chamfer_index+self.points_chamfer-1
        start_tipcirc_index  = end_chamfer_index+1  # differs from end of involute if chamfer present
        end_tipcirc_index    = start_tipcirc_index+self.points_tip-1
        upper_index          = end_tipcirc_index


        # determine boundary of half tooth segment on root circle
        rootcirc_start_point = self.data.get('d_f')/2*array([-sin(radians(self.data.get('tau')/2)),   \
                                                              cos(radians(self.data.get('tau')/2))])
        
        # determine how the root shape is defined and calculate significant points
        # root shape is circular in transverse cross-section
        if isexternal>0:        # for external gears
            if not self.data.has_key('d_Ff'):
                # root circle is tangent to involute
                if (self.data.get('d_f')**2+4*self.data.get('rho_f')*self.data.get('d_f') >= self.data.get('d_b')**2):
                    self.data.update({'d_Ff': isexternal*sqrt((sqrt((self.data.get('d_f')+2*self.data.get('rho_f'))**2-  \
                                              self.data.get('d_b')**2)-2*self.data.get('rho_f'))**2+self.data.get('d_b')**2)})
                    s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
                    fil_end_point   = array([-s_yt/2, d_yc/2])
                # no tangency possible: undercut
                elif (self.data.get('d_f')+4*self.data.get('rho_f') >= self.data.get('d_b')):  
                    self.data.update({'d_Ff': self.data.get('d_b')})
                    s_yt, d_yc = self._toothThickness(self.data.get('d_b'))
                    fil_end_point   = array([-s_yt/2, d_yc/2])          # end of involute at base circle
                    print 'Warning: undercutting occurs!'
                # in case all prior attempts to construct root fillet failed, the involute has to be extended with a straight tangential line
                else:
                    self.data.update({'d_Ff': self.data.get('d_b')})
                    d_tangent = sqrt(self.data.get('d_f')**2+4*self.data.get('rho_f')*self.data.get('d_f'))   # diameter around gear center on that tangency point of fillet curve is located
                    s_yt, d_yc = self._toothThickness(self.data.get('d_b'))
                    nu = atan(s_yt/d_yc)
                    fil_end_point   = array([-d_tangent/2*sin(nu), d_tangent/2*cos(nu)])    # tangential extension of involute beyond base circle
                    print 'Warning: involute had to be extended below base cicle to achieve root fillet tangency!'
                    inv_extension = True
            else:        
                # if root form circle diameter is supplied, it is forced strictly if possible
                if (self.data.get('d_Ff')-self.data.get('d_f'))/2 > 2*self.data.get('rho_f'):   # check if root fillet circle fits beetween root form circle and root circle
                    raise ValueError, 'root fillet radius too small: root shape cannot be determined'
                s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
                if abs(self.data.get('d_Ff')) >= abs(self.data.get('d_b')):  # fillet ends at root form circle
                    fil_end_point   = array([-s_yt/2, d_yc/2])
                else:  # base circle diameter greater than root form diameter: tangential extension of involute
                    nu = atan(s_yt/d_yc)
                    fil_end_point   = array([-self.data.get('d_Ff')*sin(nu),  \
                                              self.data.get('d_Ff')*cos(nu)])
                    print 'Warning: involute had to be extended below base cicle to enforce root form circle diameter!'
                    inv_extension = True
                    
        else:       # for internal gears
            if not self.data.has_key('d_Ff'):
                # root circle is tangent to involute
                t_b = sqrt((self.data.get('d_f')/2+self.data.get('rho_f'))**2-(self.data.get('d_b')/2)**2)
                self.data.update({'d_Ff': -2*sqrt((t_b+self.data.get('rho_f'))**2+(self.data.get('d_b')/2)**2)})
            else:
                # if root form circle diameter is supplied, it is forced strictly if possible
                if (self.data.get('d_Ff')-self.data.get('d_f'))/2 > 2*self.data.get('rho_f'):   # check if root fillet circle fits beetween root form circle and root circle
                    raise ValueError, 'root fillet radius too small: root shape cannot be determined'               
            s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
            fil_end_point = array([-s_yt/2, d_yc/2])

        # find center of root fillet circle by cutting circle around fillet end point with radius rho_f
        # with circle around center of gear wheel with radius d_f/2+rho_f
        def root_circle_center_func(phi):
            return fil_end_point+self.data.get('rho_f')*array([sin(phi[0]), cos(phi[0])])-  \
                   (self.data.get('d_f')/2+self.data.get('rho_f'))*array([sin(phi[1]), cos(phi[1])])

        phi_fil_center = fsolve(root_circle_center_func, [-pi/2, 0.0])
        fil_center_point = (self.data.get('d_f')/2+self.data.get('rho_f'))*array([sin(phi_fil_center[1]), cos(phi_fil_center[1])])

        # boundary point of root fillet and root circle
        fil_start_point = fil_center_point*self.data.get('d_f')/(self.data.get('d_f')+2*self.data.get('rho_f'))

        # if boundary point and fillet center are outside half tooth segment the shape of the root fillet
        # cannot be determined (root fillet curve is not continously differentiable and d_f is not matched)
        if abs(atan(fil_start_point[0]/fil_start_point[1])) > abs(radians(self.data.get('tau')/2)):
            raise ValueError, 'root fillet radius too large: root shape cannot be determined'

        # determine boundary points of involute
        s_yt, d_yc = self._toothThickness(self.data.get('d_Ff'))
        inv_start_point = array([-s_yt/2, d_yc/2])   # involute starts at root form circle
        s_yt, d_yc = self._toothThickness(self.data.get('d_Fa'))
        inv_end_point = array([-s_yt/2, d_yc/2])     # involute ends at tip form circle

        # determine boundary points of tip circle
        nu                  = self.data.get('s_aK')/self.data.get('d_a')
        tipcirc_start_point = array([-self.data.get('d_a')/2*sin(nu),  \
                                      self.data.get('d_a')/2*cos(nu)])   # tip circle starts at end of tip chamfer
        tipcirc_end_point = array([0.0, self.data.get('d_a')/2])         # tip circle ends at symmetry line
            
        # create array for tooth form coordinates
        formcoord_array = zeros([upper_index, 2])

        # compute points on root circle
        phi_start = -asin(2*rootcirc_start_point[0]/self.data.get('d_f'))           # starting angle of root circle
        if abs(phi_start-acos(2*rootcirc_start_point[1]/self.data.get('d_f')))>tol: # computation is not unique
            phi_start = pi-phi_start
        phi_end = -asin(2*fil_start_point[0]/self.data.get('d_f'))                  # end angle of root circle
        if abs(phi_end-acos(2*fil_start_point[1]/self.data.get('d_f')))>tol:        # computation is not unique
            phi_end = pi-phi_end
        if abs(phi_start-phi_end)>tol:                         # check if a root circle curve exists
            delta_phi = (phi_end-phi_start)/(self.points_root-1)
            n = 0
            for index in range(start_rootcirc_index, end_rootcirc_index):
                formcoord_array[index]=self.data.get('d_f')/2*array([-sin(phi_start+n*delta_phi), isexternal*cos(phi_start+n*delta_phi)])
                n += 1

        # compute points on root fillet
        print 'Warning: circular root fillet in transverse cross-section assumed!'
        phi_start = asin((fil_start_point[0]-fil_center_point[0])/self.data.get('rho_f'))              # starting angle of root fillet
        if abs(phi_start-acos(-(fil_start_point[1]-fil_center_point[1])/self.data.get('rho_f')))>tol:  # computation is not unique
            phi_start = pi-phi_start
        phi_end = asin((fil_end_point[0]-fil_center_point[0])/self.data.get('rho_f'))                  # end angle of root fillet
        if abs(phi_end-acos(-(fil_end_point[1]-fil_center_point[1])/self.data.get('rho_f')))>tol:      # computation is not unique
            phi_end = pi-phi_end
        if abs(phi_start-phi_end)>tol:                         # check if a root fillet curve exists
            delta_phi = (phi_end-phi_start)/(self.points_fillet-1)
            n = 0
            for index in range(start_fillet_index, end_fillet_index+1):
                formcoord_array[index]=fil_center_point+self.data.get('rho_f')*array([sin(phi_start+n*delta_phi), -isexternal*cos(phi_start+n*delta_phi)])
                n += 1
        if (inv_start_point-fil_end_point).any():                         # check if a root fillet circle connects directly to flank
            print 'involute was extended'                                 # placeholder for future
            
        # compute points on flank
        d_start = isexternal*norm(inv_start_point,2)*2             # start diameter of involute flank (root form diameter)
        d_end   = isexternal*norm(inv_end_point,2)*2               # end diameter of involute flank (tip form diameter)
        delta_d = (d_end-d_start)/(self.points_flank-1)
        n = 0
        for index in range(start_involute_index, end_involute_index+1):
            s_yt, d_yc = self._toothThickness(d_start+n*delta_d)
            formcoord_array[index] = array([-s_yt/2, d_yc/2])
            n += 1

        # compute points on tip chamfer
        if self.data.has_key('h_k') and (self.data.get('h_k') > 0):
            print 'Warning: straight tip chamfer assumed!'
            delta_k = 1/(self.points_chamfer-1)
            n = 0
            for index in range(end_involute_index, end_chamfer_index):
                formcoord_array[index] = inv_end_point+(tipcirc_start_point-inv_end_point)*n*delta_k
                n += 1
        
        # compute points on tip circle
        phi_start = -asin(2*tipcirc_start_point[0]/self.data.get('d_a'))                # starting angle of tip circle
        if abs(phi_start-acos(2*tipcirc_start_point[1]/self.data.get('d_a')))>tol:      # computation is not unique
            phi_start = pi-phi_start
        phi_end = -asin(2*tipcirc_end_point[0]/self.data.get('d_a'))                    # end angle of tip circle
        if abs(phi_end-acos(2*tipcirc_end_point[1]/self.data.get('d_a')))>tol:          # computation is not unique
            phi_end = pi-phi_end
        if isexternal<0:
            phi_end = phi_end+pi
        if abs(phi_start-phi_end)>tol:                         # check if a tip circle curve exists
            delta_phi = (phi_end-phi_start)/(self.points_tip-1)
            n = 1
            for index in range(end_chamfer_index+1, end_tipcirc_index):
                formcoord_array[index]=self.data.get('d_a')/2*array([-sin(phi_start+n*delta_phi), isexternal*cos(phi_start+n*delta_phi)])
                n += 1
        
        # compute points on tangential extension of involute below base circle
        if inv_extension:
            delta_k = 1/(self.points_ext-1)
            for n in range(1, self.points_ext-1):
                formcoord_array = insert(formcoord_array, start_involute_index, inv_start_point+(fil_end_point-inv_start_point)*n*delta_k, axis=0)
        
        # transform formcoords to NumPy-array
        self.formcoords = NumPyArrayToPythonOCCArray(formcoord_array)

        # remove redundant entries and set class attributes
        self.formcoords = self._makeUnique(self.formcoords)

        
        
    def makeTooth(self):
        """
        Tooth form in transverse cross-section (one tooth and one gap, tooth centered)
        points returned in 2D-cartesian coordinates, origin on wheel axis

        INPUT parameter:
        -

        OUTPUT:
        toothcoords    : coordinates of one tooth in transverse cross-section (TColgp_Array1OfPnt2d, pythonOCC)
        toothform_wire : wire of tooth form (TopoDS_Wire, pythonOCC)
        """
        
        # check if tooth form coordinates are available: calculate otherwise
        if not self.formcoords:
             self._makeFormCoords()
        if not self._formwire:     
             _self.makeFormWire()

        # set up transformation matrix
        origin_point = gp_Pnt(0.0, 0.0, 0.0)  # -2*self.data.get('b')/2)
        mirror_dir   = gp_Dir(0.0, 1.0, 0.0)
        mirror_axis  = gp_Ax1(origin_point, mirror_dir)
        mirror_trans = gp_Trsf()
        mirror_trans.SetMirror(mirror_axis)
        
        # mirror wire of tooth form to get complete tooth
        formwire_trans    = BRepBuilderAPI_Transform(self._formwire, mirror_trans)
        formwire_mirrored = TopoDS_Cast.Wire(formwire_trans.Shape())
        formwire_mirrored_reversed    = ShapeExtend_WireData()
        formwire_mirrored_reversed.AddOriented(formwire_mirrored, TopAbs_REVERSED)
        formwire_mirrored_reversed.Reverse()
        toothform_wire    = ShapeExtend_WireData()
        toothform_wire.AddOriented(self._formwire, TopAbs_REVERSED)       
        toothform_wire.AddOriented(formwire_mirrored_reversed.Wire(), TopAbs_FORWARD)
        toothcoords = TColgp_Array1OfPnt2d(1, 2*(self.formcoords.Length()-1))  # origin-point not required

        # mirroring of tooth form coordinates at center axis of tooth (x=0.0)
        lower_index = self.formcoords.Lower()
        upper_index = self.formcoords.Upper()
        for index in range(lower_index+1, upper_index+1):  # origin-point not required
             toothcoords.SetValue(index-1, self.formcoords.Value(index))
             mirror_point = gp_Pnt2d(-self.formcoords.Value(index).X(), self.formcoords.Value(index).Y())
             toothcoords.SetValue(toothcoords.Length()-index+2, mirror_point)

        return self._makeUnique(toothcoords), toothform_wire.Wire()


    def makeCrossSection(self):
        """
        Face and Sketch of transverse cross-section of gear wheel
        geometries returned with origin at [0.0, 0.0, -b/2]
        face and wire geometry entities are returned separately

        INPUT parameter:
        -

        OUTPUT:
        crosssection_face : face of the gear wheel's cross-section (TopoDS_Face, pythonOCC)
        shaft_wire        : wire of the gear wheel's cross-section without shaft contour (TopoDS_Wire, pythonOCC)
        """

        # tolerance for comparisons
        tol = self._tol_default*self.data.get('m_n')
        
        toothcoords, toothwire  = self.makeTooth()

        # indicator whether gear is external (number of teeth positive) or internal
        isexternal = sign(self.data.get('z'))

        # define plane of cross-section
        origin_cross     = gp_Pnt(0.0, 0.0, 0.0) # -self.data.get('b')/2)
        normal_cross     = gp_Dir(gp_XYZ(0.0, 0.0, -1.0))
        plane_cross      = gp_Pln(origin_cross, normal_cross)

        # set up transformation matrix
        origin_point = gp_Pnt(0.0, 0.0, 0.0)  # -self.data.get('b')/2)
        x_axis_dir   = gp_Dir(gp_XYZ(1.0, 0.0, 0.0))
        y_axis_dir   = gp_Dir(gp_XYZ(0.0, 1.0, 0.0))
        z_axis_dir   = gp_Dir(gp_XYZ(0.0, 0.0, 1.0))
        sym_axis     = gp_Ax1(origin_point, z_axis_dir)
        rot_trans    = gp_Trsf()

        # useful objects
        crosssection_wire  = ShapeExtend_WireData()

        # rotate tooth objects
        for toothnumber in range(0, abs(self.data.get('z'))):
             rotang = -toothnumber*radians(self.data.get('tau'))
             # make wire object
             rot_trans.SetRotation(sym_axis, rotang)
             crosssection_trans     = BRepBuilderAPI_Transform(toothwire, rot_trans)
             crosssection_rotated   = TopoDS_Cast.Wire(crosssection_trans.Shape())
             if isexternal>0:
                 crosssection_wire.AddOriented(crosssection_rotated, TopAbs_FORWARD)
             else:
                 crosssection_wire.AddOriented(crosssection_rotated.Reversed(), TopAbs_FORWARD)

        crosssection_face = BRepBuilderAPI_MakeFace(crosssection_wire.Wire())

        # create shaft circle
        shaftcirc_plane = gp_Ax2(origin_point, z_axis_dir, y_axis_dir)
        if abs(self.data.get('d_s'))>tol:
             shaft_curve = Geom_Circle(shaftcirc_plane, abs(self.data.get('d_s')/2))
             shaft_edge  = BRepBuilderAPI_MakeEdge(shaft_curve.GetHandle())
             shaft_wire  = BRepBuilderAPI_MakeWire(shaft_edge.Edge())
             if isexternal<0:
                  crosssection_face.Add(shaft_wire.Wire())
             else:
                  crosssection_face.Add(TopoDS_Cast.Wire(shaft_wire.Wire().Reversed()))
             has_shaft   = True
        else:
             has_shaft   = False

        return crosssection_face.Face(), crosssection_wire.WireAPIMake()



    def makeOCCSolid(self):
        """
        Returns an OpenCascade-solid (pythonOCC-implementation) of the gearwheel.
        It can be used in an OCC-context or converted to VRML-, STEP- or IGES-format
        for further use.

        INPUT parameter:
        -

        OUTPUT:
        solid_final : 3D-geometry of gear wheel (OCC.TopoDS.TopoDS_Solid-instance)
        """

        # construct center coordinate system
        dir_x_axes     = gp_Dir(1.0, 0.0, 0.0)
        dir_z_axes     = gp_Dir(0.0, 0.0, 1.0)
        dir_z_axes_neg = gp_Dir(0.0, 0.0, -1.0)
        origin_CS      = gp_Pnt(0.0, 0.0, 0.0)
        center_CS      = gp_Ax2(origin_CS, dir_z_axes_neg, dir_x_axes)

        # construct gear curve
        crosssection_face, crosssection_wire = self.makeCrossSection()

        # construct spine and auxiliary spine for extrusion
        spinepol    = BRepBuilderAPI_MakePolygon()
        auxspinepol = BRepBuilderAPI_MakePolygon()
        delta_b     = self.data.get('b')/self.points_width
        beta_a      = atan(self.data.get('d_a')/self.data.get('d')*  \
                           tan(radians(self.data.get('beta'))))   # helix angle at tip diameter
        slope       = tan(pi/2-beta_a)*self.data.get('d_a')*pi
        delta_phi   = delta_b/slope*2*pi
        for i in range(0,self.points_width+1):
             spinepol.Add(gp_Pnt(0.0, 0.0, -self.data.get('b')/2+i*delta_b))
             x_aux = cos(i*delta_phi)
             y_aux = sin(i*delta_phi)
             auxspinepol.Add(gp_Pnt(x_aux, y_aux, -self.data.get('b')/2+i*delta_b))

        # construct gear wheel from cross-section
        solid_gearwheel  = BRepOffsetAPI_MakePipeShell(spinepol.Wire())
        solid_gearwheel.SetMode(auxspinepol.Wire(), False)
        solid_gearwheel.Add(crosssection_wire)
        solid_gearwheel.Build()

        # construct faces of start and end cross-sections
        pnt_start   = gp_Pnt(0.0, 0.0, -self.data.get('b')/2)
        pnt_end     = gp_Pnt(0.0, 0.0, self.data.get('b')/2)
        plane_start = gp_Pln(pnt_start, dir_z_axes_neg)
        plane_end   = gp_Pln(pnt_end, dir_z_axes)
        if not self.data.get('d_s')==0.0:
            start_CS    = gp_Ax2(pnt_start, dir_z_axes_neg, dir_x_axes)
            end_CS      = gp_Ax2(pnt_end, dir_z_axes, dir_x_axes)
            circ_start  = gp_Circ(start_CS, abs(self.data.get('d_s')/2))
            circ_end    = gp_Circ(end_CS, abs(self.data.get('d_s')/2))
            edge_start  = BRepBuilderAPI_MakeEdge(circ_start)
            edge_end    = BRepBuilderAPI_MakeEdge(circ_end)
            wire_start  = BRepBuilderAPI_MakeWire(edge_start.Edge())
            wire_end    = BRepBuilderAPI_MakeWire(edge_end.Edge())
        wire_first  = TopoDS_Cast.Wire(solid_gearwheel.FirstShape())
        wire_last   = TopoDS_Cast.Wire(solid_gearwheel.LastShape())
        if self.data.get('z')>0:
            face_first = BRepBuilderAPI_MakeFace(plane_start, wire_first)
            if not self.data.get('d_s')==0.0:
                face_first.Add(TopoDS_Cast.Wire(wire_start.Shape().Reversed()))
            face_last  = BRepBuilderAPI_MakeFace(plane_end, TopoDS_Cast.Wire(wire_last.Reversed()))
            if not self.data.get('d_s')==0.0:
                face_last.Add(TopoDS_Cast.Wire(wire_end.Shape().Reversed()))
        else:
             
             face_first = BRepBuilderAPI_MakeFace(plane_start,  \
             TopoDS_Cast.Wire(wire_start.Shape().Reversed()))
             face_first.Add(wire_first)
             face_last  = BRepBuilderAPI_MakeFace(plane_end, TopoDS_Cast.Wire(wire_end.Shape().Reversed()))
             face_last.Add(TopoDS_Cast.Wire(wire_last.Reversed()))

        # construct shaft cylinder and final solid
        if not self.data.get('d_s')==0.0:
             vec_path      = gp_Vec(pnt_start, pnt_end)
             face_shaft    = BRepPrimAPI_MakePrism(wire_start.Wire(), vec_path)
        elif self.data.get('z')<0:
             raise ValueError, 'internal gear needs shaft diameter for geometry construction'

        # make shell of gearwheel
        shell_final = BRepBuilderAPI_Sewing()
        shell_final.Add(solid_gearwheel.Shape())
        shell_final.Add(face_first.Shape())
        shell_final.Add(face_last.Shape())
        if not self.data.get('d_s')==0.0:
            shell_final.Add(face_shaft.Shape())
        shell_final.Perform()

        # make solid of gearwheel
        shell_gear  = TopoDS_Cast.Shell(shell_final.SewedShape())
        solid_final = BRepBuilderAPI_MakeSolid(shell_gear).Solid()

        return solid_final


    def calcInertiaProperties(self, rho=None, solid=None):
        """
        Calculates the mass properties of the gear wheel. In the case no mass density
        is supplied, the volume and the quadratic moments of the volume are returned.
        These have to be multiplied with the mass density in order to obtain the mass
        and moments of inertia. If an OCC-solid of the gearwheel has been already
        generated, it is possible to supply it to this method in order to save
        computation time.
        The inertia tensor is computed relative to the center of mass of the gear wheel.
        It is a point on the axis of rotation exactly in the middle between the two
        planar faces perpendicular to the axis of rotation (distance from each face: b/2).
        The axis are the principle axis of inertia with the z-axis coincididing with the
        axis of rotation.
        Warning: - no unit consistency check is performed. The user is responsible for
                   supplying the mass density in correct units.
                 - only applicable for homogenous bodies (constant mass density)

        INPUT parameter:
        rho:       mass density (numeric)
                   optional - set to 1.0 if not supplied (see above)
        solid:     3D-geometry of gear wheel (OCC.TopoDS.TopoDS_Solid-instance)
                   optional - generated if not supplied

        OUTPUT:
        m:     mass of the solid (numeric)
        J:     inertia tensor of the solid (numeric, 3x3-matrix)
        r_g:   radius of gyration around rotation axis (numeric)
        """

        # mass density: value check
        if rho:
            if rho<=0.0:
                raise ValueError, 'non-positive mass density'
        else:
            rho = 1.0   # set to default if not supplied
            print 'Warning: mass density set to 1.0'

        # create OCC-solid if not supplied
        if not solid:
            solid = self.makeOCCSolid()
            
        eps = 1E-7      # relative error for volume computation
        reference_frame = gp_Pnt(0.0, 0.0, 0.0)    # center of gear wheel (= center of mass)

        # calculate mass properties
        inertia_prop    = GProp_GProps(reference_frame)
         # prop_calc       = brepgprop()
        brepgprop.VolumeProperties(solid, inertia_prop, eps)

        # mass
        m = inertia_prop.Mass()*rho

        # inertia tensor
        J = []
        for row_index in range(1,4):
            row = []
            for column_index in range(1,4):
                row.append(inertia_prop.MatrixOfInertia().Value(row_index, column_index)*rho)
            J.append(row)
                    
        # radius of gyration around rotation axis
        rotation_axis = gp_Ax1(reference_frame, gp_Dir(0.0, 0.0, 1.0))
        r_g = inertia_prop.RadiusOfGyration(rotation_axis)
        
        return m, J, r_g


class GearPair:
    """
    Parent Class for all gear wheel pairs.
    """

    # Attributes: gear pair data
    data = None    # dictionary containing all gear pair parameters

    # Attributes: gear wheels
    Pinion = None  # GearWheel-object
    Gear   = None  # GearWheel-object
        

    def __str__(self):
        """
        Define string conversion of GearPair objects

        INPUT parameter:
        -

        OUTPUT:
        string representation of class
        """
        
        outstr = 'gear pair data:\n'
        # output gear pair data
        for date in self.data:
             outstr += date.ljust(15)+':\t'+str(self.data.get(date))+'\n'

        # output pinion data
        if self.Pinion:    
             outstr += '\npinion data:\n'
             outstr += str(self.Pinion)

        # output gear data
        if self.Gear:    
             outstr += '\ngear data:\n'
             outstr += str(self.Gear)

        return outstr
    

    def __init__(self, pairdata, Pinion=None, Gear=None):
        """
        Initialization of GearPair-object
        Should be overwritten in derived classes

        INPUT parameters:
        pairdata   : data of gear wheel pair (dictionary)
        Pinion     : pinion (GearWheel-instance)
        Gear       : gear   (GearWheel-instance)
        """

        self.data       = deepcopy(pairdata)
        if Pinion:
             self.setPinion(Pinion)
        if Gear:
             self.setGear(Gear)


    def setPinion(self, Pinion):
        """
        Set pinion attribute

        INPUT parameter:
        Pinion     : pinion (GearWheel-instance)
        """
 
        if isinstance(Pinion, GearWheel):
             self.Pinion = Pinion
        else:
             raise TypeError, 'GearWheel instance expected'



    def setGear(self, Gear):
        """
        Set gear attribute

        INPUT parameter:
        Gear     : gear (GearWheel-instance)
        """

        if isinstance(Gear, GearWheel):
             self.Gear = Gear
        else:
             raise TypeError, 'GearWheel instance expected'


class CylindricalGearPair(GearPair):
    """
    Class representing a meshing spur or helical gear pair.
    Applicable for external and internal gears.
    Derived from GearPair-class
    """

    # Attributes: default settings for parameters
    _A_a_default       =  0.0    # default value for centre distance allowance
    
    def __init__(self, pairdata, Pinion=None, Gear=None):
        """
        Initialization of GearPair-object
        Should be overwritten in derived classes

        INPUT parameters:
        pairdata   : data of gear wheel pair (dictionary)
                     a : centre distance must be supplied
        Pinion     : pinion (GearWheel-instance)
        Gear       : gear   (GearWheel-instance)

        possible parameters (keys) in self.data (dictionary):
        epsilon_alpha : transverse contact ratio (numeric)
        epsilon_beta  : overlap ratio (numeric)
        epsilon gamma : contact ratio (numeric)
        alpha_wt      : service pressure angle (numeric) [degrees]
        alpha_t       : transverse pressure angle (numeric) [degrees]
        alpha_n       : pressure angle (numeric) [degrees]
        g_alpha       : line of action length (numeric)
        m_n           : module (numeric)
        a_d           : reference centre distance (numeric)
        a             : centre distance (numeric)
        A_a           : centre distance allowance (numeric)
        k             : tip height modification factor (numeric)
        i             : gear ratio (numeric)
        j_t           : acceptance backlash (numeric)
        """

        # value check: addendum must be supplied
        if not pairdata.has_key('a'):
             raise AttributeError, 'centre distance not found'
             
        GearPair.__init__(self, pairdata, Pinion, Gear)

        if self.Pinion and self.Gear:
            
             # compatibility check (pinion to gear)
             if not self.Pinion.data.get('m_n')==self.Gear.data.get('m_n'):
                  raise ValueError, 'gear and pinion cannot have different module'
             if sign(self.Pinion.data.get('z'))==sign(self.Gear.data.get('z')):
                  if not -self.Pinion.data.get('beta')==self.Gear.data.get('beta'):
                      raise ValueError, 'helix angles of gear and pinion not compatible'
             else:
                  if not self.Pinion.data.get('beta')==self.Gear.data.get('beta'):
                      raise ValueError, 'helix angles of gear and pinion not compatible'
                    
             # warning, because addendum will be set for GearWheel-objects
             if self.Pinion.data.has_key('a') or self.Gear.data.has_key('a'):
                  warn('addendum will be overwritten', UserWarning)
                  
             # calculate and set gear pair parameters
             # addendum
             self.data.update({'a_d':self.Pinion.data.get('m_n')/2*   \
                                     (self.Pinion.data.get('z')+self.Gear.data.get('z'))/  \
                                     cos(radians(self.Pinion.data.get('beta')))})
             
             # pitch angles
             self.data.update({'alpha_wt':degrees(acos((self.Pinion.data.get('z')+   \
                                          self.Gear.data.get('z'))*self.Pinion.data.get('m_t')*  \
                                          cos(radians(self.Pinion.data.get('alpha_t')))/   \
                                          self.data.get('a')/2))})
             self.data.update({'alpha_n':self.Pinion.data.get('alpha_n'),  \
                               'alpha_t':self.Pinion.data.get('alpha_t')})
             
             # set parameters of GearWheel-objects (old values are overwritten)
             self.Pinion.data.update({'z_2':self.Gear.data.get('z'),  \
                                      'alpha_wt':self.data.get('alpha_wt'),  \
                                      'a':self.data.get('a'),'a_d':self.data.get('a_d')})            
             self.Gear.data.update({'z_2':self.Pinion.data.get('z'),  \
                                    'alpha_wt':self.data.get('alpha_wt'),  \
                                    'a':self.data.get('a'),'a_d':self.data.get('a_d')})
             
             # addendum modification factors
             x_sum = (self.Pinion.data.get('z')+self.Gear.data.get('z'))*  \
                     (inv(self.data.get('alpha_wt'))-inv(self.data.get('alpha_t')))/  \
                      2/tan(radians(self.data.get('alpha_n')))
             if self.Pinion.data.get('x')!=0.0 and not self.Gear.data.get('x')!=0.0:
                  self.Gear.data.update({'x':x_sum-self.Pinion.data.get('x')})
             elif self.Gear.data.get('x')!=0.0 and not self.Pinion.data.get('x')!=0.0:
                  self.Pinion.data.update({'x':x_sum-self.Gear.data.get('x')})   
             elif self.Pinion.data.get('x')==0.0 and self.Gear.data.get('x')==0.0:   # ensure compatibility
                  self.Pinion.data.update({'x':x_sum})
                  self.Gear.data.update({'x':0.0})
             else:
                  self.Gear.data.update({'x':x_sum-self.Pinion.data.get('x')})
                  
             # normal modul
             self.data.update({'m_n':self.Pinion.data.get('m_n')})
             
             # tip height modification factor
             self.data.update({'k':(self.data.get('a')-self.data.get('a_d'))/  \
                               self.data.get('m_n')-x_sum})
             self.Pinion.data.update({'k':self.data.get('k')})
             self.Gear.data.update({'k':self.data.get('k')})
             
             # update gear wheel objects
             self.Pinion.formcoords = None
             self.Pinion._formwire = None
             self.Pinion.setGearData(self.Pinion.data)
             self.Gear.formcoords = None
             self.Gear._formwire = None
             self.Gear.setGearData(self.Gear.data)

             # line of action length
             self.data.update({'g_alpha':(sqrt(self.Pinion.data.get('d_a')**2-  \
                                               self.Pinion.data.get('d_b')**2)+  \
                                          sqrt(self.Gear.data.get('d_a')**2-  \
                                               self.Gear.data.get('d_b')**2)-  \
                                          (self.Pinion.data.get('d_b')+self.Gear.data.get('d_b'))*  \
                                          tan(radians(self.data.get('alpha_wt'))))/2})
             
             # transverse contact ratio
             self.data.update({'epsilon_alpha':2*self.data.get('g_alpha')/self.Pinion.data.get('d')/  \
                               radians(self.Pinion.data.get('tau'))})

             # overlap ratio
             self.data.update({'epsilon_beta':2*self.Pinion.data.get('b')*  \
                               tan(radians(abs(self.Pinion.data.get('beta'))))/  \
                               self.Pinion.data.get('d')/radians(self.Pinion.data.get('tau'))})

             # contact ratio
             self.data.update({'epsilon_gamma':self.data.get('epsilon_alpha')+  \
                               self.data.get('epsilon_beta')})

             # transmission ratio
             self.data.update({'i':-self.Pinion.data.get('z')/self.Gear.data.get('z')})

             # centre distance allowance: set to default if not supplied
             if not self.data.has_key('A_a'):
                  self.data.update({'A_a':self._A_a_default})

             # acceptance backlash (DIN 3967)
             self.data.update({'j_t':-(self.Pinion.data.get('A_s')+self.Gear.data.get('A_s'))/  \
                               cos(radians(self.Pinion.data.get('beta')))+self.data.get('A_a')*   \
                               tan(radians(self.Pinion.data.get('alpha_n')))/cos(radians(   \
                                   self.Pinion.data.get('beta')))})



class Tool:
    """
    Parent Class for all tools.
    """

    # Attributes: tool data
    data = None    # dictionary containing all tool parameters

    def __str__(self):
        """
        Define string conversion of Tool objects

        INPUT parameter:
        -

        OUTPUT:
        string representation of class
        """
        
        outstr = 'rack tool type:\n'
        # output tool type
        outstr += 'root type:   \t'.ljust(15)+self._root_type+'\n'
        outstr += 'protuberance:\t'.ljust(15)
        if self._tip_type=='standard':
             outstr += 'no\n\n'
        else:
             outstr += 'yes\n\n'
             
        # output tool data
        outstr += 'tool data:\n'
        for date in self.data:
             outstr += date.ljust(15)+':\t'+str(self.data.get(date))+'\n'

        return outstr
    

    def __init__(self, tooldata):
        """
        Initialization of Tool-object
        Should be overwritten in derived classes

        INPUT parameters:
        tooldata   : data of tool (dictionary)
        """

        self.data       = deepcopy(tooldata)


    def getToolData(self):
        """
        Return data-attribute of class

        OUTPUT:
        data attribute of class (dictionary)
        """
        
        return self.data


    def setToolData(self, tooldata):
        """
        Set data-attribute of class, overwrite current value

        INPUT parameter:
        tooldata : dictionary, containing geometric data of tool
                   for content, see method __init__
        """
        
        self.__init__(tooldata)


    def updateToolData(self, tooldata):
        """
        Set data-attribute of class, update current value

        INPUT parameter:
        tooldata : dictionary, containing geometric data of gear
                   for content, see method __init__
        """
        
        tempdata = self.data.copy()
        tempdata.update(tooldata)
        self.__init__(tooldata)


    def getY(self, X):

        def objective(param, x_value):
             curve_pnt = self.getCurvePoint(param)
             return curve_pnt[0]-x_value

        param = newton(objective, 0.0, args=(X,))
        point   = self.getCurvePoint(param)

        return point[1]

        


class ToothedRackTool(Tool):
    """
    Class representing a toothed rack tool for gear-cutting (manufacturing of involute gears).
    Derived from Tool-class

    The tool coordinate system is:
    - x-direction:   parallel to datum line of the standard basic rack tooth profile
                     positive to the right
    - y-direction:   parallel to midline of tooth of the standard basic rack tooth profile
                     positive to the top (that is to the dedendum line of the tool's profile)
    - origin:        at intersection point of midline of tooth and datum line
    """
  
    # Attributes: default settings for parameters (DIN 3972, profile II)
    _pr_P0_default      =  0.0    # default value for protuberance (no protuberance)
    _h_aP0_default      =  1.20   # default value for addendum factor of the standard basic rack tooth profile
    _h_fP0_default      =  1.25   # default value for dedendum factor of the standard basic rack tooth profile
    _alpha_P0_default   =  20.0   # default value for pressure angle of the standard basic rack tooth profile
    _alpha_prP0_default =  8.0    # default value for protuberance flank pressure angle on the standard basic rack tooth profile of the tool
    _rho_aP0_default    =  0.20   # default value for crest rounding radius factor on the standard basic rack tooth profile of the tool
    _rho_fP0_default    =  0.20   # tooth root radius factor on the standard basic rack tooth profile of the tool

    # Internal parameters of class
    _tip_type    =  None  # either "protuberance" or "standard"  (string)
    _root_type   =  None  # either "circular" or "chamfered"     (string)
    # characteristic points of half tool profile as x-y-tuples  (2x1-NumPy-arrays)
    _P_R         =  None  # reference point:       datum line cuts flank
    _P_1         =  None  # dedendum start:        right point on dedendum line
    _P_2         =  None  # dedendum end:          left point on dedendum line, start of chamfer or root radius
    _P_3         =  None  # flank start:           upper point of flank, end of chamfer or root radius
    _P_4         =  None  # flank end:             lower point of flank, start of protuberance flank or crest rounding
    _P_5         =  None  # protuberance end:      lower point of protuberance flank, start of crest rounding
                          #                        not set if self._tip_type is "standard"
    _P_6         =  None  # addendum start:        right point of addendum line, end of crest rounding
    _P_7         =  None  # addendum end:          addendum line cuts midline of tooth
    _M_1         =  None  # center root radius:    center of root radius
                          #                        set to [-1., -1.] if self._root_type is "chamfered"
    _M_2         =  None  # center crest rounding: center of crest rounding
    # curve parameters on parametric tool profile of characteristic points (numeric)
    _s_P_1       =  None  # curve parameter of self._P_1 (maximum parameter)
    _s_P_2       =  None  # curve parameter of self._P_2
    _s_P_3       =  None  # curve parameter of self._P_3
    _s_P_4       =  None  # curve parameter of self._P_4
    _s_P_5       =  None  # curve parameter of self._P_5
    _s_P_6       =  None  # curve parameter of self._P_6
    _s_P_7       =  0.0   # curve parameter of self._P_7 (minimum parameter)
    # these parameters are stored as class attributes for performance reasons only



    def __init__(self, tooldata):
        """
        Initialization of ToothedRackTool-object
        Should be overwritten in derived classes

        INPUT parameters:
        tooldata   : data of hob (dictionary)

        possible parameters (keys) in self.data (dictionary):
        m             : normal module (numeric, positive)
        h_aP0         : addendum of the standard basic rack tooth profile (numeric)
        h_fP0         : dedendum of the standard basic rack tooth profile (numeric)
        h_FaP0        : form addendum of the standard basic rack tooth profile (numeric)
        h_FfP0        : form dedendum of the standard basic rack tooth profile (numeric)
        h_prP0        : protuberance tooth depth on the standard basic rack tooth profile of the tool (numeric)
        h_P0          : tooth depth of the standard basic rack tooth profile of the tool (numeric)
        s_P0          : tooth thickness of the standard basic rack tooth profile of the tool (numeric)
        alpha_P0      : pressure angle of the standard basic rack tooth profile of the tool (numeric) [degrees]
        alpha_prP0    : protuberance flank pressure angle on the standard basic rack tooth profile of the tool (numeric) [degrees]
        alpha_KP0     : chamfered flank pressure angle on the standard basic rack tooth profile of the tool (numeric) [degrees]
        rho_aP0       : crest rounding radius on the standard basic rack tooth profile of the tool (numeric)
        rho_fP0       : tooth root radius on the standard basic rack tooth profile of the tool (numeric)
        pr_P0         : protuberance amount on the standard basic rack tooth profile of the tool (numeric)
                        do not supply if tool has no protuberance!

        All input parameters above are arranged in a dictionary. The keys are
        the names of the parameters as listed above.

        There are several possibilities for defining a complete standard basic rack tooth profile:
        1) m has to be provided plus optionally the following
        2) pr_P0, alpha_prP0 for tool with protuberance and/or
        3) alpha_KP0, h_FfP0 for chamfered flank or
        4) rho_fP0
        Some parameters can be left out, but the result might differ
        from your real tool. Missing parameters are estimated if
        possible. Some checks on the values are made, but these are far from
        being complete, thus the user is responsible for supplying valid
        parameters.
        The constructor does not check for unit consistency. The user is
        responsible for supplying all values with consistent units.
        """

        self.data = deepcopy(tooldata)

        # module must be supplied
        if not self.data.has_key('m'):
             raise AttributeError, 'no module supplied'
        # module value check
        if self.data.get('m')<=0:
             raise ValueError, 'non-positive module'
            
        # set pressure angle to default if not supplied
        if not self.data.has_key('alpha_P0'):
             self.data.update({'alpha_P0':self._alpha_P0_default})
        # pressure angle value check
        if self.data.get('alpha_P0')<=0:
             raise ValueError, 'non-positive pressure angle'

        # set type of tool
        if self.data.has_key('alpha_KP0'):
             self._root_type = 'chamfered'
        else:
             self._root_type = 'circular'
             
        if self.data.has_key('pr_P0') and self.data.get('pr_P0')>0:
             self._tip_type = 'protuberance'
        else:
             self._tip_type = 'standard'

        # calculate tip parameters with form addendum if supplied (for tool without protuberance)
        if self._tip_type=='standard' and self.data.has_key('h_FaP0'):
             if self.data.has_key('h_aP0'):
                  rho_aP0 = (self.data.get('h_aP0')-self.data.get('h_FaP0'))/    \
                            (1.0-cos(radians(self.data.get('alpha_P0'))))
                  # crest rounding radius value and consistency check
                  if self.data.has_key('rho_aP0'):
                       if not self.data.get('rho_aP0')==rho_aP0:
                            raise ValueError, 'tip definition inconsistent'
                  else:
                       self.data.update({'rho_aP0':rho_aP0})
             if self.data.has_key('rho_aP0'):
                  h_aP0 = self.data.get('h_FaP0')+self.data.get('rho_aP0')*      \
                          (1.0-cos(radians(self.data.get('alpha_P0'))))
                  # addendum value and consistency check
                  if self.data.has_key('h_aP0'):
                       if not self.data.get('h_aP0')==h_aP0:
                            raise ValueError, 'tip definition inconsistent'
                  else:
                       self.data.update({'h_aP0':h_aP0})

        # calculate tip parameters with form addendum if supplied (for tool with protuberance)
        if self._tip_type=='protuberance' and self.data.has_key('h_FaP0'):
             if self.data.has_key('h_prP0'):
                  h_aP0 = (self.data.get('h_prP0')+self.data.get('h_FaP0'))
                  # addendum value and consistency check
                  if self.data.has_key('h_aP0'):
                       if not self.data.get('h_aP0')==h_aP0:
                            raise ValueError, 'tip definition inconsistent'
                  else:
                       self.data.update({'h_aP0':h_aP0})
             if self.data.has_key('h_prP0') and self.data.has_key('alpha_prP0') and self.data.has_key('rho_aP0'):
                  pr_P0 = (self.data.get('h_prP0')-self.data.get('rho_aP0')*                    \
                          (1.0-cos(radians(self.data.get('alpha_prP0')))))*                     \
                          tan(radians(self.data.get('alpha_P0')-self.data.get('alpha_prP0')))/  \
                          cos(radians(self.data.get('alpha_prP0')))
                  # protuberance value and consistency check
                  if self.data.has_key('pr_P0'):
                       if not self.data.get('pr_P0')==pr_P0:
                            raise ValueError, 'protuberance definition inconsistent'
                  else:
                       self.data.update({'pr_P0':pr_P0})     # this case cannot occur
             if self.data.has_key('h_prP0') and self.data.has_key('alpha_prP0') and self.data.has_key('pr_P0'):
                  rho_aP0 = 1/(1.0-cos(radians(self.data.get('alpha_prP0'))))*                    \
                            (self.data.get('h_prP0')-self.data.get('pr_P0')*                      \
                            cos(radians(self.data.get('alpha_prP0')))/                            \
                            tan(radians(self.data.get('alpha_P0')-self.data.get('alpha_prP0'))))  \
                  # crest rounding radius value and consistency check
                  if self.data.has_key('rho_aP0'):
                       if not self.data.get('rho_aP0')==rho_aP0:
                            raise ValueError, 'tip definition inconsistent'
                  else:
                       self.data.update({'rho_aP0':rho_aP0})
             # the calculation of alpha_prP0 from the other parameters is omitted, as the corresponding equation
             # cannot be solved for alpha_prP0 explicitly (recursive solution necessary)

        # calculate root parameters with form dedendum if supplied (for tool without chamfered flank)
        if self._root_type=='circular' and self.data.has_key('h_FfP0'):
             if self.data.has_key('h_fP0'):
                  rho_fP0 = (self.data.get('h_fP0')-self.data.get('h_FfP0'))/    \
                            (1.0-cos(radians(self.data.get('alpha_P0'))))
                  # tooth root radius value and consistency check
                  if self.data.has_key('rho_fP0'):
                       if not self.data.get('rho_fP0')==rho_fP0:
                            raise ValueError, 'root definition inconsistent'
                  else:
                       self.data.update({'rho_fP0':rho_fP0})
             if self.data.has_key('rho_fP0'):
                  h_fP0 = self.data.get('h_FfP0')+self.data.get('rho_fP0')*      \
                          (1.0-cos(radians(self.data.get('alpha_P0'))))
                  # dedendum value and consistency check
                  if self.data.has_key('h_fP0'):
                       if not self.data.get('h_fP0')==h_fP0:
                            raise ValueError, 'root definition inconsistent'
                  else:
                       self.data.update({'h_fP0':h_fP0})

        # check and set tooth depth, addendum and dedendum
        if not self.data.has_key('h_P0'):
            
            # set addendum to default if not supplied
            if not self.data.has_key('h_aP0'):
                 self.data.update({'h_aP0':self._h_aP0_default*self.data.get('m')})

            # set dedendum to default if not supplied
            if not self.data.has_key('h_fP0'):
                 self.data.update({'h_fP0':self._h_fP0_default*self.data.get('m')})

            # calculate tooth depth
            self.data.update({'h_P0':self.data.get('h_aP0')+self.data.get('h_fP0')})

        else:
            # set addendum to default and calculate dedendum if neither addendum nor dedendum is supplied
            if not self.data.has_key('h_aP0') and not self.data.has_key('h_fP0'):
                 self.data.update({'h_aP0':self._h_aP0_default*self.data.get('m')})
                 self.data.update({'h_fP0':self.data.get('h_P0')-self.data.get('h_aP0')})

            # calculate addendum if only dedendum supplied
            if self.data.has_key('h_fP0') and not self.data.has_key('h_aP0'):
                 self.data.update({'h_aP0':self.data.get('h_P0')-self.data.get('h_fP0')})

            # calculate dedendum if only addedendum supplied
            if self.data.has_key('h_aP0') and not self.data.has_key('h_fP0'):
                 self.data.update({'h_fP0':self.data.get('h_P0')-self.data.get('h_aP0')})

        # tooth depth value check
        if self.data.get('h_P0')<=0:
             raise ValueError, 'non-positive tooth depth'
        # addendum depth value check
        if self.data.get('h_aP0')<=0:
             raise ValueError, 'non-positive addendum'
        # dedendum depth value check
        if self.data.get('h_fP0')<=0:
             raise ValueError, 'non-positive dedendum'
        # check: sum of addendum and dedendum must be tooth depth
        if not self.data.get('h_fP0')+self.data.get('h_aP0')==self.data.get('h_P0'):
             raise ValueError, 'sum of addendum and dedendum does not equal tooth depth'
            
        # set crest rounding radius to default if not supplied
        if not self.data.has_key('rho_aP0'):
             self.data.update({'rho_aP0':self._rho_aP0_default*self.data.get('m')})
        # crest rounding radius value check
        if self.data.get('rho_aP0')<0:
             raise ValueError, 'negative crest rounding radius'

        # tooth root radius and chamfered flank pressure angle are exclusive
        if self.data.has_key('rho_fP0') and self.data.has_key('alpha_KP0'):
             raise AttributeError, 'root type can either be circular or chamfered'
             
        # set tooth root radius to default if not supplied and no chamfered flank pressure angle is defined
        if not self.data.has_key('rho_fP0') and not self.data.has_key('alpha_KP0'):
             self.data.update({'rho_fP0':self._rho_fP0_default*self.data.get('m')})
        # tooth root radius value check
        if self.data.has_key('rho_fP0') and self.data.get('rho_fP0')<0:
             raise ValueError, 'negative tooth root radius'

        # check if chamfered flank definition is sufficient (pressure angle + form dedendum)
        if self.data.has_key('alpha_KP0'):
             if not self.data.has_key('h_FfP0'):
                  raise AttributeError, 'chamfered flank definition not complete'
             elif self.data.get('h_FfP0')<=0:
                  raise ValueError, 'non-positive form dedendum'

        # protuberance value check
        if self.data.has_key('pr_P0'):
             if self.data.get('pr_P0')<0:
                  raise ValueError, 'negative protuberance'
             else:
                  # set protuberance flank pressure angle to default if not supplied
                  if not self.data.has_key('alpha_prP0'):
                       self.data.update({'alpha_prP0':self._alpha_prP0_default})

        # protuberance flank pressure angle value check
        if self.data.has_key('alpha_prP0') and self.data.get('alpha_prP0')<0:
             raise ValueError, 'negative protuberance flank pressure angle'

        # calculate tooth thickness if not supplied
        if not self.data.has_key('s_P0'):
             self.data.update({'s_P0':pi*self.data.get('m')/2})
        # tooth thickness value check
        if self.data.get('s_P0')<=0:
             raise ValueError, 'tooth thickness'

        # calculate characteristic points of rack profile (details see above)
        self._P_R  = array([self.data.get('s_P0')/2.0, 0.0])
        self._P_1  = array([pi*self.data.get('m')/2.0, self.data.get('h_fP0')])
        if self._root_type=='circular':
             self._P_3 = self._P_R + array([tan(radians(self.data.get('alpha_P0'))), 1.0])*    \
                                            (self.data.get('h_fP0')-self.data.get('rho_fP0')*  \
                                            (1.0-sin(radians(self.data.get('alpha_P0')))))
             self._M_1  = self._P_3 + array([cos(radians(self.data.get('alpha_P0'))),          \
                          -sin(radians(self.data.get('alpha_P0')))])*self.data.get('rho_fP0')
             self._P_2  = array([self._M_1[0], self.data.get('h_fP0')])
        if self._root_type=='chamfered':
             self._P_3 = self._P_R + array([tan(radians(self.data.get('alpha_P0'))), 1.0])*    \
                                            self.data.get('h_FfP0')
             self._P_2  = array([self._P_3[0]+tan(radians(self.data.get('alpha_KP0')))*        \
                                 (self.data.get('h_fP0')-self.data.get('h_FfP0')),             \
                                 self.data.get('h_fP0')])
        self._P_7  = array([0.0, -self.data.get('h_aP0')])
        if self._tip_type=='standard':  # use temporary protuberance data for calculation depending on tip type
             pr_calc         = 0.0
             alpha_prP0_calc = self.data.get('alpha_P0')
        else:
             pr_calc         = self.data.get('pr_P0')
             alpha_prP0_calc = self.data.get('alpha_prP0')
        self._M_2 = array([self.data.get('s_P0')/2.0+pr_calc *                                \
                    1.0/cos(radians(self.data.get('alpha_P0')))-                              \
                    tan(radians(self.data.get('alpha_P0')))*(self.data.get('h_aP0')-          \
                    self.data.get('rho_aP0')*(1.0-sin(radians(self.data.get('alpha_P0')))))-  \
                    self.data.get('rho_aP0')*cos(radians(self.data.get('alpha_P0'))),         \
                    -self.data.get('h_aP0')+self.data.get('rho_aP0')])
        self._P_5 = self._M_2+self.data.get('rho_aP0')*                                       \
                    array([cos(radians(alpha_prP0_calc)), -sin(radians(alpha_prP0_calc))])
        if self._tip_type=='protuberance':
             self._P_4 = self._P_5+(self.data.get('pr_P0')/sin(radians(self.data.get('alpha_P0')- \
                         self.data.get('alpha_prP0')))-tan(radians(self.data.get('alpha_P0')-     \
                         self.data.get('alpha_prP0'))/2.0)*self.data.get('rho_aP0'))*             \
                         array([sin(radians(self.data.get('alpha_prP0'))),                        \
                                cos(radians(self.data.get('alpha_prP0')))])
        self._P_6  = array([self._M_2[0], -self.data.get('h_aP0')])

        # calculate form dedendum/addendum
        if self._root_type=='circular':
             self.data.update({'h_FfP0':self._P_3[1]})
        if self._tip_type=='standard':
             self.data.update({'h_FaP0':-self._P_5[1]})
        else:
             self.data.update({'h_FaP0':-self._P_4[1]})
             self.data.update({'h_prP0':self.data.get('h_aP0')-self.data.get('h_FaP0')})

        # get parameters of characteristic points on parametric rack profile curve
        self._s_P_1 = self._getParameterLimit(self._P_1)
        self._s_P_2 = self._getParameterLimit(self._P_2)
        self._s_P_3 = self._getParameterLimit(self._P_3)
        self._s_P_4 = self._getParameterLimit(self._P_4)
        self._s_P_5 = self._getParameterLimit(self._P_5)
        self._s_P_6 = self._getParameterLimit(self._P_6)
      

    def getCurvePoint(self, parameter, difforder=0):
        """
        Returns a point (Ps) on the tool's profile in Cartesian coordinates (definition of
        coordinate system see above). The profile is described as a parametric curve
        in two-dimensional space, the normal plane of the tool. The curve parameter is
        normalized with regard to curve length (--> _s).
        As an alternative return array (dP_ds), the derivative of the curve with respect to the curve
        parameter at the location defined by the supplied curve parameter is provided.
        This array defines the tangent vector on the curve.
        The third option is to return an array (d2P_ds2) that is the second derivative.
        The curve parameter is negative for profile points left of the symmetry line
        of a tooth of the tool (negative x-value) and positive for points right of
        the symmetry line (positive x-value). It is zero at the cutting point of the
        tool's profile with the symmetry line of the tooth.
        If the absolute value of the parameter is out of bounds an error is raised.

        INPUT parameter:
        parameter:  curve parameter (numeric)
        difforder:  differentiate difforder times with respect to curve parameter (0, 1, or 2)
                    optional - set to 0 if not supplied
     
        OUTPUT:
        curve_property:    the requested property of the profile (numeric, 2x1-NumPy-array)
        """ 

        # module imports
        from numpy.linalg import norm

        # same procedure for positive and negative x-values (mirroring is done later)
        s = abs(parameter)

        # point on dedendum line
        if s<=self._s_P_6:
             Ps      = self._P_7+s*(self._P_6-self._P_7)/norm(self._P_6-self._P_7)
             dP_ds   = (self._P_6-self._P_7)/norm(self._P_6-self._P_7)
             d2P_ds2 = array([0.0, 0.0])
             
        # point on crest rounding
        elif s<=self._s_P_5:
             s_      = s-self._s_P_6
             Ps      = self._M_2+self.data.get('rho_aP0')*array([sin(s_/self.data.get('rho_aP0')),   \
                                                                 -cos(s_/self.data.get('rho_aP0'))])
             dP_ds   = array([cos(s_/self.data.get('rho_aP0')), sin(s_/self.data.get('rho_aP0'))])
             d2P_ds2 = 1/self.data.get('rho_aP0')*array([-sin(s_/self.data.get('rho_aP0')),   \
                                                                 cos(s_/self.data.get('rho_aP0'))])
        # point on protuberance flank or flank
        elif s<=self._s_P_3:
             # tool with protuberance
             if self._tip_type=='protuberance':
                  # point on protuberance flank
                  if s<=self._s_P_4:
                       s_      =  s-self._s_P_5
                       Ps      =  self._P_5+s_*(self._P_4-self._P_5)/norm(self._P_4-self._P_5)
                       dP_ds   = (self._P_4-self._P_5)/norm(self._P_4-self._P_5)
                  # point on flank
                  else:
                       s_      = s-self._s_P_4
                       Ps      = self._P_4+s_*(self._P_3-self._P_4)/norm(self._P_3-self._P_4)
                       dP_ds   = (self._P_3-self._P_4)/norm(self._P_3-self._P_4)
             # tool without protuberance
             else:
                  # point on flank
                  s_      = s-self._s_P_5
                  Ps      = self._P_5+s_*(self._P_3-self._P_5)/norm(self._P_3-self._P_5)
                  dP_ds   = (self._P_3-self._P_5)/norm(self._P_3-self._P_5)
             d2P_ds2 = array([0.0, 0.0])
             
        # point on root rounding or chamfered flank
        elif s<=self._s_P_2:
             # tool with circular root rounding
             if self._root_type=='circular':
                  # point on root rounding
                  s_      = s-self._s_P_3+acos((self._P_3[0]-self._M_1[0])/-self.data.get('rho_fP0'))*   \
                            self.data.get('rho_fP0')
                  Ps      = self._M_1+self.data.get('rho_fP0')*array([-cos(s_/self.data.get('rho_fP0')), \
                                                                      sin(s_/self.data.get('rho_fP0'))])
                  dP_ds   = array([sin(s_/self.data.get('rho_fP0')), cos(s_/self.data.get('rho_fP0'))])
                  d2P_ds2 = 1/self.data.get('rho_fP0')*array([cos(s_/self.data.get('rho_fP0')), \
                                                                      -sin(s_/self.data.get('rho_fP0'))])
             # tool with chamfered root edge
             else:
                  # point on chamfered flank
                  s_      = s-self._s_P_3
                  Ps      = self._P_3+s_*(self._P_2-self._P_3)/norm(self._P_2-self._P_3)
                  dP_ds   = (self._P_2-self._P_3)/norm(self._P_2-self._P_3)
                  d2P_ds2 = array([0.0, 0.0])
        # point on addendum line (that is extended infinitely for parameters > self._s_P_1)
        else:
                  s_      = s-self._s_P_2
                  Ps      = self._P_2+s_*(self._P_1-self._P_2)/norm(self._P_1-self._P_2)
                  dP_ds   = (self._P_1-self._P_2)/norm(self._P_1-self._P_2)
                  d2P_ds2 = array([0.0, 0.0])

        # return point on profile, corresponding tangent and normal vector
        mirror = array([[-1, 0], [0, 1]])
        if parameter>=0:      # positive parameter: no transformation necassary
             if difforder==0:
                  return Ps
             elif difforder==1:
                  return dP_ds
             elif difforder==2:
                  return d2P_ds2
             else:
                  raise AttributeError, 'differentiation order maust be 0, 1 or 2'
        else:                 # negative parameter: mirror at midline of tooth and maintain continuous derivative
             if difforder==0:
                  return dot(mirror, Ps)
             elif difforder==1:
                  return -1*dot(mirror, dP_ds)
             elif difforder==2:
                  return dot(mirror, d2P_ds2)
             else:
                  raise AttributeError, 'differentiation order maust be 0, 1 or 2'


    def getY(self, X):

        s = abs(X)

        # point on dedendum line
        if s<=self._P_6[0]:
             y = self._P_6[1]
             
        # point on crest rounding
        elif s<=self._P_5[0]:
             y = self._M_2[1]-self.data.get('rho_aP0')*cos(asin((s-self._M_2[0])/self.data.get('rho_aP0')))
             
        # point on protuberance flank or flank
        elif s<=self._P_3[0]:
             # tool with protuberance
             if self._tip_type=='protuberance':
                  # point on protuberance flank
                  if s<=self._P_4[0]:
                       y = self._P_5[1]+(s-self._P_5[0])*(self._P_4[1]-self._P_5[1])/(self._P_4[0]-self._P_5[0])
                  # point on flank
                  else:
                       y = self._P_4[1]+(s-self._P_4[0])*(self._P_3[1]-self._P_4[1])/(self._P_3[0]-self._P_4[0])
             # tool without protuberance
             else:
                  # point on flank
                  y = self._P_5[1]+(s-self._P_5[0])*(self._P_3[1]-self._P_5[1])/(self._P_3[0]-self._P_5[0])
             
        # point on root rounding or chamfered flank
        elif s<=self._P_2[0]:
             # tool with circular root rounding
             if self._root_type=='circular':
                  # point on root rounding
                  y = self._M_1[1]+self.data.get('rho_fP0')*sin(acos((self._M_1[0]-s)/self.data.get('rho_fP0')))
             # tool with chamfered root edge
             else:
                  # point on chamfered flank
                  y = self._P_3[1]+(s-self._P_3[0])*(self._P_2[1]-self._P_3[1])/(self._P_2[0]-self._P_3[0])
                  
        # point on addendum line (that is extended infinitely for parameters > self._s_P_1)
        else:
             y  = self._P_2[1]

        # return y-value of point on profile
        return y


    def _getParameterLimit(self, point=NaN):
        """
        Returns the value of the curve parameter (see method 'getToolPoint') that
        is the maximum absolute value allowed. If the absolute value of the curve
        parameter exceeds this limit, the point returned by method 'getToolPoint'
        is not on the actual tooth (width equals pitch) any more.

        INPUT parameter:
        point:      curve parameter of this characteristic point is returned (2x1-NumPy-array)
                    optional: maximum value (corresponds to self._P_1) is returned
                              if not supplied
     
        OUTPUT:
        param_limit:  maximum absolute value of curve parameter (numeric)
        """

        # module imports
        from numpy.linalg import norm
        from numpy.ma import allequal

        if allequal(point,self._P_7):           # curve parameter of point self._P_7
             return 0.0
          
        s_max  = norm(self._P_6-self._P_7)      # length of line between tooth's symmetry line and begin of crest rounding
        if allequal(point,self._P_6):           # curve parameter of point self._P_6
             return s_max
          
        s_max += asin((self._P_5[0]-self._M_2[0])/self.data.get('rho_aP0'))*         \
                 self.data.get('rho_aP0')       # arc length of crest rounding (caution: expression ambiguous)
        if allequal(point,self._P_5):           # curve parameter of point self._P_5
             return s_max
          
        if self._tip_type=='protuberance':
             s_max += norm(self._P_4-self._P_5) # length of protuberance flank
             if allequal(point,self._P_4):      # curve parameter of point self._P_4
                  return s_max
             s_max += norm(self._P_3-self._P_4) # length of flank
        else:
             if allequal(point,self._P_4):      # curve parameter of point self._P_4 (point non-existent for tools without protuberance)
                  return None
             s_max += norm(self._P_3-self._P_5) # length of flank
             
        if allequal(point,self._P_3):           # curve parameter of point self._P_3
             return s_max
          
        if self._root_type=='circular':
             s_max += (acos((self._P_2[0]-self._M_1[0])/-self.data.get('rho_fP0'))-   \
                       acos((self._P_3[0]-self._M_1[0])/-self.data.get('rho_fP0')))*  \
                       self.data.get('rho_fP0')  # arc length of root rounding (caution: expression ambiguous)
        else:
             s_max += norm(self._P_2-self._P_3)  # length of chamfer flank
        if allequal(point,self._P_2):            # curve parameter of point self._P_2
             return s_max
             
        s_max += norm(self._P_1-self._P_2)       # length of line between start of root rounding and symmetry line of gap
        
        return s_max


class Machine:
    """
    Parent Class for all gear manufacturing machines.
    """

    # Attributes: machine data
    data  = None    # dictionary containing all machine parameters and settings
    tool  = None    # Tool-instance representing the cutter
    blank = None    # pre-product from which the tooth-gap is cutted out (TColgp_Array1OfPnt2d, pythonOCC)

    def __str__(self):
        """
        Define string conversion of Machine objects

        INPUT parameter:
        -

        OUTPUT:
        string representation of class
        """
        
        outstr = 'machine data:\n'

        # output machine data
        for date in self.data:
             outstr += date.ljust(15)+':\t'+str(self.data.get(date))+'\n'

        # output tool data
        outstr += '\ntool data:\n'
        outstr += str(self.tool)

        # output blank coordinates
        if self.blank:
             outstr += '\nblank data:\n'
             outstr += str(self.blank)

        return outstr
    

    def __init__(self, machinedata, tool, blank=None):
        """
        Initialization of Machine-object
        Should be overwritten in derived classes

        INPUT parameters:
        machinedata   : data of machine (dictionary)
        """

        self.data  = deepcopy(machinedata)
        self.tool  = tool
        self.blank = blank


    def getMachineData(self):
        """
        Return data-attribute of class

        OUTPUT:
        data attribute of class (dictionary)
        """
        
        return self.data


    def setMachineData(self, machinedata):
        """
        Set data-attribute of class, overwrite current value

        INPUT parameter:
        machinedata : dictionary, containing settings of machine
                      for content, see method __init__
        """
        
        self.__init__(machinedata)


    def updateMachineData(self, machinedata):
        """
        Set data-attribute of class, update current value

        INPUT parameter:
        machinedata : dictionary, containing settings of machine
                      for content, see method __init__
        """
        
        tempdata = self.data.copy()
        tempdata.update(machinedata)
        self.__init__(machinedata)


    def setTool(self, tool):
        """
        Set tool attribute

        INPUT parameters:
        tool : cutting tool of machine (Tool instance)
        """
        
        # tool: type and value check
        if not isinstance(tool, Tool):
             raise TypeError, 'instance of Tool-class expected'
        if not isinstance(tool, ToothedRackTool):
             raise TypeError, 'non rack-type cutters not implemented'

        self.tool = tool


    def getTool(self):
        """
        Get tool attribute

        INPUT parameters:
        -

        OUTPUT:
        cutting tool of machine (Tool instance)
        """

        return self.tool


    def setBlank(self, blank):
        """
        Set blank attribute

        INPUT parameters:
        blank : blank for machine (Blank)
        """
        
        # blank: type and value check
        if not isinstance(blank, Blank):
             raise TypeError, 'instance of Blank-class expected'

        self.blank = blank


    def getBlank(self):
        """
        Get blank attribute

        INPUT parameters:
        -

        OUTPUT:
        blank for machine (Blank)
        """

        return self.blank
        
            

class GearHobber(Machine):
    """
    Class representing a gear hobber for gear manufacturing (manufacturing of cylindrical gears).
    Derived from Machine-class

    The 2d machine coordinate system is:
    - x-axis     :   positive to the right
    - y-axis     :   symmetry line of tooth gap of manufactured gear
                     positive to the top (that is to the dedendum line of the tool's profile)
    - origin:        at center of manufactured gear
    """

    # Attributes: default settings for parameters
    _x_default      =  0.0   # default value for addendum modification factor
    _beta_default   =  0.0   # default value for helix angle
    _A_s_default    =  0.0   # default value for tooth thickness allowance
    _q_default      =  0.0   # default value for machining allowance
    # default settings for hobbing simulation settings
    _no_pnts_default       = [10, 10, 150, 20, 50, 10]   # default values for number of points for tool flank discretization (dedendum, root rounding/chamfer, flank, protuberance flank, tip rounding, addendum)
    # default settings for envelope check
    _no_supp_pnts_default   = 100     # default value for number of support points used to construct spline
    _zero_tolerance_default = 1.0E-8  # default value for zero tolerance (below this value is assumed to be zero)
    
   
    def __init__(self, machinedata, tool, blank=None):
        """
        Initialization of GearHobber-object

        INPUT parameters:
        machinedata   : machine settings (dictionary)

        possible parameters (keys) in self.data (dictionary):
        z        : number of teeth (numeric, integer)
        d        : pitch diameter (numeric)
                   one of the two parameters above z or d has to be supplied
        x        : addendum modification factor (numeric)
                   optional - set equal 0.0 if not supplied
        x_E      : generating addendum modification coefficient (numeric)
                   optional - calculated if not supplied
        beta     : helix angle (numeric)[degrees]
                   optional - set equal 0.0 if not supplied
        A_s      : tooth thickness allowance in normal cross-section (numeric, negative)
                   optional - set equal 0.0 if not supplied 
        q        : machining allowance (numeric)
                   optional - set equal 0.0 if not supplied
        
        All input parameters above are arranged in a dictionary. The keys are
        the names of the parameters as listed above.
        Some checks on the values are made, but these are far from
        being complete. Thus the user is responsible for supplying valid
        parameters.
        The constructor does not check for unit consistency. The user is
        responsible for supplying all values with consistent units.
        """

        self.data       = deepcopy(machinedata)
        self.setTool(tool)
        if blank:
             self.setBlank(blank)

        # check if number of teeth or pitch diameter is supplied
        if not self.data.has_key('z') and not self.data.has_key('d'):
             raise AttributeError, 'either number of teeth or pitch diameter must be supplied'

        # number of teeth: value check
        if self.data.has_key('z') and not type(self.data.get('z'))==type(1):
             raise TypeError, 'number of teeth not integer'

        # calculate pitch diameter if not supplied
        if not self.data.has_key('d'):
             self.data.update({'d':self.tool.data.get('m')*self.data.get('z')/    \
                               cos(radians(self.data.get('beta')))})

        # helix angle: set to default if not supplied   
        if not self.data.has_key('beta'):
             self.data.update({'beta':self._beta_default})

        # addendum modification factor: set to default if not supplied   
        if not self.data.has_key('x'):
             self.data.update({'x':self._x_default})

        # tooth thickness allowance: set to default if not supplied   
        if not self.data.has_key('A_s'):
             self.data.update({'A_s':self._A_s_default})
        # tooth thickness allowance: value check
        else:
             if not self.data.get('A_s')<=0:
                  raise ValueError, 'tooth thickness allowance positive'

        # machining allowance: set to default if not supplied   
        if not self.data.has_key('q'):
             self.data.update({'q':self._q_default})

        # calculate generating addendum modification coefficient:
        if not self.data.has_key('x_E'):
            self.data.update({'x_E':self.data.get('x')+(self.data.get('A_s')/2/tan(radians(self.tool.data.get('alpha_P0')))+   \
                              self.data.get('q')/sin(radians(self.tool.data.get('alpha_P0'))))/self.tool.data.get('m')})


    def _getToolPoint(self, parameter, difforder=0):
        """
        Returns a point (Ps) or the tangent (dP_ds) on the tool's profile in Cartesian coordinates (definition of
        coordinate system see above). The tool profile as seen on a cutting transverse-plane is considered.
        All length in x-direction are scaled to enable helical gear generation.
        This is basically a wrapper-function for self.tool._getCurvePoperty()

        INPUT parameter:
        parameter:  curve parameter (numeric)
        difforder:  differentiate difforder times with respect to curve parameter (0, 1, or 2)
                    optional - set to 0 if not supplied
     
        OUTPUT:
        tool_property:    the requested property of the profile (numeric, 2x1-NumPy-array)
        """

        if difforder==0:
            
            # point on flank for supplied curve parameter in  normal cross-section
            P     = self.tool.getCurvePoint(parameter, 0)

            # tool point in transverse cross-section coordinate system and considering for profile modification
            P_t     = P/array([cos(radians(self.data.get('beta'))), 1.0])+array([0.0, self.data.get('x_E')*self.tool.data.get('m')])

            return P_t

        elif difforder==1:
            
            # tangent vector for supplied curve parameter
            dP_ds = self.tool.getCurvePoint(parameter, 1)

            # projection of tool point and tangent vector in transverse cross-section
            dP_ds_t = dP_ds/array([cos(radians(self.data.get('beta'))), 1.0])

            return dP_ds_t

        elif difforder==2:
            
            # second derivative of curve for supplied curve parameter
            d2P_ds2 = self.tool.getCurvePoint(parameter, 2)

            # projection of tool point and tangent vector in transverse cross-section
            d2P_ds2_t = d2P_ds2/array([cos(radians(self.data.get('beta'))), 1.0])

            return d2P_ds2_t

        else:
            raise AttributeError, 'differentiation order must be 0, 1, or 2'


    def _getToolY(self, X):
        """

        INPUT parameter:
     
        OUTPUT:
        Y:
        """
        
        # y-value of point on flank for supplied curve parameter in  normal cross-section
        y     = self.tool.getY(X*cos(radians(self.data.get('beta'))))

        # y-value of tool point in transverse cross-section coordinate system and considering for profile modification
        y_t     = y+self.data.get('x_E')*self.tool.data.get('m')

        return y_t


    def _transformToMachineCoords(self, phi, element, elemtype='point'):
        """
        Transforms a point or vector from the tool coordinate system (with origin attached to the tool) to
        the workpiece coordinate system (origin is attached to the gear that is to be manufactured).

        INPUT parameter:
        phi:        actual position parameter (angle) of cutter (numeric) [radians]
        element:    point or vector in tool coordinate system (2x1-NumPy-array)
        elemtype:   indicator for the type of element to transform (string)
                    'point' :    transform a point, i.e. rotation and translation is applied (default)
                    'vector':    transform a vector, i.e. rotation only is applied
     
        OUTPUT:
        transelem:  point or vector in workpiece coordinate system (2x1-NumPy-array)
        """

        if isinstance(phi, ndarray):
            phi = phi[0]

        # radius of gear centrode
        rho = self.data.get('d')/2.0

        # rotation matrix from rack cutter to gear coordinate system
        L = array([[cos(phi), sin(phi)], [-sin(phi), cos(phi)]])
        
        # translation vector from rack cutter to gear coordinate system
        r = rho*array([sin(phi)-phi*cos(phi), cos(phi)+phi*sin(phi)])
       
        # point in workpiece coordinate system
        if elemtype=='point':
             return dot(L, element)+r
        # vector in workpiece coordinate system
        elif elemtype=='vector':
             return dot(L, element)            
        else:
             raise AttributeError, 'unknown option for element type'


    def _transformToToolCoords(self, phi, element, elemtype='point'):
        """
        Transforms a point or vector from the workpiece coordinate system (origin is attached to the gear that
        is to be manufactured) to the tool coordinate system (with origin attached to the tool).

        INPUT parameter:
        phi:        actual position parameter (angle) of cutter (numeric) [radians]
        element:    point or vector in workpiece coordinate system (2x1-NumPy-array)
        elemtype:   indicator for the type of element to transform (string)
                    'point' :    transform a point, i.e. rotation and translation is applied (default)
                    'vector':    transform a vector, i.e. rotation only is applied
     
        OUTPUT:
        transelem:  point or vector in tool coordinate system (2x1-NumPy-array)
        """

        if isinstance(phi, ndarray):
            phi = phi[0]

        # radius of gear centrode
        rho = self.data.get('d')/2.0

        # rotation matrix from gear to rack cutter coordinate system
        L = array([[cos(phi), -sin(phi)], [sin(phi), cos(phi)]])
        
        # translation vector from rack cutter to gear coordinate system
        r = -rho*array([sin(phi)-phi*cos(phi), cos(phi)+phi*sin(phi)])
       
        # point in tool coordinate system
        if elemtype=='point':
             return dot(L, element+r)
        # vector in tool coordinate system
        elif elemtype=='vector':
             return dot(L, element)            
        else:
             raise AttributeError, 'unknown option for element type'


    def _getVelocityVector(self, phi, parameter):
        """
        Returns the relative velocity vector between tool and workpiece coordinate system of a point
        on the tool profile. It is described in the tool's coordinate system.

        INPUT parameter:
        phi:        actual position parameter (angle) of cutter (numeric) [radians]
        parameter:  curve parameter (numeric)
     
        OUTPUT:
        v:          relative velocity vector (numeric, 2x1-NumPy-array)
        """

        # radius of gear centrode
        rho = self.data.get('d')/2.0

        # tool point in transverse cross-section coordinate system
        P_t = self._getToolPoint(parameter, 0)
        
        # relative velocity vector
        v = array([P_t[1], rho*phi-P_t[0]])

        return v
        

    def _getNormalVector(self, parameter, difforder=0):
        """
        Returns the normal vector on the tool profile curve or its first derivative with respect to
        the curve parameter. It is described in the tool's coordinate system.

        INPUT parameter:
        parameter:  curve parameter (numeric)
        difforder:  differentiate difforder times with respect to curve parameter (0 or 1)
                    optional - set to 0 if not supplied
     
        OUTPUT:
        N:          the requested vector, normal vector or its first derivative (numeric, 2x1-NumPy-array)
        """

        # tangent vector or 2nd derivative of tool profile curve in transverse cross-section coordinate system
        dnP_dsn_t = self._getToolPoint(parameter, difforder+1)

        # multiply this matrix with a vector to get orthogonal vector
        orthog = -array([[0.0, -1.0], [1.0, 0.0]])

        # normal vector on tool flank in transverse cross-section
        N = dot(orthog, dnP_dsn_t)
            
        if difforder==0 or difforder==1:
            return N
        else:
            raise AttributeError, 'differentiation order must be 0 or 1'
        
    
    def _FlankPointObjectiveFnc(self, phi, parameter):
        """
        Find rotation angle (phi) of tool for which the point with parameter s on tool is possibly a
        point of the final gear (workpiece). The necessary but not sufficient condition is that
        the vector of this tool point to the instant center of rotation is orthogonal to the tangent
        vector on the tool's surface. This is the condition for an envelope point.
        The reason why phi is searched for a fixed s (and not vice versa) is numerical stability
        (convergence to global minimum in optmization is guaranteed with gradient method).
        This is the objective function for the optimizer (in self._findFlankPoint).

        INPUT parameter:
        phi:          actual position parameter (angle) of cutter (numeric) [radians]
        parameter:    actual curve parameter of tool (numeric)
     
        OUTPUT:
        p: condition value (numeric)
        """

        # normal vector on flank and relative velocity vector
        N = self._getNormalVector(parameter, 0)
        v = self._getVelocityVector(phi, parameter)
        
        # objective function (normal vector on tool surface is perpendicular to relative velocity vector)
        # caution: this condition is necessary but not sufficient, post-processing required
        p = dot(N,v)

        return p


    def _PointInteriorMeasure(self, phi, point):
        """
        Under the assumption that the tool extends infinitely beyond the root a simple measure
        can be defined indicating whether the point is interior or exterior to the tool at the
        position determined by phi.
        For standard gear cutters the distance in y-direction from the point transformed to the
        tool coordinate system to the tool's limiting curve is unique. If this distance is less
        than zero, the point is defined to be interior. It is exterior otherwise.
        The distance computed by this method is in general not the shortest distance to the tool's
        surface but has the same roots and is preferred here as it can be calculated with very
        little effort.

        INPUT parameter:
        phi:        actual position parameter (angle) of cutter (numeric) [radians]
        point:      point to be tested represented in workpiece coordinate system (2x1-NumPy-array)
     
        OUTPUT:
        d:          relative distance of point from reference point (numeric)
        """

        # transform point to tool coordinate system
        point_tool = self._transformToToolCoords(phi, point, elemtype='point')

        # distance of point from tool in y-direction
        d = self._getToolY(point_tool[0])-point_tool[1]
        
        return d
      


    def _findFlankPoint(self, parameter):
        """
        Find rotation angle (phi) of tool for which the point with parameter s on tool is possibly a
        point of the final gear (workpiece).

        INPUT parameter:
        parameter:    actual curve parameter of tool (numeric)
     
        OUTPUT:
        P_c:      possible workpiece point in machine coordinate system (2x1-NumPy-array)
        phi_opt:  cutter angle at which condition compliant (numeric)[radians]
        """

        # find root of objective function (necessary condition for existence of envelope to family of curves)
        phi_opt = newton(self._FlankPointObjectiveFnc, 0.0, None, (parameter,), tol=1.0e-10, maxiter=500)

        # get flank point corresponding to parameter phi_opt
        P = self._getToolPoint(parameter, 0)

        # transform point on flank to machine coordinate system
        P_c = self._transformToMachineCoords(phi_opt, P, elemtype='point')

        return P_c, phi_opt

        
    def _checkEnvelope(self, point, no_supp_pnts=None):
        """
        Check if a point is member of the envelope of the family of curves defined by
        the tool surface that is transformed to the workpiece coordinate system and
        the necessary condition for envelope points (check sufficient condition). 
        The strategy is to observe the point while the tool cuts through the workpiece.
        A measure for the point's distance from the tool surface is introduced. A negative
        sign of the measure indicates that the point is interior to the tool's contour.
        The global minimum of this distance over the workpiece angle 'phi' is searched.
        If the global minimum is negative, point is not a member of the envelope as this
        means it is inside the tool at certain postions of the cutter and thus cut away.
        Finding global minima with common methods is unreliable and time-consuming, thus
        a different approach is chosen. The distance over 'phi' is a relatively smooth
        function that is firstly approximated by a spline. All roots of the slope of this
        spline can be found easily with standard methods. In a second step, the minima of
        the original distance function are searched with a local optimization algorithm
        and the roots (that are the extrema of the interpolated function) as initial values.
        If any of these minima is negative, the point is not member of the envelope.
        In tests this method proved to be very reliable if the knot-grid for the spline
        generation was not too coarse.
        Remark: it is not necassary to find the global minimum of the distance function.
                It must only be identified if there are any negative values of the distance
                function. Above algorithm serves this purpose, it doesn't perform global
                optimization!
                Alternative methods, like dividing into many intervals for minimum search
                or using global optimizers, proved to be computaionally too expensive.

        INPUT parameter:
        point:         point that is to be checked (1x2-NumPy-array)
        no_supp_pnts:  number of support points used for the spline interpolation of the
                       distance function. Accuracy and computational effort increase for
                       larger values (integer)
                       optional - set to self._no_supp_pnts_default if not supplied
     
        OUTPUT:
        condition:  True if point is envelope point, false otherwise (boolean)
        """

        # set number of support points to default if not supplied
        if no_supp_pnts==None:
             no_supp_pnts = self._no_supp_pnts_default
             
        # calculate range in that negative distances must be searched (outside of that range the tool is certainly
        # completely above tip circle of gear and thus cannot cut this tooth of workpiece any more)
        phi_lim = self._getEffectivePhiMax()

        # create linear spaced grid for knots of spline (abscissa)
        phi = linspace(-phi_lim, phi_lim, no_supp_pnts)

        # values for knots of spline (ordinate of original distance function)
        distance = map(self._PointInteriorMeasure, phi, [point]*no_supp_pnts)

        # create non-smoothing spline approximating the distance function (point from tool surface)
        distance_spline = InterpolatedUnivariateSpline(phi, distance)

        # compute approximation for derivates of distance spline
        derivatives = map(distance_spline, phi, [1]*no_supp_pnts)

        # create non-smoothing spline approximating the slope function of distance spline
        derivative_spline = InterpolatedUnivariateSpline(phi, derivatives)

        # compute roots of slope spline --> extrema of distance distance spline
        derivative_roots = derivative_spline.roots()

        # conduct local minimum search on the original distance function ith any extremum identified as initial value
        for extremum in derivative_roots:
            phi_min = fmin_cg(self._PointInteriorMeasure, extremum, args=(point,), disp=False)
            min_objective = self._PointInteriorMeasure(phi_min, point)
            if min_objective<-self.tool.data.get('m')*self._zero_tolerance_default*1.0E-3:       # if minimum is negative, point is not member of envelope
                return False

        return True    # no negative values in distance function found --> point is assumed to be member of envelope


    def _getEffectivePhiMax(self):
        """
        Find the maximum positive workpiece rotation angle that might affect generated gear shape.
        A rectangular bounding-box is drawn around one tooth of the rack-cutter. When a corner of
        the bounding-box touches the tip circle of the generated gear so that the complete box is
        outside the tip circle, contact between the corresponding cutter tooth and the workpiece
        cannot occur in the further manufacturing process.
     
        OUTPUT:
        phi_lim:        maximum positive workpiece rotation angle that might affect generated
                        gear shape (numeric)
        """

        # calculate range in that negative distances must be searched (outside of that range the tool is certainly
        # completely above tip circle of gear and thus cannot cut this tooth of workpiece any more)
        d_a = self.data.get('d')/2.0 + self.data.get('x_E')*self.tool.data.get('m') + self.tool.data.get('h_fP0')
        d_f = self.data.get('d')/2.0 + self.data.get('x_E')*self.tool.data.get('m') - self.tool.data.get('h_aP0')
        s_lim = sqrt(d_a**2 - d_f**2) + abs(self._getToolPoint((self.tool._getParameterLimit()))[0])
        phi_lim = s_lim/self.data.get('d')*2.0

        return phi_lim
            

    def _findLimitingPoint(self, param_env_false, param_env_true):
        """
        Find the exact solution (within numerical accuracy) for a point on the envelope that is
        not differentiable. Thus the envelope has a sharp bend there and undercutting occurs.
        The algorithm is searching for the tool's curve parameter where the sufficient condition
        for the point membership to the envelope to the family of curves (defined by the tool's
        shape, the kinematics of the manufacturing process and the necassary condition) switches
        from false to true.

        INPUT parameter:
        param_env_false:     tool curve parameter on that side of the singular parameter where the
                             sufficient condition evaluates to false (numeric)
        param_env_true:      tool curve parameter on that side of the singular parameter where the
                             sufficient condition evaluates to true (numeric)
     
        OUTPUT:
        P_true:              singular point in workpiece coordinate system (1x2-NumPy-array)
        """

        # initial error value so that while-loop body is evaluated
        error = 1

        # evaluate point corresponding to mean value of input parameters and compliant to
        # necessary envelope condition. This point becomes the new param_env_true or
        # param_env_false depending on the return value of the sufficient envelope condition
        # check. Repeat until error criterion (distance between the two limiting parameters)
        # is met.
        while error>1.0E-10:
            param_mean        = (param_env_false+param_env_true)/2.0
            P_mean, phi_mean  = self._findFlankPoint(param_mean)
            env_cond          = self._checkEnvelope(P_mean)

            if env_cond:
                 param_env_true = param_mean
            else:
                 param_env_false = param_mean

            error = abs(param_env_true-param_env_false)

        # return the point that is on the side of the singular point where the sufficient
        # envelope condition evaluates to true
        P_true, phi_true = self._findFlankPoint(param_env_true)

        return P_true
    

    def _subtractFromBlank(self, first_array, second_array):
        """
        Find the points that are part of the inner envelope of two curves.
        The algorithm computes the differerence in radial distance from gear center.
        Which curve is inner envelope is decided based on the sign of the distance.

        INPUT parameter:
        first_array:     points on master curve (Nx2-NumPy-array)
                         must contain point with x=0.0 as last column!
        second_array:    points on second curve (Mx2-NumPy-array)
                         must contain point with x=0.0 as last column!
     
        OUTPUT:
        envelope:        inner envelope of the two input arrays (Kx2-NumPy-array)
        """

        # module imports
        from scipy.optimize import bisect
        from scipy.interpolate import interp1d

        # check if both input arrays are valid
        if not size(first_array, 0)>1:
             raise ValueError, 'first input array contains insufficient number of points'
        if not size(second_array, 0)>1:
             raise ValueError, 'second input array contains insufficient number of points'
        if not size(first_array, 1)==2:
             raise ValueError, 'first input array dimension error'
        if not size(second_array, 1)==2:
             raise ValueError, 'second input array dimension error'

        # convert both input arrays (in reverse row order) to polar coordinates
        first_array_polar = self._removeTrailingZeros(first_array[::-1,:])  
        if self.data.get('z')<0.0:     # necessary to handle internal gears correctly
             first_array_polar = first_array_polar[::-1,:]
        for index in range(0, size(first_array_polar, 0)):
             [r, phi] = CartesianCoordinatesToPolarCoordinates(first_array_polar[index, 0], first_array_polar[index, 1])
             first_array_polar[index, 0] = sign(self.data.get('z'))*r
             first_array_polar[index, 1] = phi

        second_array_polar = self._removeTrailingZeros(second_array[::-1,:])
        for index in range(0, size(second_array_polar, 0)):
             [r, phi] = CartesianCoordinatesToPolarCoordinates(second_array_polar[index, 0], second_array_polar[index, 1])
             second_array_polar[index, 0] = sign(self.data.get('z'))*r
             second_array_polar[index, 1] = phi
             
        # linear interpolation functions for arrays in polar coordinates
        first_array_polar_interp = append(append(array([[first_array_polar[0,0], -1.0E99]]), first_array_polar, axis=0),  \
                                          array([[first_array_polar[-1,0], 1.0E99]]),axis=0)
        first_array_interpolation = interp1d(first_array_polar_interp[:,1], first_array_polar_interp[:,0])
        second_array_polar_interp = append(append(array([[second_array_polar[0,0], -1.0E99]]), second_array_polar, axis=0),  \
                                           array([[second_array_polar[-1,0], 1.0E99]]),axis=0)
        second_array_interpolation = interp1d(second_array_polar_interp[:,1], second_array_polar_interp[:,0])
        
        # radial distance between curves
        def radial_distance(phi):
             return first_array_interpolation(phi)-second_array_interpolation(phi)

        # when radial distance is positive, second_array is inner envelope, otherwise first_array
        inner_envelope = zeros([2*(size(first_array, 0)+size(second_array, 0)) ,2])
        i = 0
        phi = first_array_polar[0,1]
        distsign = 0.0
        while phi<=first_array_polar[-1,1]:
             prev_distsign = distsign
             distsign = sign(radial_distance(phi))
             if prev_distsign*distsign<0.0:         # add cutting point when envelope curve switches (linear interpolation)
                  phi_bisect = bisect(radial_distance, prev_phi, phi)
                  inner_envelope[i,:] = array([first_array_interpolation(phi_bisect), phi_bisect])
                  i += 1
             prev_phi = phi
             if distsign>0.0:
                  inner_envelope[i,:] = array([second_array_interpolation(phi), phi])
                  try:
                       phi = second_array_polar[nonzero(second_array_polar[:,1]>phi)[0][0], 1]
                  except:
                       phi = inf        # 'break'-statement doesn't work here. Why?
             else:
                  inner_envelope[i,:] = array([first_array_interpolation(phi), phi])
                  try:
                       phi = first_array_polar[nonzero(first_array_polar[:,1]>phi)[0][0], 1]
                  except:
                       phi = inf
             i += 1

        # convert envelope back in cartesian coordinates and reverse order
        polar_envelope = self._removeTrailingZeros(inner_envelope)[::-1,:]
        cartesian_envelope = zeros([size(polar_envelope, 0)+1, 2])
        for index in range(0, size(polar_envelope, 0)):
             [x, y] = PolarCoordinatesToCartesianCoordinates(polar_envelope[index,0], polar_envelope[index,1])
             cartesian_envelope[index+(sign(self.data.get('z'))+1)/2,:] = array([x, y])*sign(self.data.get('z'))
             
        return cartesian_envelope
    

    def _removeTrailingZeros(self, inmatrix):
        """
        Remove all column-vectors with both coordinate values 0.0 from the end of inputmatrix

        INPUT parameter:
        inmatrix:       matrix containing of 2d-column-vectors (Nx2-NumPy-array)
     
        OUTPUT:
        outmatrix:      matrix containing of 2d-column-vectors, trailing zero-columns removed (Mx2-NumPy-array)
        """

        # remove zero-entries at end of matrix
        out_x = trim_zeros(inmatrix[:,0], 'b')
        out_y = trim_zeros(inmatrix[:,1], 'b')
        # correct unintentionally removed entries
        # make sure both vectors (out_x and out_y) have same length (in case one has zero at end)
        size_dif = size(out_x)-size(out_y)
        if size_dif>0:
             out_y = append(out_y, zeros(abs(size_dif)))
        elif size_dif<0:
             out_x = append(out_x, zeros(abs(size_dif)))
        # create result array (join the two rows)
        outmatrix = transpose(array([out_x, out_y]))

        return outmatrix
    
    
    def createToothShape(self, number_of_points=None):
        """
        Calculate the generated tooth shape. A manufacturing (hobbing) simulation is run.
        The shape is the envelope of the tool generated by its rolling motion during
        the hobbing process.

        INPUT parameter:
        number_of_points:  number of points for tool discretization in individual sections (list of positive integers)
                           optional - set to default if not supplied
     
        OUTPUT:
        formcoords:   profile of half a tooth of generated gear in machine coordinate
                      system (Nx2-NumPy-array)
        """

        # module imports
        from scipy.linalg import norm
        from time import time

        # some user-output
        print '\nrunning manufacturing simulation'
        t0 = time()

        # number of points value check
        if number_of_points:
             if not type(number_of_points)==type([]):
                  raise TypeError, 'list expected'
             if not min(number_of_points)>1:     # having just two points is completely useless, but it's the users responsibility to do something useful with the software
                  raise ValueError, 'all number of points must be greater than one'
             if not map(type,number_of_points)==map(type,[1.0]*len(number_of_points)):
                  raise TypeError, 'number of points not integer'
        else:       # set to default if not supplied
             number_of_points = self._no_pnts_default
        
        # vector of tool parameters used for manufacturing simulation
        tool_params = []
        char_params = [self.tool._s_P_1, self.tool._s_P_2, self.tool._s_P_3, self.tool._s_P_4, self.tool._s_P_5, self.tool._s_P_6, self.tool._s_P_7]
        for param_index in range(0, size(char_params)-1):
             if not char_params[param_index]==None:
                  if not char_params[param_index+1]==None:
                       tool_params += linspace(char_params[param_index], char_params[param_index+1], number_of_points[param_index]).tolist()
                  else:
                       tool_params += linspace(char_params[param_index], char_params[param_index+2], number_of_points[param_index]).tolist()
        tool_params = list(set(tool_params))
        tool_params.sort()

        # pre-allocate array for storing tooth flank points
        gap = zeros([size(tool_params), 2])    # additional points to complete addendum circle
        phi = zeros([size(tool_params)])

        # get tooth profile points (one gap)
        gear_point_index = 0

        # get points that are compatible with necessary condition of existence of envelope
        for candidate_param in tool_params:  
             P_c, phi_c = self._findFlankPoint(candidate_param)
             gap[gear_point_index,:] = P_c
             phi[gear_point_index]   = phi_c
             gear_point_index +=1

        # create dictionary for gear data (for passing the known parameters to Gear-class constructor)
        geardata = {'m_n':self.tool.data.get('m'), 'd':self.data.get('d'), 'x':self.data.get('x'), 'beta':self.data.get('beta'),    \
                    'A_s':self.data.get('A_s'), 'q':self.data.get('q'), 'alpha_n':self.tool.data.get('alpha_P0')}

        # transform to standard representation (2D cartesian coordinates of points on the toothflank,
        # describing a half tooth) by omitting points with negative x-coordinate and rotation of half
        # pitch angle about gear's center
        # filter out points that do not comply with sufficient condition for envelope-points
        # detect singular edge-points (usually limiting the involute) and add exact point to coordinate
        # array (if edge-points are missing, spline-fitting in GearWheel-class might show unexpected results)
        half_tooth = zeros([size(gap, 0)+5, 2])                                     # allocating array
        tau = 2*pi/self.data.get('z')                                               # pitch angle (in radians)
        R = array([[cos(tau/2.0), -sin(tau/2.0)], [sin(tau/2.0), cos(tau/2.0)]])    # rotation matrix
        formcoords_index = 1                                                        # 1st point is left [0,0]
        envelope_point   = True                                                     # indicates if last point complies with sufficient condition for envelope points
        param_env_true   = None                                                     # pre-set parameter passed to singular point computation method for which the corresponding flank point evaluates to true in the sufficient envelope condition check
        for gap_index in range(0, size(gap, 0)):
             if gap[gap_index, 0]>=0.0:                         # consider only flank points of half tooth
                 if self._checkEnvelope(gap[gap_index]):        # check if point is really on envelope of family of curves (sufficient condition)
                     if mod(gap_index, round(size(gap, 0)/40.0))==0:    # more user-output (kind of a progress-bar)
                          print '.',
                     pnt_transformed = dot(R, gap[gap_index,:])         # rotate coordinates --> symmetry-line of tooth is vertical y-axis
                     if pnt_transformed[0]<=0.0:                        # exclude points on other side of symmetry-line (tooth shape is represented by half-tooth only)
                          half_tooth[formcoords_index,:] = pnt_transformed  # add point to array
                          formcoords_index += 1                             # increment array index
                          last_param_set    = [phi[gap_index], tool_params[gap_index]]  # after loop terminated, the list contains the parameter set of the last point that complied to sufficient envelope condition
                     if not envelope_point and gap_index>0:             # undercutting occurs (singular point) if two consecutive points that comply to necessesary condition for envelope points behave different regarding sufficient condition
                          param_env_false = tool_params[gap_index-1]    # store tool curve parameter of the first of two consecutive points where one does not and the next does comply to the sufficient envelope condition
                          param_env_true  = tool_params[gap_index]      # store tool curve parameter of the second of two consecutive points where one does not and the next does comply to the sufficient envelope condition
                     envelope_point = True                              # the current point complies to necessary and sufficient envelope condition 
                 else:
                     envelope_point = False                             # the current point complies to necessary but not to sufficient envelope condition 
                 if not param_env_true==None:       # if parameter to be passed to singular point computation method is set then call it
                     singular_point  = self._findLimitingPoint(param_env_false, param_env_true)   # get singular point in gear coordinate system
                     pnt_transformed = dot(R, singular_point)                                     # rotate coordinates --> symmetry-line of tooth is vertical y-axis
                     if pnt_transformed[0]<=0.0 and not (abs(half_tooth-pnt_transformed)<self._zero_tolerance_default).all(1).any():        # exclude points on wrong side of symmetry-line (tooth shape is represented by half-tooth only) and make sure that only unique points are added to result
                          half_tooth[formcoords_index,:]   = half_tooth[formcoords_index-1,:]   # copy point before current position to current position
                          half_tooth[formcoords_index-1,:] = pnt_transformed                    # add point to array at position before current position (correct edge-point position)
                          formcoords_index += 1                                                 # increment array index                 
                     param_env_true     = None                                      # reset parameter passed to singular point computation method for which the corresponding flank point evaluates to true in the sufficient envelope condition check
        
        # add border point on addendum circle if necessary
        if abs(half_tooth[formcoords_index-1,0])>self._zero_tolerance_default*self.tool.data.get('m'):
            
             def objective(p, x):
                 """
                 Compute LHS of necessary enveloping condition and distance from target-value in x-direction.
                 For envelope point with target x-coordinate both objectives are zero (necessary condition)

                 INPUT parameters:
                 p:     parameter set of point to be evaluated (tool position angle and tool curve parameter)(1x2-NumPy-array)
                 x:     target x-value of point (numeric)

                 OUTPUT parameters:
                 obj:   necessary envelope condition and x-distance from target-value (1x2-NumPy-array)
                 """
                 cond = self._FlankPointObjectiveFnc(p[0], p[1])
                 # get flank point corresponding to parameter p[1]
                 P_tool = self._getToolPoint(p[1], 0)
                 # transform point on flank to machine coordinate system
                 P_gear  = self._transformToMachineCoords(p[0], P_tool, elemtype='point')
                 # rotate coordinates --> symmetry-line of tooth is vertical y-axis
                 P_trans = dot(R, P_gear)
                 # objective: necessary envelope condition and x-distance from target-value
                 obj    = array([cond, P_trans[0]-x])
                 return obj

             # find root of objective function with last parameter set that complied to sufficient envelope condition as initial value
             p_y      = fsolve(objective, last_param_set, args=(0.0,))
             # get corresponding point on tool
             P_y_tool = self._getToolPoint(p_y[1], 0)
             # transform to workpiece coordinate system and rotate in order to get final representation
             P_y_gear = dot(R, self._transformToMachineCoords(p_y[0], P_y_tool, elemtype='point'))
             # add to result array
             half_tooth[formcoords_index,:] = array([0.0, P_y_gear[1]])

        # remove zero-entries at end of envelope
        half_tooth = self._removeTrailingZeros(half_tooth)

        # subtract from blank if supplied
        if self.blank:
             half_tooth_final = self._subtractFromBlank(half_tooth, PythonOCCArrayToNumPyArray(self.blank.getBlankCoords()))
        else:
             half_tooth_final = half_tooth

        # transform formcoords from nparray (NumPy) to TColgp_Array1OfPnt2d (pythonOCC)
        formcoords = NumPyArrayToPythonOCCArray(half_tooth_final)

        # some more user-output
        print '\nmanufacturing simulation finished.'
        print 'time elapsed: ', time()-t0, ' s'
        
        return formcoords, geardata


class Blank:
    """
    Class representing a blank for gear manufacturing.

    The 2d machine coordinate system is:
    - x-axis     :   positive to the right
    - y-axis     :   symmetry line of tooth gap of manufactured gear
                     positive to the top (that is to the dedendum line of the tool's profile)
    - origin:        at center of manufactured gear
    """

    #Attributes
    blankcoords = None   # 2d-coordinates of the blank (in transversse cross-section)
                         # (TColgp_Array1OfPnt2d, pythonOCC)

    # Attributes: default settings for parameters
    _blank_points_default  = 1000   # default value for number of points used to create blank 


    def __str__(self):
        """
        Define string conversion of Blank objects

        INPUT parameter:
        -

        OUTPUT:
        string representation of class
        """
        
        outstr = 'blank coordinates:\n'

        # output blank coordinates
        if self.blankcoords:
             # upper and lower index of point-array
             upper_index = self.blankcoords.Upper()
             lower_index = self.blankcoords.Lower()
             for index in range(lower_index, upper_index+1):
                  outstr += str(self.blankcoords.Value(index).X()).ljust(20)+'\t'+  \
                            str(self.blankcoords.Value(index).Y()).ljust(20)+'\n'

        return outstr


    
    def setBlankCoords(self, blankcoords):
        """
        Set blank coordinates

        INPUT parameters:
        blankcoords : list of 2d-coordinate points (TColgp_Array1OfPnt2d, pythonOCC)
        """
        
        # blank coordinates: type and value check (at least two points for defining a
        # tooth form (straight flanks) and two coordinates per point)
        if not isinstance(blankcoords, TColgp_Array1OfPnt2d):
             raise TypeError, 'instance of TColgp_Array1OfPnt2d expected'
        if formcoords.Length()<2:
             raise TypeError, 'too few points for blank'

        self.blankcoords = blankcoords


    def getBlankCoords(self):
        """
        Get blank coordinates

        INPUT parameters:
        -

        OUTPUT:
        blank coordinates (TColgp_Array1OfPnt2d, pythonOCC)
        """

        return self.blankcoords

     
    def __init__(self, blankcoords=None):
        """
        Initialization of Blank-object

        INPUT parameters:
        blankcoords   : machine settings (dictionary)
        """

        if blankcoords:
             self.setBlankCoords(blankcoords)


    def setCircular(self, diameter, min_angle=0.0, max_angle=2*pi, no_of_points=_blank_points_default):
        """
        create an circular blank as used with tools that don't cut the tip circle.

        INPUT parameters:
        diameter:       diameter of circular blank (numeric)
        min_angle:      start angle of blank segment (numeric)[radians]
                        optional - set to 0.0 if not supplied
        max_angle:      end angle of blank segment (numeric)[radians]
                        optional - set to 2*pi if not supplied
        no_of_points:   number of points used for blank representation (numeric, integer)
                        optional - set to default if not supplied

        OUTPUT:
        -
        """

        # create arrays of points holding coordinates
        self.blankcoords = TColgp_Array1OfPnt2d(1, no_of_points)
        # set entries
        for point_index in range(1, no_of_points+1):
             angle = -(point_index-1)/(no_of_points-1)*(max_angle-min_angle)+max_angle
             self.blankcoords.SetValue(point_index, gp_Pnt2d(diameter/2*cos(angle), diameter/2*sin(angle)))

             


""" Planned extensions """
# class Hob:
# class BevelGearWheel(GearWheel):
# class BevelGearPair(GearPair):
# class WormGearWheel(GearWheel):
# class WormGearPair(GearPair):
# class GearStage:
# class CylindricalGearStage(GearStage):
# class PlanetaryGearStage(GearStage):
# class BevelGearStage(GearStage):
# class WormGearStage(GearStage):
# class GearBox:
