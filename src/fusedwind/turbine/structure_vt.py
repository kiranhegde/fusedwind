import numpy as np
from openmdao.main.api import VariableTree
from openmdao.lib.datatypes.api import Int, Float, Array, List, Str, Enum, Bool, VarTree, Slot, Dict

from fusedwind.interface import base, implement_base
from fusedwind.turbine.geometry_vt import AirfoilShape


# too aeroelastic code specific?!
@base
class ConcentratedMass(VariableTree):

    s = Float(desc='Non-dimens the inertia is attached')
    offset = Array(desc='x, y, z offset relative to node')
    moment_of_inertia = Array(desc='Ixx, Iyy, Izz moments of inertia')
    mass = Float(desc='Concentrated mass', units='kg')


@base
class MassProperties(VariableTree):
    """mass and mass moments of inertia of a component"""

    mass = Float(units='kg', desc='mass of object')
    Ixx = Float(units='kg*m**2', desc='mass moment of inertia about x-axis') # todo: arrary or not?
    Iyy = Float(units='kg*m**2', desc='mass moment of inertia about y-axis')
    Izz = Float(units='kg*m**2', desc='mass moment of inertia about z-axis')
    Ixy = Float(units='kg*m**2', desc='mass x-y product of inertia')
    Ixz = Float(units='kg*m**2', desc='mass x-z product of inertia')
    Iyz = Float(units='kg*m**2', desc='mass y-z product of inertia')


@base
class BeamStructureVT(VariableTree):

    s = Array(desc='Running curve length of beam', units='m')
    dm = Array(desc='Mass per unit length', units='kg/m')
    x_cg = Array(desc='x-distance from blade axis to center of mass', units='m')
    y_cg = Array(desc='y-distance from blade axis to center of mass', units='m')
    ri_x = Array(desc='radius of gyration relative to elastic center.', units='m')
    ri_y = Array(desc='radius of gyration relative to elastic center', units='m')
    x_sh = Array(desc='x-distance from blade axis to shear center', units='m')
    y_sh = Array(desc='y-distance from blade axis to shear center', units='m')
    E = Array(desc='modulus of elasticity', units='N/m**2')
    G = Array(desc='shear modulus of elasticity', units='N/m**2')
    I_x = Array(desc='area moment of inertia with respect to principal bending xe axis', units='m**4')
    I_y = Array(desc='area moment of inertia with respect to principal bending ye axis', units='m**4')
    J = Array(desc='torsional stiffness constant with respect to ze axis at the shear center', units='m**4/rad')
    k_x = Array(desc='shear factor for force in principal bending xe direction', units=None)
    k_y = Array(desc='shear factor for force in principal bending ye direction', units=None)
    A = Array(desc='cross sectional area', units='m**2')
    pitch = Array(desc='structural pitch relative to main axis.', units='deg')
    x_e = Array(desc='x-distance from main axis to center of elasticity', units='m')
    y_e = Array(desc='y-distance from main axis to center of elasticity', units='m')


@base
class MaterialProps(VariableTree):
    """
    Material properties 
    """

    materialname = Str(desc='Material name')

    E1 = Float(units='GPa', desc='Youngs modulus in the radial direction')
    E2 = Float(units='GPa', desc='Youngs modulus in the tangential direction')
    E3 = Float(units='GPa', desc='Youngs modulus in the axial direction')
    nu12 = Float(desc='Poissons ratio in the radial-tangential direction')
    nu13 = Float(desc='Poissons ratio in the radial-axial direction')
    nu23 = Float(desc='Poissons ratio in the tangential-axial direction')
    G12 = Float(units='GPa', desc='Shear modulus in radial-tangential direction')
    G13 = Float(units='GPa', desc='Shear modulus in radial-axial direction')
    G23 = Float(units='GPa', desc='Shear modulus in tangential-axial direction')
    rho = Float(units='kg/m**3', desc='Mass density')

    failure_criterium = Enum('maximum_strain', ('maximum_strain', 'maximum_stress', 'tsai_wu'))
    s11_t = Float(1e6, desc='tensile strength in the 1 direction')
    s22_t = Float(1e6, desc='tensile strength in the 2 direction')
    s33_t = Float(1e6, desc='tensile strength in the 3 direction')
    s11_c = Float(1e6, desc='compressive strength in the 1 direction')
    s22_c = Float(1e6, desc='compressive strength in the 2 direction')
    s33_c = Float(1e6, desc='compressive strength in the 3 direction')
    t12 = Float(1e6, desc='shear strength in the 12 plane')
    t13 = Float(1e6, desc='shear strength in the 13 plane')
    t23 = Float(1e6, desc='shear strength in the 23 plane')
    e11_c = Float(1e6, desc='maximum tensile strain in the 1 direction')
    e22_c = Float(1e6, desc='maximum tensile strain in the 2 direction')
    e33_c = Float(1e6, desc='maximum tensile strain in the 3 direction')
    e11_t = Float(1e6, desc='maximum compressive strain in the 1 direction')
    e22_t = Float(1e6, desc='maximum compressive strain in the 2 direction')
    e33_t = Float(1e6, desc='maximum compressive strain in the 3 direction')
    g12 = Float(1e6, desc='maximum shear strain in the 12 plane ')
    g13 = Float(1e6, desc='maximum shear strain in the 13 plane')
    g23 = Float(1e6, desc='maximum shear strain in the 23 plane')

    # Reduction factors for the material safety factor
    # computed as gMa = gM0 * sum(Cia) 
    # gMa defaults to 1 so partials can be used directly in above material strengths
    gM0 = Float(0.25, desc='Material safety factor')
    C1a = Float(1., desc='influence of ageing')
    C2a = Float(1., desc='influence of temperature')
    C3a = Float(1., desc='influence of manufacturing technique')
    C4a = Float(1., desc='influence of curing technique')


@base
class Layer(VariableTree):
    """
    A layer represents a stack of multidirectional plies and is assigned
    an apparent set of material properties based on the 
    properties of the constituent materials, pre-computed using
    simple micromechanics equations and classical lamination theory.

    Layers are stacked in a region
    """

    thickness = Float()
    angle = Float()
    materialname = Str()
    plyname = Str()


@base
class Region(VariableTree):
    """
    A region covers a fraction of the curve length along the surface of
    a cross section.

    A region consists of a number of layers.
    """

    s0 = Float(0., desc='Chordwise curve fraction starting point')
    s1 = Float(1., desc='Chordwise curve fraction end point')
    thickness = Float(desc='Total thickness of region')
    width = Float(desc='Width of the region')
    layers = List()

    def add_layer(self, name):

        self.add(name, VarTree(Layer()))
        self.layers.append(name)
        return getattr(self, name)

    def compute_thickness(self):

        for lname in self.layers:
            self.thickness += np.maximum(0., getattr(self, lname).thickness)


@base
class CrossSectionStructureVT(VariableTree): # more appropriate name: CrossSectionLayupVT
    """
    Container for a cross-sectional definition of the
    internal structure of a blade.
    """
    s = Float()
    regions = List(desc='List of names of regions in the cross section')
    webs = List(desc='List of names of regions in the cross section')
    materials = Dict(desc='Dictionary of MaterialProps vartrees')
    airfoil = VarTree(AirfoilShape(), desc='Cross sectional shape')
    DPs = List(desc='Region division points (nregion + 1)')

    def add_region(self, name):

        self.add(name, VarTree(Region())) 
        self.regions.append(name)
        return getattr(self, name)

    def add_web(self, name):

        self.add(name, VarTree(Region())) 
        self.webs.append(name)
        return getattr(self, name)

    def add_material(self, name, material):

        if name in self.materials.keys():
            return
        else:
            self.materials[name] = material

        return self.materials[name]


@base
class Layer3D(VariableTree):
    """
    Same as Layer, except for being a function of span

    A layer thickness can go to zero if material disappears at
    a certain spanwise location.
    """

    thickness = Array(units='m', desc='layer thickness')
    angle = Array(units='deg', desc='layup angle')


@base
class Region3D(VariableTree):
    """
    Same as region except for being a function of span
    """
    x = Array(desc='spanwise discretization of section')
    thickness = Array(desc='Total thickness distribution as function of span')
    width = Array(desc='Width as function of span (for convenience)')
    layers = List(desc='List of names of layers')

    def add_layer(self, name):

        dubl = self.list_containers().count(name)
        if dubl > 0:
            lname = '%s%02d' % (name, dubl)
        else:
            lname = name

        self.add(lname, VarTree(Layer3D()))
        self.layers.append(lname)
        return getattr(self, lname)


@base
class BladeStructureVT3D(VariableTree): # more appropriate name:  BladeLayupVT3D
    """
    Variable tree for the structural definition of a blade.
    """

    x = Array(desc='spanwise discretization of blade')
    regions = List(desc='List of names of regions')
    webs = List(desc='List of names of webs')
    iwebs = List(desc='List of DP indices connecting webs to the surface')
    DPs = List(desc='Names of division point curves')
    materials = Dict()

    def configure_regions(self, nr, names=[]):

        for i in range(nr + 1):
            self.add('DP%02d' % i, Array(dtype=float,
                                         low=-1., high=1.,
                                         desc='Region division point curves %i' % i))
            self.DPs.append('DP%02d' % i)

        for i in range(nr):
            try:
                name = names[i]
            except:
                name = 'region%02d' % i
            self.add_region(name)

    def configure_webs(self, nw, iwebs, names=[]):

        self.iwebs = iwebs

        for i in range(nw):
            try:
                name = names[i]
            except:
                name = 'web%02d' % i
            self.add_web(name)

    def add_region(self, name):

        self.add(name, VarTree(Region3D()))
        region = getattr(self, name)
        self.regions.append(name)
        # self.nr = len(self.regions)
        return region

    def add_web(self, name):

        self.add(name, VarTree(Region3D())) 
        self.webs.append(name)
        return getattr(self, name)

    def add_material(self, name):
        """
        add a material to the blade
        
        parameters
        -----------
        name: string
            name of material to add

        returns
        --------
        MaterialProps: object
            VariableTree with material properties to set by user
        """

        mat = MaterialProps()
        mat.materialname = name
        self.add(name, VarTree(mat))
        self.materials[name] = getattr(self, name)
        return getattr(self, name)
    
    
@base
class ResultVectorArray(VariableTree):
    '''
    Container for element results
    '''
    id = Str(desc = 'Result identifier')
    comp_11 = Array(desc = 'Result component 11')
    comp_22 = Array(desc = 'Result component 22')
    comp_33 = Array(desc = 'Result component 33')
    comp_12 = Array(desc = 'Result component 12')
    comp_13 = Array(desc = 'Result component 13')
    comp_23 = Array(desc = 'Result component 23')
    
    def _fromarray(self, d):

        self.comp_11 = d[:,0]
        self.comp_22 = d[:,1]
        self.comp_33 = d[:,2]
        self.comp_12 = d[:,3]
        self.comp_13 = d[:,4]
        self.comp_23 = d[:,5]
        
    def _toarray(self):

        return np.array([self.comp_11, self.comp_22, self.comp_33, self.comp_12, self.comp_13,
                         self.comp_23]).T

@base
class KeyPointsVT(VariableTree):
    kp_name = Str('cs_xx',desc='Pointer for  key point array')
    kp_coords = Array([[0.,0.],[1.,2.]],desc='Array containing arrays of key point coordinates')
    
@base
class KeyPoint(VariableTree):
    x = Float()
    y = Float()
    

@base
class Area(VariableTree):
    '''
    Area of a cross section.
    '''
    #no_KPs,ET,MAT,ASYS,REAL,CUT_flag
    KPs = List(desc = 'List of keypoint names')
    kp_ni = Int(desc='Number of area keypoints')
    mat_name = Str(desc='Material object name')
    csys_name = Str(desc='Coordinate system object name')
    
    

    
@base
class Csys(VariableTree):
    '''
    Coordinate system used to orient the material within an area.
    '''
    #No,X,Y,Z,THETAX,THETAY,THETAZ (deg)
    # relative to origin coordinate system
    # x = pointing from left to right
    # y = pointing from bottom to top
    # z = pointing out of the screen (towards the user)
    csys_name = Str(desc='Coordinate system name')
    
    origin_x = Float(unit='m', desc='Origin of the coordinate system')
    origin_y = Float(unit='m', desc='Origin of the coordinate system')
    origin_z = Float(unit='m', desc='Origin of the coordinate system')
    theta_x = Float(unit='deg', desc='Rotation of the coordinate system')
    theta_y = Float(unit='deg', desc='Rotation of the coordinate system')
    theta_z = Float(unit='deg', desc='Rotation of the coordinate system')

@base
class CrossSectionAreasVT(VariableTree):
    """
    Container for cross-sectional definition of areas to be meshed.
    """
    
    KPs = List(desc = 'List of keypoints')
    csysts = Dict(desc='List of Csys names')
    materials = Dict(desc='Dictionary of MaterialProps vartrees')
    areas = List(desc='Dictionary of Area vartrees belonging to a cross section')
    
    #ielsets = Dict(desc='Dictionary of IndividualElset vartrees')
    
    def add_kp(self, name):
        self.add(name, VarTree(KeyPoint()))
        self.KPs.append(name)
        return getattr(self, name)
   
    def add_csys(self, name):
        self.add(name, VarTree(Csys()))
        csys = Csys()
        csys.csys_name = name
        self.add(name, VarTree(csys))
        self.csysts[name] = getattr(self, name)
        return getattr(self, name)
    
    def add_material(self, name):
        mat = MaterialProps()
        mat.materialname = name
        self.add(name, VarTree(mat))
        self.materials[name] = getattr(self, name)
        return getattr(self, name)
    
    def add_area(self, name):
        self.add(name, VarTree(Area()))
        self.areas.append(name)
        return getattr(self, name)

@base
class KeyPoint3D(VariableTree):
    x = Array()
    y = Array()

@base
class CrossSectionAreasVT3D(VariableTree):
    """
    Container for cross-sectional definition of areas as function of span.
    """
    
    z = Array(desc='spanwise discretization of blade')
    KPs = List(desc = 'List of 3D keypoints')
    csysts = Dict(desc='List of Csys names')
    materials = Dict(desc='Dictionary of MaterialProps vartrees')
    areas = List(desc='Dictionary of Area vartrees belonging to a cross section')
    
    def add_kp(self, name):
        self.add(name, VarTree(KeyPoint3D()))
        self.KPs.append(name)
        return getattr(self, name)
   
    def add_csys(self, name):
        self.add(name, VarTree(Csys()))
        csys = Csys()
        csys.csys_name = name
        self.add(name, VarTree(csys))
        self.csysts[name] = getattr(self, name)
        return getattr(self, name)
    
    def add_material(self, name):
        mat = MaterialProps()
        mat.materialname = name
        self.add(name, VarTree(mat))
        self.materials[name] = getattr(self, name)
        return getattr(self, name)
    
    def add_area(self, name):
        self.add(name, VarTree(Area()))
        self.areas.append(name)
        return getattr(self, name)
    
    
@base
class MeshProps(VariableTree):
    '''
    Mesh properties as input for a mesher
    '''
    etype = Str(desc='Mesh element type')
    

@base
class Elset(VariableTree):
    '''
    Element set container
    '''
    el_numbers = Array(desc = 'Element numbers belonging to the set' ) 
    
@base
class CrossSectionMeshVT(VariableTree):
    '''
    Container for a 2D cross sectional mesh.
    '''
    nl_2d = Array(desc='Nodal points (node nr, x, y)')
    defs = '(Element number, node 1, n2, n3, ..., n8)'
    el_2d = Array(desc='Elements %s' % defs)
    defs = '(Element nr, material nr, fiber angle, fiberplane angle)'
    emat  = Array(desc='Material per element %s' % defs)
    matprops = Array(desc='Material properties (see docs)')
    elsets = Dict(desc='List of element set names')
    
@base
class CrossSectionMeshVT3D(VariableTree):
    '''
    Container for a 3D mesh.
    '''
    nl_3d = Array(desc='Nodal points (node nr, x, y)')
    defs = '(Element number, node 1, n2, n3, ..., n8)'
    el_3d = Array(desc='Elements %s' % defs)
    defs = '(Element nr, material nr, fiber angle, fiberplane angle)'
    emat  = Array(desc='Material per element %s' % defs)
    matprops = Array(desc='Material properties (see docs)')
    elsets = List(desc='List of element set names')
    
@base
class CrossSectionElementStressRecoveryVT(VariableTree):
    '''
    Container for cross sectional element results returned from any FE code.
    '''
    el_stresses = Dict(desc='element stresses') # List of ResultVectors
    el_strains = Dict(desc='element strains') # List of ResultVectors
    
#------------------------------------------------------------------------ #@base
#---------------------------------------- #class CrossSectionMesh(VariableTree):
    #----------------------------------------------------------------------- '''
    #------------------------------------------ Container for Cross section mesh
    #----------------------------------------------------------------------- '''
    #------------------------ nl_2d = Array(desc='Nodal points (node nr, x, y)')
    #------------------------ defs = '(Element number, node 1, n2, n3, ..., n8)'
    #---------------------------------- el_2d = Array(desc='Elements %s' % defs)
    #----------------------- elsets = Dict(desc='element sets') # Lsit of Elsets
#------------------------------------------------------------------------------ 
    #------- el_stresses = Dict(desc='element stresses') # List of ResultVectors
#-------- #    el_strains = Dict(desc='element strains') # List of ResultVectors