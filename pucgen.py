#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PUCGEN - GUI for Periodic Unit Cell generator

"""

import sys
import os
import numpy as nm
from inspect import getargspec
from gen_mesh_utils import gmsh_call

FREECADPATH = '/usr/lib/freecad/lib/'
sys.path.append(FREECADPATH)
from FreeCAD import Base, newDocument
import Sketcher
import Part


class PUC(object):
    """Periodic Unit Cell object."""
    def __init__(self, cell_mat_id=1, base_cell=None):
        """Init PUC"""
        if cell_mat_id is None:
            self.components = []
        else:
            self.components = [BaseCell(mat_id=cell_mat_id)\
                if base_cell is None else base_cell] 

    def add(self, obj):
        """Add a new component to the list (inclusion, channel, layer)."""
        self.components.append(obj)

    @staticmethod
    def save_puc(filename, comps):
        """Save PUC into the text file."""
        with open(filename, 'wt') as f:
            for cls, pars, act in comps:
                aflag = '' if act else '#'
                args, _, _, _ = getargspec(cls.__init__)
                args = args[1:]
                targs = []
                for k in args:
                    val = pars.get(k)
                    if isinstance(val, nm.ndarray):
                        targs.append(k + '=' + str(tuple(val)))
                    else:
                        targs.append(k + '=' + str(val))
                f.write('%s%s;%s\n' % (aflag, cls.__name__, ';'.join(targs)))
            
    def save(self, filename):
        self.save_puc(filename, [(c.__class__, c, c.active)
            for c in self.components])

    @staticmethod
    def load_puc(filename):
        """Load PUC from a given text file."""
        from ast import literal_eval

        cls_dict = {c.__name__: c for c in pucgen_classes}
        out = []
        with open(filename, 'rt') as f:
            for line in f.readlines():
                comp = line.strip().split(';')
                pars = {}
                for arg in comp[1:]:
                    key, val = [ii.strip() for ii in arg.split('=')]
                    val = val if val.isalpha() else literal_eval(val)
                    pars[key] = val
                clsname = comp[0]
                if clsname.startswith('#'):
                    clsname = clsname[1:]
                    is_active = 0
                else:
                    is_active = 1
                out.append([cls_dict[clsname], pars, is_active])

        return out

    @staticmethod
    def from_file(filename):
        """Create the PUC object from a given text file."""
        new = PUC(cell_mat_id=None)
        for cls, pars, act in new.load_puc(filename):
            comp = cls(**pars)
            if act == 0:
                comp.deactivate()
            new.add(comp)

        return new

    def __call__(self, filename_vtk, cell_size=None, eps=1.0,
                 save_FCStd=False, centered=False):
        """Generate the finite element mesh.

        Parameters
        ----------
        filename_vtk: str
            The output VTK file name.
        cell_size: array
            The size of PUC: [sx, sy, sz].
        eps: float
            The scaling parameter.
        save_FCStd: bool
            If True, save the FreeCAD data file.
        centered: bool
            If True, the PUC is centered to the origin.
        """

        if cell_size is None:
            cell_size = self.components[0].get('size')

        element_size = self.components[0].get('el_size')

        cell_size = nm.asarray(cell_size, dtype=nm.float64)
        doc = newDocument()

        attrs = []
        mat_ids = []
        volumes = {}
        for comp in self.components:
            if not comp.active:
                continue

            mat_id = comp.params['mat_id']
            obj, attr = comp(doc, cell_size)
            if obj is None:
                continue

            if mat_id in volumes:
                volumes[mat_id].append(obj)
            else:
                volumes[mat_id] = [obj]
                mat_ids.append(mat_id)

            if attr is not None:
                attrs += attr

        bcell = volumes[mat_ids[0]]
        if len(bcell) == 2:
            bcell_ctool = bcell[1]
        elif len(bcell) >= 3:
            bcell_ctool = doc.addObject('Part::MultiFuse', 'uni_1')
            bcell_ctool.Shapes = bcell[1:]
            pass
        else:
            bcell_ctool = None

        objs = [bcell[0]]
        for mat_id in mat_ids[1:]:
            if len(volumes[mat_id]) >= 2:
                uni = doc.addObject('Part::MultiFuse', 'uni_%d' % mat_id)
                uni.Shapes = volumes[mat_id]
            else:
                uni = volumes[mat_id][0]

            if bcell_ctool is None:
                objs.append(uni)
            else:
                cut = doc.addObject('Part::Cut', 'cut_%d' % mat_id)
                cut.Base = uni
                cut.Tool = bcell_ctool
                objs.append(cut) 

        doc.recompute()

        if len(objs) > 1:
            if len(objs) > 2:
                uni = doc.addObject('Part::MultiFuse', 'uni')
                uni.Shapes = objs[1:]
            else:
                uni = objs[1]
            
            cut = doc.addObject('Part::Cut', 'cut')
            cut.Base = objs[0]
            cut.Tool = uni
            objs[0] = cut

        for io, obj in enumerate(objs):
            grp = doc.addObject('App::DocumentObjectGroup', 'grp_%d' % (io + 1))
            grp.addObject(obj)

        doc.recompute()

        physical_volumes = {}
        ivol = 1
        for io, (obj, mat_id) in enumerate(zip(objs, mat_ids)):
            if mat_id not in physical_volumes:
                physical_volumes[mat_id] = []
            for _ in range(len(obj.Shape.Solids)):
                physical_volumes[mat_id].append(ivol)
                ivol += 1

        filename_base = os.path.splitext(filename_vtk)[0]

        Part.export(objs, filename_base + '.step')
        if save_FCStd:
            doc.saveAs(filename_base + '.FCStd')

        periodicity = [('all', 0),
                       ([1, 0, 0], cell_size[0]),
                       ([0, 1, 0], cell_size[1]),
                       ([0, 0, 1], cell_size[2])]
        # periodicity = []

        if len(attrs) == 0:
            attrs = 'const'

        if not centered:
            shift = nm.asarray(cell_size) * 0.5

        gmsh_call(filename_base, filename_vtk, attrs, element_size,
                  periodicity=periodicity, merge=True,
                  physical_volumes=physical_volumes, shift=shift, scale=eps)


class BaseComponent(object):
    """The base component of the unit cell."""
    name = None
    parameters_dict = {}

    def __init__(self, mat_id=1):
        """Init parameters of the component.
    
        Parameters
        ----------
        mat_id: int
            The component material id.
        """
        print('new BaseComponent: %s' % self.name)

        self.params = {'mat_id': mat_id}
        self.active = True

    def __call__(self, cont, size):
        """Create the cell box.
    
        Parameters
        ----------
        cont: FreeCAD document
            The document to which the object is created.
        size: array
            The size of the rectangular cell: [size_x, size_y, size_z].

        Returns
        -------
        obj: FreeCAD object
            The created object.
        attrs: list
            The attractors affecting the mesh density inside and around
            the object.
        """
        pass

    def get(self, key):
        if key in self.parameters_dict:
            return(self.params['parameters'][self.parameters_dict[key]])
        else:
            return(self.params[key])

    @staticmethod
    def create_box(cont, size):
        size = nm.asarray(size)
        box = cont.addObject("Part::Box", "box")
        box.Length, box.Width, box.Height = size
        box.Placement = Base.Placement(Base.Vector(*(-size * 0.5)),
                                       Base.Rotation(0, 0, 0))

        return box


class BaseCell(BaseComponent):
    """The base cell - matrix."""
    name = 'Base Cell'

    def __init__(self, size=(1, 1, 1), el_size=0.1, mat_id=1):
        super(BaseCell, self).__init__(mat_id=mat_id)
        self.params.update({
            'size': size,
            'el_size': el_size,
        })

    def __call__(self, cont, size):
        return self.create_box(cont, size), None


class BaseEmbeddedComponent(BaseComponent):
    """The base for the inclusion and channel classes."""

    direction_tab = {'x': nm.array([1, 0, 0]),
                     'y': nm.array([0, 1, 0]),
                     'z': nm.array([0, 0, 1])}

    def __init__(self, parameters, central_point, direction,
                 es_dmin, es_dmax, es_in, mat_id):
        """Init parameters of the channel component.
    
        Parameters
        ----------
        parameters: float or array
            The geometrical parameters of the object.
        central_point: array
            The coordinates of the object center: [x, y, z].
        direction: str or array
            The object direction. If string: direction = 'x', 'y' or 'z'.
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """    
        super(BaseEmbeddedComponent, self).__init__(mat_id=mat_id)

        if isinstance(direction, list) or isinstance(direction, tuple):
            direction = nm.asarray(direction) / nm.linalg.norm(direction)

        self.params.update({
            'parameters': parameters,
            'direction': direction,
            'central_point': nm.asarray(central_point, dtype=nm.float64),
            'es_dmin': es_dmin,
            'es_dmax': es_dmax,
            'es_in': es_in,
        })

    def deactivate(self):
        self.active = False

    def activate(self):
        self.active = True


class EllipsoidalInclusion(BaseEmbeddedComponent):
    """The ellipsoidal inclusion."""
    name = 'Ellipsoidal Inclusion'
    parameters_dict = {'radius': 0}
    
    def __init__(self, radius=(0.1, 0.1, 0.1), central_point=(0, 0, 0),
                 direction=(1, 0, 0), es_dmin=1.1, es_dmax=1.3, es_in=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        radius: float or array
            The radii of the ellipsoid: r or [r1, r2, r3].
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: array
            The directional vector.
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """
        super(EllipsoidalInclusion, self).__init__(parameters=[radius],
                                                   central_point=central_point,
                                                   direction=direction,
                                                   es_dmin=es_dmin,
                                                   es_dmax=es_dmax,
                                                   es_in=es_in,
                                                   mat_id=mat_id)

    def __call__(self, cont, size):
        r = self.get('radius')
        d = self.get('direction')
        p = self.get('central_point')
        label = '%s_%d' % (self.name, self.get('mat_id'))

        attrs = []

        if nm.any(r == 0):
            return None, None

        ell = cont.addObject('Part::Ellipsoid', label)
        if isinstance(r, (int, float)):
            ell.Radius1 = ell.Radius2 = ell.Radius3 = r
            rmax = r
        else:
            ell.Radius1 = r[2]
            ell.Radius2 = r[0]
            ell.Radius3 = r[1]
            rmax = nm.max(r)

        ell.Angle1 = -90
        ell.Angle2 = 90
        ell.Angle3 = 360

        ell.Placement = Base.Placement(Base.Vector(*p), Base.Rotation())

        s = nm.cross([1, 0, 0], d)
        if nm.linalg.norm(s) > 1e-12:
            phi = nm.rad2deg(nm.arccos(nm.dot([1, 0, 0], d)))
            rot = Base.Rotation(Base.Vector(*s), phi)
            ell.Placement = Base.Placement(Base.Vector(), rot,\
                Base.Vector(*p)).multiply(ell.Placement)

        attrs = [('sphere', p, rmax * self.get('es_dmin'),
                 rmax * self.get('es_dmax'), self.get('es_in'))]

        return ell, attrs


# class SphericalInclusion(EllipsoidalInclusion):
class SphericalInclusion(BaseEmbeddedComponent):
    """The spherical inclusion."""
    name = 'Spherical Inclusion'
    parameters_dict = {'radius': 0}

    def __init__(self, radius=0.1, central_point=(0, 0, 0),
                 es_dmin=1.1, es_dmax=1.3, es_in=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        radius: float
            The radius of the sphere.
        central_point: array
            The coordinates of the center: [x, y, z].
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """
        super(SphericalInclusion, self).__init__(parameters=[radius],
                                                 central_point=central_point,
                                                #  direction=None,
                                                 direction=(1, 0.7, 0.5), # magic hack
                                                 es_dmin=es_dmin,
                                                 es_dmax=es_dmax,
                                                 es_in=es_in,
                                                 mat_id=mat_id)

    def __call__(self, cont, size):
        r = self.get('radius')
        p = self.get('central_point')
        label = '%s_%d' % (self.name, self.get('mat_id'))

        attrs = []

        if r == 0:
            return None, None

        ell = cont.addObject('Part::Sphere', label)
        ell.Radius = r
        ell.Placement = Base.Placement(Base.Vector(*p), Base.Rotation())

        d = self.get('direction')
        s = nm.cross([1, 0, 0], d)
        if nm.linalg.norm(s) > 1e-12:
            phi = nm.rad2deg(nm.arccos(nm.dot([1, 0, 0], d)))
            rot = Base.Rotation(Base.Vector(*s), phi)
            ell.Placement = Base.Placement(Base.Vector(), rot,\
                Base.Vector(*p)).multiply(ell.Placement)

        attrs = [('sphere', p, r * self.get('es_dmin'),
                 r * self.get('es_dmax'), self.get('es_in'))]

        return ell, attrs
# [(5, 12), (8, 10), (7, 11), (0, 13), (6, 9), (0, 13), (6, 9)]

class CylindricalInclusion(BaseEmbeddedComponent):
    """The cylindrical inclusion."""
    name = 'Cylindrical Inclusion'
    parameters_dict = {'radius': 0, 'length': 1}

    def __init__(self, radius=0.1, length=0.5, central_point=(0, 0, 0),
                 direction=(1, 0, 0), es_dmin=1.1, es_dmax=1.3, es_in=0.5,
                 mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        radius: float
            The cylinder radius.
        length: float
            The cylinder length.
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: array
            The directional vector.
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """
        super(CylindricalInclusion, self).__init__(parameters=[radius, length],
                                                   central_point=central_point,
                                                   direction=direction,
                                                   es_dmin=es_dmin,
                                                   es_dmax=es_dmax,
                                                   es_in=es_in,
                                                   mat_id=mat_id)

    def __call__(self, cont, size):
        r = self.get('radius')
        h = self.get('length')
        d = self.get('direction')
        p = self.get('central_point')
        attr_ext = 1.
        label = '%s_%d' % (self.name, self.get('mat_id'))

        attrs = []

        if nm.any(self.get('radius') == 0):
            return None, None

        if isinstance(d, str) or isinstance(d, unicode):
            idir = {'x': 0, 'y': 1, 'z': 2}[d]
            d = self.direction_tab[d]
            if h is None:
                h = size[idir]
                p = p.copy()
                p[idir] = 0
                attr_ext = 1.2

        cyl = cont.addObject('Part::Cylinder', label)
        cyl.Radius = r
        cyl.Height = h

        p0 = p - 0.5 * h * nm.array([0, 0, 1])
        cyl.Placement = Base.Placement(Base.Vector(*p0), Base.Rotation())

        s = nm.cross([0, 0, 1], d)
        if nm.linalg.norm(s) > 1e-12:
            phi = nm.rad2deg(nm.arccos(nm.dot([0, 0, 1], d)))
            rot = Base.Rotation(Base.Vector(*s), phi)
            cyl.Placement = Base.Placement(Base.Vector(0, 0, 0), rot,\
                Base.Vector(*p)).multiply(cyl.Placement)

        p1 = p - 0.5 * d * h * attr_ext
        p2 = p + 0.5 * d * h * attr_ext

        attrs = [('line', [p1, p2], r * self.get('es_dmin'),
                  r * self.get('es_dmax'), self.get('es_in'))]

        return cyl, attrs


class CylindricalChannel(CylindricalInclusion):
    """The cylindrical channel."""
    name = 'Cylindrical Channel'

    def __init__(self, radius=0.1, central_point=(0, 0, 0), direction='x',
                 es_dmin=1.1, es_dmax=1.3, es_in=0.5, mat_id=2):
        """Init parameters of the channel component.
    
        Parameters
        ----------
        radius: flaot
            The radius of the cylinder.
        central_point: array
            The coordinates of the cylinder center: [x, y, z].
        direction: str
            The cylinder orientation: 'x', 'y' or 'z'.
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """
        super(CylindricalChannel, self).__init__(radius=radius, length=None,
                                                 central_point=central_point,
                                                 direction=direction,
                                                 es_dmin=es_dmin,
                                                 es_dmax=es_dmax,
                                                 es_in=es_in,
                                                 mat_id=mat_id)


class BoxInclusion(BaseEmbeddedComponent):
    """The box inclusion."""
    name = 'Box Inclusion'
    parameters_dict = {'size': 0}
    
    def __init__(self, size=(0.3, 0.2, 0.1), central_point=(0, 0, 0),
                 direction=None, es_dmin=1.1, es_dmax=1.3, es_in=0.5,
                 mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        size: array
            The size of the box: [sx, sy, sz].
        central_point: array
            The coordinates of the center: [x, y, z].
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """
        super(BoxInclusion, self).__init__(parameters=[size],
                                           central_point=central_point,
                                           direction=direction,
                                           es_dmin=es_dmin,
                                           es_dmax=es_dmax,
                                           es_in=es_in,
                                           mat_id=mat_id)

    def __call__(self, cont, size):
        d = self.get('direction')
        if d is None:
            s = nm.array(self.get('size'))
            p = self.get('central_point')
        else:
            s = nm.array(size)
            idir = {'x': 0, 'y': 1, 'z': 2}[d]
            s[idir] = self.get('thickness')
            p = nm.zeros(3, dtype=nm.float64)
            p[idir] = self.get('central_point')[idir]

        label = '%s_%d' % (self.name, self.get('mat_id'))

        attrs = []

        if nm.any(s == 0):
            return None, None

        box = cont.addObject('Part::Box', label)
        box.Length, box.Width, box.Height = s
        box.Placement = Base.Placement(Base.Vector(*(p - 0.5 * s)),
                                       Base.Rotation())

        aux = s * 0.5 * self.get('es_dmin')
        attrs = [('box', p - aux, p + aux,
                  (self.get('es_dmax') - self.get('es_dmin')) * nm.mean(s),
                  self.get('es_in'))]

        return box, attrs

class SandwichLayer(BoxInclusion):
    """The sandwich layer."""
    name = 'Sandwich Layer'
    parameters_dict = {'thickness': 0}
    
    def __init__(self, thickness=0.1, central_point=(0, 0, 0),
                 direction='x', es_dmin=1.1, es_dmax=1.3, es_in=0.5,
                 mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        thickness: array
            The thicknesss of the layer.
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: str
            The orientation of the layer normal vector: 'x', 'y' or 'z'.
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """
        super(SandwichLayer, self).__init__(size=thickness,
                                            central_point=central_point,
                                            direction=direction,
                                            es_dmin=es_dmin,
                                            es_dmax=es_dmax,
                                            es_in=es_in,
                                            mat_id=mat_id)


class SweepedChannel(BaseEmbeddedComponent):
    """The sweeped channel."""
    name = 'Sweeped Channel'

    @staticmethod
    def get_sketch(path, label, cont,
                   is_closed=False, path_type='polyline', plane='xy'):
        sketch = cont.addObject('Sketcher::SketchObject', label)
        sketch.Placement = Base.Placement(Base.Vector(0, 0, 0),
                                          Base.Rotation(0.5, 0.5, 0.5, .5))

        coinc = []
        npts = len(path)
        for ii in range(npts):
            p1 = tuple(path[ii])
            jj = ii + 1
            if (ii + 1) >= npts:
                if is_closed:
                    j = 0
                else:
                    continue

            p2 = tuple(shape[jj])
            sketch.addGeometry(Part.Line(Vec(*p1), Vec(*p2)), False)
            coinc.append((ii, jj))

        for ii in coinc:
            sketch.addConstraint(Sketcher.Constraint('Coincident', ii[0], 2, ii[1], 1)) 

        return sketch

    def __init__(self, profile=[[0.1, 0], [0, 0.1], [-0.1, -0.1]],
                 path=[[0, 0], [0.3, 0.1], [1, 0]], central_point=(0, 0, 0),
                 direction='x', es_dmin=1.1, es_dmax=1.3, es_in=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        !!!

        central_point: array
            The coordinates of the center: [x, y, z].
        direction: array
            The directional vector.
        es_dmin: float
            The distance to which the "inner" element size factor is applied.
        es_dmax: float
            The distance from which the "outer" element size factor is applied.
        es_in: float
            The "inner" element size factor: in_el_size = es_in * el_size_base. 
        """
        super(SweepedChannel, self).__init__(parameters=(profile, path),
                                             central_point=central_point,
                                             direction=direction,
                                             es_dmin=es_dmin,
                                             es_dmax=es_dmax,
                                             es_in=es_in,
                                             mat_id=mat_id)

    def __call__(self, cont, size):
        profile = self.get('profile')
        path = self.get('path')
        d = self.get('direction')
        p = self.get('central_point')
        label = '%s_%d' % (self.name, self.get('mat_id'))

        attrs = []

        if isinstance(d, str):
            idir = {'x': 0, 'y': 1, 'z': 2}[d]
            d = self.direction_tab[d]
            h = size[idir]
            p = p.copy()
            p[idir] = 0

        sk_profile = self.get_sketch(profile[1], 'profile', cont,
                                     is_closed=True, path_type=profile[0],
                                     plane='yz')

        path_pts = nm.asarray(path[1])
        path_pts[:, 0] = (path_pts[:, 0] - 0.5) * h 

        sp_path = self.get_sketch(path_pts, 'path', cont,
                                  is_closed=True, path_type=path[0],
                                  plane='xy')

        ch = cont.addObject('Part::Sweep', label)
        ch.Sections = [sk_profile]
        ch.Spine = (sp_path,[])
        ch.Solid = True
        ch.Frenet = False

        attr_comm = (r * self.get('es_dmin'), r * self.get('es_dmax'),
                     self.get('es_in'))
        p1 = p - 0.5 * d * h * 1.2
        p2 = p - 0.5 * d * h * 1
        for ii in range(len(path_pts)):
            attrs = [('line', [path_pts[ii], path_pts[ii + 1]]) + attr_comm]

        return ch, attrs

pucgen_classes = [
    BaseCell,
    SphericalInclusion,
    EllipsoidalInclusion,
    CylindricalInclusion,
    BoxInclusion,
    CylindricalChannel,
    SandwichLayer,
]


def main():
    fname = sys.argv[1]
    puc = PUC.from_file(fname )
    puc(fname + '.vtk')

if __name__ == "__main__":
    main()
