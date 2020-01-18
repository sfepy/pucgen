#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PUCGEN - GUI for Periodic Unit Cell generator

"""

import sys
import os
import numpy as nm
from inspect import getargspec
from ast import literal_eval

from gen_mesh_utils import gmsh_call, repeater


def l2s(l):
    return str(list(l))[1:-1]

def b2s_delete(b):
    return 'Delete;' if b else ''

def geo_obj(obj, pars):
    return '%s(%%d) = {%s};' % (obj, l2s(pars))

def volumes_union(l1, l2, el_size, delete1=True, delete2=True):
    esize = nm.min([el_size[k] for k in l1 + l2])
    return "BooleanUnion(%%d) = {Volume{%s}; %s}{Volume{%s}; %s};"\
        % (l2s(l1), b2s_delete(delete1), l2s(l2), b2s_delete(delete2)), esize

def volumes_difference(l1, l2, el_size, delete1=True, delete2=True):
    esize = nm.min([el_size[k] for k in l1])
    return "BooleanDifference(%%d) = {Volume{%s}; %s}{Volume{%s}; %s};"\
        % (l2s(l1), b2s_delete(delete1), l2s(l2), b2s_delete(delete2),), esize

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

    def __call__(self, filename_vtk, cell_size=None, eps=1.0, centered=False):
        """Generate the finite element mesh.

        Parameters
        ----------
        filename_vtk: str
            The VTK output file name.
        cell_size: array
            The size of PUC: [sx, sy, sz].
        eps: float
            The scaling parameter.
        centered: bool
            If True, the PUC is centered to the origin.
        """
        if cell_size is None:
            cell_size = self.components[0].get('size')

        element_size = self.components[0].get('el_size')
        cell_size = nm.asarray(cell_size, dtype=nm.float64)

        geo = []

        mat_ids = []
        el_size = {}
        volumes = {}
        vid = 1
        bcell_mat_id = self.components[0].params['mat_id']
        for comp in self.components:
            if not comp.active:
                continue

            mat_id = comp.params['mat_id']
            geo_line, esize = comp(vid, cell_size)
            if geo_line is None:
                continue

            geo.append(geo_line)
            el_size[vid] = esize

            if mat_id in volumes:
                volumes[mat_id].append(vid)
            else:
                volumes[mat_id] = [vid]
                mat_ids.append(mat_id)

            vid += 1

        bcell = volumes[bcell_mat_id]
        if len(bcell) == 2:
            bcell_ctool = bcell[1]
        elif len(bcell) >= 3:
            geo_line, esize = volumes_union([bcell[1]], bcell[2:], el_size)
            geo.append(geo_line % vid)
            el_size[vid] = esize
            bcell_ctool = vid
            vid += 1
        else:
            bcell_ctool = None

        objs = [bcell[0]]
        for mat_id in mat_ids[1:]:
            if len(volumes[mat_id]) >= 2:
                geo_line, esize = volumes_union([volumes[mat_id][0]],
                                                 volumes[mat_id][1:], el_size)
                geo.append(geo_line % vid)
                el_size[vid] = esize
                uni = vid
                vid += 1
            else:
                uni = volumes[mat_id][0]

            if bcell_ctool is None:
                objs.append(uni)
            else:
                geo_line, esize = volumes_difference([uni], [bcell_ctool],
                                                     el_size)
                geo.append(geo_line % vid)
                el_size[vid] = esize
                objs.append(vid)
                vid += 1

        if len(objs) > 1:
            geo_line, esize = volumes_difference([objs[0]], objs[1:],
                                                 el_size, delete2=False)
            geo.append(geo_line % vid)
            el_size[vid] = esize
            objs[0] = vid
            vid += 1

        geo.append('Coherence;')

        geo.append('')
        for obj, mat_id in zip(objs, mat_ids):
            geo.append('Physical Volume(%d) = {%d};' % (mat_id, obj))

        geo.append('')
        for obj in objs[1:]:
            geo.append('p%d() = PointsOf{Volume{%d};};' % (obj, obj))
            geo.append('Characteristic Length{p%d()} = %e;'
                % (obj, element_size * el_size[obj]))

        peps = nm.max(cell_size) * 1e-3
        sib = 'Surface In BoundingBox'
        p1 = (cell_size + peps * nm.ones(3)) * 0.5
        p0 = -p1

        geo.append('')
        for idir, per in enumerate(['x', 'y', 'z']):
            pdir = nm.eye(3)[idir] * cell_size
            p2 = p1 - pdir
            geo.append('per%s%d() = %s{%e,%e,%e,%e,%e,%e};'
                % ((per, 1, sib) + tuple(p0) + tuple(p2)))
            geo.append('per%s%d() = %s{%e,%e,%e,%e,%e,%e};'
                % ((per, 2, sib) + tuple(p0 + pdir) + tuple(p2 + pdir)))
            geo.append('Periodic Surface{per%s%d()} = {per%s%d()} Translate{%e,%e,%e};'
                % ((per, 2, per, 1) + tuple(pdir)))

        esize_min = nm.min([k for k in el_size.values() if k is not None])
        geo_content = [
            'SetFactory("OpenCASCADE");',
            'Mesh.CharacteristicLengthMin = %e;' % element_size * esize_min,
            'Mesh.CharacteristicLengthMax = %e;' % element_size,
            ''] + geo

        filename_base = os.path.splitext(filename_vtk)[0]

        with open(filename_base + '.geo', 'wt') as f:
            f.write('\n'.join(geo_content))

        if not centered:
            shift = nm.asarray(cell_size) * 0.5

        gmsh_call(filename_base, shift=shift, scale=eps)

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



class BaseCell(BaseComponent):
    """The base cell - matrix."""
    name = 'Base Cell'

    def __init__(self, size=(1, 1, 1), el_size=0.1, mat_id=1):
        super(BaseCell, self).__init__(mat_id=mat_id)
        self.params.update({
            'size': size,
            'el_size': el_size,
        })

    def __call__(self, vid, size):
        size = nm.asarray(size)
        pars = list(-0.5 * size) + list(size)

        return geo_obj('Box', pars) % vid, None


class BaseEmbeddedComponent(BaseComponent):
    """The base for the inclusion and channel classes."""

    direction_tab = {'x': nm.array([1, 0, 0]),
                     'y': nm.array([0, 1, 0]),
                     'z': nm.array([0, 0, 1])}

    def __init__(self, parameters, central_point, direction, el_size, mat_id):
        """Init parameters of the channel component.
    
        Parameters
        ----------
        parameters: float or array
            The geometrical parameters of the object.
        central_point: array
            The coordinates of the object center: [x, y, z].
        direction: str or array
            The object direction. If string: direction = 'x', 'y' or 'z'.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """    
        super(BaseEmbeddedComponent, self).__init__(mat_id=mat_id)

        if isinstance(direction, list) or isinstance(direction, tuple):
            direction = nm.asarray(direction) / nm.linalg.norm(direction)

        self.params.update({
            'parameters': parameters,
            'direction': direction,
            'central_point': nm.asarray(central_point, dtype=nm.float64),
            'el_size': el_size,
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
                 direction=(1, 0, 0), el_size=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        radius: float or array
            The radii of the ellipsoid: r or [r1, r2, r3].
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: array
            The directional vector.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super(EllipsoidalInclusion, self).__init__(parameters=[radius],
                                                   central_point=central_point,
                                                   direction=direction,
                                                   el_size=el_size,
                                                   mat_id=mat_id)

    def __call__(self, vid, size):
        if nm.all(nm.array(self.get('radius')) > 0):
            p = list(self.get('central_point'))
            r = list(self.get('radius'))
            d = self.get('direction')
            e = nm.eye(3)[0]
            geo = geo_obj('Sphere', p + [1.]) % vid
            geo += ' Dilate {{%s}, {%s}} {Volume{%d};}' % (l2s(p), l2s(r), vid)
            ax = nm.cross(e, d)
            if nm.linalg.norm(ax) > 0.0:
                phi = nm.arccos(nm.dot(e, d) / (nm.linalg.norm(e) * nm.linalg.norm(d)))
                geo += ' Rotate {{%s}, {%s}, %e} {Volume{%d};}'\
                    % (l2s(ax), l2s(p), phi, vid)

            return geo, self.get('el_size')

class SphericalInclusion(BaseEmbeddedComponent):
    """The spherical inclusion."""
    name = 'Spherical Inclusion'
    parameters_dict = {'radius': 0}

    def __init__(self, radius=0.1, central_point=(0, 0, 0),
                 el_size=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        radius: float
            The radius of the sphere.
        central_point: array
            The coordinates of the center: [x, y, z].
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super(SphericalInclusion, self).__init__(parameters=[radius],
                                                 central_point=central_point,
                                                 direction=None,
                                                 el_size=el_size,
                                                 mat_id=mat_id)

    def __call__(self, vid, size):
        if self.get('radius') > 0:
            pars = list(self.get('central_point')) + [self.get('radius')]
            return geo_obj('Sphere', pars) % vid, self.get('el_size')

class CylindricalInclusion(BaseEmbeddedComponent):
    """The cylindrical inclusion."""
    name = 'Cylindrical Inclusion'
    parameters_dict = {'radius': 0, 'length': 1}

    def __init__(self, radius=0.1, length=0.5, central_point=(0, 0, 0),
                 direction=(1, 0, 0), el_size=0.5, mat_id=2):
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
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super(CylindricalInclusion, self).__init__(parameters=[radius, length],
                                                   central_point=central_point,
                                                   direction=direction,
                                                   el_size=el_size,
                                                   mat_id=mat_id)

    def __call__(self, vid, size):
        r = self.get('radius')
        h = self.get('length')
        d = self.get('direction')
        p = self.get('central_point')

        if nm.all(r > 0):

            if isinstance(d, str):
                idir = {'x': 0, 'y': 1, 'z': 2}[d]
                d = self.direction_tab[d]
                if h is None:
                    h = size[idir]
                    p = p.copy()
                    p[idir] = 0

            p0 = p - 0.5 * d * h
            pars = list(p0) + list(d * h) + [r]

            return geo_obj('Cylinder', pars) % vid, self.get('el_size')


class CylindricalChannel(CylindricalInclusion):
    """The cylindrical channel."""
    name = 'Cylindrical Channel'

    def __init__(self, radius=0.1, central_point=(0, 0, 0), direction='x',
                 el_size=0.5, mat_id=2):
        """Init parameters of the channel component.
    
        Parameters
        ----------
        radius: flaot
            The radius of the cylinder.
        central_point: array
            The coordinates of the cylinder center: [x, y, z].
        direction: str
            The cylinder orientation: 'x', 'y' or 'z'.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super(CylindricalChannel, self).__init__(radius=radius, length=None,
                                                 central_point=central_point,
                                                 direction=direction,
                                                 el_size=el_size,
                                                 mat_id=mat_id)


class BoxInclusion(BaseEmbeddedComponent):
    """The box inclusion."""
    name = 'Box Inclusion'
    parameters_dict = {'size': 0}
    
    def __init__(self, size=(0.3, 0.2, 0.1), central_point=(0, 0, 0),
                 el_size=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        size: array
            The size of the box: [sx, sy, sz].
        central_point: array
            The coordinates of the center: [x, y, z].
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super(BoxInclusion, self).__init__(parameters=[size],
                                           central_point=central_point,
                                           direction=None,
                                           el_size=el_size,
                                           mat_id=mat_id)

    def __call__(self, vid, size):
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

        if nm.all(s > 0):
            pars = list(p - 0.5 * s) + list(s)
            return geo_obj('Box', pars) % vid, self.get('el_size')

class SandwichLayer(BoxInclusion):
    """The sandwich layer."""
    name = 'Sandwich Layer'
    parameters_dict = {'thickness': 0}
    
    def __init__(self, thickness=0.1, central_point=(0, 0, 0),
                 direction='x', el_size=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        thickness: array
            The thicknesss of the layer.
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: str
            The orientation of the layer normal vector: 'x', 'y' or 'z'.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super(SandwichLayer, self).__init__(size=thickness,
                                            central_point=central_point,
                                            el_size=el_size,
                                            mat_id=mat_id)
        self.params['direction'] = direction


class SweepedChannel(BaseEmbeddedComponent):
    """The sweeped channel."""
    name = 'Sweeped Channel'

    def __init__(self, profile=[[0.1, 0], [0, 0.1], [-0.1, -0.1]],
                 path=[[0, 0], [0.3, 0.1], [1, 0]], central_point=(0, 0, 0),
                 direction='x', es_dmin=1.1, es_dmax=1.3, el_size=0.5, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        !!!

        central_point: array
            The coordinates of the center: [x, y, z].
        direction: array
            The directional vector.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super(SweepedChannel, self).__init__(parameters=(profile, path),
                                             central_point=central_point,
                                             direction=direction,
                                             es_dmin=es_dmin,
                                             es_dmax=es_dmax,
                                             el_size=el_size,
                                             mat_id=mat_id)

    # def __call__(self, vid, size):
    #     profile = self.get('profile')
    #     path = self.get('path')
    #     d = self.get('direction')
    #     p = self.get('central_point')
    #     label = '%s_%d' % (self.name, self.get('mat_id'))

    #     attrs = []

    #     if isinstance(d, str):
    #         idir = {'x': 0, 'y': 1, 'z': 2}[d]
    #         d = self.direction_tab[d]
    #         h = size[idir]
    #         p = p.copy()
    #         p[idir] = 0

    #     sk_profile = self.get_sketch(profile[1], 'profile', cont,
    #                                  is_closed=True, path_type=profile[0],
    #                                  plane='yz')

    #     path_pts = nm.asarray(path[1])
    #     path_pts[:, 0] = (path_pts[:, 0] - 0.5) * h

    #     sp_path = self.get_sketch(path_pts, 'path', cont,
    #                               is_closed=True, path_type=path[0],
    #                               plane='xy')

    #     ch = cont.addObject('Part::Sweep', label)
    #     ch.Sections = [sk_profile]
    #     ch.Spine = (sp_path,[])
    #     ch.Solid = True
    #     ch.Frenet = False

    #     attr_comm = (r * self.get('es_dmin'), r * self.get('es_dmax'),
    #                  self.get('el_size'))
    #     p1 = p - 0.5 * d * h * 1.2
    #     p2 = p - 0.5 * d * h * 1
    #     for ii in range(len(path_pts)):
    #         attrs = [('line', [path_pts[ii], path_pts[ii + 1]]) + attr_comm]

    #     return ch, attrs

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
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        puc = PUC.from_file(filename)
        filename_vtk = os.path.splitext(filename)[0] + '.vtk'
        puc(filename_vtk)
    elif len(sys.argv) == 4:
        filename_in = sys.argv[1]
        filename_out = sys.argv[2]
        grid = literal_eval(sys.argv[3])
        scale_x = float(sys.argv[4])
        repeater(filename_in, filename_out, grid, scale_x)
    else:
        from pucgen_gui import MainWindow
        from PyQt5.QtWidgets import QApplication

        app = QApplication(sys.argv)
        mw = MainWindow()
        sys.exit(app.exec_())

if __name__ == "__main__":
    main()
