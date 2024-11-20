#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PUCGEN - GUI for Periodic Unit Cell generator

"""

import sys
import os
from optparse import OptionParser
import numpy as nm
import gmsh
import meshio
from inspect import signature
from ast import literal_eval
from gen_mesh_utils import repeat_cell


class PUC(object):
    """Periodic Unit Cell object."""
    def __init__(self, cell_mat_id=1, base_cell=None):
        """Init PUC"""
        gmsh.initialize()
        gmsh.logger.start()
        gmsh.model.add('boolean')

        self.model = gmsh.model
        self.occ = gmsh.model.occ

        if cell_mat_id is None:
            self.components = []
        else:
            self.components = [BaseCell(mat_id=cell_mat_id)\
                               if base_cell is None else base_cell]

    def __str__(self):
        comps = [str(k) for k in self.components]
        return f'[{", ".join(comps)}]'

    def finalize(self):
        print('\n'.join(gmsh.logger.get()))
        gmsh.logger.stop()
        gmsh.finalize()

    def add(self, obj):
        """Add a new component to the list (inclusion, channel, layer)."""
        self.components.append(obj)

    @staticmethod
    def save_puc(filename, comps):
        """Save PUC into the text file."""
        with open(filename, 'wt', encoding='utf-8') as f:
            for cls, pars, act in comps:
                aflag = '' if act else '#'
                targs = []
                for key, val in pars.items():
                    if key == 'self':
                        continue
                    if isinstance(val, nm.ndarray):
                        targs.append(key + '=' + str(tuple(val)))
                    else:
                        targs.append(key + '=' + str(val))
                f.write(f'{aflag}{cls.__name__};{";".join(targs)}\n')

    def save(self, filename):
        self.save_puc(filename, [(c.__class__, c, c.active)
            for c in self.components])

    @staticmethod
    def load_puc(filename):
        """Load PUC from a given text file."""

        cls_dict = {c.__name__: c for c in pucgen_classes}
        out = []
        with open(filename, 'rt', encoding='utf-8') as f:
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

    @staticmethod
    def get_gmsh3(objs, ent=3):
        return [(ent, k) for k in objs]

    @staticmethod
    def get_obj3(objs):
        return [k for _, k in objs]

    def fuse(self, obj1, obj2):
        out, _ = self.occ.fuse(self.get_gmsh3(obj1), self.get_gmsh3(obj2))
        return self.get_obj3(out)

    def cut(self, obj1, obj2):
        out, _ = self.occ.cut(self.get_gmsh3(obj1), self.get_gmsh3(obj2),
                              removeTool=False)
        return self.get_obj3(out)

    def fragment(self, obj1, obj2):
        _, ovv = self.occ.fragment(self.get_gmsh3(obj1), self.get_gmsh3(obj2))

        out2 = [k.copy() for k in ovv]
        for g in out2[1:]:
            for k in g:
                if k in out2[0]:
                    del out2[0][out2[0].index(k)]

        return [self.get_obj3(k) for k in out2]


    def set_periodic(self, cell_size, dim, eps=1e-3):
        bbox = nm.array([cell_size * 0, cell_size])
        bbox -= cell_size * 0.5

        eye = nm.eye(4)
        eye[dim, 3] = 1

        lnb = bbox[0] - eps
        rft = bbox[1] + eps

        rft1 = rft.copy()
        rft1[dim] = lnb[dim] + 2*eps

        lnb2 = lnb.copy()
        lnb2[dim] = rft[dim] - 2*eps

        pars = list(lnb) + list(rft1) + [2]
        p1 = self.model.getEntitiesInBoundingBox(*pars)
        pars = list(lnb2) + list(rft) + [2]
        p2 = self.model.getEntitiesInBoundingBox(*pars)

        self.model.mesh.setPeriodic(2, self.get_obj3(p2), self.get_obj3(p1),
                                    eye.ravel())

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
            cell_size = self.components[0].get('dimension')

        cell_size = nm.asarray(cell_size, dtype=nm.float64)

        mat_ids = []
        volumes = {}
        esize = {}
        for comp in self.components:
            if not comp.active:
                continue

            mat_id = comp.params['mat_id']
            obj, es = comp(self.occ, cell_size)
            if obj is None:
                continue

            if mat_id in volumes:
                volumes[mat_id].append(obj)
                esize[mat_id].append(es)
            else:
                volumes[mat_id] = [obj]
                esize[mat_id] = [es]
                mat_ids.append(mat_id)

        fvolumes = []
        fesize = []
        cut_tool = None
        for mat_id in mat_ids:
            vols = volumes[mat_id]
            es = esize[mat_id]
            if len(vols) >= 2:
                if mat_id == mat_ids[0]:
                    if len(vols) >= 3:
                        cut_tool = self.fuse(vols[1:2], vols[2:])[0]
                    else:
                        cut_tool = vols[1]

                    fvolumes.append(vols[0])
                    fesize.append(es[0])
                else:
                    fvolumes.append(self.fuse(vols[:1], vols[1:])[0])
                    fesize.append(nm.average(es))
            else:
                fvolumes.append(vols[0])
                fesize.append(es[0])

        if cut_tool is not None:
            fvolumes = fvolumes[:1] + [self.cut([vol], [cut_tool])[0]
                                       for vol in fvolumes[1:]]

        out = self.fragment(fvolumes[:1], fvolumes[1:])

        self.occ.synchronize()

        model = self.model
        for k, objs_ in enumerate(out):
            mat_id = mat_ids[k]
            es = fesize[k]
            out2 = out.copy()
            del out2[k]
            others = set(sum(out2, []))
            objs = [kk for kk in objs_ if kk not in others or len(objs_) == 1]
            model.addPhysicalGroup(3, objs, mat_id)
            pts = model.getBoundary(self.get_gmsh3(objs), False, False, True)
            model.mesh.setSize(pts, es)

        for d in range(3):
            self.set_periodic(cell_size, d)

        self.model.mesh.generate(3)

        filename_base = os.path.splitext(filename_vtk)[0]
        fname = f'{filename_base}.msh'
        gmsh.write(fname)

        self.finalize()

        mesh = meshio.read(fname)
        os.remove(fname)

        if not centered:
            mesh.points += nm.asarray(cell_size) * 0.5

        mesh.point_data = {}
        mesh.cell_sets = {}
        mesh.cell_data = {'mat_id': mesh.cell_data['gmsh:physical']}

        fname = f'{filename_base}.vtk'
        mesh.write(fname, binary=False)
        print(f"VTK mesh saved to '{fname}'")


class BaseComponent(object):
    """The base component of the unit cell."""
    name = None

    def __init__(self, mat_id=1):
        """Init parameters of the component.
    
        Parameters
        ----------
        mat_id: int
            The component material id.
        """
        self.params = {'mat_id': mat_id}
        self.active = True

    def __call__(self, vid, size):
        """Create the GEO file representation of a object.
    
        Parameters
        ----------
        vid: int
            The volume indentificator.
        size: array
            The size of the cell: [size_x, size_y, size_z].

        Returns
        -------
        obj: str
            The string encoding a GEO file object.
        el_size: float
            The element size factor.
        """

    def __str__(self):
        flag = '#' if self.active is False else ''
        return flag + self.name

    def get(self, key):
        return self.params[key]


class BaseCell(BaseComponent):
    """The base cell - matrix."""
    name = 'Base Cell'

    def __init__(self, dimension=(1, 1, 1), el_size=0.1, mat_id=1):
        super().__init__(mat_id=mat_id)
        self.params.update({
            'dimension': nm.asarray(dimension),
            'el_size': el_size,
        })

    def __call__(self, occ, cell_size=None):
        size = self.get('dimension')
        pars = list(-0.5 * size) + list(size)

        return occ.addBox(*pars), self.get('el_size')


class BaseEmbeddedComponent(BaseComponent):
    """The base for the inclusion and channel classes."""

    def __init__(self, dimension, central_point, direction, el_size, mat_id):
        """Init parameters of the channel component.
    
        Parameters
        ----------
        dimension: float or array
            The dimension of the object.
        central_point: array
            The coordinates of the object center: [x, y, z].
        direction: str or array
            The object direction. If string: direction = 'x', 'y' or 'z'.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super().__init__(mat_id=mat_id)

        if isinstance(direction, (list, tuple)):
            direction = nm.asarray(direction) / nm.linalg.norm(direction)

        self.params.update({
            'dimension': dimension,
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

    def __init__(self, dimension=(0.1, 0.1, 0.1), central_point=(0, 0, 0),
                 direction=(1, 0, 0), el_size=0.05, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        dimension: float or array
            The radii of the ellipsoid: r or [r1, r2, r3].
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: array
            The directional vector.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super().__init__(dimension=dimension, central_point=central_point,
                         direction=direction, el_size=el_size, mat_id=mat_id)

    def __call__(self, occ, cell_size):
        if nm.all(nm.array(self.get('dimension')) > 0):
            p = list(self.get('central_point'))
            r = list(self.get('dimension'))
            d = self.get('direction')
            e = nm.eye(3)[0]
            pars = p + [1.]
            out = occ.addSphere(*pars)
            pars = [[(3, out)]] + p + r
            occ.dilate(*pars)
            ax = nm.cross(e, d)
            if nm.linalg.norm(ax) > 0.0:
                phi = nm.arccos(nm.dot(e, d) / (nm.linalg.norm(e) * nm.linalg.norm(d)))
                pars = [[(3, out)]] + p + list(ax) + [phi]
                occ.rotate(*pars)

            return out, self.get('el_size')

class SphericalInclusion(BaseEmbeddedComponent):
    """The spherical inclusion."""
    name = 'Spherical Inclusion'

    def __init__(self, dimension=0.1, central_point=(0, 0, 0),
                 el_size=0.05, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        dimension: float
            The radius of the sphere.
        central_point: array
            The coordinates of the center: [x, y, z].
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super().__init__(dimension=dimension, central_point=central_point,
                         direction=None, el_size=el_size, mat_id=mat_id)

    def __call__(self, occ, cell_size=None):
        r = self.get('dimension')
        if r > 0:
            pars = list(self.get('central_point')) + [r]
            return occ.addSphere(*pars), self.get('el_size')


class CylindricalInclusion(BaseEmbeddedComponent):
    """The cylindrical inclusion."""
    name = 'Cylindrical Inclusion'

    def __init__(self, dimension=(0.1, 0.5), central_point=(0, 0, 0),
                 direction=(1, 0, 0), el_size=0.05, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        dimension: (float, float)
            The cylinder radius and length.
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: array
            The directional vector.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super().__init__(dimension=dimension, central_point=central_point,
                         direction=direction, el_size=el_size, mat_id=mat_id)

    def __call__(self, occ, cell_size=None):
        r, h = self.get('dimension')
        d = self.get('direction')
        p = self.get('central_point')

        if nm.all(r > 0):
            if isinstance(d, str):
                idir = {'x': 0, 'y': 1, 'z': 2}[d]
                d = nm.eye(3)[idir]
                if h is None:
                    h = cell_size[idir]
                    p = p.copy()
                    p[idir] = 0

            p0 = p - 0.5 * d * h
            pars = list(p0) + list(d * h) + [r]

            return occ.addCylinder(*pars), self.get('el_size')


class CylindricalChannel(CylindricalInclusion):
    """The cylindrical channel."""
    name = 'Cylindrical Channel'

    def __init__(self, dimension=0.1, central_point=(0, 0, 0), direction='x',
                 el_size=0.05, mat_id=2):
        """Init parameters of the channel component.
    
        Parameters
        ----------
        dimension: float
            The radius of the cylinder.
        central_point: array
            The coordinates of the cylinder center: [x, y, z].
        direction: str
            The cylinder orientation: 'x', 'y' or 'z'.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super().__init__(dimension=(dimension, None),
                         central_point=central_point, direction=direction,
                         el_size=el_size, mat_id=mat_id)


class BoxInclusion(BaseEmbeddedComponent):
    """The box inclusion."""
    name = 'Box Inclusion'

    def __init__(self, dimension=(0.3, 0.2, 0.1), central_point=(0, 0, 0),
                 direction=None, el_size=0.05, mat_id=2):
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
        super().__init__(dimension=dimension, central_point=central_point,
                         direction=direction, el_size=el_size, mat_id=mat_id)

    def __call__(self, occ, cell_size=None):
        d = self.get('direction')
        if d is None:
            s = nm.array(self.get('dimension'))
            p = self.get('central_point')
        else:
            s = nm.array(cell_size)
            idir = {'x': 0, 'y': 1, 'z': 2}[d]
            s[idir] = self.get('dimension')
            p = nm.zeros(3, dtype=nm.float64)
            p[idir] = self.get('central_point')[idir]

        if nm.all(s > 0):
            pars = list(p - 0.5 * s) + list(s)
            return occ.addBox(*pars), self.get('el_size')

class SandwichLayer(BoxInclusion):
    """The sandwich layer."""
    name = 'Sandwich Layer'

    def __init__(self, dimension=0.1, central_point=(0, 0, 0),
                 direction='x', el_size=0.05, mat_id=2):
        """Init parameters of the component.
    
        Parameters
        ----------
        dimension: array
            The thicknesss of the layer.
        central_point: array
            The coordinates of the center: [x, y, z].
        direction: str
            The orientation of the layer normal vector: 'x', 'y' or 'z'.
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        super().__init__(dimension=dimension, central_point=central_point,
                         direction=direction, el_size=el_size, mat_id=mat_id)


# class SweepedChannel(BaseEmbeddedComponent):
#     """The sweeped channel."""
#     name = 'Sweeped Channel'

#     def __init__(self, dimension=([[0.1, 0], [0, 0.1], [-0.1, -0.1]],
#                                   [[0, 0], [0.3, 0.1], [1, 0]]),
#                  central_point=(0, 0, 0), direction='x', es_dmin=1.1,
#                  es_dmax=1.3, el_size=0.05, mat_id=2):
#         """Init parameters of the component.
    
#         Parameters
#         ----------
#         dimension:
#             ???
#         central_point: array
#             The coordinates of the center: [x, y, z].
#         direction: array
#             The directional vector.
#         el_size: float
#             The "inner" element size factor: in_el_size = el_size * el_size_base.
#         """
#         super().__init__(dimension=dimension, central_point=central_point,
#                          direction=direction,
#                          es_dmin=es_dmin, es_dmax=es_dmax, el_size=el_size,
#                          mat_id=mat_id)

#     def __call__(self, vid, size):
#         profile = self.get('profile')
#         path = self.get('path')
#         d = self.get('direction')
#         p = self.get('central_point')
#         label = '%s_%d' % (self.name, self.get('mat_id'))

#         attrs = []

#         if isinstance(d, str):
#             idir = {'x': 0, 'y': 1, 'z': 2}[d]
#             d = self.direction_tab[d]
#             h = size[idir]
#             p = p.copy()
#             p[idir] = 0

#         sk_profile = self.get_sketch(profile[1], 'profile', cont,
#                                      is_closed=True, path_type=profile[0],
#                                      plane='yz')

#         path_pts = nm.asarray(path[1])
#         path_pts[:, 0] = (path_pts[:, 0] - 0.5) * h

#         sp_path = self.get_sketch(path_pts, 'path', cont,
#                                   is_closed=True, path_type=path[0],
#                                   plane='xy')

#         ch = cont.addObject('Part::Sweep', label)
#         ch.Sections = [sk_profile]
#         ch.Spine = (sp_path,[])
#         ch.Solid = True
#         ch.Frenet = False

#         attr_comm = (r * self.get('es_dmin'), r * self.get('es_dmax'),
#                      self.get('el_size'))
#         p1 = p - 0.5 * d * h * 1.2
#         p2 = p - 0.5 * d * h * 1
#         for ii in range(len(path_pts)):
#             attrs = [('line', [path_pts[ii], path_pts[ii + 1]]) + attr_comm]

#         return ch, attrs


pucgen_classes = [
    BaseCell,
    SphericalInclusion,
    EllipsoidalInclusion,
    CylindricalInclusion,
    BoxInclusion,
    CylindricalChannel,
    SandwichLayer,
]

usage = 'Usage: %prog [[options] filename_in]'
version = '0.2'
helps = {
    'reps': 'construct grid by repeating unit cell, number of repetition defined by NX, NY, NZ',
    'sizex': 'resize geometry uniformly such that its size in x-direction is SIZE_X',
    'filename_out': 'write VTK output to FILE',
}

def main():
    parser = OptionParser(usage=usage, version='%prog ' + version)
    parser.add_option('-t', '--tile', metavar='"NX,NY,NZ"',
                      action='store', dest='reps', default=None,
                      help=helps['reps'])
    parser.add_option('-s', '--sizex', metavar='SIZE_X',
                      action='store', dest='size_x', default=None,
                      help=helps['sizex'])
    parser.add_option('-o', '--output', metavar='FILE',
                      action='store', dest='filename_out', default=None,
                      help=helps['filename_out'])

    (options, args) = parser.parse_args()

    if len(args) == 0: # run GUI
        from pucgen_gui import MainWindow
        from PyQt6.QtWidgets import QApplication

        app = QApplication(sys.argv)
        mw = MainWindow()
        mw.show()
        app.exec()
    else:
        filename_base, filename_ext = os.path.splitext(args[0])

        size_x = float(options.size_x) if options.size_x is not None else None

        if options.filename_out is not None:
            filename_out = options.filename_out
        else:
            filename_out = filename_base + '.vtk'

        if filename_ext == '.puc': # run generator
            puc = PUC.from_file(args[0])
            puc(filename_out, eps=size_x)

        if options.reps is not None or size_x is not None:
            filename_in = args[0] if filename_ext == '.vtk' else filename_out
            reps = literal_eval(options.reps)\
                if options.reps is not None else None
            repeat_cell(filename_in, filename_out, reps, size_x)


if __name__ == "__main__":
    main()
