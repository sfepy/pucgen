import os
import sys
import numpy as nm

from mshio import msh_read
from vtkio import vtk_write

def vec2list(v):
    return [v.x, v.y, v.z]


def get_periodic(shape, nvec, shift, tol=1e-9):
    def get_bb(b):
        return nm.array([b.XMin, b.XMax, b.YMin, b.YMax, b.ZMin, b.ZMax])

    Faces = shape.Faces
    e_tab = {ii.hashCode(): k for k, ii in enumerate(shape.Edges)}

    nv = nm.array([nm.abs(f.normalAt(0, 0)) for f in Faces])
    cm = nm.array([vec2list(f.CenterOfMass) for f in Faces])

    if isinstance(nvec, str) and nvec == 'all':
        nvecs = nm.unique(nv, axis=0)
    else:
        nvecs = nm.asarray([nvec])

    pfaces = []
    for nvec in nvecs:
        idxs = list(nm.where(nm.linalg.norm(nv - nvec, axis=1) < tol)[0])
        while len(idxs) > 0:
            face = idxs.pop(0)
            for pface in idxs:
                dist = nm.abs(cm[face] - cm[pface]) - nm.abs(nvec) * shift
                if nm.linalg.norm(dist) < tol:
                    darea = Faces[face].Area - Faces[pface].Area
                    if darea / Faces[face].Area < tol:
                        bb1 = get_bb(Faces[face].BoundBox)
                        bb2 = get_bb(Faces[face].BoundBox)
                        dbb = nm.linalg.norm(bb1 - bb2) / nm.linalg.norm(bb1)
                        if dbb < tol:
                            pfaces.append((face, pface))

    print('periodic faces: %s' % pfaces)
    out = []

    for pf in pfaces:
        feds = []
        for f in pf:
            eds = []
            cps = []
            for ed in shape.Faces[f].Edges:
                eds.append(e_tab[ed.hashCode()])
                cps.append(vec2list(ed.BoundBox.Center))

            feds.append((eds, nm.asarray(cps)))

        eds1, eds2 = feds
        eds = []
        for ied, ed in enumerate(eds1[0]):
            dist = nm.abs(eds2[1] - eds1[1][ied]) - nm.abs(nvec) * shift
            idx = nm.where(nm.linalg.norm(dist, axis=1) < tol)[0]
            if len(idx) == 1:
                eds.append([ed, eds2[0][idx[0]]])
            else:
                raise ValueError('Surfaces %d, %d not periodic!' % pf)
        
        out.append([pf, nm.array(eds).T])
    
    return out


def get_periodic_gmsh(in_data, periodicity, tol=1e-9):
    FREECADPATH = '/usr/lib/freecad/lib/'
    sys.path.append(FREECADPATH)
    import FreeCAD
    import Part

    if type(in_data) is not Part.Shape:
        shape = Part.Shape()
        shape.read(in_data)
    else:
        shape = in_data

    geo_out = ''
    for nvec, shift in periodicity:
        per_faces = get_periodic(shape, nvec, shift, tol=tol)
        for (f1, f2), edges in per_faces:
            es1 = ', '.join([str(ii + 1) for ii in edges[0]])
            es2 = ', '.join([str(ii + 1) for ii in edges[1]])
            geo_out += 'Periodic Surface %d {%s} = %d {%s};\n'\
                % (f1 + 1, es1, f2 + 1, es2)

    return geo_out


def mirror_mesh(nodes, elems, data, c0=0, dim=0):
    remap_tetra = nm.array([0, 2, 1, 3])
    nnodes = nodes.shape[0]
    melems = nm.vstack([elems, (elems + nnodes)[:, remap_tetra]])
    nodes2 = nodes.copy()
    nodes2[:, dim] = 2 * c0 - nodes2[:, dim]
    mnodes = nm.vstack([nodes, nodes2])

    mdata = {}
    for k, v in data.iteritems():
        if len(v[2].shape) > 1:
            mdata[k] = (v[0], v[1], nm.vstack([v[2], v[2]]))
        else:
            mdata[k] = (v[0], v[1], nm.hstack([v[2], v[2]]))

    return mnodes, melems, mdata


def find_master_slave(nodes, tol=1e-9):
    from scipy.spatial import cKDTree

    tr = cKDTree(nodes)
    mtx = tr.sparse_distance_matrix(tr, tol).tocsr()
    nrow = nm.diff(mtx.indptr)
    idxs = nm.where(nrow > 1)[0]

    npairs_max = nm.sum(nrow[idxs] - 1)

    out = nm.empty((npairs_max, 2), dtype=nm.int64)
    idx0 = 0
    for ii in idxs:
        i1, i2 = mtx.indptr[ii], mtx.indptr[ii + 1]
        cols = mtx.indices[i1:i2]
        if cols[cols < ii].shape[0] == 0:
            nc = cols.shape[0]
            if nc == 2:
                out[idx0,:] = cols
                idx0 += 1
            else:
                idx1 = idx0 + nc - 1
                out[idx0:idx1, 0] = cols[0]
                out[idx0:idx1, 1] = cols[1:]
                idx0 = idx1

    return out[:idx0, :]


def merge_nodes(nodes, elems, data=None, tol=1e-9):
    ms_tab = find_master_slave(nodes, tol=tol)
    remap = nm.ones((nodes.shape[0],), dtype=nm.int32)
    remap[ms_tab[:, 1]] = -1
    ndidxs = nm.where(remap > 0)[0]
    remap[ndidxs] = nm.arange(len(ndidxs))
    new_nodes = nodes[ndidxs, :]
    remap[ms_tab[:, 1]] = remap[ms_tab[:, 0]]

    if data is not None:
        for k in data.keys():
            v = data[k]
            if v[0] == 'p':
                data[k] = ('p', v[1], v[2][ndidxs, ...])

    return new_nodes, remap[elems], data


def mesh_size_by_fields(fields, esize):
    out = ''
    out += 'Mesh.CharacteristicLengthFromPoints = 0;\n'
    out += 'Mesh.CharacteristicLengthFromCurvature = 0;\n'
    out += 'Mesh.CharacteristicLengthExtendFromBoundary = 0;\n'

    l_id = 100000
    p_id = 100000
    f_id = 2
    f_list = []

    if fields == 'const':
        f_id, f_list = 2, [2]
        out += 'Field[%d] = MathEval;\n' % f_id
        out += 'Field[%d].F = "%f";\n' % (f_id, esize)

    if len(fields) > 0:
        for ifield in fields:
            if ifield[0] == 'box':
                out += 'Field[%d] = Box;\n' % f_id
                out += 'Field[%d].XMin = %e;\n' % (f_id, ifield[1][0])
                out += 'Field[%d].XMax = %e;\n' % (f_id, ifield[2][0])
                out += 'Field[%d].YMin = %e;\n' % (f_id, ifield[1][1])
                out += 'Field[%d].YMax = %e;\n' % (f_id, ifield[2][1])
                out += 'Field[%d].ZMin = %e;\n' % (f_id, ifield[1][2])
                out += 'Field[%d].ZMax = %e;\n' % (f_id, ifield[2][2])
                out += 'Field[%d].VIn = %e;\n' % (f_id, ifield[4] * esize)
                out += 'Field[%d].VOut = %e;\n' % (f_id, esize)
                out += 'Field[%d].Thickness = %e;\n' % (f_id, ifield[3])
                f_list.append(f_id)
                f_id += 1
            elif ifield[0] == 'line':
                pts = ifield[1]
                out += 'Point(%d) = {%e, %e, %e};\n' % ((p_id,) + tuple(pts[0]))
                out += 'Point(%d) = {%e, %e, %e};\n' % ((p_id + 1,) + tuple(pts[1]))
                out += 'Line(%d) = {%d, %d};\n' % (l_id, p_id, p_id + 1) 
                out += 'Field[%d] = Attractor;\n' % f_id
                out += 'Field[%d].EdgesList = {%d};\n' % (f_id, l_id)
                f_id += 1
                out += 'Field[%d] = Threshold;\n' % f_id
                out += 'Field[%d].IField = %d;\n' % (f_id, f_id - 1)
                out += 'Field[%d].LcMin = %e;\n' % (f_id, ifield[4] * esize)
                out += 'Field[%d].LcMax = %e;\n' % (f_id, esize)
                out += 'Field[%d].DistMin = %e;\n' % (f_id, ifield[2])
                out += 'Field[%d].DistMax = %e;\n' % (f_id, ifield[3])
                f_list.append(f_id)
                f_id += 1
                p_id += 2
                l_id += 1
            elif ifield[0] == 'sphere':
                pt = ifield[1]
                out += 'Point(%d) = {%e, %e, %e};\n' % ((p_id,) + tuple(pt))
                out += 'Field[%d] = Attractor;\n' % f_id
                out += 'Field[%d].NodesList = {%d};\n' % (f_id, p_id)
                f_id += 1
                out += 'Field[%d] = Threshold;\n' % f_id
                out += 'Field[%d].IField = %d;\n' % (f_id, f_id - 1)
                out += 'Field[%d].LcMin = %e;\n' % (f_id, ifield[4] * esize)
                out += 'Field[%d].LcMax = %e;\n' % (f_id, esize)
                out += 'Field[%d].DistMin = %e;\n' % (f_id, ifield[2])
                out += 'Field[%d].DistMax = %e;\n' % (f_id, ifield[3])
                f_list.append(f_id)
                f_id += 1
                p_id += 1

    out += 'Field[1] = Min;\n'
    out += 'Field[1].FieldsList = {%s};\n' % ', '.join([str(ii) for ii in f_list])
    out += 'Background Field = 1;\n'

    return out

def gmsh_call(filename_base, filename_out, fields, esize,
              periodicity=None, periodicity_tol=1e-6,
              node_groups=None, merge=False, dim=3,
              export_elems='3_4', only_elems_by_nodes=None,
              physical_volumes=None, shift=None, scale=None):

    def label_nodes(node_selection, nodes):
        ndgrp = nm.zeros((nodes.shape[0],), dtype=nm.int32)
        if node_selection is not None:
            for ng, sel in node_selection:
                aux = sel.strip().split(' ')
                dir = {'x': 0, 'y': 1, 'z': 2}[aux[0]]
                expr = ' '.join(aux[1:])
                idxs = eval('nm.where(nodes[:, %d] %s)[0]' % (dir, expr))
                ndgrp[idxs] = ng

        return ndgrp

    geo_file = 'Merge "%s.step";\n' % filename_base

    if fields is not None:
        geo_file += mesh_size_by_fields(fields, esize)

    if periodicity is not None:
        geo_file += get_periodic_gmsh(filename_base + '.step',
                                      periodicity, tol=periodicity_tol)
    if physical_volumes is not None:
        for k, v in physical_volumes.items():
            geo_file += 'Physical Volume(%d)={%s};\n' % (k, ', '.join(map(str, v)))

    open(filename_base + '.geo', 'wt').write(geo_file)

    os.system('gmsh -%d %s.geo -o %s.msh -format msh22'
              % (dim, filename_base, filename_base))

    nodes, elems0, _ = msh_read(filename_base + '.msh')
    elems, mats = elems0[export_elems][0], elems0[export_elems][1]
    ndgrp = label_nodes(node_groups, nodes)

    data = {
        'mat_id': ('c', 'int', mats),
        'node_groups': ('p', 'int', ndgrp)
    }

    if merge:
        # print('merging nodes')
        nodes, elems, data = merge_nodes(nodes, elems, data,
                                         tol=periodicity_tol)

    if only_elems_by_nodes is not None:
        export_elems = only_elems_by_nodes[0]
        elems, mats = elems0[export_elems][0], elems0[export_elems][1]
        ndgrp = label_nodes([(1, only_elems_by_nodes[1])], nodes)
        aux = ndgrp[elems].sum(axis=1)
        eidxs = nm.where(aux == int(export_elems[-1]))[0]
        nidxs = nm.where(ndgrp == 1)[0]
        data['node_groups'] = ('p', 'int', ndgrp[nidxs] * 0)
        data['mat_id'] = ('c', 'int', mats[eidxs])
        remap = -nm.ones((nodes.shape[0],), dtype=nm.int32)
        remap[nidxs] = nm.arange(len(nidxs))
        nodes = nodes[nidxs, :]
        elems = remap[elems[eidxs, :]]

    if shift is not None:
        nodes += shift

    if scale is not None:
        nodes *= scale

    vtk_write(filename_out, nodes, elems, export_elems, data)
