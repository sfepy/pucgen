import os
import sys
import numpy as nm

from mshio import msh_read
from vtkio import vtk_read, vtk_write


def mirror_mesh(nodes, elems, data, c0=0, dim=0):
    remap_tetra = nm.array([0, 2, 1, 3])
    nnodes = nodes.shape[0]
    melems = nm.vstack([elems, (elems + nnodes)[:, remap_tetra]])
    nodes2 = nodes.copy()
    nodes2[:, dim] = 2 * c0 - nodes2[:, dim]
    mnodes = nm.vstack([nodes, nodes2])

    mdata = {}
    for k, v in data.items():
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


def gmsh_call(filename_base, dim=3, export_elems='3_4',
              node_groups=None, shift=None, scale=None):

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

    os.system('gmsh -%d %s.geo -o %s.msh -format msh22'
              % (dim, filename_base, filename_base))

    nodes, elems0, _ = msh_read(filename_base + '.msh')
    elems, mats = elems0[export_elems][0], elems0[export_elems][1]
    ndgrp = label_nodes(node_groups, nodes)

    data = {
        'mat_id': ('c', 'int', mats),
        'node_groups': ('p', 'int', ndgrp)
    }

    # if only_elems_by_nodes is not None:
    #     export_elems = only_elems_by_nodes[0]
    #     elems, mats = elems0[export_elems][0], elems0[export_elems][1]
    #     ndgrp = label_nodes([(1, only_elems_by_nodes[1])], nodes)
    #     aux = ndgrp[elems].sum(axis=1)
    #     eidxs = nm.where(aux == int(export_elems[-1]))[0]
    #     nidxs = nm.where(ndgrp == 1)[0]
    #     data['node_groups'] = ('p', 'int', ndgrp[nidxs] * 0)
    #     data['mat_id'] = ('c', 'int', mats[eidxs])
    #     remap = -nm.ones((nodes.shape[0],), dtype=nm.int32)
    #     remap[nidxs] = nm.arange(len(nidxs))
    #     nodes = nodes[nidxs, :]
    #     elems = remap[elems[eidxs, :]]

    if shift is not None:
        nodes += shift

    if scale is not None:
        nodes *= scale

    vtk_write(filename_base + '.vtk', nodes, elems, export_elems, data)


def repeater(filename_in, filename_out, grid, size_x, tol=1e-9):
    nodes, elems, elem_type, data = vtk_read(filename_in, ret_pc_data=False)

    for idim, igrid in enumerate(grid):
        if igrid <= 0:
            raise ValueError('Incorrect numer of repetition! (%s)' % grid)

        idir = nm.eye(3)[idim]
        nnodes = nodes.shape[0]

        nodes = nm.vstack(nodes + idir * ii for ii in range(igrid))
        elems = nm.vstack(elems + nnodes * ii for ii in range(igrid))
        repdata = {}
        for k, v in data.items():
            if len(v[2].shape) > 1:
                repdata[k] = (v[0], v[1], nm.vstack([v[2]] * igrid))
            else:
                repdata[k] = (v[0], v[1], nm.hstack([v[2]] * igrid))

        nodes, elems, data = merge_nodes(nodes, elems, repdata, tol=tol)

    scale = float(size_x) / nm.max(nodes[:, 0])

    vtk_write(filename_out, nodes * scale, elems, elem_type, data)
