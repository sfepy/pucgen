import os
import sys
import numpy as nm
from scipy.spatial import cKDTree
import meshio

# def mirror_mesh(nodes, elems, data, c0=0, dim=0):
#     remap_tetra = nm.array([0, 2, 1, 3])
#     nnodes = nodes.shape[0]
#     melems = nm.vstack([elems, (elems + nnodes)[:, remap_tetra]])
#     nodes2 = nodes.copy()
#     nodes2[:, dim] = 2 * c0 - nodes2[:, dim]
#     mnodes = nm.vstack([nodes, nodes2])

#     mdata = {}
#     for k, v in data.items():
#         if len(v[2].shape) > 1:
#             mdata[k] = (v[0], v[1], nm.vstack([v[2], v[2]]))
#         else:
#             mdata[k] = (v[0], v[1], nm.hstack([v[2], v[2]]))

#     return mnodes, melems, mdata

def find_master_slave_nodes(coors, tol=1e-9):
    """
    Find matching coordinates using cKDTree.

    Parameters
    ----------
    coors: numpy.ndarray
        Nodal coordinates
    tol: float
        Tolerance, coordinates A and B match when |A - B| <= tol

    Returns
    -------
    out: numpy.ndarray
        Pairs of matching coordinates
    """
    tr = cKDTree(coors)
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
                out[idx0, :] = cols
                idx0 += 1
            else:
                idx1 = idx0 + nc - 1
                out[idx0:idx1, 0] = cols[0]
                out[idx0:idx1, 1] = cols[1:]
                idx0 = idx1

    return out[:idx0, :]


def merge_nodes(mesh, tol=1e-9):
    """
    Merge mesh nodes.

    Parameters
    ----------
    mesh: meshio.Mesh
        FE mesh
    tol: float
        Tolerance, nodes A and B are merged if |A - B| <= tol
    """
    if '_not_merge' in mesh.point_data:
        idxs = nm.where(nm.logical_not(mesh.point_data['_not_merge']))[0]
        ms_tab = find_master_slave_nodes(mesh.points[idxs], tol=tol)
        ms_tab = idxs[ms_tab]
    else:
        ms_tab = find_master_slave_nodes(mesh.points, tol=tol)

    remap = nm.ones((mesh.points.shape[0],), dtype=nm.int64)
    remap[ms_tab[:, 1]] = -1
    ndidxs = nm.where(remap > 0)[0]
    remap[ndidxs] = nm.arange(len(ndidxs))
    mesh.points = mesh.points[ndidxs, :]
    remap[ms_tab[:, 1]] = remap[ms_tab[:, 0]]

    pdata = mesh.point_data
    if pdata is not None:
        for k in pdata.keys():
            pdata[k] = pdata[k][ndidxs, ...]

    for cg in mesh.cells:
        cg.data = remap[cg.data]


def merge_cell_groups(mesh):
    """
    Merge multiple cell groups into the one.
    """
    cells = mesh.cells
    cdata = mesh.cell_data

    ctypes = list({cg.type for cg in cells})
    gcells = {ct: [cg for cg in cells if cg.type == ct] for ct in ctypes}
    ncells = [meshio.CellBlock(k, nm.vstack([c.data for c in v]))
              for k, v in gcells.items()]

    gcidxs = {ct: [k for k, cg in enumerate(cells) if cg.type == ct]
              for ct in ctypes}

    ncdata = {k: [nm.hstack([v[idx] for idx in gcidxs[ct]]) for ct in ctypes]
              for k, v in cdata.items()}

    mesh.cells = ncells
    mesh.cell_data = ncdata


def meshio_read(filename):
    mesh = meshio.read(filename)

    for k in list(mesh.point_data.keys()):
        if k == 'gmsh:dim_tags':
            key = 'node_groups'
            val = mesh.point_data[k]
            if nm.issubdtype(val.dtype, nm.floating):
                mesh.point_data[key] = nm.asarray(val, dtype=nm.float64)
            else:
                mesh.point_data[key] = nm.asarray(val, dtype=nm.int64)

        del mesh.point_data[k]

    mesh.cell_sets = {}

    for k in list(mesh.cell_data.keys()):
        if k == 'gmsh:physical':
            val = []
            for ival in mesh.cell_data[k]:
                if nm.issubdtype(ival.dtype, nm.floating):
                    val.append(nm.asarray(ival, dtype=nm.float64))
                else:
                    val.append(nm.asarray(ival, dtype=nm.int64))

            mesh.cell_data['mat_id'] = val

        del mesh.cell_data[k]

    merge_cell_groups(mesh)

    return mesh


def call_gmsh(filename_base, dim=3, export_elems='tetra',
              node_groups=None, shift=None, scale=None):

    os.system(f'gmsh -{dim} {filename_base}.geo -o {filename_base}.msh')

    mesh = meshio_read(f'{filename_base}.msh')

    if shift is not None:
        mesh.points += shift

    if scale is not None:
        mesh.points *= scale

    fname = f'{filename_base}.vtk'
    print(f'writing mesh to {fname}')
    mesh.write(fname, binary=False)


def repeat_cell(filename_in, filename_out, grid, size_x, tol=1e-9):
    mesh = meshio_read(filename_in)

    nodes = mesh.points
    pdata = mesh.point_data
    cdata = mesh.cell_data
    cell_size = nm.max(nodes, axis=0) - nm.min(nodes, axis=0)

    if grid is not None:
        for idim, igrid in enumerate(grid):
            if igrid <= 0:
                raise ValueError('Incorrect number of repetition! (%s)' % grid)

            nnodes = nodes.shape[0]
            idir = nm.eye(3)[idim] * cell_size

            # duplicate nodes
            nodes = nm.vstack([nodes + idir * ii for ii in range(igrid)])
            # duplicate cells and cell data
            for ig, cg in enumerate(mesh.cells):
                cg.data = nm.vstack([cg.data + nnodes * ii
                                     for ii in range(igrid)])

                for k, v in cdata.items():
                    v[ig] = nm.vstack([v[ig]] * igrid)

            # point data
            for k in pdata.keys():
                pdata[k] = nm.vstack([pdata[k]] * igrid)

        mesh.points = nodes
        merge_nodes(mesh, tol=tol)

    scale = float(size_x) / nm.max(nodes[:, 0]) if size_x is not None else 1.0

    mesh.points *= scale

    mesh.write(filename_out, binary=False)
