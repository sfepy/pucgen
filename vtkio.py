import numpy as nm
import sys

cell_type_table = {'2_3': 5, '2_4': 9,
                   '3_2': 3, '3_4': 10, '3_8': 12, '3_6': 13}
cell_type_table2 =\
    {v: int(k.split('_')[1]) for k, v in cell_type_table.iteritems()}
cell_type_table3 =\
    {v: k for k, v in cell_type_table.iteritems()}

_float_format = '0.6e'
float_format = '{:%s}' % _float_format
int_format = '{:d}'

dtype_table = {'float': nm.float64, 'double': nm.float64, 'int': nm.int32}
dtype_str_table = {'float': float_format, 'double': float_format,
                   'int': int_format}

vtk_head = """# vtk DataFile Version 2.6
output file
ASCII
DATASET UNSTRUCTURED_GRID

"""


def write_data(out, label, dtype, data):
    if type(data) is nm.ndarray:
        data = data.squeeze()

    dd = 1 if len(data.shape) <= 1 else nm.prod(data.shape[1:])
    dformat = dtype_str_table[dtype]
    if dd == 1:  # scalar
        out.write('\nSCALARS %s %s 1\n' % (label, dtype))
        out.write('LOOKUP_TABLE default\n')
        out.write(''.join([dformat.format(ii) + '\n' for ii in data]))

    elif dd == 3:  # vector
        dvformat = '{0:%s} {1:%s} {2:%s}\n' % tuple([_float_format] * 3)
        out.write('\nVECTORS %s %s\n' % (label, dtype))
        out.write(''.join([dvformat.format(*ii) for ii in data]))

    elif dd == 9:  # tensor
        aux = tuple([_float_format] * 3)
        dtformat = '{0:%s} {1:%s} {2:%s}\n' % aux +\
                   '{3:%s} {4:%s} {5:%s}\n' % aux +\
                   '{6:%s} {7:%s} {8:%s}\n' % aux
        out.write('\nTENSORS %s %s\n' % (label, dtype))
        out.write(''.join([dtformat.format(*ii.flatten()) for ii in data]))


def vtk_write(filename, points, cells, cell_type, data={}):
    """
    data: {'mat_id': ('c', 'int', mat_data),
           'velocity': ('p', 'float', vel_data),
           'node_group': ('p', 'int', ngrp_data)}
    """

    if type(points) is not nm.ndarray:
        points = nm.array(points)

    if type(cells) is not nm.ndarray:
        cells = nm.array(cells)

    npoints, dim = points.shape
    ncells = cells.shape[0]
    datak = data.keys()
    datak.sort()
    datac = [(k,) + data[k][1:] for k in datak if data[k][0] == 'c']
    datap = [(k,) + data[k][1:] for k in datak if data[k][0] == 'p']

    out = open(filename, 'w')
    out.write(vtk_head)

    out.write('POINTS %d float\n' % (npoints, ))
    for ii in points:
        out.write(' '.join([float_format.format(jj) for jj in ii]) + '\n')

    ntot = nm.sum([len(ii) + 1 for ii in cells])
    out.write('\nCELLS %d %d\n' % (ncells, ntot))
    for ii in cells:
        conn = ' '.join([int_format.format(jj) for jj in ii])
        out.write('%d %s\n' % (len(ii), conn))

    out.write('\nCELL_TYPES %d\n' % ncells)
    if type(cell_type) is str:
        out.write(('%d\n' % cell_type_table[cell_type]) * ncells)
    else:
        for ii in cell_type:
            out.write('%d\n' % ii)

    if len(datac) > 0:
        out.write('\nCELL_DATA %d\n' % ncells)
        for ii in datac:
            write_data(out, *ii)

    if len(datap) > 0:
        out.write('\nPOINT_DATA %d\n' % npoints)
        for ii in datap:
            write_data(out, *ii)

    out.close()


def vtk_read(filename, ret_pc_data=False):
    def read_data(f, n, dtype, dim=1):
        ntot = n * dim
        nr = 0
        out = []
        while nr < ntot:
            line = f.readline()
            if len(line) == 0:
                break
            line = line.split()
            if len(line) == 0:
                continue
            out.append(nm.array(line, dtype=dtype))
            nr += len(line)

        return nm.hstack(out).reshape((n, dim))

    f = open(filename, 'r')

    data = {}
    datakey = None
    datasize = 0

    while True:
        line = f.readline()
        if len(line) == 0:
            break

        line = line.split()
        if len(line) == 0:
            continue

        tag = line[0].strip()
        if tag == 'POINTS':
            npoints = int(line[1])
            points = read_data(f, npoints, nm.float64, dim=3)
        elif tag == 'CELLS':
            ncells = int(line[1])
            cells = read_data(f, int(line[2]), nm.int32)
        elif tag == 'CELL_TYPES':
            cell_types = read_data(f, int(line[1]), nm.int32)
        elif tag == 'CELL_DATA':
            datakey, datasize = 'c', ncells
            assert(ncells == int(line[1]))
        elif tag == 'POINT_DATA':
            datakey, datasize = 'p', npoints
            assert(npoints == int(line[1]))
        elif tag == 'SCALARS':
            label, dtype = line[1].strip(), line[2].strip()
            _ = f.readline()
            val = read_data(f, datasize, dtype_table[dtype]).squeeze()
            data[label] = (datakey, dtype, val)
        elif tag == 'VECTORS':
            label, dtype = line[1].strip(), line[2].strip()
            val = read_data(f, datasize, dtype_table[dtype], dim=3)
            data[label] = (datakey, dtype, val.reshape(datasize, 3))
        elif tag == 'TENSORS':
            label, dtype = line[1].strip(), line[2].strip()
            val = read_data(f, datasize, dtype_table[dtype], dim=9)
            data[label] = (datakey, dtype, val.reshape(datasize, 9))
        elif tag == 'FIELD':
            nfdata = int(line[2])
            for ii in range(nfdata):
                line = []
                while not len(line) == 4:
                    line = f.readline().split()
                label, dtype = line[0].strip(), line[3].strip()
                val = read_data(f, datasize, dtype_table[dtype]).squeeze()
                data[label] = (datakey, dtype, val)
        elif tag == 'POLYGONS':  # only for triangles
            ncells = int(line[1])
            cells = read_data(f, int(line[2]), nm.int32)
            cell_types = nm.ones((ncells, 1), dtype=nm.int32) * 5

    f.close()
    ctype = cell_types[0, 0]
    cells = cells.reshape((ncells, cell_type_table2[ctype] + 1))[:, 1:]

    if ret_pc_data:
        out = {'p': {}, 'c': {}}
        for label, (pc, _, val) in data.iteritems():
            out[pc][label] = val
        return points, cells, cell_type_table3[ctype], out['p'], out['c']
    else:
        return points, cells, cell_type_table3[ctype], data


def main(argv):
    filename = argv[0]
    p, c, ct, d = vtk_read(filename)
    vtk_write(filename + '.vtk', p, c, ct, d)

if __name__ == '__main__':
    main(sys.argv[1:])
