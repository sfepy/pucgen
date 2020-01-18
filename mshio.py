import sys
import os
import numpy as  nm

def msh_read(filename):
    fd = open(filename, 'r')
    phys_names = {}
    while 1:
        line = fd.readline()
        if len(line) == 0:
            break
        line = line.strip().split()
        if len(line) == 0:
            continue

        if line[0] == '$PhysicalNames':
            n = int(fd.readline())
            for ii in range(n):
                aux = fd.readline().split()
                phys_names[int(aux[1])] =\
                    aux[2].replace('"', ''), int(aux[0]), 0

        elif line[0] == '$Nodes':
            n = int(fd.readline())
            nodes = []
            for ii in range(n):
                nodes.append([float(jj) for jj in fd.readline().split()[1:]])

            nodes = nm.array(nodes, dtype=nm.float64)

        elif line[0] == '$Elements':
            n = int(fd.readline())
            elems0 = {}
            mats = {}
            for ii in range(n):
                line = [int(jj) for jj in fd.readline().split()[1:]]
                ckey, ntag, mat = line[0:3]
                if ckey in elems0:
                    elems0[ckey].append(line[(2 + ntag):])
                    mats[ckey].append(mat)
                else:
                    elems0[ckey] = [line[(2 + ntag):]]
                    mats[ckey] = [mat]

            elems = {}
            etab = {2: '2_3', 3: '2_4', 4: '3_4', 5: '3_8', 6: '3_6'}
            for ii in elems0.keys():
                if ii not in etab:
                    continue
                elems[etab[ii]] = nm.asarray(elems0[ii]) - 1, nm.asarray(mats[ii])

    fd.close()
    return nodes, elems, phys_names


def msh_write(filename, nodes, elems, phys_names):
    fd = open(filename, 'w')

    fd.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')

    fd.write('$PhysicalNames\n%d\n' % len(phys_names))
    for k, v in phys_names.items():
        fd.write('%d %d "%s"\n' % (v[1], k, v[0]))
    fd.write('$EndPhysicalNames\n')

    fd.write('$Nodes\n%d\n' % nodes.shape[0])
    for jj, ii in enumerate(nodes):
        fd.write('%d %f %f %f\n' % ((jj + 1,) + tuple(ii)))
    fd.write('$EndNodes\n')

    nel = 0
    for k, v in elems.items():
        nel += v[0].shape[0]
    fd.write('$Elements\n%d\n' % nel)
    iel = 1
    etab = {'2_3': 2, '2_4': 3, '3_4': 4, '3_8': 5, '3_6': 6}
    for k, v in elems.items():
        for conn, mat in zip(v[0], v[1]):
            sconn = ' '.join([str(ii + 1) for ii in conn])
            fd.write('%d %d %d %d %d %s\n' % (iel, etab[k], 2, mat, 1, sconn))
            iel += 1
    fd.write('$EndElements\n')

    fd.close()
