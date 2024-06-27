from pucgen import PUC, BaseCell, SphericalInclusion, CylindricalChannel

puc = PUC(cell_mat_id=1)

puc.add(BaseCell(dimension=(1, 1, 1), el_size=0.1, mat_id=5))
puc.add(SphericalInclusion(dimension=0.3, central_point=(0, 0, 0),
                           el_size=0.05, mat_id=2))
puc.add(CylindricalChannel(dimension=0.1, central_point=(0, 0, 0),
                           direction='x', el_size=0.05, mat_id=2))
puc.add(CylindricalChannel(dimension=0.15, central_point=(0, 0, 0),
                           direction='y', el_size=0.05, mat_id=2))
puc.add(CylindricalChannel(dimension=0.2, central_point=(0, 0, 0),
                           direction='z', el_size=0.05, mat_id=2))

print(puc)

puc('example1.vtk')
