from pucgen import (PUC, SphericalInclusion, CylindricalInclusion,
                    BoxInclusion, EllipsoidalInclusion)

puc = PUC()

puc.add(SphericalInclusion(dimension=0.2, central_point=(-0.2, 0.2, -0.2),
                           el_size=0.04, mat_id=2))

puc.add(CylindricalInclusion(dimension=(0.08, 0.8),
                             central_point=(0.1, -0.1, -0.2),
                             direction=(1, 1, 0.5), el_size=0.03, mat_id=3))

puc.add(BoxInclusion(dimension=(0.2, 0.25, 0.3),
                     central_point=(0.2, 0.2, 0.2), el_size=0.05, mat_id=4))

puc.add(EllipsoidalInclusion(dimension=(0.4, 0.15, 0.1),
                             central_point=(-0.15, -0.15, 0.25),
                             direction=(1, -1, 0.5), el_size=0.03, mat_id=5))

print(puc)

puc('example3.vtk')
