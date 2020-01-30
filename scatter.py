import numpy as nm
from pucgen import BaseComponent
from scipy.spatial import cKDTree


def test_particles():
    c = nm.array([[0.04703357, 0.42511719, 0.83717491],
                  [0.41861516, 0.75399372, 0.74505737],
                  [0.70565324, 0.01814634, 0.84687874],
                  [0.17751005, 0.11797473, 0.24348222],
                  [0.27154594, 0.91874972, 0.71040134],
                  [0.23774401, 0.62737545, 0.35331869],
                  [0.76373506, 0.65644267, 0.69213158],
                  [0.56160528, 0.4121705 , 0.22671179],
                  [0.06294404, 0.93708711, 0.21205263],
                  [0.85558864, 0.03838302, 0.83416308]])

    v = nm.array([[0.23604426, 0.32497936, 0.30988612], 
                  [0.46619628, 0.21399325, 0.22736252],
                  [0.08906827, 0.44527641, 0.99987815],
                  [0.24763170, 0.75352237, 0.58022644],
                  [0.89406559, 0.16519888, 0.86312528],
                  [0.32145155, 0.74658101, 0.10838125],
                  [0.26692915, 0.80306648, 0.88268946],
                  [0.88512786, 0.92227027, 0.30765224],
                  [0.48342124, 0.05867704, 0.96727921],
                  [0.35496475, 0.79872526, 0.83350329]]) - nm.array([0.5, 0.5])

    return c, v


class BaseParticles(BaseComponent):
    """The base for the particles."""

    def __init__(self, parameters, central_point, direction, grow_rate,
                 velocity, el_size, mat_id):
        """Init parameters of the particles.
    
        Parameters
        ----------
        parameter: float or array
            The geometrical parameters of the objects.
        central_point: array
            The coordinates of the object centers: [x, y, z].
        direction: array
            The object directions.
        grow_rate: array
            !!!
        velocity: array
            !!!
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """    
        super(BaseParticles, self).__init__(mat_id=mat_id)

        self.params.update({
            'parameters': parameters,
            'direction': nm.asarray(direction) / nm.linalg.norm(direction),
            'central_point': nm.asarray(central_point, dtype=nm.float64),
            'grow_rate': nm.asarray(central_point, dtype=nm.float64),
            'velocity': nm.asarray(central_point, dtype=nm.float64),
            'el_size': el_size,
        })

    def get_coors(self, dt=0):
        return self.get('central_point') + dt * self.get('velocity')

    def apply_periodicty(self, size=None):
        if size is None:
            size = self.get('size')
        
        ccoors = self.get('central_point')

        for idim in range(ccoors.shape[1]):
            d = size[idim]
            idxs1 = nm.where(ccoors[:, idim] < 0.0)[0]
            idxs2 = nm.where(ccoors[:, idim] > d)[0]
            ccoors[idxs1, idim] += d
            ccoors[idxs2, idim] -= d

    def get_ball_repr(self):
        pass

    def get_bbox(self, est_dt):
        pass

    def deactivate(self):
        self.active = False

    def activate(self):
        self.active = True


class SphericalParticles(BaseParticles):
    name = 'Sphere Particles'

    def __init__(self, central_point, grow_rate, velocity,
                 el_size=0.5, mat_id=2):
        """Init parameters of the particles.
    
        Parameters
        ----------
        central_point: array
            The coordinates of the sphere centers.
        direction: array
            The directional vectors.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        el_size: float
            The "inner" element size factor: in_el_size = el_size * el_size_base.
        """
        radius = nm.zeros_like(grow_rate)
        super(SphericalParticles, self).__init__(parameters=[radius],
                                                 central_point=central_point,
                                                 direction=None,
                                                 grow_rate=grow_rate,
                                                 velocity=velocity,
                                                 el_size=el_size,
                                                 mat_id=mat_id)

    def get_bbox(self, dt=0):
        np, dim = self.get('central_point').shape
        gr = 0 if nm.abs(dt) < 1e-16 else dt * self.get('grow')
        rad = nm.ones((1, dim)) * (self.get('radius') + gr)
        bbox = nm.zeros((np, dim, 2), dtype=nm.float64)
        coors = self.get_coors(dt=dt)
        bbox[..., 0] = coors - rad
        bbox[..., 1] = coors + rad

        return bbox


def _get_dual_particles(ccoors, bcoors, size):
    np, dim = ccoors.shape
    dual_flag = -1 * nm.ones((np,), dtype=nm.int32)

    ccoors1 = ccoors.copy()
    bcoors1 = bcoors.copy()
    dual_flag1 = dual_flag.copy()

    for idim in range(dim):
        shift = nm.eye(dim)[idim] * size[idim]
        dual_ccoors = []
        dual_bcoors = []
        dual_dual_flag = []

        idxs = nm.where(bcoors1[:, idim, 1] > size[idim])[0]
        dual_ccoors.append(ccoors1[idxs] + shift)
        dual_bcoors.append(bcoors1[idxs] + shift)
        idxs2 = nm.where(dual_flag1[idxs] >= 0)[0]
        idxs[idxs2] = dual_flag1[idxs]
        dual_dual_flag.append(idxs)

        idxs = nm.where(bcoors1[:, idim, 0] < 0)[0]
        dual_ccoors.append(ccoors1[idxs] - shift)
        dual_bcoors.append(bcoors1[idxs] - shift)
        idxs2 = nm.where(dual_flag1[idxs] >= 0)[0]
        idxs[idxs2] = dual_flag1[idxs]
        dual_dual_flag.append(idxs)

        ccoors1 = nm.vstack([ccoors1] + dual_ccoors)
        bcoors1 = nm.vstack([bcoors1] + dual_bcoors)
        dual_flag1 = nm.vstack([dual_flag1] + dual_dual_flag)

        return ccoors1, bcoors1, dual_flag1


def get_dual_particles(particles, size, est_dt=0):
    ccoors = particles.get('central_point')
    bbox = particles.get_bbox(est_dt)
    nccoors, nbbox, dual_flag = _get_dual_particles(ccoors, bbox, size)

    return nccoors, nbbox, dual_flag


def apply_pbc(coors, rve_box):
    out = coors.copy()
    for ii in range(2):
        rg = rve_box[ii]
        idxs = nm.where(coors[:,ii] < 0.0)[0]
        out[idxs,ii] += rg
        idxs = nm.where(coors[:,ii] > rg)[0]
        out[idxs,ii] -= rg

    return out


def collision_candidates(ccoors, radius):
    tr = cKDTree(ccoors)
    r_max = nm.max(radius)
    out = tr.query_ball_tree(tr, r=2*r_max)

    return out


def main(size=(1, 1, 1),
         volume_faction=0.6,
         nparticles=50,
         particle_spacing=0.05,
         grow_rate=0.02):

    genparts = 'rand'
    max_velocity = 1
    nfinal_iter = 50
    est_dt = 0.1

    dim = 3
    #  !!!!!
    volume_faction = volume_faction / (1.0 - particle_spacing)**2

    if genparts == 'rand':
        coors = nm.random.rand(nparticles, dim) * nm.array(size)
        velocity = (nm.random.rand(nparticles, dim) - nm.ones(dim) * 0.5)\
            * nm.array(size) * max_velocity

    elif genparts == 'test':
        coors, velocity = test_particles()

    grow = nm.ones((nparticles,), dtype=nm.float64) * nm.array(grow_rate)

    particles = SphericalParticles(coors, grow, velocity, el_size=0.5, mat_id=2)

    istep = 1
    while True:
        particles.apply_periodicty(size)
        coors, bbox, dual_flag = get_dual_particles(particles, size, est_dt)
        cand_idx = collision_candidates(coors, particles.get('radius'))

        if len(cand_idx) > 0:
            pass
        else:
            print('no colision, dt = %e' % est_dt)

        est_dt = nm.min(particles.get('radius')) / nm.max(particles.get('velocity'))

        cp_idxs = cps_idxs[cp,:]
        ecoors += dt * evelocity
        eradius += dt * egrow
        new_vel = collision_event(ecoors, evelocity, egrow, cp_idxs)
        if mstab.shape[0] > np:
            coors[:,:] = apply_pbc(ecoors[:np,:], size)
            radius[:] = eradius[:np]
            velocity[mstab[cp_idxs],:] = new_vel

        else:
            velocity[cp_idxs,:] = new_vel

        vf = volume_fraction(radius, size)
        print('%d: ncps = %d, dt = %f, vf = %f' % (istep, len(cps_idxs), dt, vf))

        istep += 1
        est_dt = 10 * dt

        if vf >= volume_faction:
            grow.fill(0.0)
            nfinal_iter -= 1
            if nfinal_iter <= 0:
                break

    radius *= (1.0 - particle_spacing)
    ecoors, eradius, _ = dual_particles(coors, radius, size)
    flag = nm.zeros((ecoors.shape[0],), dtype=nm.int8)
    flag[np:] = 3
    print('vf = %f' % volume_fraction(radius, size))
    # draw_state(ecoors, eradius, flag=flag, rve=size)

    savemat('particles.mat', {'coors': ecoors,
                              'radius': eradius,
                              'primar_dual': mstab})

if __name__ == "__main__":
    n = 5
    grate = nm.ones((n,), dtype=nm.float64) * 0.02
    main(nparticles=n, grow_rate=grate, volume_faction=0.2)
