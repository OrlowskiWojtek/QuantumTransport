import kwant
import numpy as np
from . import utils as utl

class Graphene:
    u = utl.Utils()

    def __init__(self, 
                 a_nm=0.25, # graphene primitive vector
                 sf=8.,    # scaling factor
                 t_eV=-3.0, # hopping pottential
                 W = 100,   # width of cell  (y dim)
                 L = 200,   # length of cell (x dim)
                 B = 0.,     # magnetic field
                 V_np = 1.,  # potential amplitude
                 d = 1.      # smothness of potential
                 ):    

        self.a0 = self.u.nm2au(a_nm) * sf
        self.t0 = self.u.eV2au(t_eV) / sf

        self.L = self.u.nm2au(L)
        self.W = self.u.nm2au(W)

        self.x_min = -self.L / 2.
        self.x_max = self.L / 2.
        self.y_min = -self.W / 2.
        self.y_max = self.W / 2.

        self.B = self.u.T2au(B)

        sin_30 = np.sin(np.pi / 6)
        cos_30 = np.cos(np.pi / 6)
        self.d = self.u.nm2au(d)
        self.V_np = self.u.eV2au(V_np)
        self.graphene = kwant.lattice.general(
            [(0, self.a0), (cos_30 * self.a0, sin_30 * self.a0)],
            [(0, 0), (self.a0 / np.sqrt(3), 0)],
            norbs=1,
        )
        self.a_sub, self.b_sub = self.graphene.sublattices

    def make_system(self):
        def pote_x(x):
            lt_term = np.tanh((x + self.L/8)/self.d)
            rt_term = np.tanh((x - self.L/8)/self.d)

            return self.V_np * ( lt_term - rt_term ) / 2

        def rect(pos):
            x, y = pos

            inside_x = self.x_min <= x <= self.x_max 
            inside_y = self.y_min <= y <= self.y_max

            return inside_x and inside_y

        def lead_shape(pos):
            x, y = pos
            return self.y_min <= y <= self.y_max

        def potential(site):
            x, y = site.pos

            return pote_x(x)

        def nn_hopping(site1, site2):
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            pereir = -self.B * (y1 + y2) * (x2 - x1) / 2

            return self.t0 * np.exp(1j * pereir)

        sys = kwant.Builder()

        graphene = self.graphene
        sys[graphene.shape(rect, (0, 0))] = potential
        sys[graphene.neighbors()] = nn_hopping

        syml = kwant.TranslationalSymmetry([-np.sqrt(3) * self.a0, 0])
        symr = kwant.TranslationalSymmetry([np.sqrt(3) * self.a0, 0])

        leadl = kwant.Builder(syml)
        leadl[graphene.shape(lead_shape, (0, 0))] = pote_x(self.x_min)
        leadl[graphene.neighbors()] = nn_hopping

        leadr = kwant.Builder(symr)
        leadr[graphene.shape(lead_shape, (0, 0))] = pote_x(self.x_max)
        leadr[graphene.neighbors()] = nn_hopping

        sys.attach_lead(leadl)
        sys.attach_lead(leadr)

        return sys.finalized()

    def get_bands(self, nr_lead=0):
        sys = self.make_system()
        dx = np.sqrt(3) * self.a0
        k_max = np.pi / dx
        momenta = np.linspace(-k_max * dx, k_max * dx, 200)
        bands = kwant.physics.Bands(sys.leads[nr_lead])
        energies = [bands(k) for k in momenta]
        return (momenta / dx), energies

        #calculates the transmission coefficient - input in eV
    def transmission(self, E):
        E = self.u.eV2au(E)
        sys= self.make_system()
        smatrix=kwant.smatrix(sys,E)
        t=smatrix.transmission(1,0) # total transmission from 0 to 1 (or reversed)
        return t

        #calculates the transmission coefficient - input in a. u.
    def _transmission(self, E):
        sys= self.make_system()
        smatrix=kwant.smatrix(sys,E)
        t=smatrix.transmission(1,0) # total transmission from 0 to 1 (or reversed)
        return t

    #calculates the conductance - the Landauer formula is used
    def conductance(self, Emax, Emin, ne):
        energies=np.linspace(Emin, Emax,ne)
        cond=[self._transmission(self.u.eV2au(E)) for E in energies]
        return energies, cond

    def wave_function(self, E, nr_lead, *, ax = None):
        E=self.u.eV2au(E)
        sys=self.make_system()
        wave=kwant.wave_function(sys, E)
        density=(abs(wave(nr_lead))**2).sum(axis=0)
        kwant.plotter.map(sys,density, ax = ax, dpi = 300)

    def current(self, E, nr_lead, nr_mod, *, ax = None):
        E = self.u.eV2au(E)
        sys = self.make_system()
        current = kwant.operator.Current(sys).bind()
        psi = kwant.wave_function(sys, E)(nr_lead)
        curr = current(psi[nr_mod])
        kwant.plotter.current(sys, curr, ax = ax)



