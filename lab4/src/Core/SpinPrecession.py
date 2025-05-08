import kwant
import numpy as np
import matplotlib.pyplot as plt
from . import Utils as utl

# class implenets simple spin orbit coupling
# It is designed to enter values in default units (nm and T)
# transformations are done inside class
class SpinSystem():
    u = utl.Utils()

    # TODO: check default alpha value
    def __init__(self, *, dx = 4, L = int(500), W = int(25), m = 0.014
                        , Bx = 0, By = 0, Bz = 0
                        , alpha = 0.000):
        self.dx = self.u.nm2au(dx)
        self.L = int(L)
        self.W = int(W)
        self.m = m
        self.Bx = self.u.T2au(Bx)
        self.By = self.u.T2au(By)
        self.Bz = self.u.T2au(Bz)
        self.alpha = self.u.nm2au(self.u.eV2au(alpha))

        self.precalc_mag_term = 0.5 * self.u.bohr_magneton_au * self.u.lande_g * \
                                       (self.Bx * self.u.sigma_x + 
                                        self.By * self.u.sigma_y + 
                                        self.Bz * self.u.sigma_z)

        self.precalc_t = 1.0/(2.0*self.m * self.dx * self.dx)
        self.precalc_tso =  self.alpha / (2 * self.dx)


    def make_system(self): 
        def onsite(site):

            return 4*self.precalc_t * np.identity(2) + self.precalc_mag_term

        def hopping_x(sitei, sitej):

            return -self.precalc_t * np.identity(2) + 1j * self.precalc_tso * self.u.sigma_y
        
        def hopping_y(sitei, sitej):

             return -self.precalc_t * np.identity(2) - 1j * self.precalc_tso * self.u.sigma_x

        sys = kwant.Builder()  
        lat = kwant.lattice.square(self.dx, norbs=2)
        sys[(lat(i,j) for i in range(self.L) for j in range(0,self.W))]= onsite
        sys[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping_x
        sys[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping_y

        sigma_law = np.matrix([[1, 0], [0, 2]])
    
        #attach the left contact to the system
        lead_left = kwant.Builder(kwant.TranslationalSymmetry((-self.dx, 0)),
                                  conservation_law = sigma_law)    
        lead_left[(lat(0,j) for j in range(0,self.W))]= onsite
        lead_left[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping_x
        lead_left[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping_y
        sys.attach_lead(lead_left)
    
        #attach the right contact to the system
        lead_right = kwant.Builder(kwant.TranslationalSymmetry((self.dx, 0))
                                   , conservation_law = sigma_law)    
        lead_right[(lat(0,j) for j in range(0,self.W))]= onsite
        lead_right[(kwant.builder.HoppingKind((-1,0), lat, lat))] = hopping_x
        lead_right[(kwant.builder.HoppingKind((0,-1), lat, lat))] = hopping_y
        sys.attach_lead(lead_right)
    
        #finalize the system
        sys = sys.finalized()
        return sys

    def dispersion(self, nr_lead, k_max, nk):
        sys = self.make_system()
        momenta = np.linspace(-k_max*self.dx,k_max*self.dx,nk)
        bands=kwant.physics.Bands(sys.leads[nr_lead])
        energies=[bands(k) for k in momenta]
        return (momenta/self.dx)*self.u.nm2au(1.0), energies

###########################################################################
##  Define various functions calculating the basic physical properties   ##
###########################################################################

#calculates the reflection and transmission coefficient
def transmission_reflection(nw, E):
    E=eV2au(E)
    sys=make_system(nw)
    smatrix=kwant.smatrix(sys,E)
    r=smatrix.transmission(0,0)
    t=smatrix.transmission(1,0)
    return r, t

#calculates the transmission coefficient
def transmission(nw, E):
    E=eV2au(E)
    sys=make_system(nw)
    smatrix=kwant.smatrix(sys,E)
    t=smatrix.transmission(1,0)
    return t

#calculates the conductance - the Landauer formula is used
def conductance(nw, Emax, ne):
    energies=np.linspace(0,Emax,ne)
    cond=[transmission(nw, E) for E in energies]
    return energies, cond

#plots the wave function of an electron with energy E incident in the contact nr_lead
def wave_function(nw, E, nr_lead, *, ax = None):
    E=eV2au(E)
    sys=make_system(nw)
    wave=kwant.wave_function(sys, E)
    density=(abs(wave(nr_lead))**2).sum(axis=0)
    kwant.plotter.map(sys,density, ax = ax, dpi = 300)

#fplots the dos of an electron with energy E
def dos(nw, E):
    E=eV2au(E)
    sys=make_system(nw)
    dos=kwant.ldos(sys, E)
    f = kwant.plotter.map(sys,dos)
    return f

#plots the current of an electron with energy E incident in the contact nr_lead in the state nr_mod
def current(nw, E, nr_lead, nr_mod, *, ax = None):
    E=eV2au(E)
    sys=make_system(nw)
    current = kwant.operator.Current(sys).bind()
    psi=kwant.wave_function(sys, E)(nr_lead)
    curr=current(psi[nr_mod])
    kwant.plotter.current(sys,curr, ax = ax)
