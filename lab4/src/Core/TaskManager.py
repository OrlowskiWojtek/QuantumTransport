from . import SpinPrecession as spin_system
import kwant as kw

class TaskManager():
    def task1(self):
        sp = spin_system.SpinSystem()
        sys = sp.make_system()
        kw.plot(sys, site_color=lambda site: sys.hamiltonian(site,site), fig_size=(10,5), colorbar=False, show=False, num_lead_cells=2);

        


