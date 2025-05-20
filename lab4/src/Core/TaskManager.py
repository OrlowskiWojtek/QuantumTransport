from . import SpinPrecession as spin_system
import kwant as kw
import numpy as np
import Core.Utils as utl

from . import GeneralPlotter as plot_utils

u = utl.Utils()
general_plotter = plot_utils.Plotter()

class TaskManager():
    def task1(self, *, do_plot = False):

        #sp = spin_system.SpinSystem()
        #sys = sp.make_system()
        #kw.plot(sys, show = do_plot, file = "../plots/ex1_system.pdf");

        #moments, enes = sp.dispersion(0, .05, 200)
        #general_plotter.plot_dispersion(moments, enes, filename = "../plots/ex1_disp.pdf", show = False)

        #sp = spin_system.SpinSystem(Bx = 1)
        #moments, enes = sp.dispersion(0, .05, 200)
        #general_plotter.plot_dispersion(moments, enes, filename = "../plots/ex1_disp_bx_1.pdf", show = False)

        #sp = spin_system.SpinSystem(By = 1)
        #moments, enes = sp.dispersion(0, .05, 200)
        #general_plotter.plot_dispersion(moments, enes, filename = "../plots/ex1_disp_by_1.pdf", show = False)

        #sp = spin_system.SpinSystem(Bz = 1)
        #moments, enes = sp.dispersion(0, .05, 200)
        #general_plotter.plot_dispersion(moments, enes, filename = "../plots/ex1_disp_bz_1.pdf", show = False)

        #enes, conds = sp.conductance(0.05, 50)
        #general_plotter.plot_conductance(conds, enes, filename = "../plots/ex1_conductance_bz_1.pdf", show = False)
        
        ## to run for night
        #t_bys_values = np.linspace(0, 1, 50) 
        #E = 0.005
        #transmittances = np.zeros((4, len(t_bys_values)), dtype=np.float64)

        #for idx, by in enumerate(t_bys_values):
        #    sp = spin_system.SpinSystem(Bz = 0.1, additional_By = by, additional_By_boundaries=[0.3, 0.8])
        #    t_upup = sp.transmission(E, [1, 0], [1, 1])
        #    t_updown = sp.transmission(E, [1, 0], [0, 1])
        #    t_downup = sp.transmission(E, [1, 0], [1, 0])
        #    t_downdown = sp.transmission(E, [1, 0], [0, 0])
        #    transmittances[0, idx] = t_upup
        #    transmittances[1, idx] = t_updown
        #    transmittances[2, idx] = t_downup
        #    transmittances[3, idx] = t_downdown

        #general_plotter.plot_transmittance(t_bys_values, transmittances, show = True, filename = "../plots/ex1_transmittances.pdf")
        sp = spin_system.SpinSystem(Bz = 0.1, additional_By = 0.6, additional_By_boundaries=[0.3, 0.8])
        sys = sp.make_system();
        E = 0.005
        up, down, both, sys = sp.wave_function(E)
        general_plotter.plot_wavefunction(sys, up, show = True, filename = "ex1_density_up.pdf")
        general_plotter.plot_wavefunction(sys, down, show = True, filename = "ex1_density_down.pdf")
        general_plotter.plot_wavefunction(sys, both, show = True, filename = "ex1_density_both.pdf")

        #up, down, both, sys = sp.current(E)
        #general_plotter.plot_current(sys, up, show = True, filename = "ex1_current_up.pdf")
        #general_plotter.plot_current(sys, down, show = True, filename = "ex1_current_down.pdf")
        #general_plotter.plot_current(sys, both, show = True, filename = "ex1_current_both.pdf")

        
        

        


        

        


