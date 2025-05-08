from . import SpinPrecession as spin_system
import kwant as kw
import numpy as np
import Core.Utils as utl
import matplotlib.pyplot as plt

u = utl.Utils()

class TaskManager():
    def task1(self, *, do_plot = False):
        sp = spin_system.SpinSystem()
        sys = sp.make_system()
        kw.plot(sys, show = do_plot);

        moments, enes = sp.dispersion(0, .1, 200)
        plt.figure();
        plt.plot(moments, np.asarray(enes) / u.eV2au(1.0),'k-')
        plt.tick_params(axis='both', which='major', labelsize=22)
        plt.ylim((0,.04))  
        plt.xlim((-0.15,0.15))  
        plt.xlabel("k [1/nm]",fontsize=22)
        plt.ylabel("E [eV]",fontsize=22)
        plt.tight_layout()
        plt.savefig("../plots/disp_ex1.pdf")
        plt.show()


        

        


