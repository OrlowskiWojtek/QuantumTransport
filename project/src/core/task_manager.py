from . import graphene
from . import utils as  utl
from . import general_plotter as plt_wrapper
import os
import kwant as kw
import matplotlib.pyplot as plt
import numpy as np

class TaskManager:
    u = utl.Utils()
    gen_plt = plt_wrapper.Plotter()

    def test_task(self):
        gr = graphene.Graphene(sf=8)
        sys = gr.make_system()
        
        save_at = os.path.join(self.u.plots_path(), 'GrapheneSystem.png')
        kw.plot(sys, site_color = lambda site: sys.hamiltonian(site,site), cmap = 'RdBu', dpi = 300, show = True, file = save_at)
        plt.show()

        moments, enes = gr.get_bands();
        save_at = os.path.join(self.u.plots_path(), 'dispersion.pdf')
        self.gen_plt.plot_dispersion(moments, enes, show = True, filename = save_at)

    def task_1(self):
        Emax    = 0.7
        Emin    = -0.7
        En      = 30
        enes = np.linspace(Emin, Emax, En)
        V_n = 100
        Vmin = -0.1
        Vmax =  0.1
        V_np_vec = np.linspace(Vmin, Vmax, V_n)
        colors = plt.cm.cool(np.linspace(0,1,V_n))

        conds = []
        for V_np in V_np_vec:
            print("Calculating...", V_np)
            gr = graphene.Graphene(V_np = V_np, B = 0, sf = 8, d = 5)
            #e, _conds = gr.conductance(Emax, Emin, En)
            _conds = gr.transmission(0.6)

            conds.append(_conds)

#        plt.imshow(np.transpose(np.diff(np.array(conds), axis = 1)), cmap = "RdBu")
        plt.imshow(np.array(conds), cmap = "RdBu")

        np.savetxt('data/conductance_v2.dat', conds)

    def task_2(self):
        # in this task map B vs n1/n2 is required
        E = -0.1;

        V_n = 100
        V_min = 0.
        V_max = 0.1
        V_np_vec = np.linspace(V_min, V_max, V_n)

        B_n = 2
        B_min = 0.
        B_max = 1.
        B_vec = np.linspace(B_min, B_max, B_n)
        B = 1
        
        colors = plt.cm.cool(np.linspace(0,1,V_n))
        B_conds = []

        for B in B_vec:
            conds = []
            for V_np in V_np_vec:
                print("Calculating: ", B, V_np)
                gr = graphene.Graphene(V_np = V_np, B = B, sf = 8, d = 10)
                conds.append(gr.transmission(E))
            #    gr = graphene.Graphene(V_np = V_np, B = B, sf = 8, d = 5)
            #    cond = gr.transmission(E)
            B_conds.append(conds)

        for conds in B_conds:
            plt.plot(conds)
        plt.show()
        #np.savetxt('data/cond_b_temp.dat', conds)
#       # plt.imshow(np.transpose(np.diff(np.array(conds), axis = 1)), cmap = "RdBu")
        #plt.imshow(np.transpose(np.array(conds)), cmap = "RdBu")

        #plt.show()
        #
        #for idx in range(len(conds)):
        #    plt.plot(B_vec, conds[idx], label = f"V_np = {V_np_vec[idx]}", color = colors[idx])

        #plt.show()

