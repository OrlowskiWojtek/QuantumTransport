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
        gr = graphene.Graphene(sf = 16)
        sys = gr.make_system()
        
        save_at = os.path.join(self.u.plots_path(), 'GrapheneSystem_sf16.png')
        kw.plot(sys, site_color = lambda site: sys.hamiltonian(site,site), cmap = 'RdBu', dpi = 300, show = True, file = save_at)
        plt.show()

        #moments, enes = gr.get_bands();
        #save_at = os.path.join(self.u.plots_path(), 'dispersion.pdf')
        #self.gen_plt.plot_dispersion(moments, enes, show = True, filename = save_at)

    def task_1(self):
        Emax    = 0.5
        Emin    = -0.5
        En      = 100
        enes = np.linspace(Emin, Emax, En)
        V_n = 100
        Vmin = -0.2
        Vmax =  0.2
        V_np_vec = np.linspace(Vmin, Vmax, V_n)
        colors = plt.cm.cool(np.linspace(0,1,V_n))

        conds = []
        for V_np in V_np_vec:
            print("Calculating...", V_np)
            gr = graphene.Graphene(V_np = V_np, B = 0, sf = 8, d = 5)
            e, _conds = gr.conductance(Emax, Emin, En)
            #_conds = gr.transmission(0.6)

            conds.append(_conds)

#        plt.imshow(np.transpose(np.diff(np.array(conds), axis = 1)), cmap = "RdBu")
        plt.imshow(np.array(conds), cmap = "RdBu")
        plt.show()

        np.savetxt('data/conductance_temp.dat', conds)

    def task_2(self):
        # in this task map B vs n1/n2 is required
        E = 0.01;

        V_n = 50
        V_min = 0.02
        V_max = 0.08
        V_np_vec = np.linspace(V_min, V_max, V_n)

        B_n = 50
        B_min = -0.8
        B_max = 0.8
        B_vec = np.linspace(B_min, B_max, B_n)
        B = 1
        
        colors = plt.cm.cool(np.linspace(0,1,V_n))
        B_conds = []

        for B in B_vec:
            conds = []
            for V_np in V_np_vec:
                print("Calculating: ", B, V_np)
                gr = graphene.Graphene(V_np = V_np, B = B, sf = 8, d = 5)
                cond = gr.transmission(E)
                conds.append(cond)
            B_conds.append(conds)

        for (idx, conds) in enumerate(B_conds):
            plt.plot(conds, label = f"B = {B_vec[idx]}")
        plt.legend()
        plt.show()

        np.savetxt('data/cond_b_v2.dat', B_conds)
#       # plt.imshow(np.transpose(np.diff(np.array(conds), axis = 1)), cmap = "RdBu")
        plt.imshow(np.transpose(np.array(B_conds)), cmap = "RdBu")
        plt.show()

        #plt.show()
        #
        #for idx in range(len(conds)):
        #    plt.plot(B_vec, conds[idx], label = f"V_np = {V_np_vec[idx]}", color = colors[idx])

        #plt.show()

    def task_3(self):
        E = 0.
        fig, axs = plt.subplots(2, 2, figsize = (10,6));

        V_np = 0.0
        gr = graphene.Graphene(V_np = V_np, B = 0, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[0, 0])
        axs[0, 0].set_title(f"E = {E} eV, V_np = {V_np} eV")

        V_np = 0.05
        gr = graphene.Graphene(V_np = V_np, B = 0, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[0, 1])
        axs[0, 1].set_title(f"E = {E} eV, V_np = {V_np} eV")

        E = 0.01
        V_np = 0.0
        gr = graphene.Graphene(V_np = V_np, B = 0, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[1, 0])
        axs[1, 0].set_title(f"E = {E} eV, V_np = {V_np} eV")

        V_np = 0.05
        gr = graphene.Graphene(V_np = V_np, B = 0, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[1, 1])
        axs[1, 1].set_title(f"E = {E} eV, V_np = {V_np} eV")

        plt.tight_layout()
        plt.savefig("../figures/currents.pdf")
        plt.show()

    def task_4(self):
        E = 0.01
        fig, axs = plt.subplots(2, 2, figsize = (10, 6));

        V_np = 0.05
        B = 0.875
        V_np = 0.0205
        gr = graphene.Graphene(V_np = V_np, B = B, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[0, 0])
        axs[0, 0].set_title(f"B = {B} T, V_np = {V_np} eV")

        V_np = 0.03
        gr = graphene.Graphene(V_np = V_np, B = B, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[0, 1])
        axs[0, 1].set_title(f"B = {B} T, V_np = {V_np} eV")

        B = 0.875
        V_np = 0.0395
        gr = graphene.Graphene(V_np = V_np, B = B, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[1, 0])
        axs[1, 0].set_title(f"B = {B} T, V_np = {V_np} eV")

        V_np = 0.055
        gr = graphene.Graphene(V_np = V_np, B = B, sf = 8, d = 5)
        gr.current(E, 0, 0, ax = axs[1, 1])
        axs[1, 1].set_title(f"B = {B} T, V_np = {V_np} eV")

        plt.tight_layout()
        plt.savefig("../figures/currents_changing_b.pdf")
        plt.show()

