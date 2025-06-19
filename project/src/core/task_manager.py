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
        gr = graphene.Graphene(sf=15)
        sys = gr.make_system()
        
        save_at = os.path.join(self.u.plots_path(), 'GrapheneSystem.png')
        kw.plot(sys, site_color = lambda site: sys.hamiltonian(site,site), cmap = 'RdBu', file = save_at, dpi = 300)

        
        moments, enes = gr.get_bands();
        self.gen_plt.plot_dispersion(moments, enes, show = True)
        
