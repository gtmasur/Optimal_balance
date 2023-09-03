# Update plot objects if saving
import numpy as np
#import matplotlib.pyplot as plt
from pylab import *

rc('text', usetex = True)
rc('font', **{'family' : "sans-serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
rcParams.update(params)

def update_save_2D(sim):

    Nx, Ny = sim.Nx, sim.Ny
    for var_cnt in range(len(sim.plot_vars)):

        var = sim.plot_vars[var_cnt]

        for L in range(sim.Nz):
        
            # attain data for system
            if sim.system == 'Real':
                u = sim.soln.u[0:sim.Nx,0:sim.Ny,L]
                v = sim.soln.v[0:sim.Nx,0:sim.Ny,L]
                h = sim.soln.h[0:sim.Nx,0:sim.Ny,L]
                time = sim.time
                proc = 'real'
                if sim.prognostic:
                    q = sim.solnp.q[:,:,L] 
                    delta = sim.solnp.delta[:,:,L] 
                    gamma = sim.solnp.gamma[:,:,L] 
                else: 
                    (q, delta, gamma) = sim.classic_to_prognostic(sim)
                
            if sim.system == 'Ramped':
                u = sim.ramp_soln.u[0:sim.Nx,0:sim.Ny,L]
                v = sim.ramp_soln.v[0:sim.Nx,0:sim.Ny,L]
                h = sim.ramp_soln.h[0:sim.Nx,0:sim.Ny,L]
                if sim.prognostic:
                    q = sim.ramp_solnp.q[:,:,L] 
                    delta = sim.ramp_solnp.delta[:,:,L] 
                    gamma = sim.ramp_solnp.gamma[:,:,L] 
                else: 
                    (q, delta, gamma) = sim.classic_to_prognostic(sim)
                
                time = sim.ramp_time
                if sim.ramp_direction == 'Forward':
                    proc = 'forward'
                    ramp_dir = 'Forward'
                elif sim.ramp_direction == 'Backward':
                    proc = 'backward' 
                    ramp_dir = 'Backward'

            if var == 'u':
                sim.ttls[var_cnt][L].set_text(r'Zonal Velocity $u$ ({0:s}): {1:.4f}'.format(proc,time))
                to_plot = u
                
            elif var == 'v':
                sim.ttls[var_cnt][L].set_text(r'Meridional Velocity $v$ ({0:s}): {1:.4f}'.format(proc,time))
                to_plot = v
                
            elif var == 'h':
                sim.ttls[var_cnt][L].set_text(r'Free Surface $h$ ({0:s}): {1:.4f}'.format(proc,time))
                to_plot = h
                
            elif var == 'zeta':
                # zeta = relative vorticity
                to_plot =  sim.ddx(v,sim) - sim.ddy(u,sim)
                sim.ttls[var_cnt][L].set_text(r'Rel. vorticity $zeta$ ({0:s}): {1:.4f}'.format(proc,time))
                
            elif var == 'q':
                # q = potential vorticity
                to_plot = q - sim.qbar
                #print to_plot.min(), to_plot.max()
                sim.ttls[var_cnt][L].set_text(r'Potential Vorticity $(q-\bar{{q}})$ ({0:s}): {1:.4f}'.format(proc,time))
                
            elif var == 'delta':
                # delta = velocity divergence 
                to_plot = delta
                sim.ttls[var_cnt][L].set_text(r'Velocity divergence $\delta$ ({0:s}): {1:.4f}'.format(proc,time))

            elif var == 'gamma':
                # gamma = ageostrophic vorticity
                to_plot = gamma
                sim.ttls[var_cnt][L].set_text(r'Agest. vorticity $\gamma$ ({0:s}): {1:.4f}'.format(proc,time))

            sim.Qs[var_cnt][L].set_array(to_plot.ravel())
            sim.Qs[var_cnt][L].changed()

        plt.draw() 
        if sim.system == 'Real':
            sim.figs[var_cnt].savefig(sim.path + '/Realt_Data/Frames/{0:s}_{1:04d}.png'.format(var,sim.frame_counter))
        elif sim.system == 'Ramped':
            """if sim.ramp_direction == 'Forward':
                ramp_dir = 'Forward'
            elif sim.ramp_direction == 'Backward':
                ramp_dir = 'Backward'"""
            sim.figs[var_cnt].savefig(sim.path + '/Artit_Data/' + ramp_dir
                                      +'/Frames/{0:s}_{1:04d}_cycle{2:d}_iter{3:d}.png'.format(var,
                                      sim.ramp_frame_counter,sim.nramping,sim.iter))
        
    if sim.system == 'Real':    
        sim.frame_counter += 1
    elif sim.system == 'Ramped':
        sim.ramp_frame_counter += 1

