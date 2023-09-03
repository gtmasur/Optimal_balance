# Initialize plot objects for anim or save
# Assume that the field is 2-dimensional

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.ticker as mtick
from pylab import *
from numpy.fft import *
from .update_anim_2D import update_anim_2D
from .update_save_2D import update_save_2D

rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


def initialize_plots_animsave_2D(sim):

    figs = []
    Qs   = []
    ttls = []

    # Loop through each element of plot_vars
    # Each field will be given its own figure
    for var_cnt in range(len(sim.plot_vars)):

        Qs   += [[]]
        ttls += [[]]

        var = sim.plot_vars[var_cnt]
        
        # GTM: set plot size
        fig = plt.figure(figsize = (6,4))
        #ax = fig.add_subplot(111)
        #ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        #plt.xlabel("X-Position")
        #plt.ylabel("Y-Position")
        
        #fig = plt.figure()
        figs += [fig]

        # GTM2 : In python xy axes are switched, be sure that you put right
        #       data to right place after labeling axes !!!!!!        

        # Plot the data
        for L in range(sim.Nz):
        
            # Plotting tools are adjusted only for "Spectral" method

            plt.subplot(1,sim.Nz,L+1)
            
            # attain data for system
            if sim.system == 'Real':
                u = sim.soln.u[:,:,L]
                v = sim.soln.v[:,:,L]
                h = sim.soln.h[:,:,L]
                proc = 'real'
                if sim.prognostic:
                    q = sim.solnp.q[:,:,L] 
                    delta = sim.solnp.delta[:,:,L] 
                    gamma = sim.solnp.gamma[:,:,L] 
                else: 
                    (q, delta, gamma) = sim.classic_to_prognostic(sim)
                
            if sim.system == 'Ramped':
                u = sim.ramp_soln.u[:,:,L]
                v = sim.ramp_soln.v[:,:,L]
                h = sim.ramp_soln.h[:,:,L]
                if sim.prognostic:
                    q = sim.ramp_solnp.q[:,:,L] 
                    delta = sim.ramp_solnp.delta[:,:,L] 
                    gamma = sim.ramp_solnp.gamma[:,:,L] 
                else: 
                    (q, delta, gamma) = sim.classic_to_prognostic(sim)
                    
                if sim.ramp_direction == 'Forward':
                    proc = 'forward'
                elif sim.ramp_direction == 'Backward':
                    proc = 'backward'
            
            if var == 'u':
                ttl = fig.suptitle(r'Zonal Velocity $u$ ({0:s}): t=0'.format(proc))
                to_plot = u 
                X = sim.grid_x.u
                Y = sim.grid_y.u
                
            elif var == 'v':
                ttl = fig.suptitle(r'Meridional Velocity $v$ ({0:s}): t=0'.format(proc))
                to_plot = v 
                X = sim.grid_x.v
                Y = sim.grid_y.v
                
            elif var == 'h':
                ttl = fig.suptitle(r'Free Surface $h$ ({0:s}): t=0'.format(proc))
                to_plot = h 
                X = sim.grid_x.h
                Y = sim.grid_y.h
                    
            elif var == 'zeta':
                # zeta = relative vorticity
                to_plot =  sim.ddx(v,sim) - sim.ddy(u,sim)
                X = sim.grid_x.h
                Y = sim.grid_y.h
                ttl = fig.suptitle(r'Rel. vorticity $\zeta$ ({0:s}): t=0'.format(proc))
            
            elif var == 'q':
                # q = potential vorticity
                to_plot = q - sim.qbar
                X = sim.grid_x.h
                Y = sim.grid_y.h
                ttl = fig.suptitle(r'Potential Vorticity $(q-\bar{{q}})$ ({0:s}): t=0'.format(proc))    

            elif var == 'delta':
                # delta = velocity divergence 
                to_plot = delta
                X = sim.grid_x.h
                Y = sim.grid_y.h
                ttl = fig.suptitle(r'Velocity divergence $\delta$ ({0:s}): t=0'.format(proc))

            elif var == 'gamma':
                # gamma = ageostrophic vorticity
                to_plot = gamma
                #print to_plot.min(), to_plot.max()
                X = sim.grid_x.h
                Y = sim.grid_y.h
                ttl = fig.suptitle(r'Agest. vorticity $\gamma$ ({0:s}): t=0'.format(proc))

            # Has the user specified plot limits?
            if len(sim.clims[var_cnt]) == 2:
                vmin = sim.clims[var_cnt][0]
                vmax = sim.clims[var_cnt][1]
            else:
                cv = np.max(np.abs(to_plot.ravel()))
                vmin = -cv
                vmax =  cv

            # Extend the grid to account for pcolor peculiarities
            Nx = X.shape[0]
            Ny = Y.shape[1]
            X_plot = np.zeros((Nx+1,Ny+1))
            Y_plot = np.zeros((Nx+1,Ny+1))

            X_plot[1:,1:] = X + sim.dx[0]/2.
            X_plot[1:,0]  = X[:,0] + sim.dx[0]/2.
            X_plot[0,:]   = X[0,0] - sim.dx[0]/2.

            Y_plot[1:,1:] = Y + sim.dx[1]/2.
            Y_plot[0,1:]  = Y[0,:] + sim.dx[1]/2.
            Y_plot[:,0]   = Y[0,0] - sim.dx[1]/2.

            #Q = plt.pcolormesh(X_plot/1e3, Y_plot/1e3, to_plot, cmap=sim.cmap, 
            #            vmin = vmin, vmax = vmax)
            Q = plt.pcolormesh(X_plot, Y_plot, to_plot, cmap=sim.cmap, 
                        vmin = vmin, vmax = vmax)
            Qs[var_cnt] += [Q]
            ttls[var_cnt] += [ttl]

            plt.colorbar()   # GTM : Change the size of the colorbar!

            plt.axis('tight')

            if 1./1.1 <= sim.Ly/sim.Lx <= 1.1:
                plt.gca().set_aspect('equal')

    if sim.animate == 'Anim':
        # update frames in animation only for real time evolution
        #if sim.system == 'Real':
        sim.update_plots = update_anim_2D
    elif sim.animate == 'Save':
        sim.update_plots = update_save_2D
        plt.ioff()
        #plt.pause(0.01)

    if sim.animate == 'Anim':
        plt.ion()
        plt.pause(0.01)
        plt.draw()

    sim.figs = figs
    sim.Qs = Qs
    sim.ttls = ttls    

