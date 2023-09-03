# Modified version of Francis Poulin's open souce PyRsw code 
"""
GTM: This file includes the main classes for the nondimensional
     rotating shallow water model. 
     The model is evolved by the classical variables but
     the prognostic varibles can be also computed by kinematic
     PV-inversion equations.
     
    classic : primitive variables (u,v,h)
    prognostic: geo-ageo variables (q,delta,gamma)     
"""

from pylab import *
import Plot_tools
import Steppers
import Fluxes as Fluxes
import Ramping as Ramp 
import Simdata
import Siminfo
import sys

    
# possible to add one more layer at Nz for topography
class Solution():
    def __init__(self,Nx,Ny,Nz):
        self.u = zeros((Nx,Ny,Nz)) 
        self.v = zeros((Nx,Ny,Nz))
        self.h = zeros((Nx,Ny,Nz))
        
class Solutionprog():
    def __init__(self,Nx,Ny,Nz):
        self.q = zeros((Nx,Ny,Nz)) 
        self.delta = zeros((Nx,Ny,Nz))
        self.gamma = zeros((Nx,Ny,Nz))

class Flux():
    def __init__(self):
        self.u = []
        self.v = []
        self.h = []

class Simulation:

    # First-level initialization, default values
    def __init__(self):
    
        # Note: not all attributions are needed to be stated here!
    
        self.nfluxes = 0            # Length of flux history
        self.fluxes = Flux()   
        self.num_steps = 0
    
        self.time = 0.              # initial time
        self.min_dt = 1e-7          # minimum timestep
        
        self.system = 'Real'        # Real or Ramped
        self.cmap = 'jet'       # Default colour map
        

    def build_grid(self):
        
        dxs = [1,1]
        #dx = self.Lx/(self.Nx-1)
        #self.x = arange(dx/2,self.Lx,dx)- self.Lx/2.
        dx = self.Lx/(self.Nx-1)  
        self.x = linspace(0,self.Lx,self.Nx)
        dxs[0] = dx
#        dy = self.Ly/self.Ny
#        self.y = arange(dy/2,self.Ly,dy) - self.Ly/2.
        dy = self.Ly/(self.Ny-1)  
        self.y = linspace(0,self.Ly,self.Ny)
        dxs[1] = dy
        [self.X, self.Y]  = meshgrid(self.x, self.y, indexing='ij')  
        self.dx = dxs
        
        [self.grid_x.u, self.grid_y.u] = meshgrid(self.x,self.y, indexing='ij')  
        [self.grid_x.v, self.grid_y.v] = meshgrid(self.x,self.y, indexing='ij')
        [self.grid_x.h, self.grid_y.h] = meshgrid(self.x,self.y, indexing='ij')

    
    def set_scaling(self):

        # 1. arrange h0 and time scalings depending on scaling limits
        # 2. set mean PV values for plotting and initial conditions
        # this is needed to derive initial conditions for rand PVs
        if self.scalinglimit == 'quasi-geostrophic':
            print( "" )
            print( "Quasi-geostrophic scaling limit....")
            self.h0 = 1./self.eps              # scaled mean height
            self.T_ramp    = self.T            # ramping time
            self.t_forward = self.tprime       # forward time
            self.qbar = self.eps
        elif self.scalinglimit == 'semi-geostrophic':
            print( "")
            print( "Semi-geostrophic scaling limit....")
            self.h0 = 1.
            self.T_ramp    = self.T/self.eps
            self.t_forward = self.tprime/self.eps
            self.qbar = 1.
        
        self.end_time = self.t_forward 
        
    
    # Full initialization for once the user has specified parameters
    def initialize(self):
            
        # Initial conditions
        self.soln      = Solution(self.Nx,self.Ny,self.Nz)   # primitive solution
        self.curr_flux = Solution(self.Nx,self.Ny,self.Nz)
        self.grid_x    = Solution(self.Nx,self.Ny,self.Nz)
        self.grid_y    = Solution(self.Nx,self.Ny,self.Nz)
        # Prognostic solution
        if self.prognostic:
            self.solnp = Solutionprog(self.Nx,self.Ny,self.Nz)  
       
        # Fluxes are taken as Fluxes since we have Flux class.
        self.flux_method = Fluxes.spectral_sw
        # Define (q,delta,gamma) in the case of plotting
        self.classic_to_prognostic = Ramp.classic_to_prognostic   
        # Create random initial data wrt specific choices
        self.initialize_fields = Ramp.initialize_fields 
        # Set tools to use in initial or boundary conditions
        self.set_conditions = Ramp.set_conditions
        
        # Initialize grid, and cell centers
        self.build_grid()
        # Initialize parameters depending on scaling limits
        self.set_scaling()
        # Initialize differentiation and operators
        self.flux_method(self)             
        self.spectral_operators(self)     # set inside Flux.spectral_sw
        self.set_conditions(self)
          
          
    def prepare_for_run(self):   
    
        # If prognostic, initialize prognostic solution
        if self.prognostic:
            self.classic_to_prognostic(self)
        
        # If saving, initialize those too        
        if self.output:
            self.next_save_time = self.savet 
            self.output_counter = 0 
            self.Specs = []
        
        # If saving, initialize the directory
        if self.output or (self.animate == 'Save'):  
            Simdata.initialize_saving(self)
            
        # If plotting, then initialize the plots  
        if self.animate != 'None':
            self.frame_counter = 0       
            # add empty array for unspecified varlims
            self.clims += [[]]*(len(self.plot_vars) - len(self.clims)) 
            self.initialize_plots = Plot_tools.initialize_plots_animsave_2D
            self.initialize_plots(self)
            self.update_plots(self) 
            
            # GTM: no need to keep things complicated    
            self.frame_counter = 1
            self.next_plot_time = self.plott 
                 
                 
    # Compute the current flux              
    def flux(self):
        return self.spectral_sw_flux(self)    
        
        
    # Compute time-step using CFL condition.
    def compute_dt(self):
    
        # The fixed dt is set by Dritchlel's code but it fails.
        # Update Poulin's dt with c = 1./sim.eps
        if self.adaptive:
        
            def dt(max_val):
                if self.Nx != self.Ny:
                    print( "Error! Reset compute_dt function...."); sys.exit()
                return self.cfl*(self.dx[0]/(max_val+2./self.eps)) 
                
            t = self.end_time
            max_u = abs(self.soln.u).max()
            max_v = abs(self.soln.v).max()
            dt_gw = self.fixed_dt
            dt_u  = min(dt(max_u),dt(max_v))
            if self.prognostic: 
                max_delta = abs(self.solnp.delta).max()
                max_gamma = abs(self.solnp.gamma).max()  
                dt_gamma = dt(max_gamma)
                dt_delta = dt(max_delta)
                
            if self.tstep == 'dtu':  
                self.dt = dt_u
            elif self.tstep == 'dtmin':
                if not self.prognostic: 
                    print( "Error! Set prognostic true."); sys.exit()
                self.dt = min([dt_gw, dt_u, dt_gamma, dt_delta]) 
             
        elif not (self.adaptive):
            # This fails....
            self.dt = self.fixed_dt
         
        # Slowly ramp-up the dt for the first 20 steps
        if self.num_steps <= 20:
            self.dt *= 1./(5*(21-self.num_steps))
    
       
    # Adjust time-step when necessary    
    def adjust_dt(self):
        
        t = self.time + self.dt
        # next time to save or plot
        nt = self.end_time    
        if self.animate != 'None':
            nt = min([self.next_plot_time, nt])
        if self.output:
            nt = min([self.next_save_time, nt])
        # needed for different end_time and t_forward
        if self.optimal_balance:
            nt = min([self.end_time, nt])
            
        # to stop at desired time
        if (nt < t):   
            self.dt = nt - self.time
        t = self.time + self.dt

        # GTM: plot and save at the end of time
        do_plot = False
        if self.animate != 'None':
            if t >= self.next_plot_time or nt==t: 
                # to plot at given plott intervals
                do_plot = True
        do_save = False
        if self.output:
            if t >= self.next_save_time or nt==t:
                # to save at given savet intervals
                do_save = True
        return do_plot, do_save
      

    # Evolve the model one real time-step.
    def step(self):
    
        self.compute_dt() 
        # to match an output time  
        do_plot, do_save = self.adjust_dt()
        self.stepper(self)    # flux() inside stepper
        self.apply_filter(self)
        if self.prognostic:
            self.classic_to_prognostic(self)  # update (q,d,g)
        self.time += self.dt
        if do_plot:
            self.update_plots(self)   
            self.next_plot_time += self.plott
        if do_save:
            Simdata.save_state(self)
            Simdata.save_specs(self)
            self.next_save_time += self.savet 
        # update the records
        self.num_steps += 1
          
   
    # NOTE: after (re)balancing, no need to compute prognostic 
    # variables since new ones already attained in Nonlinear_BC
    def run_algorithm1(self):
        # The main diagnosed imbalance algorithm
    
        # initialize optimal balance states and parameters
        Ramp.initialize_ramp(self)
        # TODO: replace set_kappa.....
        self.set_kappa(self)
        
        self.prepare_for_run() 
      
        # balance the system
        print( ""  )
        print( "===== balancing ======" )
        self.system = 'Ramped' 
        self.balance_wave(self)
        print( self.optbal_iterate)  # added to check
        
        # real time evolution
        print( ""  )
        print( "======== real =======")
        self.system = 'Real'
        # save balanced states    
        if self.output:
            Simdata.save_state(self)     # save balanced state for output_counter=0
            Simdata.save_specs(self)     # being saved starting balanced state
        while self.time < self.end_time:
            self.step()
        
        # re-balance the system
        print( ""  )
        print( "==== rebalancing ====")
        self.system = 'Ramped'   
        self.balance_wave(self)
        print( self.optbal_iterate)
        
        # Return the system to Real
        self.system = 'Real'
        # save rebalanced data for output_counter=end and imbalances
        # write and print( imbalances array)
        self.imbalances_nonlinear(self)
        if self.output:
            Simdata.save_additional_state(self) 
            Simdata.save_imbalances(self)     
          
        return
        
    def run_algorithm2(self):
        # Balance the wave and after physical time evolution,
        # quantify imbalances at linear-end.
        
        # initialize optimal balance states and parameters
        Ramp.initialize_ramp(self)   
        self.prepare_for_run() 
        
        # balance the system  
        print( ""  )
        print( "===== balancing ======" )
        self.system = 'Ramped' 
        self.balance_wave(self)
        
        # real time evolution
        print( ""  )
        print( "======== real =======")
        self.system = 'Real'
        # save balanced states    
        if self.output:
            Simdata.save_state(self)     # save balanced data for output_counter=0
            Simdata.save_specs(self)     # being saved starting balanced state
        while self.time < self.end_time:
            self.step()
        
        # re-balance the system
        print( ""  )
        print( "==== reaching linear end ====")
        self.system = 'Ramped'   
        self.reach_linearend(self)
         
        # Return the system to Real
        self.system = 'Real'
        # save rebalanced data for output_counter=end and imbalances
        # write and print( imbalances array)
        self.imbalances_linear(self)
        if self.output:
            Simdata.save_additional_state(self) 
            Simdata.save_imbalances(self)           
        
        return
    
    
    def run_algorithm3(self):
        # For given ramping time T, balance the wave for T and 2T
        # rebalance it with T or 2T: decide with sim.rebalT
        
        # initialize optimal balance states and parameters
        Ramp.initialize_ramp(self)
        self.prepare_for_run() 
        
        # balance the system for T ramp time 
        print( ""  )
        print( "===== balancing for T ======" )
        self.system = 'Ramped' 
        print( self.T_ramp)
        self.balance_wave(self)
        
        # balance the system for 2T ramp time 
        print( ""  )
        print( "===== balancing for 2T ======" )
        self.system = 'Ramped' 
        self.T_ramp *= 2  # double Tramp value
        print( self.T_ramp)
        self.balance_wave(self)
        self.T_ramp /= 2  # go back to original value
        
        # real time evolution
        print( ""  )
        print( "======== real =======")
        self.system = 'Real'
        # save balanced states    
        if self.output:
            Simdata.save_state(self)     # save balanced data for output_counter=0
            Simdata.save_specs(self)     # being saved starting balanced state
        while self.time < self.end_time:
            self.step()
        
        # re-balance the system with ramp time T or 2T
        print( ""  )
        print( "==== rebalancing for {0:s} ====".format(self.rebaltime))
        if self.rebaltime == 'T' :  self.T_ramp = self.T_ramp
        if self.rebaltime == '2T':  self.T_ramp *= 2
        self.system = 'Ramped'   
        print( self.T_ramp)
        self.balance_wave(self)
        
        # Return the system to Real
        self.system = 'Real'
        # save rebalanced data for output_counter=end and imbalances
        # write and print( imbalances array)
        self.imbalances_nonlinear(self)
        if self.output:
            Simdata.save_additional_state(self) 
            Simdata.save_imbalances(self)            
        return
         
    
    def run_algorithm_L1orGeo(self):  
        
        # by default...
        self.base_point == "h"
        # we use ramped and real to compute relative norm
        Ramp.initialize_ramp(self)
        self.prepare_for_run() 
        
        print( ""  )
        if self.optbal_algorithm == 'algorithmL1':
            print( "===== L1 balancing ======" )
            self.set_L1u(self)
        elif self.optbal_algorithm == 'algorithmGeo':
            print( "===== geo balancing ======" )
            self.set_geostrophy(self)
        print( ""  )
        print( "======== real =======")
        self.system = 'Real'   
        if self.output:
            Simdata.save_state(self)    
            Simdata.save_specs(self)     
        while self.time < self.end_time:
            self.step()
        
        print( "")
        # rebalancing in balance_wave_L1orGeo...
        self.system = 'Ramped'  
        self.balance_wave_L1orGeo(self)
        
        self.system = 'Real'
        self.imbalances_nonlinear(self)
        if self.output:
            Simdata.save_additional_state(self) 
            Simdata.save_imbalances(self) 
        return                 
           
                       
    def run(self):
        if self.optimal_balance:
            
            if self.optbal_algorithm == 'algorithm1':
                print( "")
                print( "Optimal balance: algorithm 1")
                self.run_algorithm1()  
                return self.imb
            elif self.optbal_algorithm == 'algorithm2':
                print( "")
                print( "Optimal balance: algorithm 2")
                self.run_algorithm2()  
                return self.imb
            elif self.optbal_algorithm == 'algorithm3':
                print( "")
                print( "Optimal balance: algorithm 3")
                self.run_algorithm3()  
                return self.imb  
            elif self.optbal_algorithm == 'algorithmL1':
                print( "")
                print( "Balance: L1-balanced model")
                self.run_algorithm_L1orGeo()  
                return self.imb 
            elif self.optbal_algorithm == 'algorithmGeo':
                print( "")
                print( "Balance: Geo-balanced model")
                self.run_algorithm_L1orGeo()  
                return self.imb       
        else:
            # Run SWEs without optimal balance
            self.run_swe()             
            return    
    
    
    def run_swe(self):   
    
        self.prepare_for_run()
        while self.time < self.end_time:
            self.step()
        return      
        
        
