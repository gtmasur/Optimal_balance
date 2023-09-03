"""
GTM : This file includes main optimal balance procedures.
"""

from pylab import *
import Steppers
import Fluxes as Fluxes
import Simdata
from scipy.fftpack import fftfreq, fft2, ifft2
import sys

from .Linear_BC import *
from .Nonlinear_BC import *


class Solution():
    def __init__(sim,Nx,Ny,Nz):
        sim.u = zeros((Nx,Ny,Nz)) 
        sim.v = zeros((Nx,Ny,Nz))
        sim.h = zeros((Nx,Ny,Nz))
        
class Solutionprog():
    def __init__(sim,Nx,Ny,Nz):
        sim.q = zeros((Nx,Ny,Nz)) 
        sim.delta = zeros((Nx,Ny,Nz))
        sim.gamma = zeros((Nx,Ny,Nz))

class Flux():
    def __init__(sim):
        sim.u = []
        sim.v = []
        sim.h = []
     
def initialize_ramp(sim):   

    # set parameters
    # 1st cycle: balancing loop
    # 2nd cycle: rebalancing loop
    sim.nramping = 0           # needed to separate 1st and 2nd cycle in optbal 
    sim.iters = []             # empty array to save iteration numbers
    sim.ramp_direction = 'Backward' # Backward: backward time in ramping
    sim.wave = 'Rossby'

    # initialize ramping conditions
    sim.ramp_soln        = Solution(sim.Nx,sim.Ny,sim.Nz)
    sim.ramp_curr_flux   = Solution(sim.Nx,sim.Ny,sim.Nz)
    if sim.prognostic:                         
        sim.ramp_solnp   = Solutionprog(sim.Nx,sim.Ny,sim.Nz)

    # TODO: new...
    # Decide kappa value
    sim.set_kappa = set_kappa

    # introduce ramping procedure operators
    sim.ramp_flux_method = Fluxes.spectral_ramp_sw    
    sim.balance_wave     = balance_wave  
    sim.reach_linearend  = reach_linearend       
    sim.build_proj_mats  = build_proj_mats
    if sim.optbal_algorithm in ['algorithmL1', 'algorithmGeo'] :
        sim.balance_wave_L1orGeo   = balance_wave_L1orGeo
    sim.ramp_flux = ramp_flux
    sim.imbalances_nonlinear = imbalances_nonlinear  
    sim.imbalances_linear = imbalances_linear    
    # define a function to calculate norms
    sim.compute_norm = compute_norm
    sim.project_on_physicalspace = project_on_physicalspace
    
    # initialize ramping procedure operators
    sim.ramp_flux_method(sim)  
    
    # build projection matrices if needed
    if sim.linear_bc == "obliqpr" or sim.linear_bc == "orthpr":
        sim.build_proj_mats(sim)


def prepare_for_ramp_run(sim):   
    
    # If prognostic, initialize prognostic solution
    if sim.prognostic:
        sim.classic_to_prognostic(sim)
    
    # goal: same restriction on ramp_dt as in dt
    # and make initializations for ramping
    sim.ramp_time = 0.      
    sim.ramp_num_steps = 0  
    sim.ramp_nfluxes = 0   
    sim.ramp_fluxes = Flux()
       
    # If saving, initialize those too  
    # If not specificially required, Specs not saved      
    if sim.ramp_output:
        # adjust ramp_savet in test files
        sim.ramp_next_save_time = sim.ramp_savet 
        sim.ramp_output_counter = 0   
        # initiliaze energy spectrum
        if sim.ramp_spec_output:
            sim.ramp_Specs = []   
        # sim.initialize_saving is already run 
        Simdata.ramp_save_state(sim)     # first and end data is saved
        if sim.ramp_spec_output:
            Simdata.ramp_save_specs(sim)     
        
    # plotting already initialized...
    if sim.animate != 'None':
        sim.ramp_frame_counter = 0
        # update plots
        sim.update_plots(sim)          
        sim.ramp_frame_counter = 1  # GTM: needed, dont delete
        sim.ramp_next_plot_time = sim.plott
    
     
# GTM: Compute the ramped flux_function
def ramp_flux(sim):
    return sim.spectral_ramp_sw_flux(sim)    

        
def compute_ramp_dt(sim):

    if sim.adaptive:
        def dt(max_val):
            if sim.Nx != sim.Ny:
                print ("Error! Reset compute_dt function...."); sys.exit()
            return sim.cfl*(sim.dx[0]/(max_val+2./sim.eps)) 
            
        t = sim.end_time
        max_u = abs(sim.ramp_soln.u).max()
        max_v = abs(sim.ramp_soln.v).max()
        dt_x     = dt(max_u)  
        dt_y     = dt(max_v)
        dt_u     = min(dt_x,dt_y)
        if sim.prognostic:
            max_delta = abs(sim.ramp_solnp.delta).max()
            max_gamma = abs(sim.ramp_solnp.gamma).max()
            dt_gamma = dt(max_gamma)
            dt_delta = dt(max_delta)
           
        if sim.tstep == 'dtmin':
            sim.ramp_dt = min([dt_u, dt_gamma, dt_delta]) 
        elif sim.tstep == 'dtu':    
            sim.ramp_dt = dt_u     
               
    elif not (sim.adaptive):
        sim.ramp_dt = sim.fixed_dt
    
    # Slowly ramp-up the dt for the first 20 steps
    if sim.ramp_num_steps <= 20:
        sim.ramp_dt *= 1./(5*(21-sim.ramp_num_steps))  
            

# Adjust ramp_time-step when necessary         
def adjust_ramp_dt(sim):
    
    t = sim.ramp_time + sim.ramp_dt    
    # next time to save or plot
    nt = sim.T_ramp    
    if sim.animate != 'None':
        nt = min([sim.ramp_next_plot_time, nt])
    if sim.ramp_output:
        nt = min([sim.ramp_next_save_time, nt])
    
    # finish ramping at exactly T_ramp
    if (nt < t):
        sim.ramp_dt = nt - sim.ramp_time
    t = sim.ramp_time + sim.ramp_dt

    # GTM: plot and save at the end of time
    do_plot = False
    if sim.animate != 'None':
        if t >= sim.ramp_next_plot_time or nt==t: 
            # to plot at given plott intervals
            do_plot = True
    do_save = False
    if sim.ramp_output:
        if t >= sim.ramp_next_save_time or nt==t:
            # to save at given savet intervals
            do_save = True
    return (do_plot, do_save)
    

# Evolve the ramped system regarding ramp_direction
def step_ramp(sim):
    
    compute_ramp_dt(sim)
    # Check ramp_dt to match T_ramp
    do_plot, do_save = adjust_ramp_dt(sim)
    sim.stepper(sim)  # ramp_flux() inside stepper
    sim.apply_ramp_filter(sim)  
    if sim.prognostic:
        classic_to_prognostic(sim)  # update (q,d,g)
    sim.ramp_time += sim.ramp_dt
    if do_plot:
        sim.update_plots(sim)   
        sim.ramp_next_plot_time += sim.plott
    if do_save:
        #print "Eroor! Check step_ramp in Optbal."; sys.exit()
        Simdata.ramp_save_state(sim)
        if sim.ramp_spec_output:
            Simdata.ramp_save_specs(sim)     
        sim.ramp_next_save_time += sim.ramp_savet
    # Update the records
    sim.ramp_num_steps += 1
   
   
def run_backward_scheme(sim):   
    sim.ramp_direction = 'Backward'
    prepare_for_ramp_run(sim)
    while sim.ramp_time < sim.T_ramp:
        step_ramp(sim)   

def run_forward_scheme(sim):
    sim.ramp_direction = 'Forward'
    prepare_for_ramp_run(sim)
    while sim.ramp_time < sim.T_ramp:
        step_ramp(sim)
    

def set_kappa(sim):        
    # Reset kappa to default value again
    sim.kappa = float("1e-{0:g}".format(sim.n))
    
    if sim.base_point == "q":
        print ("kappa takes the default value.")
        print ("kappa:", sim.kappa)
        
    elif sim.base_point == "h":
        
        # TODO....
        # This is arranged to test stopping criteria of basec h
        
        # TODO: renew he ames depending on below values....!!!!
#        if sim.eps >= 0.4:
#            sim.kappa = 2e-2
#        elif 0.4 > sim.eps >= 0.25:
#            sim.kappa = 1e-2
#        elif 0.25 > sim.eps >= 0.1:
#            # modified to have different kappa values for eps=0.1...
#            # only for frames!!!
#            if sim.n == 3:
#                sim.kappa = 5e-3
#            elif sim.n == 4: 
#                sim.kappa = 1e-3
#        elif 0.1 > sim.eps >= 0.07:
#            sim.kappa = 1e-3
#        elif 0.07 > sim.eps :
#            if sim.n == 3:
#                sim.kappa = 1e-3
#            elif sim.n == 4: 
#                sim.kappa = 5e-4
#                """no difference between 5e-4 and 1e-4"""
                
        if sim.eps >= 0.4:
            sim.kappa = 2e-2
        elif 0.4 > sim.eps >= 0.25:
            sim.kappa = 1e-2
        elif 0.25 > sim.eps >= 0.1:
            # modified to have different kappa values for eps=0.1...
            # only for frames!!!
            if sim.n == 25:
                sim.kappa = 5e-3
            elif sim.n == 3: 
                sim.kappa = 1e-3
            elif sim.n == 35: 
                sim.kappa = 5e-4
            elif sim.n == 4: 
                # this is added only for this test case
                sim.kappa = 1e-4
                
        elif 0.1 > sim.eps >= 0.07:
            sim.kappa = 1e-3
        elif 0.07 > sim.eps :
            if sim.n == 3:
                sim.kappa = 1e-3
            elif sim.n == 35: 
                sim.kappa = 5e-4
            elif sim.n == 4: 
                sim.kappa = 1e-4
                """no difference between 5e-4 and 1e-4"""        
                
        print ("kappa:", sim.kappa )
        
        
def balance_wave(sim):
    # Goal: Balance the wave by optimal balance method. 
    # Run the balancing loop over T_ramp until converging base_point.
    # The whole optimal balance is summarized below: 
    
    # Start with u,v,h from soln
    sim.ramp_soln.h[:,:,:] = sim.soln.h[:,:,:]
    sim.ramp_soln.u[:,:,:] = sim.soln.u[:,:,:]
    sim.ramp_soln.v[:,:,:] = sim.soln.v[:,:,:]
    
    # TODO: commented out...
#    # Decide kappa value
#    set_kappa(sim)
    #print "kappa takes the default value."
    #print "kappa:", sim.kappa 
    
    # Calculate norm = (xTnn - xTn)/xTnn in nonlinear_bc  
    if sim.base_point == "q":
        sim.qTn = zeros((sim.Nx,sim.Ny))  # qTn = qT(n) 
    elif sim.base_point == "h": 
        sim.hTn = zeros((sim.Nx,sim.Ny))   # hTn = hT(n)
     
    sim.iter = 0
    sim.optbal_iterate = True
    while sim.optbal_iterate:
        sim.iter += 1
        print ("")
        print ("........ iteration .........", sim.iter)
        print (sim.optbal_iterate)
        # 1. Start at nonlinear swe (appr. slow manifold)
        #    Evolve the ramped system back in time 
        run_backward_scheme(sim)
        # 2. Reaches at linear swe (trivial slow manifold)
        #    Set BC on linear swe end by projection or geostrophy
        set_linear_bc(sim)
        # 3. Evolve the ramped system forward in time.
        run_forward_scheme(sim)
        # 4. Reaches at nonlinear swe (appr. slow manifold)   
        #    Set BC at nonlinear swe end depending on base_point
        set_nonlinear_bc(sim) # outcome is in ramp_sol   
        
    # increase for second cycle   
    # 0: balancing cycle
    # 1: rebalancing cycle
    sim.nramping += 1  
    sim.iters += [sim.iter]
    
    return
    

def reach_linearend(sim):    
    # Compute diagnosed imbalances at linear end
    # For run_algorithm2, compute relative norms of Pg(u,v,h)
    
    # start with u,v,h from soln
    sim.ramp_soln.h[:,:,:] = sim.soln.h[:,:,:]
    sim.ramp_soln.u[:,:,:] = sim.soln.u[:,:,:]
    sim.ramp_soln.v[:,:,:] = sim.soln.v[:,:,:]
  
    # evolve nonliner system to linear-end
    run_backward_scheme(sim)
    # apply linear BC for relative comparison
    set_linear_bc(sim)              
    

def balance_wave_L1orGeo(sim):
	# Set L1-balanced model and also calculate norms
	
	# initialize ramp function
    sim.ramp_soln = Solution(sim.Nx,sim.Ny,sim.Nz)
    if sim.prognostic:                         
        sim.ramp_solnp = Solutionprog(sim.Nx,sim.Ny,sim.Nz)
    
    # save L1- or geo-balanced fields in ramped solution
    sim.ramp_soln.h[:,:,:] = sim.soln.h[:,:,:]
    if sim.optbal_algorithm == 'algorithmL1':
        print ("===== L1 rebalancing ======" )
        sim.set_L1u(sim)
    elif sim.optbal_algorithm == 'algorithmGeo':
        print ("===== geo rebalancing ======" )
        sim.set_geostrophy(sim)
        
    # update (q,d,g)
    if sim.prognostic:
        sim.classic_to_prognostic(sim)  
    # calculate norms
    sim.compute_norm(sim)    
    
    
# TODO: norm of u is updated for jet flow
def imbalances_nonlinear(sim):                        
    if sim.base_point == "q":      
        # compute imb array including all norms 
        sim.imb = array([sim.norm_q, sim.norm_delta, sim.norm_gamma, 
                         sim.norm_u, sim.norm_v, 
                         sim.norm_uv, sim.norm_h ])
        sim.iters = ravel(sim.iters)                    
        print ("")
        if sim.norm_type == 'absolute':
            print ("Diagnosed imbalances (absolute):")
        elif sim.norm_type == 'relative':
            print ("Diagnosed imbalances (relative):")
        print ("      (u,v,h)    : ", sim.norm_u, sim.norm_v, sim.norm_uv, sim.norm_h)
        print (" (q,delta,gamma) : ", sim.norm_q, sim.norm_delta, sim.norm_gamma )
        print ("iterations :", sim.iters)
        
    elif sim.base_point == "h":
        # compute imb array including all norms 
        sim.imb = array([sim.norm_q, sim.norm_delta, sim.norm_gamma, 
                         sim.norm_u, sim.norm_v, 
                         sim.norm_uv, sim.norm_h ])
        print ("")
        if sim.norm_type == 'absolute':
            print ("Diagnosed imbalances (absolute):")
        elif sim.norm_type == 'relative':
            print ("Diagnosed imbalances (relative):")
        print ("     (u,v,h)     : ", sim.norm_u, sim.norm_v, sim.norm_uv, sim.norm_h)
        print (" (q,delta,gamma) : ", sim.norm_q, sim.norm_delta, sim.norm_gamma )
        print ("iterations :", sim.iters)
    
        
def imbalances_linear(sim):
    # prognostic quantities are added....
    # since norm is computed at linear_end, no effect of balance state.
    sim.imb = array([sim.norm_q, sim.norm_delta, sim.norm_gamma, 
                     sim.norm_u, sim.norm_v, sim.norm_h, 
                     sim.norm_qg, sim.norm_deltag, sim.norm_gammag, 
                     sim.norm_ug, sim.norm_vg, sim.norm_hg,
                     sim.norm_qr, sim.norm_deltar, sim.norm_gammar, 
                     sim.norm_ur, sim.norm_vr, sim.norm_hr]) 
    print ("")
    if sim.norm_type == 'absolute':
        print ("Imbalances at linear-end (absolute):")
    elif sim.norm_type == 'relative':
        print ("Imbalances at linear-end (relative):")
    print ("     (u,v,h)     : ", sim.norm_u, sim.norm_v, sim.norm_h)
    print (" (q,delta,gamma) : ", sim.norm_q, sim.norm_delta, sim.norm_gamma )
    print ("    (u,v,h)_g    : ", sim.norm_ug, sim.norm_vg, sim.norm_hg   )
    print ("(q,delta,gamma)_g: ", sim.norm_qg, sim.norm_deltag, sim.norm_gammag  )
    print ("    (u,v,h)_r    : ", sim.norm_ur, sim.norm_vr, sim.norm_hr   )
    print ("(q,delta,gamma)_r: ", sim.norm_qr, sim.norm_deltar, sim.norm_gammar    )
    
    
