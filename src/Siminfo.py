import numpy as np
import os,sys
               
               
def print_parameter(sim):
    print('Parameters:')
    print('-----------')
    print('stepper      = {0:s}'.format(sim.stepper.__name__))
    print('CFL          = {0:g}'.format(sim.cfl))
    print('Nx           = {0:d}'.format(sim.Nx))
    print('Ny           = {0:d}'.format(sim.Ny))
    print('Nz           = {0:d}'.format(sim.Nz))
    print('Coriolis     = f-plane')    


def run_ramping_print(sim):
    
    # The whole optimal balance method is summarized: 
    # Goal: run the balancing over T_ramp until converging to a balance state.
    
    sim.location = "before assignment:"
    print_solns(sim)
    
    # start with u,v,h from soln
    sim.ramp_soln.h = sim.soln.h.copy()
    sim.ramp_soln.u = sim.soln.u.copy()
    sim.ramp_soln.v = sim.soln.v.copy() 
  
    sim.norm = 1. 
    sim.iteration = 0
    sim.thrs = 1e-4
    while sim.norm > sim.thrs:
        
        print ("")
        print ("........ iteration .........", sim.iteration)
        
        # 1. the sytem starts at nonlinear swe (appr. slow manifold)
        
        sim.location = "before ramping starts:"
        print_solns(sim)
        
        # initialize the settings for backward in time
        sim.ramp_direction = 'Backward'
        sim.ramp_prepare_for_run()
        
        # 2. step the ramped system back in time 
        while sim.ramp_time < sim.T_ramp:
           
            sim.step_ramp()
            #print((sim.ramp_time)     )
  
        # 3. the system reaches at linear swe (slow manifold)
        sim.location = "before setting linear_bc:"
        print_solns(sim)
        
        # 4. Set BC on linear swe end by projection or geostrophy
        sim.set_linear_bc(sim)
        
        sim.location = "after setting linear_bc:"
        print_solns(sim)
       
        
        # initialize the settings for forward in time
        sim.ramp_direction = 'Forward'
        sim.ramp_prepare_for_run()
              
        # 5. step the ramped system forward in time.
        while sim.ramp_time < sim.T_ramp:
          
            sim.step_ramp()
            #print((sim.ramp_time))
         
        # 6. the system reaches at nonlinear swe (appr. slow manifold)           
       
        sim.location = "after forward ramping:"
        print_solns(sim)
        
        # 7. set BC at nonlinear swe end depending on balance_state 
        sim.set_nonlinear_bc(sim)
            
        sim.location = "after nonliner_bc:"    
        print_solns(sim)
        
        sim.iteration += 1
        
        n = 20
        if sim.iteration > n:
            print ("The method does not converge!")
            #for i in range(n):
            #    print((imb[i]))
            sim.norm = 1e-16
            
    sim.location = "after ramping procedure:"
    print_solns(sim)
            
    # increase for second direction
    # 0: backward ramping
    # 1: forward ramping
    sim.nramping += 1
   
    return

def print_solns(sim):
    # to check each step in run_ramping()
    
    if sim.location == "before assignment:":
        print ("")
        print ("before assignment:")
        print ("max_u:", sim.soln.u.max())
        print ("max_v:", sim.soln.v.max())
        print ("max_h:", sim.soln.h.max())
        """print ("max_q:", sim.solnp.q.max())
        print ("max_d:", sim.solnp.delta.max())
        print ("max_g:", sim.solnp.gamma.max())"""
    
    elif sim.location == "before ramping starts:":
        print ("")
        print ("before ramping starts:")
        print ("max_ramp_u:", sim.ramp_soln.u.max())
        print ("max_ramp_v:", sim.ramp_soln.v.max())
        print ("max_ramp_h:", sim.ramp_soln.h.max())
        """print ("max_ramp_q:", sim.ramp_solnp.q.max())
        print ("max_ramp_d:", sim.ramp_solnp.delta.max())
        print ("max_ramp_g:", sim.ramp_solnp.gamma.max())"""
    elif sim.location == "before setting linear_bc:":
        print ("")
        print ("before setting linear_bc:")
        print ("max_ramp_u:", sim.ramp_soln.u.max())
        print ("max_ramp_v:", sim.ramp_soln.v.max())
        print ("max_ramp_h:", sim.ramp_soln.h.max())
        """print ("max_ramp_q:", sim.ramp_solnp.q.max())
        print ("max_ramp_d:", sim.ramp_solnp.delta.max())
        print ("max_ramp_g:", sim.ramp_solnp.gamma.max())"""
    elif sim.location == "after setting linear_bc:":
        print ("")
        print ("after setting linear_bc:")
        print ("max_ramp_u:", sim.ramp_soln.u.max())
        print ("max_ramp_v:", sim.ramp_soln.v.max())
        print ("max_ramp_h:", sim.ramp_soln.h.max())
        """print ("max_ramp_q:", sim.ramp_solnp.q.max())
        print ("max_ramp_d:", sim.ramp_solnp.delta.max())
        print ("max_ramp_g:", sim.ramp_solnp.gamma.max()  )"""
    elif sim.location == "after forward ramping:":
        print ("")
        print ("after forward ramping:")
        print ("max_ramp_soln_u:", sim.ramp_soln.u.max())
        print ("max_ramp_soln_v:", sim.ramp_soln.v.max())
        print ("max_ramp_soln_h:", sim.ramp_soln.h.max())
        """print ("max_ramp_q:", sim.ramp_solnp.q.max())
        print ("max_ramp_d:", sim.ramp_solnp.delta.max())
        print ("max_ramp_g:", sim.ramp_solnp.gamma.max())"""
    elif sim.location == "after nonliner_bc:":
        print ("")
        print ("after nonliner_bc:")
        print ("max_ramp_soln_u:", sim.ramp_soln.u.max())
        print ("max_ramp_soln_v:", sim.ramp_soln.v.max())
        print ("max_ramp_soln_h:", sim.ramp_soln.h.max())
        """print ("max_ramp_q:", sim.ramp_solnp.q.max())
        print ("max_ramp_d:", sim.ramp_solnp.delta.max())
        print ("max_ramp_g:", sim.ramp_solnp.gamma.max())"""
    elif sim.location == "after ramping procedure:":
        print ("")
        print ("after ramping procedure:")
        print ("max_u:", sim.soln.u.max())
        print ("max_v:", sim.soln.v.max())
        print ("max_h:", sim.soln.h.max())
        """print ("max_ramp_q:", sim.ramp_solnp.q.max())
        print ("max_ramp_d:", sim.ramp_solnp.delta.max())
        print ("max_ramp_g:", sim.ramp_solnp.gamma.max())"""

    
def print_step_information(sim):

    h0 = sim.soln.h[:,:,0] - sim.soln.h[:,:,1]
    minh = np.min(np.ravel(h0))
    maxh = np.max(np.ravel(h0))

    maxu = np.max(np.ravel(sim.soln.u[:,:,0]))
    minu = np.min(np.ravel(sim.soln.u[:,:,0]))
    maxv = np.max(np.ravel(sim.soln.v[:,:,0]))
    minv = np.min(np.ravel(sim.soln.v[:,:,0]))
    
    if sim.Nz == 1:             
        pstr  = '{0:s}'.format(Plot_tools.smart_time(sim.time))    
        pstr += ' '*(13-len(pstr))
        try:
            pstr += ',  dt = {0:0<7.1e}'.format(mean(sim.dts))
        except:
            pstr += ',  dt = {0:0<7.1e}'.format(sim.dt)
        L = len(pstr) - 13
        pstr += ', min(u,v,h) = ({0: < 8.4e},{1: < 8.4e},{2: < 8.4e})'.format(minu,minv,minh)
        #pstr += ', del_mass = {0: .2g}'.format(mass/sim.Ms[0]-1)
        pstr += '\n'
        tmp = '  = {0:.3%}'.format(sim.time/sim.end_time)
        pstr += tmp
        pstr += ' '*(L - len(tmp))            
        pstr += 'avg = {0:0<7.1e}'.format(sim.mean_dt)
        pstr += ', max(u,v,h) = ({0: < 8.4e},{1: < 8.4e},{2: < 8.4e})'.format(maxu,maxv,maxh)
        #pstr += ', del_enrg = {0: .2g}'.format(enrg/(sim.KEs[0]+sim.PEs[0])-1)
    else:
        pstr  = 't = {0:.4g}s'.format(sim.time)
        try:
            pstr += ', dt = {0:0<7.1e}'.format(np.mean(sim.dts))
        except:
            pstr += ', dt = {0:0<7.1e}'.format(sim.dt)
        #pstr += ', del_mass = {0:+.2g}'.format(mass/sim.Ms[0]-1)
        #pstr += ', del_enrg = {0:+.2g}'.format(enrg/(sim.KEs[0]+sim.PEs[0])-1)
        pstr += '\n'
        pstr += '  = {0:.3%}'.format(sim.time/sim.end_time)

    head_str = ('\n{0:s}'.format(sim.run_name))
    if sim.animate != 'None':
        head_str += ': frame {0:d}'.format(sim.frame_counter-1)
    if sim.output:
        head_str += ': output {0:d}'.format(sim.output_counter-1)
    head_str += ': step {0:d}'.format(sim.num_steps)
    #print((head_str))
    #print((pstr))
