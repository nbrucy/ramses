import gc
import osyris
import numpy as np
import os


def check_res(h_over_r=0.1,alpha=0.1):

    inner_boundary = 0.5
    outer_boundary = 3.
    r0             = 1.

    mydata = osyris.RamsesDataset(path='.',nout=1)
    mydata.load()

    disk = np.where( np.logical_and( mydata.get("r",only_leafs=True)/r0 >= inner_boundary , \
                                         mydata.get("r",only_leafs=True)/r0 <= outer_boundary  ) )

    max_dx_over_r = max(mydata.get('dx',only_leafs=True)[disk]/mydata.get('r',only_leafs=True)[disk])
    dt_visc       = min(mydata.r.values**(-3./2.)*mydata.dx.values**2/(mydata.info["gamma"]*mydata.pressure.values/mydata.density.values*4.*alpha))

    print('spatial res. : max(dx/r) = ')
    print(max_dx_over_r)
    print('must be smaller than')
    print(4*h_over_r**2.)
    print('timestep in sim =')
    print(mydata.info['dtold'])
    print('much be smaller than viscous timescale: ')
    print(dt_visc)

def print_1D_data(nout=1,nr=100,disk_radius=1.,H_over_R=0.1,alpha=0.1,inner_boundary=0.,outer_boundary=3., write=False):
#plot_lisa.print_1D_data(nout=1,nr=100,disk_radius=1.,H_over_R=0.1,alpha=0.1)
    mydata = osyris.RamsesDataset(path='.',nout=nout)
    mydata.load()

    disk_density   = 1.0
    if inner_boundary==0.:
        inner_boundary = min(mydata.dx.values)

    t = mydata.info['time'] # c.u.
    t_Porb = int(t / (2.*np.pi))
    t_nu = outer_boundary**1.5  / (3./2. * alpha * H_over_R**2. )
    print(f"t_nu[outer_boundary]= {t_nu}")

    dr = (outer_boundary-inner_boundary)/nr
    range_r = np.arange(inner_boundary,outer_boundary,dr)
    density = np.zeros(nr) ; vr = np.zeros(nr) ; vphi = np.zeros(nr)

    mydata.new_field(name="volume",operation="dx**2.", unit="[L]^2",label=r"V_{cell}") #c.u  
    mydata.new_field(name="vr",operation="(x*velocity_x+y*velocity_y)/r",unit="",label="vr")
    mydata.new_field(name="vphi",operation="(x*velocity_y-y*velocity_x)/r",unit="",label="vphi")

    for i in range(nr):
        r_min = range_r[i]
        if i < nr - 1 :
            r_max = range_r[i+1]
        else :
            r_max = range_r[i] + dr

        ring = np.where( np.logical_and( mydata.get("r",only_leafs=True) >= r_min , \
                                         mydata.get("r",only_leafs=True) <= r_max  ) )
        density_ring = mydata.get('density',only_leafs=True)[ring]
        vr_ring      = mydata.get('vr',only_leafs=True)[ring]
        vphi_ring    = mydata.get('vphi',only_leafs=True)[ring]
        volume_ring  = mydata.get('volume',only_leafs=True)[ring]

        density[i] = np.nansum(density_ring*volume_ring)/np.nansum(volume_ring)
        vr[i] = np.nansum(vr_ring*volume_ring)/np.nansum(volume_ring)
        vphi[i] = np.nansum(vphi_ring*volume_ring)/np.nansum(volume_ring)

        density_t0_r0 = disk_density
        density[i] = density[i] / density_t0_r0
        
    result = {'r':range_r,'density':density,'vr':vr,'vphi':vphi}

    #write results to file 
    if write:
        file = "t"+str(t_Porb)+"Porb_H_over_R"+str(H_over_R)+"_"+"alpha"+str(alpha)+"_nr"+str(nr)+"_1D.dat" #automatize name 
        f = open(file,"w")
        #f.write('#r, density, vr, vphi \n')
        for i in range(nr):
            f.write(str(range_r[i]) + ","  )
            f.write(str(density[i]) + ","  )
            f.write(str(vr[i])      + ","  )
            f.write(str(vphi[i])    + "\n" )
        f.close()
    
    return result


def timeseries(start=1,nend=1,nstep=1,disk_radius=1.,gravity_params=[0.9999,0.0,4.,4.,4.,1e-4,1.],soft_secondary=0.06):
#plot_lisa.timeseries(start=1,nend=401,nstep=10) #1em1
#plot_lisa.timeseries(start=1,nend=401,nstep=10,disk_radius=1.,gravity_params=[0.9999,0.0,4.,4.,4.,1e-4,1.],soft_secondary=0.018) #3em2

    # Time arrays
    time           = np.zeros((int(nend/nstep)))
    x_secondary    = np.zeros((int(nend/nstep)))
    y_secondary    = np.zeros((int(nend/nstep)))
    torque         = np.zeros((int(nend/nstep)))
    torque_cu      = np.zeros((int(nend/nstep)))
    torque_hill    = np.zeros((int(nend/nstep)))
    torque_primary = np.zeros((int(nend/nstep)))
    peak_density   = np.zeros((int(nend/nstep)))
    diskmass_cu    = np.zeros((int(nend/nstep)))

    torque_sec_andrea = np.zeros((int(nend/nstep)))

    disk_density   = 1.0
    r0             = disk_radius # same notation as in RAMSES : in units of half the boxsize
    inner_boundary = 0.5
    outer_boundary = 3

    # gravity params
    gmass      = gravity_params[0] #gmass1 actually
    gmass2     = gravity_params[5]
    emass1     = gravity_params[1] # code units
    separation = gravity_params[6] # c.u.
    emass2     = soft_secondary

    omega      = np.sqrt((gmass+gmass2) / separation**3) # Keplerian rotation speed, in c.u.
    Porb       = 2.*np.pi/omega
    
    fact1      = gmass2 / (gmass + gmass2)
    fact2      = gmass  / (gmass + gmass2)

    # loop over outputs
    for nout in range(start,nend,nstep):

        mydata = osyris.RamsesData(path='.',nout=nout,center=[0.5,0.5,0.],verbose=True)
    
        t = mydata.info['time'] # c.u. 
        dxmin = np.nanmin(mydata.get('dx',only_leafs=True))

        # position of BH1 and BH2
        xmass1 = - fact1 * separation * np.cos(omega * t)
        ymass1 = - fact1 * separation * np.sin(omega * t)
        rmass1 =   np.sqrt(xmass1**2. + ymass1**2.)
        
        xmass2 =   fact2 * separation * np.cos(omega * t)
        ymass2 =   fact2 * separation * np.sin(omega * t)
        rmass2 =   np.sqrt(xmass2**2. + ymass2**2.)

        #--- torque on the secondary and on primary in code units ---#
        mydata.new_field(name="mass2_to_x",operation="x-"+str(xmass2)) # c.u
        mydata.new_field(name="mass2_to_y",operation="y-"+str(ymass2)) # c.u
        mydata.new_field(name="mass1_to_x",operation="x-"+str(xmass1)) # c.u
        mydata.new_field(name="mass1_to_y",operation="y-"+str(ymass1)) # c.u

        mydata.new_field(name="phi_cell",operation="np.arctan2(y,x)",unit="",label=r"")

        mydata.new_field(name="distance_to_mass2",operation="(mass2_to_x**2. + mass2_to_y**2.)**0.5") #c.u
        mydata.new_field(name="distance_to_mass1",operation="(mass1_to_x**2. + mass1_to_y**2.)**0.5") #c.u 

        mydata.new_field(name="acc_from_mass2", operation=str(gmass2)+ \
                         "*distance_to_mass2 / (distance_to_mass2**2. + ("+str(emass2)+")**2.)**(3./2.) ",\
                         unit="[L][T]^-2",label=r"a_{secondary}" ) #c.u
        mydata.new_field(name="acc_from_mass1", operation=str(gmass)+ \
                         "*distance_to_mass1 / (distance_to_mass1**2. + ("+str(emass1)+")**2.)**(3./2.) ",\
                         unit="[L][T]^-2",label=r"a_{primary}" ) #c.u

        mydata.new_field(name="force_from_mass2", operation="acc_from_mass2*density*dx**2.", \
                         unit="[M][L][T]^-2",label=r"F_{secondary}" ) #c.u
        mydata.new_field(name="force_from_mass1", operation="acc_from_mass1*density*dx**2.", \
                         unit="[M][L][T]^-2",label=r"F_{primary}" ) #c.u

        mydata.new_field(name="force_phi", operation="force_from_mass2*(mass2_to_y*"\
                         "np.cos("+str(omega * t)+") - mass2_to_x*np.sin("+str(omega * t)+")  )/distance_to_mass2", \
                         unit="[M][L][T]^-2",label=r"F_{2,phi}" ) #c.u                                 
        mydata.new_field(name="force_1_phi", operation="force_from_mass1*(mass1_to_y*"\
                         "-np.cos("+str(omega * t)+") - mass1_to_x*-np.sin("+str(omega * t)+")  )/distance_to_mass1", \
                         unit="[M][L][T]^-2",label=r"F_{1,phi}" ) #c.u

        mydata.new_field(name="torque_sec", operation=str(rmass2)+"*force_phi",unit="[M][L]^2[T]^-2",label=r"T") #c.u
        mydata.new_field(name="torque_primary", operation=str(rmass1)+"*force_1_phi",unit="[M][L]^2[T]^-2",label=r"T") #c.u

        disk = np.where( np.logical_and( mydata.get("r",only_leafs=True) >= inner_boundary , \
                                         mydata.get("r",only_leafs=True) <= outer_boundary  ) ) 
        rH   = (gmass2/3.)**(1./3.)*separation
        hill = np.where( np.logical_and( np.logical_and( mydata.get("r",only_leafs=True) >= inner_boundary , \
                                                         mydata.get("r",only_leafs=True) <= outer_boundary  ),
                                         mydata.get("distance_to_mass2",only_leafs=True) > rH) )
        ring_secondary = np.where( np.logical_and( mydata.get("r",only_leafs=True) >= rmass2 - dxmin/2. , \
                                         mydata.get("r",only_leafs=True) <= rmass2 + dxmin/2.  ) )

        mydata.new_field(name="mass", operation="density*dx**2.", \
                         unit="[M]",label=r"M_{cell}") #c.u

        torque_cu_nout   = np.nansum(mydata.get('torque_sec',only_leafs=True)[disk])
        torque_hill_nout = np.nansum(mydata.get('torque_sec',only_leafs=True)[hill])
        torque_primary_nout = np.nansum(mydata.get('torque_primary',only_leafs=True)[disk])
        #--- torque on the secondary and on primary in code units ---#

        diskmass_cu_nout = np.nansum(mydata.get('mass',only_leafs=True)[disk])
        #azim_avg_density = np.nanmean(mydata.get('density',only_leafs=True)[ring_secondary])
        peak_density_cu_nout = np.nanmax(mydata.get('density',only_leafs=True)[disk])/disk_density#azim_avg_density
        

        time[int(nout/nstep)]            = t/Porb
        torque_cu[int(nout/nstep)]       = torque_cu_nout
        torque_hill[int(nout/nstep)]     = torque_hill_nout
        torque_primary[int(nout/nstep)]  = torque_primary_nout
        diskmass_cu[int(nout/nstep)]     = diskmass_cu_nout
        peak_density[int(nout/nstep)]    = peak_density_cu_nout


        #--- torques using Andrea's code ---#
        fg_phi_sec = 0.
        
        disk = np.where( np.logical_and( mydata.get("r",only_leafs=True) >= inner_boundary , \
                                         mydata.get("r",only_leafs=True) <= outer_boundary  ) )
        x = mydata.get("x",only_leafs=True)[disk]
        y = mydata.get("y",only_leafs=True)[disk]

        r_arr     = np.sqrt(x**2 + y**2)
        phi_arr   = np.arctan2(y,x)
        dens_arr  = mydata.get("density",only_leafs=True)[disk]
        dV_ramses = (mydata.get("dx",only_leafs=True)[disk])**2.
        mgas      = dens_arr * dV_ramses

        # update for location of secondary in different files - phi = 0 or pi?, r = 1 or 1/(1+q)
        q       = 1.e-4
        r_pri   =  q/(1.+q)
        r_sec   = 1./(1.+q)
        phi_sec = omega * t
        phi_pri = phi_sec + np.pi

        # update for other codes that smooth potentials differently?
        eps_pri = 0.
        eps_sec = emass2#0.6 * r_sec / mach

        dx = r_arr*np.cos(phi_arr) - r_sec*np.cos(phi_sec)
        dy = r_arr*np.sin(phi_arr) - r_sec*np.sin(phi_sec)
            
        # distance from fluid element to planet
        script_r = np.sqrt(dx**2 + dy**2)
            
        f = gmass2*script_r/(script_r**2 + eps_sec**2)**(3./2.)
            
        cosa = dx/script_r ; sina = dy/script_r
        
        #cosa_p = cosa*np.cos(phi_sec) + sina*np.sin(phi_sec)
        sina_p = sina*np.cos(phi_sec) - cosa*np.sin(phi_sec)
        
        fg_phi_sec_i = sina_p * f #Eq. 7 overleaf

        fg_phi_sec = np.sum(mgas*fg_phi_sec_i)
        torque_sec_andrea[int(nout/nstep)] = fg_phi_sec

        print('torque [raph formula]=',torque_cu_nout)
        print('torque [andrea formula]=',fg_phi_sec)
        #--- torques using Andrea's code ---#

        del mydata
        gc.collect()
    
    #write results to files
    #timeseries file obtained with my method
    file = "mignonrisse_ramses_l12_q1e-4_hr01.dat"
    f = open(file,"w")
    f.write('#time[orb secondary],torque_secondary,torque_hill,torque_primary,peak_density[azim avg] \n')
    for nout in range(start,nend,nstep):
        f.write(str(time[nout/nstep])           + "," )
        f.write(str(torque_cu[nout/nstep])      + "," )
        f.write(str(torque_hill[nout/nstep])    + "," )
        f.write(str(torque_primary[nout/nstep]) + "," )
        f.write(str(peak_density[nout/nstep])   + "\n")
    f.close()
    
    #timeseries file obtained with andrea's method
    file = "mignonrisse_torquecheck_ramses_l12_q1e-4_hr01.dat"
    f = open(file,"w")
    f.write('#time[Porb],torque_sec[c.u.] \n')
    for nout in range(start,nend,nstep):
        f.write(str(time[nout/nstep])               + "," )
        f.write(str(torque_sec_andrea[nout/nstep])  + "\n")
    f.close()

    #file with some basic diagnostics (e.g. disk mass)
    file = "datatests_l12_q1e-4_hr01.dat"
    f = open(file,"w")
    f.write('#time[Porb],x2,y2,mass \n')
    for nout in range(start,nend,nstep):
        f.write(str(time[nout/nstep])           + "," )
        f.write(str(x_secondary[nout/nstep])    + "," )
        f.write(str(y_secondary[nout/nstep])    + "," ) 
        f.write(str(diskmass_cu[nout/nstep])       + "\n")
    f.close()


#returns a data file with r, \Sigma (azimuthally-averaged), torque (azim. avg and divided by ring volume)
def data_1Dprofiles(nout=1,nr=400,disk_radius=1.,gravity_params=[0.9999,0.0,4.,4.,4.,1e-4,1.],soft_secondary=0.06,H_over_R_suffix="1em1",q_suffix="1em4"):
#plot_lisa.data_1Dprofiles(nout=401,nr=400,disk_radius=1.,gravity_params=[0.9999,0.0,4.,4.,4.,1e-4,1.],soft_secondary=0.06,H_over_R_suffix="1em1")
#plot_lisa.data_1Dprofiles(nout=401,nr=400,disk_radius=1.,gravity_params=[0.9999,0.0,4.,4.,4.,1e-4,1.],soft_secondary=0.018,H_over_R_suffix="3em2")


    disk_density   = 1.0
    r0             = disk_radius # same notation as in RAMSES : in units of half the boxsize
    inner_boundary = 0.5
    outer_boundary = 3

    # gravity params     
    gmass      = gravity_params[0] #gmass1 actually
    gmass2     = gravity_params[5]
    emass1     = gravity_params[1] # code units
    separation = gravity_params[6] # c.u.
    emass2     = soft_secondary

    omega      = np.sqrt((gmass+gmass2) / separation**3) # Keplerian rotation speed, in c.u.
    Porb       = 2.*np.pi/omega
    
    fact1      = gmass2 / (gmass + gmass2)
    fact2      = gmass  / (gmass + gmass2)

    mydata = osyris.RamsesData(path='.',nout=nout,center=[0.5,0.5,0.],verbose=True)
    
    t = mydata.info['time'] # c.u.
    #print "Time[Porb]= ", t/Porb, "Porb=", Porb
    dxmin = np.nanmin(mydata.get('dx',only_leafs=True))
    
    # position of BH1 and BH2
    xmass1 = - fact1 * separation * np.cos(omega * t)
    ymass1 = - fact1 * separation * np.sin(omega * t)
    rmass1 =   np.sqrt(xmass1**2. + ymass1**2.)

    xmass2 =   fact2 * separation * np.cos(omega * t)
    ymass2 =   fact2 * separation * np.sin(omega * t)
    rmass2 =   np.sqrt(xmass2**2. + ymass2**2.)
    
    #--- torque from disk on the secondary and on primary in code units ---#                                              
    mydata.new_field(name="mass2_to_x",operation="x-"+str(xmass2)) # c.u                                                      
    mydata.new_field(name="mass2_to_y",operation="y-"+str(ymass2)) # c.u                                                      
    mydata.new_field(name="mass1_to_x",operation="x-"+str(xmass1)) # c.u                                                      
    mydata.new_field(name="mass1_to_y",operation="y-"+str(ymass1)) # c.u                                                      
    
    
    mydata.new_field(name="distance_to_mass2",operation="(mass2_to_x**2. + mass2_to_y**2.)**0.5") #c.u                        
    mydata.new_field(name="distance_to_mass1",operation="(mass1_to_x**2. + mass1_to_y**2.)**0.5") #c.u                        
    
    mydata.new_field(name="acc_from_mass2", operation=str(gmass2)+ \
                     "*distance_to_mass2 / (distance_to_mass2**2. + ("+str(emass2)+")**2.)**(3./2.) ",\
                     unit="[L][T]^-2",label=r"a_{secondary}" ) #c.u                                                           
    mydata.new_field(name="acc_from_mass1", operation=str(gmass)+ \
                     "*distance_to_mass1 / (distance_to_mass1**2. + ("+str(emass1)+")**2.)**(3./2.) ",\
                     unit="[L][T]^-2",label=r"a_{primary}" ) #c.u

    mydata.new_field(name="force_from_mass2", operation="acc_from_mass2*density*dx**2.", \
                     unit="[M][L][T]^-2",label=r"F_{secondary}" ) #c.u                                                        
    mydata.new_field(name="force_from_mass1", operation="acc_from_mass1*density*dx**2.", \
                     unit="[M][L][T]^-2",label=r"F_{primary}" ) #c.u 

    mydata.new_field(name="force_2_phi", operation="force_from_mass2*(mass2_to_y*"\
                         "np.cos("+str(omega * t)+") - mass2_to_x*np.sin("+str(omega * t)+")  )/distance_to_mass2", \
                         unit="[M][L][T]^-2",label=r"F_{2,phi}" ) #c.u 
    mydata.new_field(name="force_1_phi", operation="force_from_mass1*(mass1_to_y*"\
                         "-np.cos("+str(omega * t)+") - mass1_to_x*-np.sin("+str(omega * t)+")  )/distance_to_mass1", \
                         unit="[M][L][T]^-2",label=r"F_{1,phi}" ) #c.u 

    mydata.new_field(name="torque_sec",     operation=str(rmass2)+"*force_2_phi",unit="[M][L]^2[T]^-2",label=r"T") #c.u
    mydata.new_field(name="torque_primary", operation=str(rmass1)+"*force_1_phi",unit="[M][L]^2[T]^-2",label=r"T") #c.u       
    mydata.new_field(name="mass", operation="density*dx**2.", unit="[M]",label=r"M_{cell}") #c.u 
    mydata.new_field(name="volume",operation="dx**2.", unit="[L]^2",label=r"V_{cell}") #c.u 

    dr = (outer_boundary-inner_boundary)/nr
    range_r = np.arange(inner_boundary,outer_boundary,dr)
    torque = np.zeros(nr) ; density = np.zeros(nr)
    for i in range(nr):
        r_min = range_r[i]
        if i < nr - 1 :
            r_max = range_r[i+1]
        else :
            r_max = range_r[i] + dr

        ring = np.where( np.logical_and( mydata.get("r",only_leafs=True) >= r_min , \
                                         mydata.get("r",only_leafs=True) <= r_max  ) )

        # torque onto BBH for a given ring
        torque_ring  = (mydata.get('torque_sec',only_leafs=True)[ring] + mydata.get('torque_primary',only_leafs=True)[ring])
        volume_ring  = mydata.get('volume',only_leafs=True)[ring]
        density_ring = mydata.get('density',only_leafs=True)[ring]

        torque[i]  = np.nansum(torque_ring) /np.nansum(volume_ring)
        density[i] = np.nanmean(density_ring)
        density_t0_r0 = disk_density
        density[i] = density[i] / density_t0_r0

    #write results to file
    file = "mignonrisse_ramses_l12_"+q_suffix+"_"+H_over_R_suffix+"_"+str(int(t/Porb))+"Porb_1D.dat"#automatize name
    f = open(file,"w")
    f.write('#r, density, torque \n')
    for i in range(nr):
        f.write(str(range_r[i]) + ","  )
        f.write(str(density[i]) + ","  )
        f.write(str(torque[i])  + "\n" )
    f.close()
