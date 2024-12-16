import sys
import numpy as np

def basecoeff(ibase):
    
    nbc = np.copy(0)
    A_in=np.zeros((4,5), dtype=np.float)
    
    if ibase[6] > 0:                          # if we have a CNM
        A_in[:,0]=[0.00,0.0001,-2.,2.1]       # set the starting value and limits of the 0th-order start at 0
    else:
        A_in[:,0]=[1.,0.001,0.,2.1]        # no CNM: set the starting value and limits of the 0th-order start at 1        
    
    nbc = nbc+1
    
    if ibase[0] > 0:
        A_in[:,1]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the first-order A_in
        nbc = nbc+1
        
    if ibase[0] > 1:
        A_in[:,2]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the second-order A_in
        nbc = nbc+1
         
    if ibase[0] > 2:
        A_in[:,3]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the second-order A_in
        nbc = nbc+1       

    if ibase[0] > 3:
        A_in[:,4]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the second-order A_in
        nbc = nbc+1      
        
    # B coeff => AM:    B[0]*AM + B[1]*AM^2  
    B_in=np.zeros((4,2), dtype=np.float)
    if ibase[1] > 0:
        B_in[:,0]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the first-order B_in
        nbc = nbc+1
        
    if ibase[1] > 1:
        B_in[:,1]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the second-order B_in
        nbc = nbc+1
     
    #  C coeff => XY:    C[0]*x + C[1]*x^2 + C[2]*y + C[3]*y^2 + C[4]*xy 
    C_in=np.zeros((4,5), dtype=np.float)
    if ibase[2] > 0:
        C_in[:,0]=[0.,0.0001,-1.e7,1.e7]  # set the starting value and limits of the first-order C_in
        C_in[:,2]=[0.,0.0001,-1.e7,1.e7]  # set the starting value and limits of the first-order C_in
        nbc = nbc+2
        
    if ibase[2] > 1:
        C_in[:,1]=[0.,0.00001,-1.e7,1.e7]  # set the starting value and limits of the second-order C_in    
        C_in[:,3]=[0.,0.00001,-1.e7,1.e7]  # set the starting value and limits of the second-order C_in    
        C_in[:,4]=[0.,0.0001,-1.e7,1.e7]   # set the starting value and limits of the second-order C_in   
        nbc = nbc+3
    
    # D coeff => fwhm:  D[0]*fwhm + D[1]*fwhm^2
    D_in=np.zeros((4,2), dtype=np.float)
    if ibase[3] > 0:
        D_in[:,0]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the first-order D_in
        nbc = nbc+1
        
    if ibase[3] > 1:
        D_in[:,1]=[0.,0.001,-1.e7,1.e7]  # set the starting value and limits of the second-order D_in
        nbc = nbc+1
        
    # E coeff => sky:   E[0]*sky + E[1]*sky^2    
    E_in=np.zeros((4,2), dtype=np.float)
    if ibase[4] > 0:
        E_in[:,0]=[0.,0.00001,-1.e8,1.e8]    # set the starting value and limits of the first-order E_in
        nbc = nbc+1
        
    if ibase[4] > 1:
        E_in[:,1]=[0.,0.0000001,-1.e8,1.e8]  # set the starting value and limits of the second-order E_in
        nbc = nbc+1
    # E coeff => sky:   E[0]*sky + E[1]*sky^2    

    #       G coeff => sin:   G[0]*np.sin(G[1]*ts+G[2])
    G_in=np.zeros((4,3), dtype=np.float)
    if ibase[5] > 0:
        G_in[:,0]=[0.0001,0.0001,0,1]  # set the starting value and limits of the sinus amplitude
        G_in[:,1]=[50.,0.1,2.5,333]  # set the starting value and limits of the sinus frequency (between 30 min and 10h)
        G_in[:,2]=[np.pi,np.pi/40.,0,2.*np.pi]  # set the starting value and limits of the sinus offset (between 0 min and 2pi)
        nbc = nbc+3
                    
    # H coeff => CNM                
    H_in=np.zeros((4,2), dtype=np.float)
    if ibase[6] > 0:
        H_in[:,0]=[1.,0.001,0,1.e8]  # set the starting value and limits of the first-order H_in
        nbc = nbc+1
        
    if ibase[6] > 1:
        H_in[:,1]=[0.,0.0001,-1.e8,1.e8]  # set the starting value and limits of the second-order H_in       
        nbc = nbc+1
                    
    return A_in, B_in, C_in, D_in, E_in, G_in, H_in, nbc
        


