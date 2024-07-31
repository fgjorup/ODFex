import numpy as np

def vector_euler_angles(v):
    """return the euler angles for a vector v with respect to z=(0,0,1)"""
    v = v/np.sqrt(np.sum(v**2))
    omega = np.arctan2(v[0],v[2])*180/np.pi
    chi = np.arctan2(-v[1],np.sqrt(np.sum((v**2)[[0,2]])))*180/np.pi
    return omega, chi, 0.
        
def rgba_0_1_to_0_255(rgba):
    """convert a tuple of float [0-1] to a tuple of int [0-255]"""
    rgba = tuple([round(x*255) for x in rgba])
    return rgba

def R_yxy(omega=0,chi=0,phi=0):
    """
    YXY rotation matrix - Angles in degrees
    Relates the sample coordinates A and laboratory coordinates A_L, such that RA_L = A
    The transpose of R is equal to the inverse of R: R_T = R^-1 and RR_T=E, where E is the identity matrix.
    """
    o = omega*np.pi/180
    c = chi*np.pi/180
    p = phi*np.pi/180
    
    
    co,cc,cp = np.cos(o),np.cos(c),np.cos(p)
    so,sc,sp = np.sin(o),np.sin(c),np.sin(p)
    
    R = np.array([[ co*cp-cc*so*sp , sc*sp , -cp*so-co*cc*sp],
                  [ so*sc          , cc    , co*sc          ],
                  [ co*sp+cc*cp*so , -cp*sc, co*cc*cp-so*sp ]])
    return(R)

def R_zyz(omega=0,chi=0,phi=0):
    """
    ZYZ rotation matrix - Angles in degrees
    The transpose of R is equal to the inverse of R: R_T = R^-1 and RR_T=E, where E is the identity matrix.
    """
    o = omega*np.pi/180
    c = chi*np.pi/180
    p = phi*np.pi/180
    
    
    co,cc,cp = np.cos(o),np.cos(c),np.cos(p)
    so,sc,sp = np.sin(o),np.sin(c),np.sin(p)
    

    R_z1 = np.array([[ co,-so,  0],
                    [ so, co,  0],
                    [  0,  0,  1]])
    
    R_y1 = np.array([[ cc,  0, sc],
                    [  0,  1,  0],
                    [-sc,  0, cc]])
    
    R_z2 = np.array([[cp,-sp, 0],
                    [sp, cp, 0],
                    [ 0,  0, 1]])

    R = np.dot(R_z1,np.dot(R_y1,R_z2))
    return(R)