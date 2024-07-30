import numpy as np
from crystal import Crystal
   
class FibreSample():
    def __init__(self,xyz=np.zeros((1,3))):
        # initialize variables
        self.omega = 0
        self.chi = 0
        self.phi = 0

        # azimuthal angle
        etas = np.arange(0,360,1,dtype=np.float32) # origin at 3 o'clock,cw
        deltas = np.linspace(0,360,180,dtype=np.float32)
        self.etas, self.deltas = np.meshgrid(etas,deltas)

        self.S0 = np.array([0.,0.,1.],dtype=np.float32) # incident beam direction
        self.POI = (0,0,1) # preferred orientation hkl index
        self.hkl = (0,1,1) # hkl of interest
        #KNN (Amm2)a = 4.009, b = 5.591, c = 5.632, α = 90.00° β = 90.00° γ = 90.00°
        self.UC = np.array([4.009, 5.591, 5.632, 90, 90, 90]) #a,b,c,alpha,beta,gamma
        self.SG_centering = 'A'

        self.hkl_sets = {'1 0 0':np.array([1,0,0])}

        self.hkl_max = 12

        # uniform 1D ODF
        self._xyz = xyz
        self.xyz = xyz
        self.MD_r = 1. # March-Dollase r
        self.ODF = np.ones(xyz.shape)

        self.energy = 35 #keV

        self.crystal = Crystal()
        self.crystal.setUC(*self.UC)

        self._updateD()
        self._updateTth()
        self._updateQs()
        self._updateReciprocalBasis()
        self._updateHKLsets()

        self._updateRot()
        self._updatePOA()
        self._updatePODangles()
        self._updateXyz()
        self._updateAlphas()
        self._updateODF()

        self._updatePODs()

    def setEnergy(self,E):
        """set the X-ray energy in keV"""
        self.energy = E
        self._updateEVariables()

    def setUC(self,UC):
        """set the unit cell parameters a,b,c,alpha,beta,gamma (angstrom and degrees)"""
        self.UC = np.array(UC) #a,b,c,alpha,beta,gamma
        self.crystal.setUC(*UC)
        self._updateUCVariables()

    def setPOI(self,hkl):
        """set the preferred orientation hkl index"""
        self.POI = hkl
        self._updateUCVariables()

    def setCentering(self,cen):
        """set the space group centering 'P', 'F', 'I', 'A', 'B', or 'C'"""
        self.SG_centering = cen
        self._updateUCVariables()

    def setHKLmax(self,hkl_max):
        """set the maximum h,k,l index"""
        self.hkl_max = hkl_max
        self._updateUCVariables()

    def setEulerAngles(self,euler):
        """set the Euler angles euler = (omega, chi, phi) in degrees"""
        self.omega, self.chi, self.phi = euler
        self._updateOVariables()

    def setMDparam(self,r):
        """set March-Dollase parameter(s)"""
        self.MD_r = r
        self._updateODFVariables()

    def setHKL(self,hkl):
        self.hkl = hkl
        self._updateUCVariables()

    def getHKLsets(self):
        return self.hkl_sets
    
    def getPODangles(self):
        return self.POD_angles

    def getQs(self):
        return self.Qs
    
    def getEtas(self):
        return self.etas[0]

    def getEulerAngles(self):
        """return the current Euler angles (omega, chi, phi) in degrees """
        euler = (self.omega, self.chi, self.phi)
        return euler
    
    def getPOA(self):
        return self.POA
    
    def getAlphas(self):
        return self.alphas
    
    def getODF(self):
        return self.ODF
    
    def getI_azis(self):
        return self.I_azis 
    
    def importCIF(self,fname):
        self.crystal = Crystal(fname,energy=self.energy)
        return self.crystal.getUC(), self.crystal.getCentering()

    def _updateEVariables(self):
        """_update energy-dependent variables""" 
        self._updateTth()
        self._updateQs()
        self._updatePODs()
        self._updateI_azis()

    def _updateUCVariables(self):
        """update unit cell-dependent variables"""
        self._updateD()
        self._updateTth()
        self._updateQs()
        self._updateHKLsets()
        self._updateReciprocalBasis()
        self._updatePODangles()
        self._updatePODs()
        self._updateI_azis()

    def _updateOVariables(self):
         """update orientation-dependent variables"""
         self._updateRot()
         self._updatePOA()
         self._updateXyz()
         self._updateAlphas()
         self._updateODF()
         self._updateI_azis()

    def _updateODFVariables(self):
        """update ODF-dependent variables"""
        self._updateODF()
        self._updateI_azis()

    def _updateD(self):
        """return d-spacing based on the current hkl and UC"""
        d = self._hkl2d(self.hkl,*self.UC)
        self.d = d
        return d
    
    def _updateTth(self):
        """return tth based on the current d-spacing"""
        d = self.d
        tth = self._d2tth(d,self.energy)
        self.tth = tth
        return tth

    def _updateQs(self):
        """return (n,3) Q vectors for the n azimutal angles"""
        th = self.tth*np.pi/360
        etas = self.getEtas()*np.pi/180
        # calculate Q magnitude
        Q = self._tth2Q(self.tth,self.energy)
        # calculate vectors around the scherrer ring
        # starting at (-x,0,-z)
        Qs = Q * np.array([-np.cos(th)*np.cos(etas),
                           -np.cos(th)*np.sin(etas),
                           -np.sin(th)*np.ones(etas.shape)]).T
        self.Qs = Qs
        return Qs
    
    def _updateHKLsets(self):
        """return a dictionary of hkl indices with similar d-spacing"""
        # hkl_sets = self._getHKLPODsets(self.UC,self.SG_centering,6)
        hkl_sets = self._getHKLsets(self.UC,self.SG_centering,self.hkl_max)
        self.hkl_sets = hkl_sets
        return hkl_sets
    
    def _updateReciprocalBasis(self):
        """return the reciprocal basis transform matrix"""
        a,b,c,alpha,beta,gamma = self.UC
        alpha,beta,gamma = alpha*np.pi/180 ,beta*np.pi/180 ,gamma*np.pi/180
        B_r = self._reciprocalBasis(a,b,c,alpha,beta,gamma)
        self.B_r = B_r
        return B_r

    def _updateRot(self):
        """return the rotation matrix for the current euler angles"""
        omega, chi, phi = self.omega, self.chi, self.phi
        rot = self._R_yxy(omega, chi, phi ).T
        self.rot = rot
        return rot

    def _updatePOA(self):
        """return the preferred oriantation axis for the current euler angles"""
        POA = self.rot.T@self.S0
        self.POA = POA
        return POA
    
    def _updatePODangles(self):
        """
        return the angle between the current hkl and the preferred orientation hkl
        Angle in degrees
        """
        # convert to string
        hkl = ' '.join([str(h) for h in self.hkl])
        if not hkl in self.hkl_sets:
            hkl = list(self.hkl_sets.keys())[0]
        
        hkl_pod_sets = self._getHKLPODsets(self.hkl_sets[hkl]) # {hkl: [np.array(hkl, ...),POD_angle]}
        POD_angles = {key: [hkl_pod_sets[key][1], hkl_pod_sets[key][0].shape[0]*2] for key in hkl_pod_sets}
        self.POD_angles = POD_angles # {hkl:[POD_angle],m}
        return POD_angles
    
    def _updatePODs(self):
        """return a dictionary of (n,m,3) POD vectors expressed in laboratory coordinates"""
        theta = self.tth*np.pi/360
        etas = self.etas*np.pi/180
        deltas = self.deltas*np.pi/180
        
        st, ct = np.sin(theta), np.cos(theta)
        se, ce = np.sin(etas), np.cos(etas)
        sd, cd = np.sin(deltas), np.cos(deltas)
        POD = {}
        for key in self.POD_angles:
            psi, m = self.POD_angles[key]
            psi = psi*np.pi/180
            sps, cps = np.sin(psi), np.cos(psi)

            POD[key] = [np.array([se*sps*cd - st*ce*sps*sd - ct*ce*cps,
                            -ce*sps*cd - st*se*sps*sd - ct*se*cps,
                            ct*sps*sd - st*cps,
                            ]).T,
                        m]
        self.PODs = POD
        return POD
    
    def _updateI_azis(self):
        """return a dictionary of azimuthal intensity variation for the current reflections"""
        S0 = self.S0
        omega, chi, phi = self.omega, self.chi, self.phi
        rot = self._R_yxy(omega, chi, phi ).T
        POA = rot.T@S0
        I_azis = {}
        # hkls =  np.array([key.split() for key in self.PODs_],dtype=int)
        # F2s = self.crystal.getF2(hkls)
        # print(hkls)
        # print(F2s)
        for key in self.PODs:
            PODs,m = self.PODs[key]
            alpha = np.arccos(np.round(PODs@POA,6))*180/np.pi
            mdp = self.MDP(alpha,self.MD_r)
            I_azis[key]=[np.trapz(mdp,x=self.deltas[:,0]*np.pi/180,axis=1)/(2*np.pi),
                         m,
                         self.crystal.getF2(np.array(key.split(),dtype=int))]
        self.I_azis = I_azis
        return I_azis

    def _updateXyz(self):
        """return the xyz coordinates (n,3) based on the current rotation matrix"""
        # get the original xyz coordinate (n,3)
        xyz = self._xyz
        # rotate
        xyz = (self.rot@xyz.T).T
        self.xyz = xyz.astype(np.float32)
        return xyz
        
    def _updateAlphas(self):
        """
        return the angle away from the POA in degrees for each xyz coordinate
        """
        xyz = self.xyz
        alphas = np.arccos(xyz[:,2])*180/np.pi
        self.alphas = alphas
        return alphas
        
    def _updateODF(self):
        """
        return the ODF for the current March-Dollase parameter and
        current alpha values
        """
        r = self.MD_r
        alphas = self.alphas
        ODF = self.MDP(alphas,r)
        self.ODF = np.round(ODF,6)
        return ODF

    def MDP(self,a,r):
        """March-Dollase preferred orientation - a in degrees"""
        a = a*np.pi/180
        if r > 1:
            a = np.pi/2-a
        return (1/r*np.sin(a)**2 + r**2*np.cos(a)**2)**(-3/2)

    def getAngleBetweenHKLs(self,hkl_1,hkl_2):
        """return the angle between two hkl-reflections (degrees)"""
        # get the reciprocal basis matrix
        B_r = self.B_r
        # ensure the hkls are numpy arrays
        hkl_1,hkl_2 = np.array(hkl_1),np.array(hkl_2)
        # transform to reciprocal basis and normalize
        v_1, v_2 = np.matmul(hkl_1,B_r) , np.matmul(hkl_2,B_r)
        v_1 /= np.linalg.norm(v_1)
        v_2 /= np.linalg.norm(v_2)
        # calculate the angle
        angle = np.arccos(np.dot(v_1,v_2))*180/np.pi
        if angle-90>0:
            angle =180-angle
        return angle
    
    def hkl2d(self,hkl):
        """return d-spacings for the provided hkl reflections"""""
        d = self._hkl2d(hkl,*self.UC)
        return d
    
    def d2tth(self,d):
        """return 2theta values for the provided d-spacings"""
        tth = self._d2tth(d,self.energy)
        return tth
    
    def hkl2tth(self,hkl):
        """"return 2thete values for the provided hkl reflections"""
        tth = self._d2tth(self._hkl2d(hkl,*self.UC),self.energy)
        return tth

    def _hkl2d(self,hkl,a,b=None,c=None,alpha=90.,beta=90.,gamma=90.):
        """calculate d-spacing based on hkl and UC parameters"""""
        #Giacovazzo, p. 71
        h,k,l = hkl
        if b is None:
            b=a
        if c is None:
            c=a
        alpha, beta, gamma = [x*np.pi/180 for x in [alpha, beta, gamma]]
        ca,cb,cg = np.cos(alpha), np.cos(beta), np.cos(gamma)
        sa,sb,sg = np.sin(alpha), np.sin(beta), np.sin(gamma)
        d_2 = (1-ca**2-cb**2-cg**2+2*ca*cb*cg)**-1 *\
            (h**2/a**2*sa**2 + k**2/b**2*sb**2 + l**2/c**2*sg**2 + 2*k*l/(b*c)*(cb*cg-ca) \
            +2*l*h/(c*a)*(cg*ca-cb) + 2*h*k/(a*b)*(ca*cb-cg) )
        return np.sqrt(1/d_2)
    
    def _keV2A(self,E):
        """Convert keV to Angstrom"""
        try:
            hc = 1.23984198e-06 # planck constant * speed of light (eV m)
            lambd = hc/(E*10**3)*10**10
        except ZeroDivisionError:
            lambd = np.full(E.shape,0.0)
        return lambd

    def _d2tth(self,d,E):
        """
        Convert d-spacing to 2theta. Provide the energy E in keV
        Returns NaN where d<(lambda/2)
        """
        d = np.atleast_1d(d)
        lambd = self._keV2A(E)
        d[d<(lambd/2)] = np.nan
        tth = 2*np.arcsin(lambd/(2*d))*180/np.pi
        if d.shape[0]==1 and len(d.shape)==1:
            return tth[0]
        return tth
        
    def _tth2Q(self,tth,E):
        """Convert 2theta to Q. Provide the energy E in keV"""
        try:
            if len(E)>1:
                E = np.mean(E)
                print(f'More than one energy provided. Using average: {E:.3f} keV')
        except TypeError:
            pass
        try:
            return 4*np.pi*np.sin(tth/2*np.pi/180)/self._keV2A(E)
        except ZeroDivisionError:
            return np.full(tth.shape,0.0)
        
    def _R_yxy(self,omega=0,chi=0,phi=0):
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
    
    def _reciprocal_lattice(self,a,b,c,alpha,beta,gamma):
        """
        return the reciprocal lattice paramters a*, b*, c*, α*, β*, γ*, and V*
        angles in radians
        """
        # Giacovazzo, p. 67
        # calculate reciprocal angles
        ca, cb, cg = np.cos(alpha),np.cos(beta),np.cos(gamma)
        sa, sb, sg = np.sin(alpha),np.sin(beta),np.sin(gamma)

        alpha_r = np.arccos((cb*cg-ca)/(sb*sg))
        beta_r  = np.arccos((ca*cg-cb)/(sa*sg))
        gamma_r = np.arccos((ca*cb-cg)/(sa*sb))
        # calculate reciprocal volume
        V_r = 1/(a*b*c*np.sin(alpha_r)*sb*sg)
        # calculate reciprocal unit cell lengths
        a_r = b*c*sa*V_r
        b_r = a*c*sb*V_r
        c_r = a*b*sg*V_r
        return a_r, b_r, c_r, alpha_r, beta_r, gamma_r, V_r
    
    def _triclinic2orthogonal(self,a,b,c,alpha,beta,gamma):
        """
        Return an orthonormal basis for the lattice vectors
        Angles in radians
        """
        # Giacovazzo, p. 75
        a_r, b_r, c_r, alpha_r, beta_r, gamma_r, V_r = self._reciprocal_lattice(a,b,c,alpha,beta,gamma)
        
        # triclinic basis
        A = np.array([[a , 0. , 0. ],
                      [0. , b , 0. ],
                      [0. , 0. , c ]])
        
        # orthogonal basis
        B = np.array([[1/a_r, -1/np.tan(gamma_r)/a_r , a*np.cos(beta)  ],
                      [0    , 1/(b_r*np.sin(gamma_r)), b*np.cos(alpha) ],
                      [0    ,          0             ,     c           ]])
        
        # orthonormal basis
        #A_inv = np.linalg.inv(A) #np.true_divide(1,A,where=A!=0)
        #B = np.matmul(A_inv,B)
        #np.round(B,12)
        return B
    
    def _metricMatrix(self,a,b,c,alpha,beta,gamma):
        """return the metric matrix. Angles in radians"""
        # get the orthonormal basis
        B = self._triclinic2orthogonal(a,b,c,alpha,beta,gamma)
        # calculate the metric matrix
        G = np.matmul(B,B.T)
        return G
    
    def _reciprocalBasis(self,a,b,c,alpha,beta,gamma):
        """return the reciprocal basis. Angles in radians"""
        B = self._triclinic2orthogonal(a,b,c,alpha,beta,gamma)
        # calculate the metric matrix
        G = np.matmul(B,B.T)
        B_r = np.matmul(np.linalg.inv(G),B)
        return B_r

    def _applyExtinctionRules(self,hkl,centering='P'):
        hkl = np.atleast_2d(hkl)
        if centering == 'P':
            hkl = hkl
        elif centering == 'I':
            # h+k+l = even
            hkl = hkl[np.sum(hkl,axis=1)%2==0]
        elif centering == 'A':
            # k + l = even
            hkl = hkl[np.sum(hkl[:,1:],axis=1)%2==0]
        elif centering == 'B':
            # h + l = even
            hkl = hkl[np.sum(hkl[:,[0,2]],axis=1)%2==0]
        elif centering == 'C':
            # h + k = even
            hkl = hkl[np.sum(hkl[:,:2],axis=1)%2==0]
        elif centering == 'F':
            # h, k, l all odd or all even
            hkl = hkl[np.sum(hkl%2==0,axis=1)%3==0]
        return hkl

    def _getHKLsets(self,UC,SG_centering='P',h_max=10,k_max=None,l_max=None,deci=1):
        """
        return a dictionary of hkl indices with similar 2theta for the given unit cell
            parameters
                UC                - list of unit cell parameters [a,b,c,alpha,beta,gamma] (angles in degrees)
                SG_centering      - space group centering (default 'P')
                h_max,k_max,l_max - maximum h,k, and l indices. (default h=k=l=10)
            return
                hkl_sets          - dictionary of hkl indices with similar tth {'1 0 0': np.array([[1,0,0],[0,1,0], ...] )}
        """
        if k_max is None:
            k_max = h_max
        if l_max is None:
            l_max = h_max
        h = np.arange(0,h_max+1,1)
        k = np.append(np.arange(0,k_max+1,1),-np.arange(0,k_max+1,1)[1:])
        l = np.append(np.arange(0,l_max+1,1),-np.arange(0,l_max+1,1)[1:])

        # generate all hkl combinations where h>=0
        hkls = np.array(np.meshgrid(k,l,h)).T.reshape(-1,3)[1:,[2,0,1]]
        # remove combinations where all hkl<=0
        hkls = hkls[np.any(hkls>0,axis=1)]
        # apply space group centering extinction rules
        hkls = self._applyExtinctionRules(hkls,SG_centering)
        # calculate d-spacings
        d = self._hkl2d(hkls.T,*UC)
        d_allowed = d>(self._keV2A(self.energy)/2)
        hkls = hkls[d_allowed]
        d = d[d_allowed]
        # calculate rounded tth values
        tths = np.round(self._d2tth(d,self.energy),deci)

        unique_tth, unique_i = np.unique(tths,return_index=True)
        # find unique d-spacings and indices
        # unique_d, unique_i = np.unique(d,return_index=True)
        # reverse order (descending d)
        # unique_d, unique_i = unique_d[::-1], unique_i[::-1]
        # find unique hkls an convert to string
        unique_hkls = hkls[unique_i].astype(str)
        # create a dictionary of similar hkls
        # hkl_sets = {' '.join(hkl):hkls[d==unique_d[i]] for i,hkl in enumerate(unique_hkls)}
        hkl_sets = {' '.join(hkl):hkls[tths==unique_tth[i]] for i,hkl in enumerate(unique_hkls)}
        return hkl_sets
    
    def _getHKLPODsets(self,hkls):
        
        PODs = np.round(np.array([self.getAngleBetweenHKLs(self.POI,hkl) for hkl in hkls]),2)
        unique_PODs, unique_i = np.unique(PODs,return_index=True)
        unique_hkls = hkls[unique_i].astype(str)
        hkl_pod_sets = {' '.join(hkl):[hkls[PODs==unique_PODs[i]],unique_PODs[i]] for i,hkl in enumerate(unique_hkls)}
        
        # PRINT
        # for key in hkl_pod_sets:
        #     print(key)
        #     for hkl in hkl_pod_sets[key]:
        #         print('  ',hkl)
        #     print()
        return hkl_pod_sets
