import Dans_Diffraction as dif
import numpy as np

class Crystal(dif.Crystal):
    """crystal object. Inherited from Dans_Diffraction.Crystal"""
    def __init__(self, parent=None,energy=35.):
        super().__init__(parent)
        self.Scatter.setup_scatter(scattering_type='xray',
                                   energy_kev=energy,
                                   #scattering_factors='waaskirf',
                                   output=False)
    def getUC(self):
        """get the crystal unit cell parameters: [a,b,c,alpha,beta,gamma]"""
        UC = np.array([self.Cell.a,
                       self.Cell.b,
                       self.Cell.c,
                       self.Cell.alpha,
                       self.Cell.beta,
                       self.Cell.gamma])
        return UC

    def getF2(self,hkl,E=None):
        """return F^2 for a given set of hkls (n,3) at a given X-ray energy (keV)"""
        if E is None:
            E = self.Scatter.get_energy()
        return self.Scatter.xray_dispersion(hkl,E)
    
    def setUC(self,a,b,c,alpha,beta,gamma):
        """set the crystal unit cell parameters: [a,b,c,alpha,beta,gamma]"""
        self.Cell.latt(a,b,c,alpha,beta,gamma)

    def getCentering(self):
        """get the crystal space group centering"""
        return self.Symmetry.spacegroup_name()[0]

