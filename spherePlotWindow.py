import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtWidgets import QApplication
from PyQt5 import QtGui, QtWidgets
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from misc import vector_euler_angles, rgba_0_1_to_0_255, R_yxy

class SpherePlotWindow(gl.GLViewWidget):
    """
    OpenGL view widget window plotting a sphere with colorsurface, used for visualizing the
    orientation distrubution function (ODF) along with other graphical elements, such as 
    scattering cones, coordinate system axes, and vectors.
    """
    def __init__(self, parent=None,custom_geometry=None,icon=QtGui.QIcon()):
        super().__init__(parent)

        if not custom_geometry is None:
            self.setGeometry(*custom_geometry)     
        self.setWindowTitle('ODF Explorer')
        self.setWindowIcon(icon) #QtGui.QIcon('ODFex_logo.png'))
        self.setCameraPosition(distance=5)

        # azi - z
        # elevation - y
        # yz-plane
        # self.setCameraPosition(elevation=30,azimuth=45)# first elevation, then azimuth
        # yx-plane
        # self.setCameraPosition(elevation=90,azimuth=-90)
        
        ## Add a grid to the view
        g = gl.GLGridItem()
        g.scale(.1,.1,.1)
        g.setDepthValue(1)  # draw grid after surfaces since they may be translucent
        self.addItem(g)

        # add axes
        axis = gl.GLAxisItem()
        axis.setSize(2,2,2)
        self.addItem(axis)

        # add labels to the xyz axes
        font = QtGui.QFont('Helvetica', 16)
        font.setBold(False)
        # font.setPointSize(12)
        axis_labels = {'x': {'text':'x',
                            'color':(0,0,255,153),
                            'pos':(2,0,0),
                            'font': font,
                            },
                       'y': {'text':'y',
                            'color':(255,255,0,153),
                            'pos':(0,2,0),
                            'font': font,
                            },
                       'z': {'text':'z',
                            'color':(0,255,0,153),
                            'pos':(0,0,2),
                            'font': font,
                            },
                    }
        for label in axis_labels:
            text = gl.GLTextItem(**axis_labels[label])
            self.addItem(text)
    
        self.setCmap('viridis')
        self.addSphere()

        # add POA
        self.POA_line = gl.GLLinePlotItem(pos=(0,0,1),color=(1,1,1,.6))
        self.POA_label = gl.GLTextItem(text='POA',
                                       color=(255,255,255,153),
                                       pos=(0,0,1),
                                       font= QtGui.QFont('Helvetica', 10),
                                       )
        self.addItem(self.POA_line)
        self.addItem(self.POA_label)
        self.setPOA(pos=(0,0,1),scale=1.2)

        self.cones = {}


    def addSphere(self):

        # create a meshdata sphere
        self.md = gl.MeshData.sphere(rows=20,cols=20,radius=1)
        self.xyz = self.md.vertexes() # (n,3)

        # create a mesh item
        self.mi = gl.GLMeshItem(meshdata=self.md, smooth=False, shader=None)
        self.addItem(self.mi)

    def getCmap(self):
        return self.cmap
    
    def setCmap(self,cmap):
        self.cmap = plt.get_cmap(cmap)

    def getColorsFromScalar(self,values):
        smap = cm.ScalarMappable(cmap=self.cmap)
        colors = smap.to_rgba(values,norm=True)
        return colors

    def setSphereValues(self,values):
        colors = self.getColorsFromScalar(values)
        if len(colors.shape) == 3:
            flat_shape = self.surface_shape[0]*self.surface_shape[1]
            colors = colors.reshape((flat_shape,4))
        self.md.setVertexColors(colors)
        self.mi.setMeshData(meshdata=self.md)

    def addCone(self,label,psi=0.,v=np.array([0,0,1]),color=(1.,1.,1.,1.)):
        """Add a cone ring to the 3D plot"""
        if label in self.cones:
            self.removeCone(label)

        color = tuple(color)
        # calculate a cone ring around z
        x,y,z = self._coneRing(psi*np.pi/180)
        pos=np.stack([x,y,z]).T

        v = v/np.sqrt(np.sum(v**2))
        # get the euler angles for the cone axis vector
        omega, chi, phi = vector_euler_angles(v)
        # get the euler rotation matrix
        rot = R_yxy(omega, chi, phi).T
        # rotate the cone from the z-axiz to the vector position
        pos = (rot@pos.T).T

        center_line = gl.GLLinePlotItem(pos=np.stack([np.array([0,0,0]),v]),color=color)
        self.addItem(center_line)
        ring_line = gl.GLLinePlotItem(pos=pos,color=color)
        self.addItem(ring_line)

        font = QtGui.QFont('Helvetica', 10)
        text_item = gl.GLTextItem(**{'text':label,
                                    'color':rgba_0_1_to_0_255(color), # (R,G,B,A) 0-255
                                    'pos':v,
                                    'font': font,
                                    })
        self.addItem(text_item)
        
        self.cones[label]={'psi':psi,
                          'pos'    : v ,
                          'color':color,
                          'label':label,
                          'center_line':center_line,
                          'ring_line':ring_line,
                          'text_item':text_item,
                          }
        
    def closeEvent(self,event):
        app = QtWidgets.QApplication.instance()
        app.closeAllWindows()
        
    def removeCone(self,label):
        if not label in self.cones:
            return
        cone = self.cones.pop(label)
        for item in ['center_line','ring_line','text_item']:
            self.removeItem(cone[item])
            
    def setPOA(self,pos,scale=1.2):
        """set preferred orientation axis"""
        pos = np.array(pos)
        pos = pos/np.sqrt(np.sum(pos**2))*scale
        self.POA_line.setData(pos=np.stack([-pos,pos]))
        self.POA_label.setData(pos=pos)

    def _coneRing(self,psi):
        h = 1
        r = np.tan(psi)*h
        mag = np.sqrt(h**2+r**2)
        x = np.cos(np.linspace(0,2*np.pi,100))*r
        y = np.sin(np.linspace(0,2*np.pi,100))*r
        z = np.ones(x.shape)
        # mag = np.sqrt(x**2+y**2+z**2)
        x /= mag
        y /= mag
        z /= mag
        return x,y,z

def main():
    app = QApplication(sys.argv)
    window = SpherePlotWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()


