import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt6.QtWidgets import QApplication
from PyQt6 import QtGui, QtWidgets
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from misc import vector_euler_angles, rgba_0_1_to_0_255, R_yxy, R_zyz

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

        self.lw = 1.5
        self.global_font = QtGui.QFont('Helvetica', 12)
        self.global_font.setWeight(75)
        # azi - z
        # elevation - y
        # yz-plane
        # self.setCameraPosition(elevation=30,azimuth=45)# first elevation, then azimuth
        # yx-plane
        # self.setCameraPosition(elevation=90,azimuth=-90)
        
        ## Add a grid to the view
        self._initGrid()
        self.addGrid()

        self.setCmap('viridis')
        self.addSphere()

        self._initAxis()
        self.addAxis()     

        self._initPOALine()
        self.addPOALine()
        self.setPOA(pos=(0,0,1),scale=1.2)

        self.cones = {}
        self._initGoniometerLines()
        self.addGoniometerLines()
        # self.addGoniometerCircle()

        self._initScatterCone()
        self.addScatterCone()
        self.updateScatterCone(12,0)


    def addSphere(self):
        # create a meshdata sphere
        self.md = gl.MeshData.sphere(rows=20,cols=20,radius=1)
        self.xyz = self.md.vertexes() # (n,3)

        # create a mesh item
        self.mi = gl.GLMeshItem(meshdata=self.md, smooth=False, shader=None)
        self.addItem(self.mi)

    def getXyz(self):
        """return the sphere xyz coordinates in DanMAX convention (y vertical, z along beam)"""
        return self.transform_GL2DM(self.xyz)

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

        # v = v/np.sqrt(np.sum(v**2))
        v = self.transform_DM2GL(v)/np.sqrt(np.sum(v**2))
        # get the euler angles for the cone axis vector
        omega, chi, phi = vector_euler_angles(v)
        # get the euler rotation matrix
        rot = R_yxy(omega, chi, phi).T
        # rotate the cone from the z-axiz to the vector position
        pos = (rot@pos.T).T # (n,3)

        center_line = gl.GLLinePlotItem(pos=np.stack([np.array([0,0,0]),v]),color=color,width=self.lw)
        self.addItem(center_line)
        ring_line = gl.GLLinePlotItem(pos=pos,color=color,width=self.lw)
        self.addItem(ring_line)

        # font = QtGui.QFont('Helvetica', 10)
        text_item = gl.GLTextItem(**{'text':label,
                                    'color':rgba_0_1_to_0_255(color), # (R,G,B,A) 0-255
                                    'pos':v,
                                    'font': self.global_font,
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
        
    # def addGoniometerCircle(self):
    #     # chi axis (in GL coordinates)
    #     # start at x (1,0,0), rotate towards -z (0,0,-1)
    #     radius = 1.5
    #     color = (255,255,0,153)
    #     _chi = np.linspace(0,2*np.pi,100)
    #     pos = np.zeros((_chi.shape[0],3))
    #     pos[:,0] = np.cos(_chi)*radius
    #     pos[:,2] = -np.sin(_chi)*radius
    #     self.goniometer_circle_chi = gl.GLLinePlotItem(pos=pos,color=color)
    #     self.addItem(self.goniometer_circle_chi)

    # def updateGoniometerCircle(self,omega,chi,phi):
    #     radius = 1.5
    #     color = (255,255,0,153)
    #     _chi = np.linspace(0,2*np.pi,100)
    #     pos = np.zeros((_chi.shape[0],3))
    #     pos[:,0] = np.cos(omega*np.pi/180)*np.cos(_chi)*radius
    #     pos[:,1] = np.sin(omega*np.pi/180)*np.cos(_chi)*radius
    #     pos[:,2] = -np.sin(_chi)*radius

    #     self.goniometer_circle_chi.setData(pos=pos)

    def _initGrid(self):
        self.grid = gl.GLGridItem()
        self.grid.scale(.1,.1,.1)
        self.grid.setDepthValue(1)# draw grid after surfaces since they may be translucent
    
    def addGrid(self):
        self.addItem(self.grid)

    def removeGrid(self):
        self.removeItem(self.grid)

    def _initAxis(self):
        self.axis = gl.GLAxisItem()
        self.axis.setSize(2,2,2)

        # add labels to the xyz axes
        # GL: x,y,z -> DM: z,x,y
        font = QtGui.QFont('Helvetica', 16)
        font.setBold(False)
        # font.setPointSize(12)
        axis_labels = {'x': {'text':'z',
                            'color':(0,0,255,153),
                            'pos':(2,0,0),
                            'font': font,
                            },
                       'y': {'text':'x',
                            'color':(255,255,0,153),
                            'pos':(0,2,0),
                            'font': font,
                            },
                       'z': {'text':'y',
                            'color':(0,255,0,153),
                            'pos':(0,0,2),
                            'font': font,
                            },
                    }
        self.axis_text = {}
        for label in axis_labels:
            self.axis_text[label] = gl.GLTextItem(**axis_labels[label])

    def addAxis(self):
        self.addItem(self.axis)
        for item in self.axis_text.values():
            self.addItem(item)

    def removeAxis(self):
        self.removeItem(self.axis)
        for item in self.axis_text.values():
            self.removeItem(item)

    def _initPOALine(self):
        # add POA
        pos = np.zeros((50,3))
        pos[:,0] = np.linspace(0,1,pos.shape[0])
        self.POA_line = gl.GLLinePlotItem(pos=pos,color=(1,1,1,.6),width=self.lw,mode='lines')
        self.POA_label = gl.GLTextItem(text='POA',
                                       color=(255,255,255,153),
                                       pos=(0,0,1),
                                       font= self.global_font,#QtGui.QFont('Helvetica', 10),
                                       )

    def addPOALine(self):
        self.addItem(self.POA_line)
        self.addItem(self.POA_label)

    def removePOALine(self):
        self.removeItem(self.POA_line)
        self.removeItem(self.POA_label)

    def setPOA(self,pos,scale=1.2):
        """set preferred orientation axis (DanMAX coordinates)"""
        pos = np.array(pos)
        pos = self.transform_DM2GL(pos)/np.sqrt(np.sum(pos**2))*scale
        pos=np.array([np.linspace(-x,x,100) for x in pos]).T
        self.POA_line.setData(pos=pos)#np.stack([-pos,pos]))
        self.POA_label.setData(pos=pos[-1,:])

    def _initGoniometerLines(self):
        pos = np.array([[0,0,0],[0,0,0]])
        self.goniometer_lines = {'omega':gl.GLLinePlotItem(pos=pos,color=(1. ,1. ,0. ,0.4),width=self.lw),
                                 'chi':gl.GLLinePlotItem(pos=pos  ,color=(0. ,1. ,0. ,0.4),width=self.lw),
                                 'phi':gl.GLLinePlotItem(pos=pos  ,color=(0. ,0. ,1. ,0.4),width=self.lw),
                                }

    def addGoniometerLines(self):
        for line in self.goniometer_lines.values():
            self.addItem(line)

    def removeGoniometerLines(self):
        for line in self.goniometer_lines.values():
            self.removeItem(line)

    def updateGoniometerLines(self,omega,chi,phi):
        radius = 1.3
        radius_increment = 0.05
        # omega
        _omega = np.append(0,np.arange(0,omega+1,1))*np.pi/180
        pos = np.zeros((_omega.shape[0],3))
        pos[:,0] = np.cos(_omega)
        pos[:,1] = np.sin(_omega)
        pos[0] *= radius+radius_increment
        pos[1:] *= radius
        self.goniometer_lines['omega'].setData(pos=pos)

        # chi
        radius -= radius_increment
        _chi = np.append(0,np.arange(0,chi+1,1)*np.pi/180)
        pos = np.zeros((_chi.shape[0],3))
        pos[:,0] = np.cos(omega*np.pi/180)*np.cos(_chi)
        pos[:,1] = np.sin(omega*np.pi/180)*np.cos(_chi)
        pos[:,2] = -np.sin(_chi)
        pos[0] *= radius+radius_increment
        pos[1:] *= radius
        self.goniometer_lines['chi'].setData(pos=pos)
        
        # phi
        radius -= radius_increment
        rot = R_zyz(omega,chi,0)
        _phi = np.append(0,np.arange(0,phi+1,1)*np.pi/180)
        pos = np.zeros((_phi.shape[0],3))
        pos[:,0] = np.cos(_phi)
        pos[:,1] = np.sin(_phi)
        pos = (rot@pos.T).T
        pos[0] *= radius+radius_increment
        pos[1:] *= radius
        self.goniometer_lines['phi'].setData(pos=pos)

    def _initScatterCone(self):
        pos = np.array([[0,0,0],[1,0,0]])
        self.scatter_line = gl.GLLinePlotItem(pos=pos,color=(1,0,0,.6),width=self.lw)
        self.scatter_ring = gl.GLLinePlotItem(pos=pos,color=(1,0,0,.6),mode='lines',width=self.lw)
        pos = np.zeros((50,3))
        pos[:,0] = np.linspace(0,1,pos.shape[0])
        self.scatter_transmitted = gl.GLLinePlotItem(pos=pos,color=(1,0,0,.6),mode='lines',width=self.lw)
        self.scatter_S0_label = gl.GLTextItem(text='S'+'\u2080',  #'S₀',
                                      color=(255,0,0,153),
                                      pos=(-1.3,0,0),
                                      font=self.global_font,# QtGui.QFont('Helvetica', 10),
                                      )
        
        self.scatter_S_label = gl.GLTextItem(text='S',
                                     color=(255,0,0,153),
                                     pos=(1.,0,0),
                                     font=self.global_font,# QtGui.QFont('Helvetica', 10),
                                     )
        

    def addScatterCone(self):
        self.addItem(self.scatter_line)
        self.addItem(self.scatter_transmitted)
        self.addItem(self.scatter_ring)
        self.addItem(self.scatter_S0_label)
        self.addItem(self.scatter_S_label)

    def removeScatterCone(self):
        self.removeItem(self.scatter_line)
        self.removeItem(self.scatter_transmitted)
        self.removeItem(self.scatter_ring)
        self.removeItem(self.scatter_S0_label)
        self.removeItem(self.scatter_S_label)

    def updateScatterCone(self,tth,eta):
        tth, eta = tth*np.pi/180, eta*np.pi/180
        pos = np.stack(self._coneRing(tth)).T
        pos=pos[:,[2,0,1]]
        pos[:,0]=pos[:,0]*np.sign(np.cos(tth))
        self.scatter_ring.setData(pos=pos)
        S = [np.cos(tth),
             -np.cos(eta)*np.sin(tth),
             -np.sin(eta)*np.sin(tth),
            ]
        pos = np.array([[-1.3,0,0],
                        [0,0,0],
                        S,
                        ])
        self.scatter_line.setData(pos=pos)
        self.scatter_S_label.setData(pos=S)



    def closeEvent(self,event):
        app = QtWidgets.QApplication.instance()
        app.closeAllWindows()
        
    def removeCone(self,label):
        if not label in self.cones:
            return
        cone = self.cones.pop(label)
        for item in ['center_line','ring_line','text_item']:
            self.removeItem(cone[item])

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
    
    def transform_GL2DM(self,a):
        """transform an array from OpenGL coordinates to DanMAX coordinates"""
        axes = [1,2,0]
        return self._transform(a,axes)
    
    def transform_DM2GL(self,a):
        """transform an array from DanMAX coordinates to OpenGL coordinates"""
        axes = [2,0,1]
        return self._transform(a,axes)

    def _transform(self,a,axes):
        if len(a.shape)==1:
            return a[axes]
        if len(a.shape)==2:
            return a[:,axes]
        if len(a.shape)==3:
            return a[:,:,axes]


def main():
    app = QApplication(sys.argv)
    window = SpherePlotWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()


