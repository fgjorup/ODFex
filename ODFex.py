import sys
import os
import numpy as np
import pyqtgraph as pg
from PyQt6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QSlider
from PyQt6 import QtCore, QtGui, QtWidgets
import Dans_Diffraction as dif
from text_editor import TextEditor
from textureSample import FibreSample
from spherePlotWindow import SpherePlotWindow
from palettes import get_palette

class ODFexplorer(QMainWindow):
    def __init__(self,icon=QtGui.QIcon()):
        super().__init__()
        self.setWindowTitle('ODF Explorer')
        self.setWindowIcon(icon) #QtGui.QIcon('ODFex_logo.png'))

        # set the default geometry of the main window and sphere plot window
        # based on the current screen geometry, such that a margin of 1/6
        # of the screen height is left free around the two windows
        screen_geometry = QtWidgets.QApplication.instance().primaryScreen().geometry()
        h,w = screen_geometry.height(), screen_geometry.width()
        self.setGeometry(h//6+h//2,    # x pos
                         h//6,         # y pos
                         int(w-5/6*h), # width
                         2*h//3)       # height
        
        sphere_window_geometry = (h//6,h//6,h//2,h//2) # x,y,width,height

        self.setAcceptDrops(True)

        self.cif_file_name = 'default.cif'
        
        # Create a central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Create a vertical layout for the central widget
        layout = QVBoxLayout(central_widget)

        # Create a 3D surface plot
        self.sphere_window = SpherePlotWindow(custom_geometry=sphere_window_geometry,icon=icon)
        # layout.addWidget(self.sphere_widget)
        self.sphere_window.show()
        
        # xyz = self.sphere_window.xyz
        xyz = self.sphere_window.getXyz()
        
        self.sample = FibreSample(xyz)
        self.sample.setMDparam(0.6)
        self.sample.setHKL([1,1,1])
        self.sample.setEnergy(35)

        POA = self.sample.getPOA()
        ODF = self.sample.getODF()

        tth = self.sample.tth

        # self.sphere_window.addCone('S0',psi=90-(tth/2),v=-self.sample.S0,color=(1,1,1,.5))

        self.sphere_window.setSphereValues(ODF)
        self.sphere_window.setPOA(POA)
        
        self.text_edit = TextEditor(fname=self.cif_file_name,icon=QtGui.QIcon('ODFex_logo.png'))
        self.text_edit.fileChanged.connect(self.editor_file_changed)
        

        self.invalid_cif_message = QtWidgets.QErrorMessage()
        
        # Create a 1D line plot widget
        self.line_widget = pg.PlotWidget()
        self.line_widget.setLabel('left', 'intensity')
        self.line_widget.setLabel('bottom', 'eta (°)')
        self.line_widget.setWindowTitle('1D Line Plot')
        self.line_widget.showGrid(x=True, y=True)
        self.line_widget.setLimits(minYRange=1,
                                   xMin=-10,
                                   xMax=370,
                                   yMin=-0.1)
        layout.addWidget(self.line_widget,1)
        
        self.line_widget.getAxis('bottom').setTickSpacing(30,1)

        # Create a line plot item
        self.line_plot = self.line_widget.plot([0], [0], pen='r')
        self.line_plot_dot = self.line_widget.plot([0], [0], brush='w',symbol='o')

        self.use_dark = True
        self.apply_palette()

        self._init_layout(layout)
        self.updatePOI()
        self.updateUC()
        self.updateSGcentering(self.combo_centering.currentText())
        self.initI_azi_plots()
        self.updateI_azis()
        
        self.updateEuler()
        self.updateCone()

        self.importCIF()

    def _init_layout(self,layout):
        
        # Create menu bar
        menubar = self.menuBar()
        edit_action = QtGui.QAction("Edit CIF", self)
        edit_action.triggered.connect(self.open_editor)
        menubar.addAction(edit_action)

        # create layouts
        sub_panel_layout = QtWidgets.QHBoxLayout()
        input_layout = QtWidgets.QVBoxLayout()
        spin_1_layout = QtWidgets.QHBoxLayout()
        spin_2_layout = QtWidgets.QHBoxLayout()
        dial_label_layout = QtWidgets.QHBoxLayout()
        dial_layout = QtWidgets.QHBoxLayout()

        # add layouts to main layout
        # layout -> sub_panel
        layout.addLayout(sub_panel_layout)
        # layout -> sub_panel -> input
        sub_panel_layout.addLayout(input_layout,stretch=1)
        # layout -> sub_panel -> input -> spin_1
        input_layout.addLayout(spin_1_layout,stretch=1)
        # layout -> sub_panel -> input -> spin_2
        input_layout.addLayout(spin_2_layout,stretch=1)
        # layout -> sub_panel -> input -> dial
        input_layout.addLayout(dial_label_layout)
        input_layout.addLayout(dial_layout)

        ### layout -> sub_panel -> input -> spin_1
        #import cif
        self.toolbutton_cif = QtWidgets.QToolButton()
        self.toolbutton_cif.setText('Import CIF')
        self.toolbutton_cif.setToolTip('Browse for a CIF file\nAlternatively drag-n-drop a CIF file to import it')
        self.toolbutton_cif.clicked.connect(self.importCIFdialog)
        spin_1_layout.addWidget(self.toolbutton_cif)

        # energy
        label = QtWidgets.QLabel(text='energy')
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        label.setToolTip('X-ray energy (keV)')
        spin_1_layout.addWidget(label,1)
        self.spinbox_E = QtWidgets.QDoubleSpinBox()
        self.spinbox_E.setDecimals(1)
        self.spinbox_E.setMinimum(5)
        self.spinbox_E.setMaximum(120)
        self.spinbox_E.setValue(35) 
        self.spinbox_E.setSingleStep(1)
        self.spinbox_E.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
        self.spinbox_E.valueChanged.connect(self.updateEnergy)
        spin_1_layout.addWidget(self.spinbox_E,1)

        # centering
        label = QtWidgets.QLabel(text='centering')
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        label.setToolTip('Space group centering\nSet to primitive (P) for the most general case')
        spin_1_layout.addWidget(label,1)
        self.combo_centering = QtWidgets.QComboBox()
        self.combo_centering.addItems(['P','F','I','A','B','C'])
        self.combo_centering.setCurrentIndex(2)
        self.combo_centering.currentTextChanged.connect(self.updateSGcentering)
        spin_1_layout.addWidget(self.combo_centering,1)

        # hklmax
        label = QtWidgets.QLabel(text='<html><head/><body><p>hkl<span style=" vertical-align:sub;">max</span></p></body></html>')
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
        label.setToolTip('Maximum h, k, and l index to calculate\nAdjust with caution!')
        spin_1_layout.addWidget(label,1)
        self.spinbox_HKLmax = QtWidgets.QSpinBox()
        self.spinbox_HKLmax.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
        self.spinbox_HKLmax.setMinimum(1)
        self.spinbox_HKLmax.setMaximum(100)
        self.spinbox_HKLmax.setValue(12)
        self.spinbox_HKLmax.valueChanged.connect(self.updateHKLmax)
        spin_1_layout.addWidget(self.spinbox_HKLmax,1)

        # preferred orientation index
        self.spinbox_POI = {key:QtWidgets.QSpinBox() for key in ['h','k','l']}
        self.spinbox_POI['h'].setValue(1)
        self.spinbox_POI['k'].setValue(0)
        self.spinbox_POI['l'].setValue(0)
        
        for key in self.spinbox_POI:
            label = QtWidgets.QLabel(text=key)
            self.spinbox_POI[key].valueChanged.connect(self.updatePOI)
            self.spinbox_POI[key].setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
            self.spinbox_POI[key].setMinimum(-100)
            self.spinbox_POI[key].setMaximum(100)
            label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            label.setToolTip('Preferred orientation '+key+' index')
            spin_1_layout.addWidget(label,1)
            spin_1_layout.addWidget(self.spinbox_POI[key],1)


        spin_1_layout.addStretch(8)

        ### layout -> sub_panel -> input -> spin_2
        # UC input 
        self.spinbox_UC = {key: QtWidgets.QDoubleSpinBox() for key in ['a','b','c','alpha','beta','gamma']}
        for key in self.spinbox_UC:
            label = QtWidgets.QLabel(text=key)
            label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            if key in ['a','b','c']:
                self.spinbox_UC[key].setDecimals(3)
                self.spinbox_UC[key].setMinimum(1)
                self.spinbox_UC[key].setValue(2.866) # alpha-iron  (BCC) 2.866
                self.spinbox_UC[key].setSingleStep(0.1)
                label.setToolTip('Unit cell parameter (Å)')
            else:
                self.spinbox_UC[key].setDecimals(2)
                self.spinbox_UC[key].setMinimum(60)
                self.spinbox_UC[key].setMaximum(120)
                self.spinbox_UC[key].setValue(90)
                label.setToolTip('Unit cell parameter (°)')
            self.spinbox_UC[key].valueChanged.connect(self.updateUC)
            self.spinbox_UC[key].setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
            spin_2_layout.addWidget(label,1)
            spin_2_layout.addWidget(self.spinbox_UC[key],1)

        spin_2_layout.addStretch(8)

        ### layout -> sub_panel -> input -> dial
        # create euler dials
        self.euler_dial_labels = {}
        self.euler_dials = {}
  
        start_val = {'omega':270,'chi':90,'phi':0}
        tool_tips = {'omega':'First Euler angle (Y)',
                    'chi'  :"Second Euler angle (X')",
                    'phi'  :"Third Euler angle (Y'')",
                    }
        for dial in ['omega','chi','phi']:
            self.euler_dials[dial] = QtWidgets.QDial()
            self.euler_dials[dial].setMinimum(0)
            self.euler_dials[dial].setMaximum(359)
            self.euler_dials[dial].setValue(start_val[dial])
            self.euler_dials[dial].setPageStep(1)
            self.euler_dials[dial].setWrapping(True)
            self.euler_dials[dial].setNotchTarget(45)
            self.euler_dials[dial].setNotchesVisible(True)
            self.euler_dials[dial].valueChanged.connect(self.updateEuler)

            self.euler_dial_labels[dial] = QtWidgets.QLabel()
            self.euler_dial_labels[dial].setText(f'{dial}: {self.euler_dials[dial].value():>3d} °')
            self.euler_dial_labels[dial].setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
            self.euler_dial_labels[dial].setToolTip(tool_tips[dial])

            subdial_layout = QtWidgets.QVBoxLayout()
            subdial_layout.addWidget(self.euler_dial_labels[dial])
            subdial_layout.addWidget(self.euler_dials[dial])
            subdial_layout.addStretch()
            dial_layout.addLayout(subdial_layout,stretch=1)

        self.eta_dial = QtWidgets.QDial()
        self.eta_dial.setMinimum(0)
        self.eta_dial.setMaximum(359)
        self.eta_dial.setValue(270)
        self.eta_dial.setWrapping(True)
        self.eta_dial.setNotchTarget(45)
        self.eta_dial.setNotchesVisible(True)
        self.eta_dial.valueChanged.connect(self.updateCone)

        self.eta_dial_label = QtWidgets.QLabel()
        self.eta_dial_label.setText(f'eta: {self.getEtaValue():>3d} °')
        self.eta_dial_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.eta_dial_label.setToolTip('Azimuthal scattering angle')
        
        subdial_layout = QtWidgets.QVBoxLayout()
        subdial_layout.addWidget(self.eta_dial_label)
        subdial_layout.addWidget(self.eta_dial)
        subdial_layout.addStretch()
        dial_layout.addLayout(subdial_layout,stretch=1)
        
        self.slider_MDparam = QSlider(QtCore.Qt.Orientation.Vertical)
        self.slider_MDparam.setMinimum(1)
        self.slider_MDparam.setMaximum(100)
        self.slider_MDparam.setPageStep(1)
        self.slider_MDparam.setValue(90)
        self.slider_MDparam.valueChanged.connect(self.updateMDparam)
        
        self.MDparam_label = QtWidgets.QLabel()
        self.MDparam_label.setText(f'r: {self.slider_MDparam.value()/100:>3.2f}')
        self.MDparam_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.MDparam_label.setToolTip('March-Dollase texture parameter')

        subdial_layout = QtWidgets.QVBoxLayout()
        subdial_layout.addWidget(self.MDparam_label)
        subdial_layout.addWidget(self.slider_MDparam)
        dial_layout.addLayout(subdial_layout,stretch=1)

        ### layout -> sub_panel
        self.treewidget = QtWidgets.QTreeWidget()
        self.treewidget.setHeaderLabels(['h', 'k', 'l', 'd (Å)','2θ (°)','POD (°)'])
        self.treewidget.header().setStretchLastSection(False)
        self.treewidget.header().setSectionResizeMode(self.treewidget.header().ResizeMode.ResizeToContents)
        
        [self.treewidget.headerItem().setTextAlignment(i,QtCore.Qt.AlignmentFlag.AlignHCenter) for i in range(5)]
        self.treewidget.headerItem().setToolTip(0,'h-index')
        self.treewidget.headerItem().setToolTip(1,'k-index')
        self.treewidget.headerItem().setToolTip(2,'l-index')
        self.treewidget.headerItem().setToolTip(3,'d-spacing')
        self.treewidget.headerItem().setToolTip(4,'2θ angle')
        self.treewidget.headerItem().setToolTip(5,'Preferred orientation direction angle')
        self.treewidget.setSortingEnabled(True)
        self.treewidget.setIndentation(10)

        self.treewidget.currentItemChanged.connect(self.setHKL)
        self.treewidget.itemExpanded.connect(self.treeItemExpended)

        self._setHKLtree()
        self.treewidget.setCurrentItem(self.treewidget.topLevelItem(0))

        sub_panel_layout.addWidget(self.treewidget,stretch=3)
        
        # graphic settings check boxes
        checkbox_layout = QtWidgets.QVBoxLayout()
        label = QtWidgets.QLabel()
        label.setText('Display settings')
        checkbox_layout.addWidget(label)
        # goniometer
        self.checkbox_goniometer = QtWidgets.QCheckBox('goniometer')
        self.checkbox_goniometer.setChecked(True)
        self.checkbox_goniometer.stateChanged.connect(self.showGoniometer)
        self.checkbox_goniometer.setToolTip('Show goniometer guide lines')
        # axis
        self.checkbox_axis = QtWidgets.QCheckBox('axis')
        self.checkbox_axis.setChecked(True)
        self.checkbox_axis.stateChanged.connect(self.showAxis)
        self.checkbox_axis.setToolTip('Show xyz axes')
        # beam
        self.checkbox_beam = QtWidgets.QCheckBox('scattered beam')
        self.checkbox_beam.setChecked(True)
        self.checkbox_beam.stateChanged.connect(self.showBeam)
        self.checkbox_beam.setToolTip('Show incident (S'+'\u2080'+') and scattered beam (S)')
        # grid
        self.checkbox_grid = QtWidgets.QCheckBox('grid')
        self.checkbox_grid.setChecked(True)
        self.checkbox_grid.stateChanged.connect(self.showGrid)
        self.checkbox_grid.setToolTip('Show horizontal plane surface grid')
        # add to layout
        checkbox_layout.addWidget(self.checkbox_goniometer)
        checkbox_layout.addWidget(self.checkbox_axis)
        checkbox_layout.addWidget(self.checkbox_beam)
        checkbox_layout.addWidget(self.checkbox_grid)
        checkbox_layout.addStretch(1)
        sub_panel_layout.addLayout(checkbox_layout,stretch=1)

    def showGoniometer(self,val):
        if self.checkbox_goniometer.isChecked():
            self.sphere_window.addGoniometerLines()
        else:
            self.sphere_window.removeGoniometerLines()
    
    def showAxis(self,val):
        if self.checkbox_axis.isChecked():
            self.sphere_window.addAxis()
        else:
            self.sphere_window.removeAxis()

    def showBeam(self,val):
        if self.checkbox_beam.isChecked():
            self.sphere_window.addScatterCone()
        else:
            self.sphere_window.removeScatterCone()

    def showGrid(self,val):
        if self.checkbox_grid.isChecked():
            self.sphere_window.addGrid()
        else:
            self.sphere_window.removeGrid()


    def open_editor(self):
        self.text_edit._load_file(self.cif_file_name)
        self.text_edit.show()
    
    def editor_file_changed(self,fname):
        if dif.functions_crystallography.cif_check(dif.functions_crystallography.readcif(fname)):
            self.cif_file_name = fname
            self.importCIF()
        else:
            self.invalid_cif_message.showMessage('Invalid CIF file')

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls() and len(event.mimeData().urls())==1:
            if event.mimeData().urls()[0].toLocalFile().lower().endswith('.cif'):
                event.accept()
            else:
                event.ignore()
        else:
            event.ignore()

    def dropEvent(self, event):
        fname = event.mimeData().urls()[0].toLocalFile()
        self.cif_file_name = fname
        self.importCIF()

    def importCIFdialog(self):
        defualt_path = ''
        if not self.cif_file_name is None:
            defualt_path = os.path.dirname(self.cif_file_name)
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self,'Load CIF',defualt_path,'CIF file (*.cif)')
        if fname.lower().endswith('.cif') and os.path.isfile(fname):
            self.cif_file_name = fname
        else:
            self.cif_file_name = None
        self.importCIF()
        
    def importCIF(self):
        if self.cif_file_name is None:
            print('No valid CIF file found')
            return
        UC, centering = self.sample.importCIF(self.cif_file_name)
        [self.spinbox_UC[key].blockSignals(True) for key in self.spinbox_UC]
        [self.spinbox_UC[key].setValue(UC[i]) for i,key in enumerate(self.spinbox_UC)]
        [self.spinbox_UC[key].blockSignals(False) for key in self.spinbox_UC]
        self.combo_centering.setCurrentText(centering)
        self.updateUC()
        self.setWindowTitle(f'ODF Explorer - {self.cif_file_name}')

    def treeItemExpended(self,item):
        self._setPODtreeItems(item)
        self.treewidget.setCurrentItem(item)

    def updateUC(self):
        UC = [self.spinbox_UC[key].value() for key in self.spinbox_UC]
        self.sample.setUC(UC)
        self._setHKLtree()
        self.updateI_azis()
        self.updateCone()
        if not self.cif_file_name is None:
            self.setWindowTitle(f'ODF Explorer - {self.cif_file_name}*')
    
    def updateEnergy(self):
        E = self.spinbox_E.value()
        self.sample.setEnergy(E)
        self._setHKLtree()
        self.updateI_azis()
        self.updateCone()

    def updatePOI(self):
        hkl = [self.spinbox_POI[key].value() for key in self.spinbox_POI]
        if hkl == [0,0,0]:
            self.spinbox_POI['h'].setValue(1)
            return
        self.sample.setPOI(hkl)
        self._setHKLtree()
        self.updateI_azis()
        self.updateCone()

    def updateSGcentering(self,value):
        self.sample.setCentering(value)
        self._setHKLtree()
        if not self.cif_file_name is None:
            self.setWindowTitle(f'ODF Explorer - {self.cif_file_name}*')

    def updateHKLmax(self,value):
        self.sample.setHKLmax(value)
        self._setHKLtree()


    def updateMDparam(self,value):
        MDparam = value/100 # convert to float
        self.MDparam_label.setText(f'r: {MDparam:>3.2f}')
        self.sample.setMDparam(MDparam)
        self.updateEuler()
        self.updateCone()

    def initI_azi_plots(self):
        self.line_widget.clear()
        pen = pg.mkPen(color='w',
                       width=2)
        self.line_plots = {'sum':self.line_widget.plot([0],[0],name='sum',pen=pen)}
        x = self.sample.getEtas()
        I_azis = self.sample.getI_azis()
        for i,hkl in enumerate(I_azis):
            pen = pg.mkPen(color=(i,len(self.sample.getI_azis())),
                           style=QtCore.Qt.PenStyle.DashLine,
                           width=2,
                           )
            name = f'{hkl} ({I_azis[hkl][1]})'
            self.line_plots[hkl]=self.line_widget.plot(x,np.ones(x.shape[0]),name=name,pen=pen)
        self.line_widget.addLegend()
        self.line_plot_dot = self.line_widget.plot([0], [0], brush='w',symbol='o')


    def updateI_azis(self):
        x = self.sample.getEtas()
        I_azis = self.sample.getI_azis()
        y_sum, m_sum = 0., 0
        for hkl in I_azis:
            y, m, F2 = I_azis[hkl]
            # self.line_plots[hkl].setData(x,y*m/2)
            self.line_plots[hkl].setData(x,y*m*F2)
            y_sum += y*m*F2
            # m_sum += m
        self.line_plots['sum'].setData(x,y_sum)#/m_sum)
        
        i = self.getEtaValue()
        self.line_plot_dot.setData([i],(y_sum)[i:i+1])
        # self.line_plot_dot.setData([i],(y_sum/m_sum)[i:i+1])

    def updateEuler(self):
        omega,chi,phi = [self.euler_dials[dial].value() for dial in ['omega','chi','phi']]
        for label in self.euler_dial_labels:
            self.euler_dial_labels[label].setText(f'{label} : {self.euler_dials[label].value():>3d} °')
        self.sample.setEulerAngles((omega,chi,phi))

        POA = self.sample.getPOA()
        ODF = self.sample.getODF()
        
        self.sphere_window.setSphereValues(ODF)
        self.sphere_window.setPOA(POA)
        if self.checkbox_goniometer.isChecked():
            self.sphere_window.updateGoniometerLines(omega,chi,phi)
        # self.sphere_window.updateGoniometerCircle(omega,chi,phi)

        self.updateI_azis()


    def updateCone(self):
        # get dial value
        eta = self.getEtaValue()
        self.eta_dial_label.setText(f'eta: {eta:>3d} °')
        Q = self.sample.getQs()[eta]
        pod_hkl, psi = self.getCurrentPOD()
        self.sphere_window.addCone('Q',psi=psi,v=Q)
        y = self.line_plots['sum'].getData()[1]
        if y.shape[0]>1:
            self.line_plot_dot.setData([eta],[y[eta]])
        
        self.sphere_window.updateScatterCone(self.sample.tth,eta)


    def getCurrentPOD(self):
        current = self.treewidget.currentItem()
        hkl, pod_angle = None, 0.
        if not current is None:
            hkl = ' '.join([current.text(i) for i in range(3)]).strip()
            pod_angle = current.text(5)
            if not pod_angle:
                pod_angle = 0.
        return hkl, float(pod_angle)
        
    def getEtaValue(self):
        dial_offset = 90
        value = self.eta_dial.value()
        return (360+value+dial_offset)%360
    

    def _setHKLtree(self):
        self.treewidget.clear()
        hkl_sets = self.sample.getHKLsets() # {hkl:np.array([hkl1,hkl2,...])}
        for hkl in hkl_sets.keys():
            h = [int(h) for h in hkl.split()]
            d = self.sample.hkl2d(h)
            tth = self.sample.d2tth(d)
            if np.isnan(tth):
                continue
            item = QtWidgets.QTreeWidgetItem(self.treewidget)
            for i,val in enumerate(h):
                item.setText(i,f'{val:2d}')
                item.setTextAlignment(i,QtCore.Qt.AlignmentFlag.AlignTrailing)
            item.setText(3,f'{d:6.3f}')
            item.setText(4,f'{tth:6.2f}')
            item.setTextAlignment(3,QtCore.Qt.AlignmentFlag.AlignTrailing)
            item.setTextAlignment(4,QtCore.Qt.AlignmentFlag.AlignTrailing)
            for _hkl in hkl_sets[hkl]:
                if not np.all(_hkl == np.array(hkl.split(),dtype=int)):
                    _item = QtWidgets.QTreeWidgetItem(item)
                    _hkld = list(_hkl)+[]
                    for i,val in enumerate(_hkld):
                        _item.setText(i,f'{val}')
                        _item.setTextAlignment(i,QtCore.Qt.AlignmentFlag.AlignTrailing)
        
        self.treewidget.sortItems(3,QtCore.Qt.SortOrder.DescendingOrder)
        current = self.treewidget.findItems(f'{self.sample.d:.4f}',QtCore.Qt.MatchFlag.MatchContains,3)
        if current:
            self.treewidget.setCurrentItem(current[0])
        self.initI_azi_plots()

    def _setPODtreeItems(self,current):
        items = [current]+[current.child(i) for i in range(current.childCount())]
        for item in items:
            hkl = [int(item.text(i)) for i in range(3)]
            pod_angle = self.sample.getAngleBetweenHKLs(self.sample.POI,hkl)
            item.setText(5,f'{pod_angle:.0f}')
            item.setTextAlignment(5,QtCore.Qt.AlignmentFlag.AlignTrailing)
        current.sortChildren(5,QtCore.Qt.SortOrder.AscendingOrder)
    
    def setHKL(self,current,previous):
        if not current is None:
            if not current.parent() is None:
                current = current.parent()
            hkl = [int(current.text(i)) for i in range(3)]
            self.sample.setHKL(hkl)
            self.initI_azi_plots()
            self.updateCone()
            self.updateI_azis()

    def closeEvent(self,event):
        app = QtWidgets.QApplication.instance()
        app.closeAllWindows()

    def apply_palette(self):
        palette = get_palette(self.use_dark)
        app = QtWidgets.QApplication.instance()
        app.setPalette(palette)
        app.setStyle('Fusion')


def main():
    app = QApplication(sys.argv)
    icon = QtGui.QIcon('ODFex_logo.png')
    odf_ex = ODFexplorer(icon=icon)
    odf_ex.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()