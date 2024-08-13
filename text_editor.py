import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QTextEdit, QFileDialog,QMessageBox
from PyQt6.QtGui import QAction, QKeySequence
from PyQt6.QtCore import pyqtSignal

class TextEditor(QMainWindow):

    fileChanged = pyqtSignal(str)

    def __init__(self,fname='',icon=None):
        super().__init__()
        self.init_ui()
        if not icon is None:
            self.setWindowIcon(icon)
        self.current_file = fname
        self._load_file(self.current_file)

    def init_ui(self):
        self.text_edit = QTextEdit(self)
        self.setCentralWidget(self.text_edit)
        self.setWindowTitle('Text editor')

        self.text_edit.textChanged.connect(self.contentChanged)

        # Create actions
        load_action = QAction("Load", self)
        load_action.triggered.connect(self.load_file)
        

        save_action = QAction("Save", self)
        save_action.triggered.connect(self.save_file)
        save_action.setShortcut(QKeySequence('Ctrl+s'))
        
        save_as_action = QAction("Save As", self)
        save_as_action.triggered.connect(self.save_as_file)

        # Create menu bar
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")
        file_menu.addAction(load_action)
        file_menu.addAction(save_action)
        file_menu.addAction(save_as_action)

        self.setGeometry(100, 100, 800, 600)
        self.setWindowTitle("Simple Text Editor")
        # self.show()

    def load_file(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Load File", "", "CIF Files (*.cif)", options=options)
        if file_name:
            self._load_file(file_name)
        self.fileChanged.emit(self.current_file)

    def _load_file(self,file_name):
        if file_name:
            self.current_file = file_name
            with open(file_name, "r") as file:
                self.text_edit.setPlainText(file.read())
            self.setWindowTitle('Text editor - '+self.current_file)

    def save_file(self):
        if self.current_file:
            with open(self.current_file, "w") as file:
                file.write(self.text_edit.toPlainText())
            self.setWindowTitle('Text editor - '+self.current_file)
            self.fileChanged.emit(self.current_file)
            self.text_edit.document().setModified(False)
        else:
            self.save_as_file()

    def save_as_file(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getSaveFileName(self, "Save As", "", "CIF Files (*.cif)", options=options)
        if file_name:
            self.current_file = file_name
            self.save_file()
    
    def contentChanged(self):
        if not self.current_file is None:
            self.setWindowTitle('Text editor - '+self.current_file+'*')

    def closeEvent(self,event):
        if self.text_edit.document().isModified():
            reply = QMessageBox.question(self,
                                         "Unsaved Changes",
                                         "Do you want to save your changes?",
                                          QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            if reply == QMessageBox.Save:
                self.save_file()
            elif reply == QMessageBox.Cancel:
                event.ignore()  # Cancel the close event
            else:
                event.accept()  # No unsaved changes, close the window
        

if __name__ == "__main__":
    app = QApplication(sys.argv)
    editor = TextEditor()
    editor.show()
    sys.exit(app.exec_())