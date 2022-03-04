import sys
from os.path import expanduser
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QStatusBar
from PyQt5.QtWidgets import QToolBar
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QDirModel, QTreeView,  QFileSystemModel, QFileDialog, QDesktopWidget

import misc_functions

class OpenMenuWindow(QWidget):
    """
        This "window" is a QWidget. If it has no parent, it
        will appear as a free-floating window as we want.
        """

    def __init__(self):
        super().__init__()
        home_directory = expanduser('~')
        layout = QVBoxLayout()

        self.setGeometry(400, 400, 640, 480)
        self.label = QLabel("Another Window")
        layout.addWidget(self.label)

        windowLayout = QVBoxLayout()
        self.setLayout(windowLayout)

        self.show()


class Window2(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Window22222")

class MainWindow(QMainWindow):
    """Main Window."""
    def __init__(self, parent=None):
        """Initializer for the main window."""
        super().__init__(parent)
        self.title = 'PhotometryPipeline'
        centerPoint = QDesktopWidget().availableGeometry().center()
        self.width = 640
        self.height = 480
        # left, top, width, height
        self.openwindow = None  # No external window yet
        self.initUI()

    def initUI(self):
        """initialize user interface of the main window"""
        self.setWindowTitle(self.title)
        self.resize(self.width, self.height)
        self.center()
        self.setCentralWidget(QLabel("I'm the Central Widget"))
        self._createMainMenu()
        self._createConfigMenu()
        self._createToolBar()
        self._createStatusBar()

    def center(self):
        """locate the main window at the center of the screen"""
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def _createMainMenu(self):
        """creates the first pop-up of the main menu"""
        self.menu = self.menuBar().addMenu("&Menu")
        self.menu.addAction('New', self._createNewMainWindow)
        self.menu.addAction('Open Dir', self._selectWorkingDirectory)
        self.menu.addAction('Open',  self._openFileNamesDialog)
        self.menu.addAction('Save', self._file_save)
        self.menu.addAction('Exit', self.close)

    def _createConfigMenu(self):
        """creates the second pop-up of the main menu"""
        self.menu = self.menuBar().addMenu("&Config")
        self.menu.addAction('New', self.close)
        self.menu.addAction('Open', self._openConfigDialog)
        self.menu.addAction('Save', self._file_save)
        self.menu.addAction('Exit', self.close)

    def _createToolBar(self):
        tools = QToolBar()
        self.addToolBar(tools)
        tools.addAction('Exit', self.close)

    def _createStatusBar(self):
        status = QStatusBar()
        status.showMessage("I'm the Status Bar")
        self.setStatusBar(status)

    def _createOpenWindow(self):
        if self.openwindow is None:
            self.openwindow = OpenMenuWindow()
        self.openwindow.show()

    def _createNewMainWindow(self):
        self.newWindow = Window2()
        self.newWindow.show()

    def _openFileNameDialog(self, filetype):
        """"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        if filetype:
            fileName, _ = QFileDialog.getOpenFileName(self, "Open file", "",
                                                      "{};;All Files (*)".format(filetype),
                                                      options=options)
        else:
            fileName, _ = QFileDialog.getOpenFileName(self, "Open file", "",
                                                    "All Files (*);;Python Files (*.py)",
                                                    options=options)
        if fileName:
            print(fileName)

    def _openFileNamesDialog(self):
        """"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        files, _ = QFileDialog.getOpenFileNames(self, "Open file(s)", "",
                                                "All Files (*);;Python Files (*.py)",
                                                options=options)
        if files:
            print(files)

    def _file_save(self):
        """"""
        name = QFileDialog.getSaveFileName(self, 'Save File')
        file = open(name, 'w')
        text = self.textEdit.toPlainText()
        file.write(text)
        file.close()

    def _selectWorkingDirectory(self):
        folderpath = QFileDialog.getExistingDirectory(self, 'Select Project Folder')
        if folderpath:
            print(folderpath)

    def _openConfigDialog(self):
        """"""
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        configName, _ = QFileDialog.getOpenFileName(self, "Open config", "",
                                                      "Config file (*.yml);;All Files (*)",
                                                      options=options)

        if configName:
            print(configName)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())

