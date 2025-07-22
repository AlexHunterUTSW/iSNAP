###################
# Import Packages #
###################
print('- Importing os...')
from os.path import (
    join
)
print('- Importing PyQt...')
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import ( 
    QWidget, 
    QVBoxLayout, 
    QHBoxLayout, 
    QGridLayout,
    QPushButton, 
    QFileDialog, 
    QListWidget, 
    QMessageBox, 
    QLineEdit, 
    QLabel,
    QCheckBox,
    QComboBox,
    QSpacerItem,
    QSizePolicy,
    QTreeView
)


class SetModality(QWidget):
    def __init__(self):
        super().__init__()

        # Initiate Layouts
        self.layoutMain = QVBoxLayout()
        self.layoutParam = QGridLayout()
        self.layoutDevMode = QHBoxLayout()

        # Create Widgets
        icon = QPixmap('iSNAP_Icon.png')
        resizedIcon = icon.scaled(200, 200, Qt.KeepAspectRatio)
        labelIcon = QLabel('iSNAP Icon')
        labelIcon.setPixmap(resizedIcon)
        labelWelcome = QLabel('Welcome to iSNAP!')
        labelWelcome.setStyleSheet('font-size: 24px;')


        modalities = ['Xenium', 'Cell Ranger', 'Single .h5ad file']
        self.comboModality = QComboBox()
        self.comboModality.addItems(modalities)
        self.labelModality = QLabel('Select Sequencing Software')
        
        self.inSeed = QLineEdit('12345')
        self.labelSeed = QLabel('Random State Seed')

        self.checkDev = QCheckBox('Dev Mode')

        hspacer = QSpacerItem(1, 1, QSizePolicy.Expanding, QSizePolicy.Minimum)
        vspacer = QSpacerItem(1, 1, QSizePolicy.Minimum, QSizePolicy.Expanding)
        
        # Set Layouts
        self.layoutParam.addWidget(self.inSeed, 0, 0, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelSeed, 0, 1)
        self.layoutParam.addWidget(self.comboModality, 1, 0, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelModality, 1, 1)
        self.layoutParam.addItem(hspacer,0,2)
        self.layoutParam.addItem(hspacer,1,2)


        self.layoutMain.addItem(vspacer)
        self.layoutMain.addWidget(labelIcon)
        self.layoutMain.addWidget(labelWelcome)
        self.layoutMain.addItem(vspacer)
        self.layoutMain.addLayout(self.layoutParam)
        self.layoutMain.addItem(vspacer)
        self.layoutMain.addWidget(self.checkDev)

        self.setLayout(self.layoutMain)

        
class FileFinder(QWidget):
    def __init__(self):
        super().__init__()

        self.folderIn_paths = []
        self.folderOut_path = ""

        # Initiate Layouts
        self.layoutMain = QVBoxLayout()
        self.layoutIn = QHBoxLayout()
        self.layoutInBtns = QVBoxLayout()
        self.layoutOut = QHBoxLayout()
        self.layoutConfirm = QHBoxLayout()

        # Create Widgets

        # Displays list of folder names
        self.labelIn = QLabel()
        self.labelIn.setText('Input Directories')
        self.listIn = QListWidget()

        # Add Input Folder Button
        self.selectIn_btn = QPushButton("Add Folder")
        self.selectIn_btn.clicked.connect(self.add_folderIn)

        self.removeIn_btn = QPushButton('Remove Selected')
        self.removeIn_btn.clicked.connect(self.remove_folderIn)

        self.isSpatial = QCheckBox('Spatial Transcriptomics')
        self.isSpatial.setChecked(True)

        #
        self.lableOut = QLabel('Output Folder')
        self.lineOut = QLineEdit()

        # Add Output Folder Button
        self.selectOut_button = QPushButton("Add Folder")
        self.selectOut_button.clicked.connect(self.add_folderOut)

        # Confirm Selection Button
        self.confirm_button = QPushButton("Confirm")
        self.confirm_button.clicked.connect(self.confirm_folder)

        # Add widgets to layout, organize layout
        self.layoutInBtns.addWidget(self.selectIn_btn, alignment = Qt.AlignTop)
        self.layoutInBtns.addWidget(self.removeIn_btn, alignment=Qt.AlignTop)
        self.layoutIn.addWidget(self.listIn) 
        self.layoutIn.addLayout(self.layoutInBtns) # Add Select button on same row as list
        
        self.layoutOut.addWidget(self.lineOut)
        self.layoutOut.addWidget(self.selectOut_button)

        self.setLayout(self.layoutMain)
        
    def setPage(self, modality):
        
        self.clearLayout(self.layoutMain)
        self.layoutMain.addWidget(self.labelIn)
        self.layoutMain.addLayout(self.layoutIn)
        if modality in ['Xenium', 'Single .h5ad file']:
            self.layoutMain.addWidget(self.isSpatial)
        self.layoutMain.addWidget(self.lableOut)
        self.layoutMain.addLayout(self.layoutOut)

    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)

    def add_folderIn(self):
        fileDialog = FileDialog()
        fileDialog.exec()
        files = fileDialog.selectedFiles
        try:
            for folder in files:
                if folder not in self.folderIn_paths and folder != "":
                    self.folderIn_paths.append(folder)
                    self.listIn.addItem(folder)
        except:
            pass
        
    def remove_folderIn(self):
        listItems=self.listIn.selectedItems()
        if not listItems: return        
        for item in listItems:
            self.listIn.takeItem(self.listIn.row(item))
            self.folderIn_paths.remove(item.text())
    
    def add_folderOut(self):
        fileDialog = QFileDialog()
        folder = fileDialog.getExistingDirectory(self, "Select Folder", options=QFileDialog.Option.DontUseNativeDialog)
        self.folderOut_path = folder
        self.lineOut.setText(folder)

    def confirm_folder(self):
        if self.folderIn_paths == [''] or self.folderIn_paths == [] or self.folderOut_path == "":
            error_dialog = QMessageBox()
            error_dialog.setWindowTitle('Error!')
            error_dialog.setText('Please select folders')
            error_dialog.exec()
        else:
            self.close()


class FileDialog(QFileDialog):
    def __init__(self, *args):
        QFileDialog.__init__(self, *args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setFileMode(self.ExistingFiles)
        self.setGeometry(100, 100, 1280, 720)
        btns = self.findChildren(QPushButton)
        self.openBtn = [x for x in btns if 'open' in str(x.text()).lower()][0]
        self.openBtn.clicked.disconnect()
        self.openBtn.clicked.connect(self.openClicked)
        self.tree = self.findChild(QTreeView)
        self.show()

    def openClicked(self):
        inds = self.tree.selectionModel().selectedIndexes()
        files = []
        for i in inds:
            if i.column() == 0:
                files.append(join(str(self.directory().absolutePath()),str(i.data())))
        self.selectedFiles = files
        self.close()