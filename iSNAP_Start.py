print('Welcome to iSNAP!\n')
print('Importing Essential Packages...')

print('- Importing sys...')
from sys import (
    exit,
    argv
)

print('- Importing os...')

from os import (
    makedirs
)
from os.path import (
    exists as pathexists, 
    join
)

print('- Importing pandas...')
from pandas import (
    concat,
)

print('- Importing Numpy...')
from numpy import (
    average,
    argmax,
    ndarray,
    zeros
)

print('- Importing matplotlib...')
from matplotlib.pyplot import (
    close as pltclose
)
from matplotlib.figure import Figure
from matplotlib.image import imread
from matplotlib.cm import Reds

from copy import copy

print('- Importing Scanpy...')
from scanpy.plotting import (
    umap as pltumap,
    rank_genes_groups_dotplot,
    violin
)
from scanpy.tools import (
    rank_genes_groups,
    dendrogram
)
from scanpy.plotting.palettes import default_102


print('- Importing PyQt...')
from PyQt6.QtCore import QCoreApplication, Qt, QObject, QThread, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication,
    QVBoxLayout,
    QHBoxLayout,
    QStackedWidget,
    QWidget,
    QPushButton,
    QScrollArea,
    QGroupBox,
    QSizePolicy,
    QSpacerItem,
    QLabel,
    QCheckBox
)
from qt_material import apply_stylesheet

from PyQt6.QtGui import QFont, QIcon

print('\nImporting Input Packages')
from iSNAP_Input import SetModality, FileFinder
print('\nImporting Reading Packages')
from iSNAP_Read import read_data, err_NoFiles, ReadQC, save_adatas
print('\nImporting Preprocessing Packages')
from iSNAP_PP import (
    calc_qc_metrics, 
    quality_control, 
    normalization, 
    pcaVarianceRatio, 
    integration, 
    cluster,
    silhouetteLeiden, 
    TableNormalize, 
    PCAInt, 
    TestLeiden, 
    SelectLeiden
)
print('\nImporting Analysis Packages')
from iSNAP_Analysis import (
    deg, 
    degtoCSV,
    specificClusterUMAP,
    leidenToCellType, 
    typeToCSV, 
    DEGDotplot, 
    AnnotateCluster, 
    SpecificCluster, 
    UMAPCellType,
    FeatureMap
    )
print('\nImporting Post Analysis Packages')
from iSNAP_Post import (
    TypeSampleTable,
    FinishedWidget, 
    DGEWindow, 
    VolcanoPlot, 
    dgetoCSV, 
    plotVolcano
    ) 

from iSNAP_CellSorter import CellSortGame

###########
# Objects #
###########

# Main Window that displays operation. Acts as communication link between MainWorker object that does the work with the data and the page for the stacked widget.
class MainWindow(QWidget):
        toWorkerNextPG = pyqtSignal(dict, int)          # Signal that outputs dictionary of parameters and the current page for MainWorker's processData function
        toWorkerCluster = pyqtSignal(bool, int, int)    # Signal goes to MainWorker's toCluster Function to set up a specific cluster's page during annotation. Outputs a bool to tell whether to use the main or subcluster data and pages, an integer for how many top DEG to display, and the cluster index to find the cluster page.
        toWorkerAnnotate = pyqtSignal(int, str, bool)   # Signal goes to MainWorker's toAnnotate function to update and return to the annotation page. Outputs the cluster's index, new name, and if it is from the main or subcluster data.
        
        def __init__(self):
            super(MainWindow, self).__init__()
            
            self.setWindowTitle("iSNAP")
            self.setGeometry(100, 100, 1280, 720)
            self.setWindowIcon(QIcon('iSNAP_Icon.png'))
            
            self.devMode = False

            self.h5adImport = False

            self.thread = QThread()
            self.worker = MainWorker()
            self.worker.moveToThread(self.thread)

            self.toWorkerNextPG.connect(self.worker.processData)
            self.toWorkerCluster.connect(self.worker.toCluster)
            self.toWorkerAnnotate.connect(self.worker.toAnnotate)
            self.worker.toNextpgSet.connect(self.setNextpg)
            self.worker.skipSignal.connect(self.createSkipWindow)
            self.worker.toClusterSignal.connect(self.toClusterSetPage)
            self.worker.toAnnotateSignal.connect(self.toAnnotateSetPage)
            self.worker.initiateFeatureSignal.connect(self.createFeatureMap)
            self.worker.updateFeatureSignal.connect(self.updateFeatureMap)
            self.worker.toSubClusterSignal.connect(self.setSubCluster)
            self.worker.toDGESignal.connect(self.createDGE)
            self.worker.toVolSignal.connect(self.createVol)
            self.worker.toVolFailSignal.connect(self.enableDGEbtn)
            
            self.mainbox = QHBoxLayout()
            self.subbox = QVBoxLayout()

            # Assign pages
            self.page0 = SetModality()
            self.page1 = FileFinder()
            self.page2 = ReadQC()
            self.page3 = TableNormalize()
            self.page4 = PCAInt()
            self.page5 = TestLeiden()
            self.page6 = SelectLeiden()
            self.page7 = DEGDotplot()
            self.page8 = AnnotateCluster(True)
            self.page8sub = []
            self.page9 = UMAPCellType()
            self.page10 = [PCAInt(), TestLeiden(), SelectLeiden(), DEGDotplot(), AnnotateCluster(False)]
            self.page15 = TypeSampleTable()
            self.page16 = FinishedWidget()

            # Page to Main Window Signals
            self.page8.toClusterSignal.connect(self.toCluster)
            self.page8.typeToCSVSignal.connect(self.worker.typeToCSVWork)
            self.page8.toFeatureSignal.connect(self.worker.initiateFeatureMap)

            self.page9.typeToCSVSignal.connect(self.worker.typeToCSVWork)
            self.page9.combTypesSignal.connect(self.worker.combTypes)
            self.page9.toFeatureSignal.connect(self.worker.initiateFeatureMap)
            self.page9.toSubClusterSignal.connect(self.worker.toSubCluster)

            self.page10[4].toClusterSignal.connect(self.toCluster)
            self.page10[4].typeToCSVSignal.connect(self.worker.typeToCSVWork)
            self.page10[4].toFeatureSignal.connect(self.worker.initiateFeatureMap)

            self.page15.toDGESignal.connect(self.worker.toDGE)

            # Indices for masked cluster plots during annotation
            self.sub8Index = []
            self.sub10Index = []

            self.err_NoFiles = err_NoFiles()

            # Create Stacked Widget
            self.stackedWidget = QStackedWidget()
            self.stackedWidget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            
            # Add Pages
            self.stackedWidget.addWidget(self.page0)
            self.stackedWidget.addWidget(self.page1)
            self.stackedWidget.addWidget(self.page2)
            self.stackedWidget.addWidget(self.page3)
            self.stackedWidget.addWidget(self.page4)
            self.stackedWidget.addWidget(self.page5)
            self.stackedWidget.addWidget(self.page6)
            self.stackedWidget.addWidget(self.page7)
            self.stackedWidget.addWidget(self.page8)
            self.stackedWidget.addWidget(self.page9)
            for page in self.page10:
                self.stackedWidget.addWidget(page)
            self.stackedWidget.addWidget(self.page15)
            self.stackedWidget.addWidget(self.page16)

            # Set Main Features
            self.layoutPages = QVBoxLayout()
            self.prepPageNames = [
                'Modality',
                'Folders',
                'QC',
                'Normalization',
                'Integration',
                'Test Leiden',
                'Select Leiden',
                'DEG Dotplot',
                'Annotation',
                'Subclustering'
            ]
            
            self.subPageNames = [
                'Integration',
                'Test Leiden',
                'Select Leiden',
                'DEG Dotplot',
                'Annotation',
            ]
            self.endPageNames = [
                'DGE Analysis',
                'Finished'
            ]

            self.tabPage = []
            for i in range(len(self.prepPageNames)):
                self.tabPage.append(TabButton(self, self.prepPageNames[i], i))
                spacer = QSpacerItem(20, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)
                layoutLine=QHBoxLayout()
                layoutLine.addWidget(self.tabPage[i])
                layoutLine.addItem(spacer)

                self.layoutPages.addLayout(layoutLine)
            for i in range(len(self.subPageNames)):
                self.tabPage.append(TabButton(self, self.subPageNames[i], i+len(self.prepPageNames)))
                spacer = QSpacerItem(20, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)
                layoutLine=QHBoxLayout()
                layoutLine.addItem(spacer)
                layoutLine.addWidget(self.tabPage[i+len(self.prepPageNames)])
                layoutLine.addItem(spacer)

                self.layoutPages.addLayout(layoutLine)

            for i in range(len(self.endPageNames)):
                self.tabPage.append(TabButton(self, self.endPageNames[i], i+len(self.prepPageNames)+len(self.subPageNames)))
                spacer = QSpacerItem(20, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)
                layoutLine=QHBoxLayout()
                layoutLine.addWidget(self.tabPage[i+len(self.prepPageNames)+len(self.subPageNames)])
                layoutLine.addItem(spacer)

                self.layoutPages.addLayout(layoutLine)
            

            self.groupNav = QGroupBox('Go to Step')
            self.groupNav.setLayout(self.layoutPages)
            self.scrollNav = QScrollArea()
            self.scrollNav.setWidget(self.groupNav)
            self.scrollNav.setWidgetResizable(True)
            self.scrollNav.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            self.scrollNav.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)

            self.mainPageCount = self.stackedWidget.count()
            self.setProgress()

            self.btn_box = QHBoxLayout()

            self.back_btn = QPushButton('Close')
            self.back_btn.clicked.connect(self.backpg)

            self.save_btn = QPushButton('Save .h5ad of Analysis')
            self.save_btn.clicked.connect(self.worker.save)

            self.next_btn =  QPushButton('Next')
            self.next_btn.clicked.connect(self.nextpgStart)

            self.checkDarkmode = QCheckBox('Dark Mode')
            self.checkDarkmode.toggled.connect(self.setDarkmode)

            hspacer =  QSpacerItem(1, 1, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
            
            # Layout
            self.btn_box.addWidget(self.back_btn, alignment=Qt.AlignmentFlag.AlignLeft) 
            self.btn_box.addWidget(self.checkDarkmode, alignment=Qt.AlignmentFlag.AlignLeft)
            self.btn_box.addItem(hspacer)
            self.btn_box.addWidget(self.save_btn)
            self.btn_box.addWidget(self.next_btn, alignment=Qt.AlignmentFlag.AlignRight)

            self.save_btn.hide()
            
            self.layoutPages.setContentsMargins(0,0,0,0)
            
            self.subbox.addWidget(self.stackedWidget)
            self.subbox.addLayout(self.btn_box)
            self.mainbox.addWidget(self.scrollNav)
            self.mainbox.addLayout(self.subbox)
            self.mainbox.setSpacing(0)
                 
            self.setLayout(self.mainbox)
            self.setContentsMargins(0,0,0,0)

            self.thread.start()
        
        def createSkipWindow(self):
            self.skipWindow = skipToCelltype()
        
        def setDarkmode(self):
            if self.checkDarkmode.isChecked():
                apply_stylesheet(app, theme='dark_blue.xml')
            else:
                apply_stylesheet(app, theme='light_blue.xml', invert_secondary=True)

        def createDGE(self, bySamples, mainItems, groupItems):
            try:
                if hasattr(self, 'DGEWindow'):
                    delattr(self, 'DGEWindow')
                self.DGEWindow = DGEWindow(bySamples, mainItems, groupItems)
                self.DGEWindow.toVolcanoSignal.connect(self.worker.toVolcano)
                print('Done!')
            except Exception as e:
                print(e)
        
        def createVol(self, figVol):
            try:
                if hasattr(self, 'volWindow'):
                    delattr(self, 'volWindow')
                self.volWindow = VolcanoPlot(figVol)
                self.DGEWindow.btnDGE.setEnabled(True)
            except Exception as e:
                self.DGEWindow.btnDGE.setEnabled(True)
                print(e)

        def enableDGEbtn(self):
            self.DGEWindow.btnDGE.setEnabled(True)

        def createFeatureMap(self, figUMAP, geneList, isMain):
            self.featureMap = FeatureMap(figUMAP, geneList, isMain)
            self.featureMap.updateFeatureSignal.connect(self.worker.toFeatureMapUpdate)
            self.featureMap.show()

        def updateFeatureMap(self, figUMAP):
            self.featureMap.updateUMAP(figUMAP)

        def toCluster(self, clusterIndex, nTop):
            try:
                if self.stackedWidget.currentIndex() == self.stackedWidget.indexOf(self.page8):         
                    isMain = True                    

                elif self.stackedWidget.currentIndex() == self.stackedWidget.indexOf(self.page10[4]):
                    isMain = False

                self.toWorkerCluster.emit(isMain, nTop, clusterIndex)
            except Exception as e:
                print(e)
        
        def toClusterSetPage(self, topNames, clusterIndex, isMain):
            try:
                if isMain:
                    self.page8sub[clusterIndex].setPage(topNames)
                    self.stackedWidget.setCurrentIndex(self.stackedWidget.indexOf(self.page8sub[clusterIndex]))
                    self.page8sub[clusterIndex].listDEG.verticalScrollBar().setValue(0)
                else:
                    self.page10sub[clusterIndex].setPage(topNames)
                    self.stackedWidget.setCurrentIndex(self.stackedWidget.indexOf(self.page10sub[clusterIndex]))
                    self.page10sub[clusterIndex].listDEG.verticalScrollBar().setValue(0)

                self.setProgress()
                self.setBackNextButton()
            except Exception as e:
                print(e)

        def toAnnotate(self, groupIndex):
            try:
                if self.stackedWidget.currentIndex() in self.sub8Index:
                    isMain = True
                    newName = self.page8sub[groupIndex].lineAnnotate.text()
                    self.page8sub[groupIndex].setName(self.page8sub[groupIndex].lineAnnotate.text())

                elif self.stackedWidget.currentIndex() in self.sub10Index:
                    isMain = False
                    newName = self.page10sub[groupIndex].lineAnnotate.text()
                    self.page10sub[groupIndex].setName(self.page10sub[groupIndex].lineAnnotate.text())

                self.toWorkerAnnotate.emit(groupIndex, newName, isMain)
            
            except Exception as e:
                print(e)

        def toAnnotateSetPage(self, groupIndex, figWholeUMAP, figSpecificUMAP, tempTypes, isSpatial, isMain):
            try:
                if isMain:
                    self.page8sub[groupIndex].setUMAP(figSpecificUMAP)
                    
                    self.page8.setUMAP(figWholeUMAP)
                    self.page8.setPage(tempTypes.values(), isSpatial)

                    self.stackedWidget.setCurrentIndex(8)

                else:
                    self.page10sub[groupIndex].setUMAP(figSpecificUMAP)
                    
                    self.page10[4].setUMAP(figWholeUMAP)
                    self.page10[4].setPage(tempTypes.values(), isSpatial)

                    self.stackedWidget.setCurrentIndex(14)

                self.setProgress()
                self.setBackNextButton()
            except Exception as e:
                print(e)
        
        def setSubCluster(self, figPCA, recPCs, multiSample):
            try:
                self.page10[0].setPage(figPCA, recPCs, multiSample)
                self.stackedWidget.setCurrentIndex(self.stackedWidget.currentIndex()+1)

                self.setProgress()
                self.setBackNextButton()

            except Exception as e:
                print(e)

        def gotopg(self, pageIndex):
            self.stackedWidget.setCurrentIndex(pageIndex)

            self.setProgress()
            self.setBackNextButton()
            
        
        def nextpgStart(self, currentPage = None, h5adDetected=False):
            try:
                self.next_btn.setEnabled(False)
                QApplication.processEvents()
                QApplication.processEvents()
                if self.devMode:
                    self.loadScreen = CellSortGame()
                else:
                    self.loadScreen = LoadingScreen()
                self.loadScreen.show()

                QApplication.processEvents()

                paramDict = {}
                if not currentPage:
                    currentPage = self.stackedWidget.currentIndex()
                    
                if currentPage == 8 and h5adDetected:
                    self.h5adImport = True
                    paramDict['h5adDetected'] = True
                
                if currentPage == 16:
                    self.close()

                if currentPage == 14:
                    print('Page 15: Setting Cell Types...')
                    
                    cellTypes = []

                    paramDict['h5adImport'] = self.h5adImport

                    for i in range(len(self.page10sub)):
                        cellTypes.append(self.page10sub[i].lineAnnotate.text())
                    paramDict['cellTypesSub'] = cellTypes

                elif currentPage == 13:
                    print('Page 14: Setting Annotation Page...')

                        
                # Page 12 Logic: Subcluster, select Leiden Resolution
                elif currentPage == 12:
                    print('Page 13: Performing Differentially Expressed Gene Analysis...')

                    paramDict['DEGMethodSub'] = self.page6.comboDEGMeth.currentText()
                    if self.page10[1].singleRes:
                        paramDict['chosenResSub'] = 0
                    elif self.page10[2].rbRes1.isChecked():
                        paramDict['chosenResSub'] = 0
                    elif self.page10[2].rbRes2.isChecked():
                        paramDict['chosenResSub'] = 1
                    elif self.page10[2].rbRes3.isChecked():
                        paramDict['chosenResSub'] = 2
                    else:
                        if self.devMode:
                            self.loadScreen._load_finished()
                        else:
                            self.loadScreen.close()
                        self.setBackNextButton()
                        raise ValueError('No Option Selected.')
                        

                # page 11 Logic: Subcluster, get test Leiden resolutions, Plot with Leiden
                elif currentPage == 11:
                    print('Page 12: Testing Resolutions...')
                    if self.page10[1].singleRes:
                        self.resList = [float(self.page10[1].inRes1.text())]
                    else:
                        self.resList = [
                            float(self.page10[1].inRes1.text()),
                            float(self.page10[1].inRes2.text()),
                            float(self.page10[1].inRes3.text())
                        ]
                    paramDict['resListSub'] = self.resList
                    paramDict['singleResSub'] = self.page10[1].singleRes

                # Page 10 Logic: Subcluster, Get NPC and Integration method, Create UMAP
                elif currentPage == 10:
                    print('Page 11: Integrating...')
                    paramDict['NPCSub'] = int(self.page10[0].inNPC.text())
                    paramDict['nNeighSub'] = int(self.page10[0].inNeigh.text())
                    paramDict['intMethod'] = self.page10[0].inInt.currentText()

                # Page 9 Logic: Go to Post Analysis
                elif currentPage == 9:
                    print('Page 10: Going to DGE Analysis...')
                        

                # Page 8 Logic: Annotate Cluster, Go to Cell type page
                elif currentPage == 8:
                    print('Page 9: Setting Cell Types...')
                    cellTypes = []

                    paramDict['h5adImport'] = self.h5adImport

                    if not self.h5adImport:
                        for i in range(len(self.page8sub)):
                            cellTypes.append(self.page8sub[i].lineAnnotate.text())
                        paramDict['cellTypes'] = cellTypes
                                        
                
                # Page 7 Logic: Differentially expressed genes Heatmap, Annotation Set up
                elif currentPage == 7:
                    print('Page 8: Setting Annotation Page...')


                # Page 6 Logic: Differentially expressed genes
                elif currentPage == 6: 
                    print('Page 7: Performing Differentially Expressed Gene Analysis...')
                    paramDict['DEGMethod'] = self.page6.comboDEGMeth.currentText()
                    if self.page5.singleRes:
                        paramDict['chosenRes'] = 0
                    elif self.page6.rbRes1.isChecked():
                        paramDict['chosenRes'] = 0
                    elif self.page6.rbRes2.isChecked():
                        paramDict['chosenRes'] = 1
                    elif self.page6.rbRes3.isChecked():
                        paramDict['chosenRes'] = 2
                    else:
                        if self.devMode:
                            self.loadScreen._load_finished()
                        else:
                            self.loadScreen.close()
                        self.setBackNextButton()
                        raise ValueError('No Option Selected.')

                # Page 5 Logic: Input Resolution, Clustering
                elif currentPage == 5: 
                    print('Page 6: Testing Resolutions...')
                    try:
                        if self.page5.singleRes:
                            self.resList = [float(self.page5.inRes1.text())]
                        else:
                            self.resList = [
                                float(self.page5.inRes1.text()),
                                float(self.page5.inRes2.text()),
                                float(self.page5.inRes3.text())
                            ]
                        paramDict['resList'] = self.resList
                        paramDict['singleRes'] = self.page5.singleRes
                    except Exception as e:
                        print(e)

                # Page 4 Logic: Inputs PC and Integration, Integration & UMAP
                elif currentPage == 4: 
                    print('Page 5: Integrating...')
                    paramDict['NPC'] = int(self.page4.inNPC.text())
                    paramDict['nNeigh'] = int(self.page4.inNeigh.text())
                    paramDict['intMethod'] = self.page4.inInt.currentText()

                # Page 3 Logic: Normalization
                elif currentPage == 3:
                    print('Page 4: Normalizing...')
                    paramDict['boolHVG'] = self.page3.checkHVG.isChecked()
                    paramDict['nHVG'] = int(self.page3.inHVG.text())
                                        
                # Page 2 Logic: Quality Control
                elif currentPage == 2:
                    print('Page 3: Filtering with Quality Control Thresholds...')
                    paramDict['scrublet'] = self.page2.checkScrublet.isChecked()
                    paramDict['perMt'] = float(self.page2.inPerMt.text())
                    paramDict['perRib'] = float(self.page2.inPerRib.text())
                    paramDict['minCells'] = int(self.page2.inMinCells.text())
                    paramDict['minGenes'] = int(self.page2.inMinGenes.text())

                # Page 1 Logic: Input Folders & Read
                elif currentPage == 1:
                    print('Page 2: Processing input files...')
                    if self.page1.folderIn_paths == [''] or self.page1.folderIn_paths == [] or self.page1.folderOut_path == "":
                        self.err_NoFiles.exec()
            
                    else:
                        # Try read to check if paths are valid
                        try:
                            paramDict = {}
                            paramDict['input_paths'] = self.page1.folderIn_paths
                            paramDict['output_path'] = self.page1.folderOut_path
                            paramDict['h5adImport'] = False
                            self.h5adImport = False
                            
                            i=1
                            if pathexists(join(paramDict['output_path'], 'iSNAP Output')):
                                while pathexists(join(paramDict['output_path'], f'iSNAP Output ({i})')):
                                    i += 1
                                
                                paramDict['output_path'] = join(paramDict['output_path'], f'iSNAP Output ({i})')
                            else:
                                paramDict['output_path'] = join(paramDict['output_path'], 'iSNAP Output')
                            

                            if self.modality in ['Xenium', 'Single .h5ad file']:
                                self.isSpatial = self.page1.isSpatial.isChecked()
                            else:
                                self.isSpatial = False

                            paramDict['isSpatial'] = self.isSpatial
                            paramDict['AvailableDEGMethods'] = [self.page6.comboDEGMeth.itemText(i) for i in range(self.page6.comboDEGMeth.count())]
                            

                        except Exception as e:
                            print(e)
                            self.err_NoFiles.exec()
                            raise ValueError(e)

                # Page 0 Logic: Set Data Source
                elif currentPage == 0:
                    print('Page 1: Setting File Folder...')
                    paramDict['seed'] = int(self.page0.inSeed.text())
                    paramDict['modality'] = self.page0.comboModality.currentText()
                    self.modality = paramDict['modality']


                self.toWorkerNextPG.emit(paramDict, currentPage)
            except Exception as e:
                print(e)
            

        def setNextpg(self, inputSetPage, currentPage):
            if currentPage == 15:
                self.stackedWidget.setCurrentIndex(currentPage+1)

            if currentPage == 14:
                try:
                    self.page9.setPage(inputSetPage['figUMAP'], inputSetPage['cellTypes'], inputSetPage['isSpatial'])

                    # Clear added clusters to save RAM
                    for widget in self.page10sub:
                        self.stackedWidget.removeWidget(widget)

                    self.stackedWidget.setCurrentIndex(self.stackedWidget.indexOf(self.page9))

                except Exception as e:
                    print(e)

            elif currentPage == 13:
                try: 
                    self.sub10Index = [] # Get sub page indexes for sub pages so that next and back is easier

                    self.page10sub = [None] * len(inputSetPage['clustersSub'])                    

                    for i in range(len(inputSetPage['clustersSub'])):
                        self.page10sub[i] = SpecificCluster(parent=self, group=inputSetPage['tempTypes10'][inputSetPage['clustersSub'][i]], group_index = str(i))
                        self.page10sub[i].setUMAP(inputSetPage['figSpecificUMAPs'][i])
                        self.page10sub[i].toAnnotateSignal.connect(self.toAnnotate)
                        self.stackedWidget.addWidget(self.page10sub[i])
                        self.sub10Index.append(self.stackedWidget.indexOf(self.page10sub[i])) # Get sub page indexes for sub pages so that next and back is easier
                    
                    self.page10[4].setUMAP(inputSetPage['figWholeUMAP'])
                    self.page10[4].setPage(inputSetPage['tempTypes10'].values(), inputSetPage['isSpatial'])

                    self.stackedWidget.setCurrentIndex(currentPage+1)

                except Exception as e:
                    print(e)
                    
            # Page 12 Logic: Subcluster, select Leiden Resolution
            elif currentPage == 12:
                try:
                    self.page10[3].setPage(inputSetPage['figDEG'])
                    self.stackedWidget.setCurrentIndex(currentPage+1)
                except Exception as e:
                    print(e)

            # page 11 Logic: Subcluster, get test Leiden resolutions, Plot with Leiden
            elif currentPage == 11:
                try:
                    self.page10[2].setPage(
                        inputSetPage['figListLeiden'],
                        inputSetPage['silhouetteScore'],
                        inputSetPage['resListLeiden'],
                        inputSetPage['singleRes']
                        )
                    self.stackedWidget.setCurrentIndex(currentPage+1)

                except Exception as e:
                    print(e)

            # Page 10 Logic: Subcluster, Get NPC and Integration method, Create UMAP
            elif currentPage == 10:
                try:
                    self.page10[1].setPage(inputSetPage['figUMAP'])
                    self.stackedWidget.setCurrentIndex(currentPage+1)
                except Exception as e:
                    print(e)

            # Page 9 Logic: Go to Post Analysis
            elif currentPage == 9:
                try:
                    self.page15.setPage(inputSetPage['sampleTypeMtx'], inputSetPage['samples'], inputSetPage['types'])
                    self.stackedWidget.setCurrentIndex(15)
                except Exception as e:
                    print(e)

            # Page 8 Logic: Annotate Cluster, Go to Cell type page
            elif currentPage == 8: 
                try:
                    self.page9.setPage(inputSetPage['figUMAP'], inputSetPage['cellTypes'], inputSetPage['isSpatial'])

                    self.stackedWidget.setCurrentIndex(self.stackedWidget.indexOf(self.page8)+1)
                    
                except Exception as e:
                    print(e)
            
            
            # Page 7 Logic: Differentially expressed genes Heatmap, Annotation Set up
            elif currentPage == 7:
                try: 
                    self.sub8Index = [] # Get sub page indexes for sub pages so that next and back is easier

                    self.page8sub = [None] * len(inputSetPage['clusters'])                    

                    for i in range(len(inputSetPage['clusters'])):
                        self.page8sub[i] = SpecificCluster(parent=self, group=inputSetPage['tempTypes8'][inputSetPage['clusters'][i]], group_index = str(i))
                        self.page8sub[i].setUMAP(inputSetPage['figSpecificUMAPs'][i])
                        self.page8sub[i].toAnnotateSignal.connect(self.toAnnotate)
                        self.stackedWidget.addWidget(self.page8sub[i])
                        self.sub8Index.append(self.stackedWidget.indexOf(self.page8sub[i])) # Get sub page indexes for sub pages so that next and back is easier
                    
                    self.page8.setUMAP(inputSetPage['figWholeUMAP'])
                    self.page8.setPage(inputSetPage['tempTypes8'].values(), inputSetPage['isSpatial'])

                    self.stackedWidget.setCurrentIndex(currentPage+1)
                except Exception as e:
                    print(e)

            # Page 6 Logic: Differentially expressed genes
            elif currentPage == 6: 
                try:
                    self.page7.setPage(inputSetPage['figDEG'])
                    
                    self.stackedWidget.setCurrentIndex(currentPage+1)

                except Exception as e:
                    print(e)

            # Page 5 Logic: Input Resolution, Clustering
            elif currentPage == 5: 
                try:
                    self.page6.setPage(
                        inputSetPage['figListLeiden'],
                        inputSetPage['silhouetteScore'],
                        inputSetPage['resListLeiden'],
                        inputSetPage['singleRes']
                        )
                    self.stackedWidget.setCurrentIndex(currentPage+1)
                
                except Exception as e:
                    print(e)


            # Page 4 Logic: Inputs PC and Integration, Integration & UMAP
            elif currentPage == 4: 
                try:
                    self.page5.setPage(inputSetPage['figUMAP'])
                    self.stackedWidget.setCurrentIndex(currentPage+1)

                except Exception as e:
                    print(e)
            # Page 3 Logic: Normalization
            elif currentPage == 3:
                try:
                    self.page4.setPage(inputSetPage['figPCA'], inputSetPage['recPCs'], inputSetPage['multiSample'])
                    self.stackedWidget.setCurrentIndex(currentPage+1)
                
                except Exception as e:
                    print(e)

            # Page 2 Logic: Quality Control
            elif currentPage == 2:
                try:
                    self.page3.createTable(
                        inputSetPage['uniqueSamples'],
                        inputSetPage['numCellsRaw'],
                        inputSetPage['avgGenesRaw'],
                        inputSetPage['numCellsQC'],
                        inputSetPage['avgGenesQC']
                        )

                    self.stackedWidget.setCurrentIndex(currentPage+1)

                except Exception as e:
                    print(e)
                    # *** Err invalid number

            # Page 1 Logic: Input Folders & Read
            elif currentPage == 1:
                # Try read to check if paths are valid
                try:                                            
                        self.page2.setPage(
                            inputSetPage['figCells'], 
                            inputSetPage['figGenes'], 
                            inputSetPage['figMT'], 
                            inputSetPage['figRP'], 
                            inputSetPage['input_paths'], 
                            inputSetPage['modality']
                            )
                        self.stackedWidget.setCurrentIndex(currentPage+1)

                except Exception as e:
                    print(e)
                    self.err_NoFiles.exec()

            # Page 0 Logic: Set Data Source
            elif currentPage == 0:
                self.page1.setPage(inputSetPage['modality'])
                
                self.stackedWidget.setCurrentIndex(currentPage+1)


            QApplication.processEvents()
            QApplication.processEvents()
            
            self.next_btn.setEnabled(True)             
            
            self.setProgress()
            self.setBackNextButton()

            if self.devMode:
                self.loadScreen._load_finished()
            else:
                self.loadScreen.close()
                
            self.devMode = self.page0.checkDev.isChecked()
            
            print('Done! \n')
            

        def backpg(self):
            if self.stackedWidget.currentIndex() == 0:
                exit()

            elif self.stackedWidget.currentIndex() in self.sub8Index:
                self.stackedWidget.setCurrentIndex(8)
            
            elif self.stackedWidget.currentIndex() == 9 and self.h5adImport:
                self.stackedWidget.setCurrentIndex(1)

            elif self.stackedWidget.currentIndex() in self.sub10Index:
                self.stackedWidget.setCurrentIndex(14)

            elif self.stackedWidget.currentIndex() == 10:
                self.stackedWidget.setCurrentIndex(9)
            
            elif self.stackedWidget.currentIndex() == 15:
                self.stackedWidget.setCurrentIndex(9)

            else:
                try:
                    self.stackedWidget.setCurrentIndex(self.stackedWidget.currentIndex()-1)
                except Exception as e:
                    print(e)

            self.setProgress()
            self.setBackNextButton()

        def setProgress(self):
            currentIndex = self.stackedWidget.currentIndex()
            if currentIndex in self.sub8Index:
                currentIndex = 8
            if currentIndex in self.sub10Index:
                currentIndex = 10

            for i in range(currentIndex):
                self.tabPage[i].btnTab.setEnabled(True)
                self.tabPage[i].setStyleSheet('background-color:#90B6F7')
            
            self.tabPage[currentIndex].btnTab.setEnabled(True)
            self.tabPage[currentIndex].setStyleSheet('background-color:None')

            if self.h5adImport:
                for i in range(2, 9):
                    self.tabPage[i].btnTab.setEnabled(False)
            
            for i in range(currentIndex+1, self.mainPageCount):
                self.tabPage[i].btnTab.setEnabled(False)
                self.tabPage[i].setStyleSheet('background-color:None')
        
        def setBackNextButton(self):
            if self.stackedWidget.currentIndex() == 0:
                self.back_btn.setText('Close')
                self.save_btn.hide()
                self.next_btn.setText('Next')
                self.next_btn.setEnabled(True)
                self.back_btn.setEnabled(True)
            elif self.stackedWidget.currentIndex() == 8 or self.stackedWidget.currentIndex() == 14:
                self.back_btn.setText('Back')
                self.save_btn.hide()
                self.next_btn.setText('Finish Annotation')
                self.next_btn.setEnabled(True)
                self.back_btn.setEnabled(True)
            elif self.stackedWidget.currentIndex() == 9:
                self.back_btn.setText('Back')
                self.save_btn.show()
                self.next_btn.setText('Next')
                self.next_btn.setEnabled(True)
                self.back_btn.setEnabled(True)
            elif self.stackedWidget.currentIndex() in self.sub8Index or self.stackedWidget.currentIndex() in self.sub10Index:
                self.back_btn.setText('Back')
                self.save_btn.hide()
                self.next_btn.setText('Finish Annotation')
                self.next_btn.setEnabled(False)
                self.back_btn.setEnabled(True)
            elif self.stackedWidget.currentIndex() == 15:
                self.back_btn.setText('Back')
                self.save_btn.show()
                self.next_btn.setText('Finished')
                self.next_btn.setEnabled(True)
                self.back_btn.setEnabled(True)
            elif self.stackedWidget.currentIndex() == 16:
                self.back_btn.setText('Back')
                self.save_btn.show()
                self.next_btn.setText('Close (Unsaved Data will be lost)')
                self.next_btn.setEnabled(True)
                self.next_btn.setEnabled(True)
            else:
                self.back_btn.setText('Back')
                self.save_btn.hide()
                self.next_btn.setText('Next')
                self.next_btn.setEnabled(True)
                self.back_btn.setEnabled(True)


class TabButton(QWidget):
    def __init__(self, parent, pageName, pageIndex):
        super(TabButton, self).__init__(parent)
        self.layoutMain = QVBoxLayout()

        self.btnTab = QPushButton(pageName)
        self.btnTab.clicked.connect(lambda: parent.gotopg(pageIndex))
        self.btnTab.setEnabled(False)

        self.layoutMain.addWidget(self.btnTab)

        self.setLayout(self.layoutMain)     


class MainWorker(QObject):
    toNextpgSet = pyqtSignal(dict, int)
    skipSignal = pyqtSignal()
    toAnnotateSignal = pyqtSignal(int, Figure , Figure, dict, bool, bool)
    toClusterSignal = pyqtSignal(ndarray, int, bool)
    initiateFeatureSignal = pyqtSignal(Figure, list, bool)
    updateFeatureSignal = pyqtSignal(Figure)
    toSubClusterSignal = pyqtSignal(Figure, int, bool)
    toDGESignal = pyqtSignal(bool, list, list)
    toVolSignal = pyqtSignal(Figure)
    toVolFailSignal = pyqtSignal()

    def __init__(self):
        super().__init__()

        self.paramDict = {}
    
    def save(self):
            try:
                save_adatas(self.adatasWhole.copy(), self.paramDict['output_path'], suffix='All')

            except Exception as e:
                print(e)
    
    def toDGE(self, bySamples):
        print('Creating DGE Window...')
        try:
            if bySamples:
                mainItems = self.adatasWhole.obs['celltype'].cat.categories.tolist()
                groupItems = self.adatasWhole.obs['sample'].cat.categories.tolist()
            else:
                mainItems = self.adatasWhole.obs['sample'].cat.categories.tolist()
                groupItems = self.adatasWhole.obs['celltype'].cat.categories.tolist()
        
            self.toDGESignal.emit(bySamples, mainItems, groupItems)

        except Exception as e:
            print(e)

    def toVolcano(self, mainList, groupAList, groupBList, bySamples, DGEMethod, nTopGenes, listGenes, pvalThresh, logFCThresh):
            try:
                print('Creating Volcano Plot...')
                if bySamples:
                    adatas = self.adatasWhole[self.adatasWhole.obs['celltype'].isin(mainList)].copy()
                    adatas = adatas[(adatas.obs['sample'].isin(groupAList)) | (adatas.obs['sample'].isin(groupBList))].copy()
                    mapDGE = {}
                    for item in groupAList:
                        mapDGE[item] = 'Group A'
                    for item in groupBList:
                        mapDGE[item] = 'Group B'

                    adatas.obs['groupDGE'] = adatas.obs['sample'].copy().map(mapDGE)
                else:
                    adatas = self.adatasWhole[self.adatasWhole.obs['sample'].isin(mainList)].copy()
                    adatas = adatas[(adatas.obs['celltype'].isin(groupAList)) | (adatas.obs['celltype'].isin(groupBList))].copy()
                    
                    mapDGE = {}
                    for item in groupAList:
                        mapDGE[item] = 'Group A'
                    for item in groupBList:
                        mapDGE[item] = 'Group B'

                    adatas.obs['groupDGE'] = adatas.obs['celltype'].copy().map(mapDGE)

                rank_genes_groups(adatas, groupby='groupDGE', key_added=f'DGE {DGEMethod}', reference='Group B', method=DGEMethod, use_raw=False)

                dgetoCSV(adatas, mainList, groupAList, groupBList, DGEMethod, self.paramDict['output_path'])

                figVol = plotVolcano(adatas, DGEMethod, nTopGenes, listGenes, pvalThresh, logFCThresh)

                outpath = join(self.paramDict['output_path'], 'Figures', 'DGE Volcano Plot')
                if not pathexists(outpath):
                    makedirs(outpath)

                i=1
                if pathexists(join(outpath, 'VolcanoPlot.png')):
                    while pathexists(join(outpath, f'VolcanoPlot({i}).png')):
                        i += 1
                    figVol.savefig(join(outpath, f'VolcanoPlot({i}).png'))
                    
                else:
                    figVol.savefig(join(outpath, 'VolcanoPlot.png'))
                    
                print(f'Volcano Plot Saved Under: {outpath}')

                self.toVolSignal.emit(figVol)

            except Exception as e:
                self.toVolFailSignal.emit()
                print(e)

            

    
    def toSubCluster(self, chosenType):
        
        try:
            self.paramDict['chosenType'] = chosenType
            self.adatasSub = self.adatasWhole[self.adatasWhole.obs['celltype'] == chosenType].copy()

            self.adatasSub = normalization(
                        self.adatasSub.copy(), 
                        self.paramDict['boolHVG'], 
                        self.paramDict['nHVG'], # Keeps only highly variable genes if hvg parameter is set to true.
                        self.paramDict['seed']
                        )
            
            suffix = self.paramDict['chosenType']
            # Create a PCA Variance Ratio to help PC number choice, save so that it can display and for future paper use
            outpath = join(self.paramDict['output_path'], 'Figures', 'PCA Variance Ratio')
            if not pathexists(outpath):
                makedirs(outpath)

            figPCA = pcaVarianceRatio(self.adatasSub.copy(), n_pcs=50)
            figPCA.savefig(join(outpath, f'PCAVarianceRatio_{suffix}.png'))

            # Find optimal PCs by using maximum distance from endpoint line
            var = self.adatasSub.uns['pca']['variance_ratio'].copy()
            pcList = range(1, len(var)+1)
            
            slope = (var[-1]-var[0])/(pcList[-1]-pcList[0])
            distsqPC = []
            for i in range(len(var)):
                x_closest = (pcList[0]*slope + pcList[i]/slope + var[i] - var[0])/(slope + 1/slope)
                y_closest = (x_closest-pcList[0]) * slope + var[0]
                distsqPC.append((x_closest-pcList[i])**2 + (y_closest-var[i])**2)
            
            recPCs = argmax(distsqPC)+1

            self.toSubClusterSignal.emit(figPCA, recPCs, len(self.adatasSub.obs['sample'].cat.categories)>1)
        except Exception as e:
            print(e)
    
    def initiateFeatureMap(self, isAnnotated, isMain):
        if isAnnotated:
            size = 240000/self.adatasWhole.n_obs
            figUMAP = pltumap(self.adatasWhole.copy(), show=False, color="celltype", frameon = False, palette=default_102, size=size, title='Cell Type', return_fig=True)

            genes = self.adatasWhole.var_names.tolist()
            genesList = []
            for gene in genes:
                if self.adatasWhole.var['n_counts'][gene] > 0:
                    genesList.append(gene)
            
        elif isMain:
            suffix = 'All'
            size = 240000/self.adatasWhole.n_obs
            figUMAP = pltumap(self.adatasWhole, show=False, color='tempType', frameon=False, palette=default_102, size=size, title=f'Annotations UMAP, {suffix}', return_fig=True)

            genes = self.adatasWhole.var_names.tolist()
            genesList = []
            for gene in genes:
                if self.adatasWhole.var['n_counts'][gene] > 0:
                    genesList.append(gene)
            
        else:
            suffix = self.paramDict['chosenType']
            genes = self.adatasSub.var_names.tolist()
            size = 240000/self.adatasSub.n_obs
            figUMAP = pltumap(self.adatasSub, show=False, color='tempType', frameon=False, palette=default_102, size=size, title=f'Annotations UMAP, {suffix}', return_fig=True)

            genesList = []
            for gene in genes:
                if self.adatasSub.var['n_counts'][gene] > 0:
                    genesList.append(gene)

        figUMAP.set_constrained_layout(True)
        self.initiateFeatureSignal.emit(figUMAP, genesList, isMain)

    def toFeatureMapUpdate(self, isMain, marker, vmax, vmin):
        outpath = join(self.paramDict['output_path'], 'Figures', 'Feature Maps')
        reds = copy(Reds)
        reds.set_under('lightgray')

        if vmax=='':
            vmax = 'p95'
        else:
            try:
                vmax = float(vmax)
            except:
                pass
        
        if vmin=='':
            vmin='p1.5'
        else:
            try:
                vmin = float(vmin)
            except:
                pass

        try:
            if isMain:
                suffix = 'All'
                size = 240000/self.adatasWhole.n_obs
                figUMAP = pltumap(self.adatasWhole, color=marker, frameon=False, show=False, vmin=vmin, vmax=vmax, cmap=reds, size=size, return_fig=True)               
            else:
                suffix = self.paramDict['chosenType']
                size = 240000/self.adatasSub.n_obs
                figUMAP = pltumap(self.adatasSub, color=marker, frameon=False, show=False, vmin=vmin, vmax=vmax, cmap=reds, size=size, return_fig=True)
            
            # figUMAP.set_constrained_layout(True)
            outpath = join(self.paramDict['output_path'], 'Figures', 'Feature Maps', f'{suffix}')

            if not pathexists(outpath):
                makedirs(outpath)
            
            figUMAP.savefig(join(outpath, f'{marker} Feature Map.png'))

            self.updateFeatureSignal.emit(figUMAP)

        except Exception as e:
            print(e)


    def typeToCSVWork(self, isMain, groupKey, suffix):
        try:
            if isMain:
                typeToCSV(self.adatasWhole, groupKey, suffix+' All', self.paramDict['output_path'])
            else:
                typeToCSV(self.adatasSub, groupKey, suffix+f' {self.paramDict["chosenType"]}', self.paramDict['output_path'])
        except Exception as e:
            print(e)
        
    

    def combTypes(self, type1, type2, newName):
        celltypes = self.adatasWhole.obs['celltype'].cat.categories
        typemap = {}
        for celltype in celltypes:
            typemap[celltype] = celltype
        
        typemap[type1] = newName
        typemap[type2] = newName

        self.adatasWhole.obs['celltype'] = self.adatasWhole.obs['celltype'].copy().map(typemap).astype('category')
        size = 240000/self.adatasWhole.n_obs
        figUMAP = pltumap(self.adatasWhole, show=False, color="celltype", frameon = False, palette=default_102, size=size, title='Cell Type', return_fig=True)
        figUMAP.set_constrained_layout(True)

        outpath = join(self.paramDict['output_path'], 'Figures' ,'UMAP', 'Cell Types')
        if not pathexists(outpath):
            makedirs(outpath)
        i=1
        if pathexists(join(outpath, 'UMAP CellType.png')):
            while pathexists(join(outpath, f'UMAP CellType({i})')):
                i+=1
            figUMAP.savefig(join(outpath, f'UMAP CellType({i})'))
        else:
            figUMAP.savefig(join(outpath, 'UMAP CellType.png'))
        
        inputSetPage = {}
        inputSetPage['figUMAP'] = figUMAP
        inputSetPage['cellTypes'] = self.adatasWhole.obs['celltype'].cat.categories
        inputSetPage['isSpatial'] = self.paramDict['isSpatial']

        with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:
            lines = [
                'Combine Types:\n',
                f'Type 1 = {type1}\n',
                f'Type 2 = {type2}\n',
                f'New Name = {newName}\n\n',
            ]
            log.writelines(lines)

        self.toNextpgSet.emit(inputSetPage, 8)
    
    def toAnnotate(self, groupIndex, newName, isMain):
        try:
            if isMain:
                self.tempTypes8[self.clusters[groupIndex]] = newName
                self.adatasWhole.obs['tempType'] = self.adatasWhole.obs['leiden'].copy().map(self.tempTypes8)

                suffix = 'All'

                # Copy to output
                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Annotation', suffix)
                if not pathexists(outpath):
                    makedirs(outpath)

                # Specific UMAP
                figUMAPSpecific = specificClusterUMAP(self.adatasWhole.copy(), self.tempTypes8[self.clusters[groupIndex]], key='tempType', suffix=suffix)
                figUMAPSpecific.set_constrained_layout(True)
                figUMAPSpecific.savefig(join(outpath, f'UMAP Specific Cluster {self.tempTypes8[self.clusters[groupIndex]]}.png'))

                # Whole Unannotated UMAP
                size = 240000/self.adatasWhole.n_obs
                figUMAP = pltumap(self.adatasWhole, show=False, color='tempType', frameon=False, palette=default_102, size=size, title=f'Unannotated, {suffix}', return_fig=True)
                figUMAP.set_constrained_layout(True)
                figUMAP.savefig(join(outpath, f'UMAP Unannotated {suffix}.png'))

                self.toAnnotateSignal.emit(groupIndex, figUMAP, figUMAPSpecific, self.tempTypes8, self.paramDict['isSpatial'], isMain)

            else:
                self.tempTypes10[self.clustersSub[groupIndex]] = newName
                self.adatasSub.obs['tempType'] = self.adatasSub.obs['leiden'].copy().map(self.tempTypes10)

                suffix = self.paramDict['chosenType']

                # Copy to output
                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Annotation', suffix)
                if not pathexists(outpath):
                    makedirs(outpath)

                # Specific UMAP
                figUMAPSpecific = specificClusterUMAP(self.adatasSub.copy(), self.tempTypes10[self.clustersSub[groupIndex]], key='tempType', suffix=suffix)
                figUMAPSpecific.set_constrained_layout(True)
                figUMAPSpecific.savefig(join(outpath, f'UMAP Specific Cluster {self.tempTypes10[self.clustersSub[groupIndex]]}.png'))

                # Whole Unannotated UMAP
                size = 240000/self.adatasSub.n_obs
                figUMAP = pltumap(self.adatasSub.copy(), show=False, color='tempType', frameon=False, palette=default_102, size=size, title=f'Unannotated, {suffix}', return_fig=True)
                figUMAP.set_constrained_layout(True)
                figUMAP.savefig(join(outpath, f'UMAP Unannotated {suffix}.png'))

                self.toAnnotateSignal.emit(groupIndex, figUMAP, figUMAPSpecific, self.tempTypes10, self.paramDict['isSpatial'], isMain)

        except Exception as e:
            print(e)
    
    def toCluster(self, isMain, nTop, clusterIndex):
        try:
            if isMain:
                topNames = self.adatasWhole.uns['DEG']['names'][str(clusterIndex)]
            else:
                topNames = self.adatasSub.uns['DEG']['names'][str(clusterIndex)]
             
            if len(topNames) > nTop:
                topNames = topNames[:nTop]

            
            self.toClusterSignal.emit(topNames, clusterIndex, isMain)
        except Exception as e:
            print(e)

    
    def processData(self, inParamDict, currentPage):
        inputSetPage = {}
        for param in inParamDict.keys():
            self.paramDict[param] = inParamDict[param] 

        if 'h5adDetected' in inParamDict.keys():
            if currentPage == 8:
                self.h5adImport = True
                inputSetPage['h5adDetected'] = True
                self.paramDict['boolHVG'] = False
                self.paramDict['nHVG'] = 0
                self.adatasWhole = self.adatasRaw.copy()

        if currentPage == 14:
            try:
                print('')
                self.adatasSub = leidenToCellType(self.adatasSub.copy(), self.paramDict['cellTypesSub'])
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:
                    lines = [
                        f'Page 14: Annotations ({self.paramDict["chosenType"]})\n',
                        'Cell Types:\n'
                    ]
                    for i in range(len(self.paramDict['cellTypesSub'])):
                        text = f'\t Cluster {i} = {self.paramDict["cellTypesSub"][i]}\n'
                        lines.append(text)
                    lines.append('\n')

                    log.writelines(lines)

                # save_adatas(self.adatasWhole.copy(), self.paramDict['output_path'], suffix='Initial')

                newcelltype = self.adatasSub.obs['celltype'].copy()
                oldcelltype = self.adatasWhole[self.adatasWhole.obs['celltype']!=self.paramDict['chosenType']].obs['celltype'].copy()

                combcelltype = concat([newcelltype, oldcelltype])    

                self.adatasWhole.obs['celltype'] = combcelltype.astype('category')
                
                size = 240000/self.adatasWhole.n_obs
                figUMAP = pltumap(self.adatasWhole.copy(), show=False, color="celltype", frameon = False, palette=default_102, size=size, title='Cell Type', return_fig=True)
                figUMAP.set_constrained_layout(True)

                outpath = join(self.paramDict['output_path'], 'Figures' ,'UMAP', 'Cell Types')
                if not pathexists(outpath):
                    makedirs(outpath)
                i=1
                if pathexists(join(outpath, 'UMAP CellType.png')):
                    while pathexists(join(outpath, f'UMAP CellType({i})')):
                        i+=1
                    figUMAP.savefig(join(outpath, f'UMAP CellType({i})'))
                else:
                    figUMAP.savefig(join(outpath, 'UMAP CellType.png'))
                
                inputSetPage['figUMAP'] = figUMAP
                inputSetPage['cellTypes'] = self.adatasWhole.obs['celltype'].cat.categories
                inputSetPage['isSpatial'] = self.paramDict['isSpatial']
                
                save_adatas(self.adatasSub, self.paramDict['output_path'], suffix=self.paramDict['chosenType'])
                
                del self.adatasSub
                del self.leidenListSub

            except Exception as e:
                print(e)

        elif currentPage == 13:
            try: 
                self.tempTypes10 = {}

                degResults = self.adatasSub.uns['DEG'].copy()
                self.clustersSub = degResults['names'].dtype.names

                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines = [
                        f'Page 13: DEG Dotplot ({self.paramDict["chosenType"]})\n\n',
                    ]
                    log.writelines(lines)

                suffix = self.paramDict['chosenType']

                for i in range(len(self.clustersSub)):
                    self.tempTypes10[self.clustersSub[i]] = f'{suffix} {i}'

                self.adatasSub.obs['tempType'] = self.adatasSub.obs['leiden'].copy().map(self.tempTypes10)

                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Annotation', suffix)
                if not pathexists(outpath):
                    makedirs(outpath)

                figSpecificUMAPs = []
                for i in range(len(self.clustersSub)):

                    figUMAP = specificClusterUMAP(self.adatasSub.copy(), self.tempTypes10[self.clustersSub[i]], key='tempType', suffix=suffix)
                    figUMAP.set_constrained_layout(True)
                    figUMAP.savefig(join(outpath, f'UMAP Specific Cluster {self.tempTypes10[self.clustersSub[i]]}.png'))
                    figSpecificUMAPs.append(figUMAP)
                    
                
                # Whole Unannotated UMAP
                size = 240000/self.adatasSub.n_obs
                figUMAP = pltumap(self.adatasSub, show=False, color='tempType', frameon=False, palette=default_102, size=size, title=f'Unannotated, {suffix}', return_fig=True)
                figUMAP.set_constrained_layout(True)

                figUMAP.savefig(join(outpath, f'UMAP Unannotated {suffix}.png'))
                
                inputSetPage['clustersSub'] = self.clustersSub
                inputSetPage['tempTypes10'] = self.tempTypes10
                inputSetPage['figSpecificUMAPs'] = figSpecificUMAPs
                inputSetPage['figWholeUMAP'] = figUMAP
                inputSetPage['isSpatial'] = self.paramDict['isSpatial']

                pltclose('all')
            except Exception as e:
                print(e)
                
        # Page 12 Logic: Subcluster, select Leiden Resolution
        elif currentPage == 12:
            try:
                self.adatasSub.obs['leiden'] = self.leidenListSub[self.paramDict['chosenResSub']][0]
                self.adatasSub.uns['leiden'] = self.leidenListSub[self.paramDict['chosenResSub']][1]
                
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines = [
                        f'Page 12: Choose Leiden Clustering ({self.paramDict["chosenType"]})\n',
                        f'Chosen Resolution = {self.paramDict['resListSub'][self.paramDict["chosenResSub"]]}\n',
                        f'Number of Clusters = {len(set(self.adatasSub.obs["leiden"]))}\n'
                        f'DEG Statistical Method = {self.paramDict["DEGMethodSub"]}\n\n'
                    ]
                    log.writelines(lines)

                self.adatasSub = deg(self.adatasSub.copy(), self.paramDict['DEGMethodSub'])

                suffix = self.paramDict['chosenType']

                degtoCSV(self.adatasSub.copy(), self.paramDict['output_path'], self.paramDict['DEGMethodSub'], suffix=suffix)

                # Set up Dot Plot
                dendrogram(self.adatasSub, groupby='leiden', n_pcs=self.paramDict['NPCSub'])
                figDEG = rank_genes_groups_dotplot(self.adatasSub.copy(), n_genes=5, key='DEG', groupby="leiden", show=False, return_fig=True)
                # *** Make own for clean labeling


                # Copy to output
                outpath = join(self.paramDict['output_path'], 'Figures', 'Dotplot')
                if not pathexists(outpath):
                    makedirs(outpath)
                
                figDEG.savefig(join(outpath, f'Dotplot {suffix}.png'))
                
                # Temp Fix
                figDEG = Figure()
                axDEG = figDEG.add_subplot(111)
                pltDEG = imread(join(outpath, f'Dotplot {suffix}.png'))
                axDEG.imshow(pltDEG)
                axDEG.axis('off')

                inputSetPage['figDEG'] = figDEG
                
                pltclose('all')
                
            except Exception as e:
                print(e)

        # page 11 Logic: Subcluster, get test Leiden resolutions, Plot with Leiden
        elif currentPage == 11:
            try:
                if self.paramDict['singleResSub']:
                    self.leidenListSub = [
                    cluster(self.adatasSub.copy(), self.paramDict['resListSub'][0], self.paramDict['seed']), 
                    ]
                    lines = [
                        f'Page 11: Leiden Parameters ({self.paramDict["chosenType"]})\n',
                        f'Test Resolution = {self.paramDict['resListSub'][0]}\n\n',
                    ]
                else:
                    self.leidenListSub = [
                        cluster(self.adatasSub.copy(), self.paramDict['resListSub'][0], self.paramDict['seed']), 
                        cluster(self.adatasSub.copy(), self.paramDict['resListSub'][1], self.paramDict['seed']), 
                        cluster(self.adatasSub.copy(), self.paramDict['resListSub'][2], self.paramDict['seed'])
                        ]
                    lines = [
                        f'Page 11: Leiden Parameters ({self.paramDict["chosenType"]})\n',
                        f'Test Resolution 1 = {self.paramDict['resListSub'][0]}\n',
                        f'Test Resolution 2 = {self.paramDict['resListSub'][1]}\n',
                        f'Test Resolution 3 = {self.paramDict['resListSub'][2]}\n\n',
                    ]

                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:
                    log.writelines(lines)

                suffix = self.paramDict['chosenType']               

                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Leiden Resolution', f'{suffix}')
                if not pathexists(outpath):
                    makedirs(outpath)  
                
                figUMAPList = []
                silhouetteScore = []
                for i in range(len(self.paramDict['resListSub'])):
                    adatasLeiden = self.adatasSub.copy()
                    adatasLeiden.obs['leiden'] = self.leidenListSub[i][0]
                    adatasLeiden.uns['leiden'] = self.leidenListSub[i][1]

                    silhouetteScore.append(silhouetteLeiden(adatasLeiden))

                    size = 240000/adatasLeiden.n_obs
                    figUMAPList.append(pltumap(adatasLeiden, show=False, color="leiden", frameon = False, palette=default_102, title=f'Resolution {self.paramDict['resListSub'][i]}, {suffix}', size=size, return_fig=True))
                    figUMAPList[i].set_constrained_layout(True)

                    figUMAPList[i].savefig(join(outpath, f'Leiden Resolution {self.paramDict["resListSub"][i]}.png'))
                
                inputSetPage['figListLeiden'] = figUMAPList
                inputSetPage['silhouetteScore'] = silhouetteScore
                inputSetPage['resListLeiden'] = self.paramDict['resListSub']
                inputSetPage['singleRes'] = self.paramDict['singleResSub']
            
                pltclose('all')
                                
            except Exception as e:
                print(e)

        # Page 10 Logic: Subcluster, Get NPC and Integration method, Create UMAP
        elif currentPage == 10:
            try:
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines = [
                        f'Page 10: Integration & UMAP Projection ({self.paramDict['chosenType']})\n',
                        f'Number of PCs = {self.paramDict["NPCSub"]}\n',
                        f'Number of nNeigh = {self.paramDict["nNeighSub"]}\n',
                        f'Integration Method = {self.paramDict["intMethod"]}\n\n',
                    ]
                    log.writelines(lines)

                self.adatasSub = integration(self.adatasSub.copy(), self.paramDict["NPCSub"], self.paramDict["nNeighSub"], self.paramDict["intMethod"], self.paramDict['seed'])
                suffix = self.paramDict['chosenType']

                size = 240000/self.adatasSub.n_obs
                figUMAP = pltumap(self.adatasSub.copy(), show=False, frameon=False, color='sample', size=size, title=f'Samples, {suffix}', return_fig=True) 
                figUMAP.set_constrained_layout(True)
                # *** Make own for clean labelling

                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Samples', f'{suffix}')
                if not pathexists(outpath):
                    makedirs(outpath)
                figUMAP.savefig(join(outpath, f'SampleUMAP_{suffix}.png'))

                inputSetPage['figUMAP'] = figUMAP

                pltclose('all')
            except Exception as e:
                print(e)

        # Page 9 Logic: Go to Post Analysis
        elif currentPage == 9:
            print('Setting cell count table...')
            try:
                samples = list(self.adatasWhole.obs['sample'].cat.categories)
                samples.append('All Samples')

                types = list(self.adatasWhole.obs['celltype'].cat.categories)
                types.append('All Types')
                
                sampleTypeMtx = zeros((len(types), len(samples)))

                for i in range(len(samples)-1):
                    sampleTypeMtx[len(types)-1, i] = self.adatasWhole[self.adatasWhole.obs['sample']==samples[i]].n_obs
                    for j in range(len(types)-1):
                        sampleTypeMtx[j, i] = self.adatasWhole[(self.adatasWhole.obs['sample']==samples[i]) & (self.adatasWhole.obs['celltype']==types[j])].n_obs
                
                for i in range(len(types)-1):
                    sampleTypeMtx[i, len(samples)-1] = self.adatasWhole[self.adatasWhole.obs['celltype']==types[i]].n_obs
                
                sampleTypeMtx[len(types)-1, len(samples)-1] = self.adatasWhole.n_obs

                inputSetPage['samples'] = samples
                inputSetPage['types'] = types
                inputSetPage['sampleTypeMtx'] = sampleTypeMtx

            except Exception as e:
                print(e)
            

        # Page 8 Logic: Annotate Cluster, Go to Cell type page
        elif currentPage == 8: 
            try:
                if not self.paramDict['h5adImport']:
                    self.adatasWhole = leidenToCellType(self.adatasWhole.copy(), self.paramDict['cellTypes'])
                    with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                        lines = [
                            'Page 8: Annotations\n',
                            'Cell Types:\n'
                        ]
                        for i in range(len(self.paramDict['cellTypes'])):
                            text = f'\t Cluster {i} = {self.paramDict["cellTypes"][i]}\n'
                            lines.append(text)
                        lines.append('\n')

                        log.writelines(lines)

                    # save_adatas(self.adatasWhole.copy(), self.paramDict['output_path'], suffix='Initial')
                
                size = 240000/self.adatasWhole.n_obs
                figUMAP = pltumap(self.adatasWhole.copy(), show=False, color="celltype", frameon = False, palette=default_102, size=size, title='Cell Type', return_fig=True)
                figUMAP.set_constrained_layout(True)

                outpath = join(self.paramDict['output_path'], 'Figures' ,'UMAP', 'Cell Types')
                if not pathexists(outpath):
                    makedirs(outpath)
                i=1
                if pathexists(join(outpath, 'UMAP CellType.png')):
                    while pathexists(join(outpath, f'UMAP CellType({i})')):
                        i+=1
                    figUMAP.savefig(join(outpath, f'UMAP CellType({i})'))
                else:
                    figUMAP.savefig(join(outpath, 'UMAP CellType.png'))
                
                inputSetPage['figUMAP'] = figUMAP
                inputSetPage['cellTypes'] = self.adatasWhole.obs['celltype'].cat.categories
                inputSetPage['isSpatial'] = self.paramDict['isSpatial']
                
            except Exception as e:
                print(e)
        
        
        # Page 7 Logic: Differentially expressed genes Heatmap, Annotation Set up
        elif currentPage == 7:
            try: 
                if not hasattr(self, 'tempTypes8'):
                    self.tempTypes8 = {}

                degResults = self.adatasWhole.uns['DEG'].copy()
                self.clusters = degResults['names'].dtype.names

                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines = [
                        'Page 7: DEG Dotplot\n\n',
                    ]
                    log.writelines(lines)
                
                for i in range(len(self.clusters)):
                    self.tempTypes8[self.clusters[i]] = f'Cluster {i}'

                self.adatasWhole.obs['tempType'] = self.adatasWhole.obs['leiden'].copy().map(self.tempTypes8)

                suffix = 'All'
                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Annotation', suffix)
                if not pathexists(outpath):
                    makedirs(outpath)

                figSpecificUMAPs = []
                for i in range(len(self.clusters)):

                    figUMAP = specificClusterUMAP(self.adatasWhole.copy(), self.tempTypes8[self.clusters[i]], key='tempType', suffix=suffix)
                    figUMAP.set_constrained_layout(True)
                    figUMAP.savefig(join(outpath, f'UMAP Specific Cluster {self.tempTypes8[self.clusters[i]]}.png'))
                    figSpecificUMAPs.append(figUMAP)
                    
                
                # Whole Unannotated UMAP
                size = 240000/self.adatasWhole.n_obs
                figUMAP = pltumap(self.adatasWhole, show=False, color='tempType', frameon=False, palette=default_102, size=size, title=f'Unannotated, {suffix}', return_fig=True)
                figUMAP.set_constrained_layout(True)

                figUMAP.savefig(join(outpath, f'UMAP Unannotated {suffix}.png'))
                
                inputSetPage['clusters'] = self.clusters
                inputSetPage['tempTypes8'] = self.tempTypes8
                inputSetPage['figSpecificUMAPs'] = figSpecificUMAPs
                inputSetPage['figWholeUMAP'] = figUMAP
                inputSetPage['isSpatial'] = self.paramDict['isSpatial']

                pltclose('all')

            except Exception as e:
                print(e)

        # Page 6 Logic: Differentially expressed genes
        elif currentPage == 6: 
            try:
                self.adatasWhole.obs['leiden'] = self.leidenList[self.paramDict['chosenRes']][0]
                self.adatasWhole.uns['leiden'] = self.leidenList[self.paramDict['chosenRes']][1]
                
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines = [
                        'Page 6: Choose Leiden Clustering\n',
                        f'Chosen Resolution = {self.paramDict['resList'][self.paramDict["chosenRes"]]}\n',
                        f'Number of Clusters = {len(set(self.adatasWhole.obs["leiden"]))}\n'
                        f'DEG Statistical Method = {self.paramDict["DEGMethod"]}\n\n'
                    ]
                    log.writelines(lines)

                self.adatasWhole = deg(self.adatasWhole.copy(), self.paramDict['DEGMethod'])

                degtoCSV(self.adatasWhole.copy(), self.paramDict['output_path'], self.paramDict['DEGMethod'])

                suffix = 'All'

                # Set up Dot Plot
                dendrogram(self.adatasWhole, groupby='leiden', n_pcs=self.paramDict['NPC'])
                figDEG = rank_genes_groups_dotplot(self.adatasWhole.copy(), n_genes=5, key='DEG', groupby="leiden", show=False, return_fig=True)
                # *** Make own for clean labeling


                # Copy to output
                outpath = join(self.paramDict['output_path'], 'Figures', 'Dotplot')
                if not pathexists(outpath):
                    makedirs(outpath)
                
                figDEG.savefig(join(outpath, f'Dotplot {suffix}.png'))
                
                # Temp Fix
                figDEG = Figure()
                axDEG = figDEG.add_subplot(111)
                pltDEG = imread(join(outpath, f'Dotplot {suffix}.png'))
                axDEG.imshow(pltDEG)
                axDEG.axis('off')

                inputSetPage['figDEG'] = figDEG
                
                pltclose('all')
            except Exception as e:
                print(e)

        # Page 5 Logic: Input Resolution, Clustering
        elif currentPage == 5: 
            try:
                if self.paramDict['singleRes']:
                    self.leidenList = [
                    cluster(self.adatasWhole.copy(), self.paramDict['resList'][0], self.paramDict['seed']), 
                    ]
                    lines = [
                        'Page 5: Leiden Parameters\n',
                        f'Test Resolution = {self.paramDict['resList'][0]}\n\n',
                    ]
                else:
                    self.leidenList = [
                        cluster(self.adatasWhole.copy(), self.paramDict['resList'][0], self.paramDict['seed']), 
                        cluster(self.adatasWhole.copy(), self.paramDict['resList'][1], self.paramDict['seed']), 
                        cluster(self.adatasWhole.copy(), self.paramDict['resList'][2], self.paramDict['seed'])
                        ]

                    lines = [
                        'Page 5: Leiden Parameters\n',
                        f'Test Resolution 1 = {self.paramDict['resList'][0]}\n',
                        f'Test Resolution 2 = {self.paramDict['resList'][1]}\n',
                        f'Test Resolution 3 = {self.paramDict['resList'][2]}\n\n',
                    ]

                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:
                    log.writelines(lines)

                suffix = 'All'               

                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Leiden Resolution', f'{suffix}')
                if not pathexists(outpath):
                    makedirs(outpath)  
                
                figUMAPList = []
                silhouetteScore = []
                for i in range(len(self.paramDict['resList'])):
                    adatasLeiden = self.adatasWhole.copy()
                    adatasLeiden.obs['leiden'] = self.leidenList[i][0]
                    adatasLeiden.uns['leiden'] = self.leidenList[i][1]

                    silhouetteScore.append(silhouetteLeiden(adatasLeiden))

                    size = 240000/adatasLeiden.n_obs
                    figUMAPList.append(pltumap(adatasLeiden, show=False, color="leiden", frameon = False, palette=default_102, title=f'Resolution {self.paramDict['resList'][i]}, {suffix}', return_fig=True))
                    figUMAPList[i].set_constrained_layout(True)

                    figUMAPList[i].savefig(join(outpath, f'Leiden Resolution {self.paramDict["resList"][i]}.png'))
                
                inputSetPage['figListLeiden'] = figUMAPList
                inputSetPage['silhouetteScore'] = silhouetteScore
                inputSetPage['resListLeiden'] = self.paramDict['resList']
                inputSetPage['singleRes'] = self.paramDict['singleRes']
            
                pltclose('all')
                                
            except Exception as e:
                print(e)


        # Page 4 Logic: Inputs PC and Integration, Integration & UMAP
        elif currentPage == 4: 
            try:
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines = [
                        'Page 4: Integration & UMAP Projection\n',
                        f'Number of PCs = {self.paramDict['NPC']}\n',
                        f'Number of nNeigh = {self.paramDict['nNeigh']}\n',
                        f'Integration Method = {self.paramDict['intMethod']}\n\n',
                    ]
                    log.writelines(lines)
            
                self.adatasWhole = integration(self.adatasWhole.copy(), self.paramDict['NPC'], self.paramDict['nNeigh'], self.paramDict['intMethod'], self.paramDict['seed'])

                suffix = 'All'

                size = 240000/self.adatasWhole.n_obs
                figUMAP = pltumap(self.adatasWhole.copy(), show=False, frameon=False, color='sample', size=size, title=f'Samples, {suffix}', return_fig=True) 
                figUMAP.set_constrained_layout(True)
                # *** Make own for clean labelling

                outpath = join(self.paramDict['output_path'], 'Figures', 'UMAP Samples', f'{suffix}')
                if not pathexists(outpath):
                    makedirs(outpath)
                figUMAP.savefig(join(outpath, f'SampleUMAP_{suffix}.png'))

                inputSetPage['figUMAP'] = figUMAP

                pltclose('all')

            except Exception as e:
                print(e)
        # Page 3 Logic: Normalization
        elif currentPage == 3:
            try:
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                        lines = [
                            'Page 3: Normalization & HVG\n',
                            f'Do HVG: {self.paramDict["boolHVG"]}\n',
                            f'Number of HVG: {self.paramDict["nHVG"]}\n\n'
                        ]
                        log.writelines(lines)

                # Normalize Total, log1p, HVG, and scaling
                self.adatasWhole = normalization(
                    self.adatasBeforeHVG.copy(), 
                    self.paramDict['boolHVG'], 
                    self.paramDict['nHVG'], # Keeps only highly variable genes if hvg parameter is set to true.
                    self.paramDict['seed']
                    ) 

                suffix = 'All'
                # Create a PCA Variance Ratio to help PC number choice, save so that it can display and for future paper use
                outpath = join(self.paramDict['output_path'], 'Figures', 'PCA Variance Ratio')
                if not pathexists(outpath):
                    makedirs(outpath)

                figPCA = pcaVarianceRatio(self.adatasWhole, n_pcs=50)
                figPCA.savefig(join(outpath, f'PCAVarianceRatio_{suffix}.png'))

                # Find optimal PCs by using maximum distance from endpoint line
                var = self.adatasWhole.uns['pca']['variance_ratio'].copy()
                pcList = range(1, len(var)+1)
                
                slope = (var[-1]-var[0])/(pcList[-1]-pcList[0])
                distsqPC = []
                for i in range(len(var)):
                    x_closest = (pcList[0]*slope + pcList[i]/slope + var[i] - var[0])/(slope + 1/slope)
                    y_closest = (x_closest-pcList[0]) * slope + var[0]
                    distsqPC.append((x_closest-pcList[i])**2 + (y_closest-var[i])**2)
                
                recPCs = argmax(distsqPC)+1
                
                inputSetPage['figPCA'] = figPCA
                inputSetPage['recPCs'] = recPCs # Recommended PCs for integration
                inputSetPage['multiSample'] = len(self.adatasWhole.obs['sample'].cat.categories)>1

                pltclose('all')

            except Exception as e:
                print(e)

        # Page 2 Logic: Quality Control
        elif currentPage == 2:
            try:
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines = [
                        'Page 2: Quality Control\n',
                        f'Use Scrublet = {self.paramDict["scrublet"]}\n',
                        f'Enriched Mitochondrial Gene, Maximum % = {self.paramDict["perMt"]}\n',
                        f'Enriched Ribosomal Gene, Maximum % = {self.paramDict["perRib"]}\n',
                        f'Filter Genes, Minimum Cells = {self.paramDict["minCells"]}\n',
                        f'Filter Cells, Minimmum Genes = {self.paramDict["minGenes"]}\n\n'
                    ]
                    log.writelines(lines)

                self.adatasRaw, self.adatasBeforeHVG = quality_control(
                    self.adatasRaw.copy(), 
                    self.paramDict['minCells'], 
                    self.paramDict['minGenes'], 
                    self.paramDict['perMt'], 
                    self.paramDict['perRib'], 
                    self.paramDict['seed'], 
                    self.paramDict['scrublet'], 
                    self.paramDict['output_path']
                    )
                
                uniqueSamples = self.adatasBeforeHVG.obs['sample'].unique()
                numCellsRaw = []
                avgGenesRaw = []
                numCellsQC = []
                avgGenesQC = []
                for i in range(len(uniqueSamples)):
                    # Metrics before QC
                    numCellsRaw.append(len(self.adatasRaw[self.adatasRaw.obs['sample'] == uniqueSamples[i]].obs_names))
                    avgGenesRaw.append(average(self.adatasRaw[self.adatasRaw.obs['sample'] == uniqueSamples[i]].obs['n_genes_by_counts']))
                    # Metrics after QC
                    numCellsQC.append(len(self.adatasBeforeHVG[self.adatasBeforeHVG.obs['sample'] == uniqueSamples[i]].obs_names))
                    avgGenesQC.append(average(self.adatasBeforeHVG[self.adatasBeforeHVG.obs['sample'] == uniqueSamples[i]].obs['n_genes_by_counts']))

                inputSetPage['uniqueSamples'] = uniqueSamples
                inputSetPage['numCellsRaw'] = numCellsRaw
                inputSetPage['avgGenesRaw'] = avgGenesRaw
                inputSetPage['numCellsQC'] = numCellsQC
                inputSetPage['avgGenesQC'] = avgGenesQC

            except Exception as e:
                print(e)
                # *** Err invalid number

        # Page 1 Logic: Input Folders & Read
        elif currentPage == 1:
            # Try read to check if paths are valid
            try:
                if not pathexists(self.paramDict['output_path']):
                    makedirs(self.paramDict['output_path'])

                # Modality Params, Check Spatial
                self.h5adImport = False
                
                with open(join(self.paramDict['output_path'], 'iSNAP User Log.txt'), 'a+') as log:

                    lines1 = [
                        'iSNAP User Log\n\n',
                        'Start of new log.\n\n',
                        f'Modality = {self.paramDict["modality"]}\n',
                        f'Random State = {self.paramDict["seed"]}\n\n',
                        'Page 1: File Input\n',
                        'Input Paths =\n'
                    ]
                    for path in self.paramDict['input_paths']:
                        lines1.append(f'\t {path}\n')
                    
                    lines2 = [
                        f'Output Path = {self.paramDict["output_path"]}\n',
                        f'Use Spatial Features = {self.paramDict["isSpatial"]}\n\n'
                    ]
                    lines = lines1+lines2
                    log.writelines(lines)

                # Read data and save as adatas and adatasRaw
                self.adatasRaw = read_data(self.paramDict['input_paths'], self.paramDict['modality'])

                self.adatasRaw = calc_qc_metrics(self.adatasRaw)

                figCells = Figure(constrained_layout=True)
                axCells = figCells.add_subplot(111)
                figGenes = Figure(constrained_layout=True)
                axGenes = figGenes.add_subplot(111)
                figMT = Figure(constrained_layout=True)
                axMT = figMT.add_subplot(111)
                figRP = Figure(constrained_layout=True)
                axRP = figRP.add_subplot(111)

                violin(self.adatasRaw, 'total_counts', log=True, ylabel='Cell Counts', xlabel='Gene Counts', ax=axCells, stripplot=False, show=False)
                violin(self.adatasRaw, 'n_genes_by_counts', log=True, ylabel='Gene Counts', xlabel='Cell Counts', ax=axGenes, stripplot=False, show=False)
                violin(self.adatasRaw, 'pct_counts_mito', log=True, ylabel='Percent Mitochondrial Genes', xlabel='Cell Counts', ax=axMT, stripplot=False, show=False)
                violin(self.adatasRaw, 'pct_counts_RP', log=True, ylabel='Percent Ribosomal Genes', xlabel='Cell Counts', ax=axRP, stripplot=False, show=False)

                axCells.set(xticklabels=[])
                axGenes.set(xticklabels=[])
                axMT.set(xticklabels=[])
                axRP.set(xticklabels=[])

                outpath = join(self.paramDict['output_path'], 'Figures', 'Quality Control Violins')
                if not pathexists(outpath):
                    makedirs(outpath)

                figCells.savefig(join(outpath, 'Cell Counts.png'))
                figGenes.savefig(join(outpath, 'Gene Counts.png'))
                figMT.savefig(join(outpath, 'Percent Mitochondrial.png'))
                figRP.savefig(join(outpath, 'Percent Ribosomal.png'))

                inputSetPage['figCells'] = figCells
                inputSetPage['figGenes'] = figGenes
                inputSetPage['figMT'] = figMT
                inputSetPage['figRP'] = figRP
                inputSetPage['input_paths'] = self.paramDict['input_paths'] 
                inputSetPage['modality'] = self.paramDict['modality']

                if self.paramDict['modality'] == 'Single .h5ad file':
                    if 'celltype' in self.adatasRaw.obs.columns:
                        self.skipSignal.emit()

                    else:
                        print('Unprocessed Anndata detected, starting new analysis.')
            except Exception as e:
                print(e)

        # Page 0 Logic: Set Data Source
        elif currentPage == 0:
            inputSetPage['modality'] = self.paramDict['modality']

        self.toNextpgSet.emit(inputSetPage, currentPage)


class LoadingScreen(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("iSNAP - Loading...")
        self.setGeometry(420, 280, 640, 360)

        self.layout = QVBoxLayout()
        self.label = QLabel("Loading, please wait...")

        self.layout.addWidget(self.label)

        self.setLayout(self.layout)


class skipToCelltype(QWidget):
    def __init__(self):
        super().__init__()    
        try: 
            self.setWindowTitle("iSNAP - Processed .h5ad Detected")
            self.setGeometry(100, 100, 1280, 720)
            
            self.layoutMain = QVBoxLayout()
            self.layoutBtn = QHBoxLayout()

            self.labelChoice = QLabel('You have read a iSNAP processed .h5ad file. \nWould you like to continue to post-annotation pages?')
            self.btnInitiate = QPushButton('Start Over')
            self.btnInitiate.clicked.connect(lambda: self.skipTo(False))
            self.btnGoToCellType = QPushButton('Continue')
            self.btnGoToCellType.clicked.connect(lambda: self.skipTo(True))

            self.layoutBtn.addWidget(self.btnInitiate, alignment = Qt.AlignmentFlag.AlignLeft)
            self.layoutBtn.addWidget(self.btnGoToCellType, alignment = Qt.AlignmentFlag.AlignRight)

            self.layoutMain.addWidget(self.labelChoice)
            self.layoutMain.addLayout(self.layoutBtn)

            self.setLayout(self.layoutMain)

            self.show()
        except Exception as e:
            print(e)

    def skipTo(self, skip):
        if skip:
            window.nextpgStart(currentPage=8, h5adDetected=True)
        else:
            pass
        self.close()


app = QCoreApplication.instance()
if app is None:
    app = QApplication(argv)
apply_stylesheet(app, theme='light_blue.xml', invert_secondary=True)
qss_for_dropdown = """
        QComboBox {
            /* Example: Restore some desired padding for the selected text in the main box */
            padding-right: 32px; /* Make space for the dropdown arrow if it's too close */
            /* Add any other desired styling for the main combo box */
        }
        QComboBox::drop-down{
        border-left: 0px;
        padding-left: 0px;
        }
        
        """
app.setStyleSheet(app.styleSheet()+qss_for_dropdown)
font = QFont()
font.setPointSize(12)
app.setFont(font)
window = MainWindow()
window.show()
print('\niSNAP window was created.\n')
app.exec()

print('\nGood bye, iSNAP.')