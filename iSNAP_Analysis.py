###################
# Import Packages #
###################
print('- Importing os...')
from os import (
    makedirs
)
from os.path import (
    exists as pathexists, 
    join
)

print('- Importing Scanpy...')
from scanpy.tools import (
    rank_genes_groups,
)

print('- Importing pandas...')
from pandas import DataFrame


print('- Importing matplotlib...')
from matplotlib.pyplot import (
    close as pltclose
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

print('- Importing PyQt...')
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QWidget, 
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QComboBox,
    QPushButton,
    QPlainTextEdit,
    QSizePolicy,
    QSpacerItem
)

print('- Importing Scanpy...')
from scanpy.plotting import (
    umap as pltumap
)
from scanpy.plotting.palettes import default_102

#############
# Functions #
#############

def deg(adatas, DEGmethod):
    print('Ranking DEG for each cluster...')
    rank_genes_groups(adatas, groupby='leiden', method=DEGmethod, key_added='DEG', use_raw=False)
    adatas.uns['DEGMethod'] = DEGmethod

def degtoCSV(adatas, output_folder, DEGmethod, suffix='Initial'):
    # Extract results
    print('Saving DEG to CSV...')
    result = adatas.uns['DEG']
    groups = result['names'].dtype.names

    # Create an empty list to collect rows
    records = []

    # Loop through each cluster and build rows
    if DEGmethod == 'logreg':
        for group in groups:
            for i in range(len(result['names'][group])):
                records.append({
                    'cluster': group,
                    'gene': result['names'][group][i],
                    'score': result['scores'][group][i],
                })
    else:
        for group in groups:
            for i in range(len(result['names'][group])):
                records.append({
                    'cluster': group,
                    'gene': result['names'][group][i],
                    'logfoldchange': result['logfoldchanges'][group][i],
                    'score': result['scores'][group][i],
                    'pval': result['pvals'][group][i],
                    'pval_adj': result['pvals_adj'][group][i],
                })
    
    # Convert to dataframe
    deg_table = DataFrame(records)
    
    # Save to CSV
    outpath = join(output_folder, 'DEG Table')

    if not pathexists(outpath):
        makedirs(outpath)
    deg_table.to_csv(join(outpath, f'DEGTable {suffix}.csv'), index=False)

def specificClusterUMAP(adatas, focusType, key='cell_type', suffix = ''):
    clusterMask = adatas.obs[key] == focusType
    size = 240000/adatas.n_obs
    figUMAP = pltumap(adatas, show=False, color=key, mask_obs=clusterMask, title=f'{focusType}, {suffix}', legend_loc='none', frameon=False, palette=default_102, size=size, return_fig=True)

    pltclose('all')

    return figUMAP
    
def leidenToCellType(adatas, cellTypes, column='cell_type'):
    clusters = adatas.uns['DEG']['names'].dtype.names
    if len(clusters) == len(cellTypes):
        cellTypeMap = {}
        for i in range(len(clusters)):
            cellTypeMap[clusters[i]] = cellTypes[i]
            
        if 'cell_type' in adatas.obs.columns and column=='cell_type':
            j=1
            while f'celltypeBackup{j}' in adatas.obs.columns:
                j += 1
            adatas.obs[f'celltypeBackup{j}'] = adatas.obs['cell_type']
        
        adatas.obs[column] = adatas.obs['leiden'].copy().map(cellTypeMap).astype('category')

    else:
        raise ValueError('Leiden and Annotations Misaligned.')

def typeToCSV(adatas, groupKey, suffix, outfolder):
    savePathOut = join(outfolder, 'Xenium Labels', suffix)
    if not pathexists(savePathOut):
        makedirs(savePathOut)
    for sample in adatas.obs['sample'].cat.categories:
        print(f'Exporting {sample} labels...')
        adatasSample = adatas[adatas.obs['sample']==sample]
        type_df = DataFrame({
            'cell_id': adatasSample.obs['cell_id'], # Cell Barcode
            'group': adatasSample.obs[groupKey] # Cell Type Annotations
            })
        i=1
        # Save csv
        if pathexists(join(savePathOut, f'{sample} Xenium Labels {suffix}.csv')):
            while pathexists(join(savePathOut, f'{sample} Xenium Labels {suffix}({i}).csv')):
                i += 1
            
            type_df.to_csv(join(savePathOut + f'{sample} Xenium Labels {suffix}({i}).csv'), index=False)
        else:
            type_df.to_csv(join(savePathOut, f'{sample} Xenium Labels {suffix}.csv'))
        
    print('Done!')


###########
# Widgets #
###########

class DEGDotplot(QWidget):
    def __init__(self):
        super().__init__()

        # Assign Layouts
        self.layoutMain = QVBoxLayout()
        
        # Create Widget
        self.labelHeatmap = QLabel('Differentially Expressed Genes Heatmap')

        #Initiate Layout
        self.setLayout(self.layoutMain)

    
    def setPage(self, figDEG):
        self.DEGDotplot = FigureCanvasQTAgg(figDEG)
        self.DEGDotplot.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self.clearLayout(self.layoutMain)

        self.layoutMain.addWidget(self.labelHeatmap)
        self.layoutMain.addWidget(self.DEGDotplot)

        pltclose('all')
        
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)
                

class AnnotateCluster(QWidget):
    typeToCSVSignal = pyqtSignal(bool, str, str)
    toClusterSignal = pyqtSignal(int, int)
    toFeatureSignal = pyqtSignal(bool, bool)

    def __init__(self, isMain):
        super().__init__()
        self.isMain = isMain

        # Assign Layouts
        self.layoutUMAP = QVBoxLayout()
        self.layoutMain = QHBoxLayout()
        self.layoutSub = QVBoxLayout()
        
        # Create Widget
        self.spacer = QSpacerItem(1, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)
        self.labelExplore = QLabel('Exploration:')
        self.btnSpatial = QPushButton('Save Xenium Labels')
        self.btnFeature = QPushButton('Explore Gene Markers')
        self.btnFeature.clicked.connect(lambda: self.toFeatureSignal.emit(False, self.isMain))
        listNTop = ['50', '100', '200']
        self.comboNTop = QComboBox()
        self.comboNTop.addItems(listNTop)
        self.labelNTop = QLabel('Number of top ranked genes to show ')
        self.btnAnnotate = QPushButton('Annotate')
        self.btnAnnotate.clicked.connect(self.toClusterSignalFunc)

        #Initiate Layout
        self.setLayout(self.layoutMain)

    def toClusterSignalFunc(self):
        self.toClusterSignal.emit(self.comboCluster.currentIndex(), int(self.comboNTop.currentText()))

    def setUMAP(self, figUMAP):

        self.UMAPAnnotate = FigureCanvasQTAgg(figUMAP)
        self.toolbar = NavigationToolbar2QT(self.UMAPAnnotate, self)
        self.UMAPAnnotate.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        pltclose('all')

    def setPage(self, groups, isSpatial):        
        self.comboCluster = QComboBox()

        for i in groups:
            self.comboCluster.addItem(f'{i}')
        
        self.btnSpatial = QPushButton('Save Xenium Labels')
        self.btnSpatial.clicked.connect(lambda: self.typeToCSVSignal.emit(self.isMain, 'tempType', 'Annotation'))
        self.clearLayout(self.layoutSub)
        self.layoutSub.addItem(self.spacer)
        self.layoutSub.addWidget(self.labelExplore)
        self.layoutSub.addWidget(self.btnFeature)
        if isSpatial:
            self.layoutSub.addWidget(self.btnSpatial)
        self.layoutSub.addItem(self.spacer)
        self.layoutSub.addWidget(self.labelNTop, alignment=Qt.AlignmentFlag.AlignBottom)
        self.layoutSub.addWidget(self.comboNTop)
        self.layoutSub.addWidget(self.comboCluster)
        self.layoutSub.addWidget(self.btnAnnotate)
        self.layoutSub.addItem(self.spacer)

        self.comboCluster.setMaximumWidth(300)

        self.clearLayout(self.layoutUMAP)
        self.layoutUMAP.addWidget(self.UMAPAnnotate)
        self.layoutUMAP.addWidget(self.toolbar)

        self.clearLayout(self.layoutMain)
        self.layoutMain.addLayout(self.layoutUMAP)
        self.layoutMain.addLayout(self.layoutSub)      
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)

                
class SpecificCluster(QWidget):
    toAnnotateSignal = pyqtSignal(int)
    def __init__(self, parent, group, group_index):
        super(SpecificCluster, self).__init__(parent)
        # What group index is it
        self.group = group
        self.group_index = group_index
        # Assign Layouts
        self.layoutMain = QHBoxLayout()
        self.layoutUMAP = QVBoxLayout()
        self.layoutAnnotate = QHBoxLayout()
        self.layoutSub = QVBoxLayout()
        
        # Create Widget
        self.labelDEG = QLabel('Genes Ranked by Significance and Up-Regulation')
        self.btnConfirm = QPushButton('Confirm')
        self.btnConfirm.clicked.connect(self.toAnnotateSignalFunc)
        
        self.labelAnnotate = QLabel('Enter Cell Type')
        self.lineAnnotate = QLineEdit(f'{group}')

        #Initiate Layout
        self.setLayout(self.layoutMain)
    
    def toAnnotateSignalFunc(self):
        self.toAnnotateSignal.emit(int(self.group_index))

    def setName(self, newName):
        self.group = newName

    def setUMAP(self, figUMAP):
        self.UMAPSpecific = FigureCanvasQTAgg(figUMAP)
        self.toolbar = NavigationToolbar2QT(self.UMAPSpecific, self)
        self.UMAPSpecific.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        pltclose('all')

    def setPage(self, topNames):
        self.listDEG = QPlainTextEdit(readOnly=True)
        for gene in topNames:
            self.listDEG.appendPlainText(gene)

        self.clearLayout(self.layoutAnnotate)
        self.layoutAnnotate.addWidget(self.labelAnnotate)
        self.layoutAnnotate.addWidget(self.lineAnnotate)

        self.clearLayout(self.layoutSub)
        self.layoutSub.addWidget(self.labelDEG)
        self.layoutSub.addWidget(self.listDEG)
        self.layoutSub.addLayout(self.layoutAnnotate)
        self.layoutSub.addWidget(self.btnConfirm)
        self.layoutSub.setSizeConstraint
        
        self.clearLayout(self.layoutUMAP)
        self.layoutUMAP.addWidget(self.UMAPSpecific)
        self.layoutUMAP.addWidget(self.toolbar)

        self.clearLayout(self.layoutMain)
        self.layoutMain.addLayout(self.layoutUMAP)
        self.layoutMain.addLayout(self.layoutSub)
        
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)

class UMAPCellType(QWidget):
    typeToCSVSignal = pyqtSignal(bool, str, str)
    combTypesSignal = pyqtSignal(str, str, str)
    toFeatureSignal = pyqtSignal(bool, bool)
    toSubClusterSignal = pyqtSignal(str)

    def __init__(self):
        super().__init__()

        # Assign Layouts
        self.layoutMain = QHBoxLayout()
        self.layoutUMAP = QVBoxLayout()
        self.layoutSelCluster = QVBoxLayout()
        self.layoutComb = QVBoxLayout()
        self.layoutCombsub = QHBoxLayout()
        
        # Create Widget
        self.spacer = QSpacerItem(1, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)
        
        self.labelSubClust = QLabel('Select a SubCluster to annotate or \nClick Finish to save')
        self.btnConfirm = QPushButton('Subcluster')
        self.btnConfirm.clicked.connect(lambda: self.toSubClusterSignal.emit(self.comboType.currentText()))
        self.labelComb = QLabel('Combine Types')
        self.btnComb = QPushButton('Combine')
        self.btnComb.clicked.connect(self.combTypes)
        self.labelExplore = QLabel('Exploration:')
        self.btnFeature = QPushButton('Explore Gene Markers')
        self.btnFeature.clicked.connect(lambda: self.toFeatureSignal.emit(True, True))
        
        #Initiate Layout
        self.setLayout(self.layoutMain)

    def setPage(self, figUMAP, celltypes, isSpatial):
        self.UMAPCellType = FigureCanvasQTAgg(figUMAP)
        self.toolbar = NavigationToolbar2QT(self.UMAPCellType, self)

        self.UMAPCellType.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self.comboType = QComboBox()
        self.comboType.setMaximumWidth(300)
        self.comboComb1 = QComboBox()
        self.comboComb1.setMaximumWidth(150)
        self.comboComb2 = QComboBox()
        self.comboComb2.setMaximumWidth(150)
        self.lineComb = QLineEdit()
        

    
        for i in range(len(celltypes)):
            self.comboType.addItem(celltypes[i])     
            self.comboComb1.addItem(celltypes[i])    
            self.comboComb2.addItem(celltypes[i])    
        
        self.btnSpatial = QPushButton('Save Xenium Labels')
        self.btnSpatial.clicked.connect(lambda: self.typeToCSVSignal.emit(True, 'cell_type', 'Cell Types'))

        pltclose('all')
        self.clearLayout(self.layoutCombsub)
        self.layoutCombsub.addWidget(self.comboComb1)
        self.layoutCombsub.addWidget(self.comboComb2)

        self.clearLayout(self.layoutComb)
        self.layoutComb.addWidget(self.labelComb)
        self.layoutComb.addLayout(self.layoutCombsub)
        self.layoutComb.addWidget(self.lineComb)
        self.layoutComb.addWidget(self.btnComb)

        self.clearLayout(self.layoutSelCluster)
        self.layoutSelCluster.addItem(self.spacer)
        self.layoutSelCluster.addWidget(self.labelSubClust, alignment=Qt.AlignmentFlag.AlignBottom)
        self.layoutSelCluster.addWidget(self.comboType)
        self.layoutSelCluster.addWidget(self.btnConfirm)
        self.layoutSelCluster.addItem(self.spacer)
        self.layoutSelCluster.addWidget(self.labelExplore)
        self.layoutSelCluster.addWidget(self.btnFeature)
        if isSpatial:
            self.layoutSelCluster.addWidget(self.btnSpatial)
        self.layoutSelCluster.addItem(self.spacer)
        self.layoutSelCluster.addLayout(self.layoutComb)
        self.layoutSelCluster.addItem(self.spacer)
        
        self.clearLayout(self.layoutUMAP)
        self.layoutUMAP.addWidget(self.UMAPCellType)
        self.layoutUMAP.addWidget(self.toolbar)

        self.clearLayout(self.layoutMain)
        self.layoutMain.addLayout(self.layoutUMAP, stretch=5)
        self.layoutMain.addLayout(self.layoutSelCluster, stretch=1)
    
    def combTypes(self):
        if not self.lineComb.text() == '':
            type1 = self.comboComb1.currentText()
            type2 = self.comboComb2.currentText()
            newName = self.lineComb.text()

            self.combTypesSignal.emit(type1, type2, newName)  
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)

                
class FeatureMap(QWidget):
    updateFeatureSignal = pyqtSignal(bool, str, str, str)
    def __init__(self, UMAP, genesList, isMain):
        super().__init__()
        self.setWindowTitle("iSNAP - Feature Maps")
        self.setGeometry(200, 200, 1280, 720)

        # Add to Figure for each leiden resolution & Page
        figUMAP = FigureCanvasQTAgg(UMAP)
        toolbar = NavigationToolbar2QT(figUMAP, self)
        
        self.isMain = isMain
        self.genesList = genesList

        # Assign Layouts
        self.layoutMain = QHBoxLayout()
        self.layoutUMAP = QVBoxLayout()
        self.layoutSub = QVBoxLayout()
        
        # Create Widget
        self.spacer = QSpacerItem(1, 1, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)
        self.labelGene = QLabel('Type Gene Marker')
        self.inGene = QLineEdit('')
        self.btnTyped = QPushButton('Plot Typed Marker')
        self.btnTyped.clicked.connect(lambda: self.updateFeatureSignal.emit(self.isMain, self.inGene.text(), self.inVmax.text(), self.inVmin.text()))
        self.labelOR = QLabel('Or Select Gene Marker')
        self.comboGene = QComboBox()
        self.comboGene.addItems(self.genesList)
        self.btnCombo = QPushButton('Plot Selected Marker')
        self.btnCombo.clicked.connect(lambda: self.updateFeatureSignal.emit(self.isMain, self.comboGene.currentText(), self.inVmax.text(), self.inVmin.text()))
        self.labelVmax = QLabel('Set Vmax (Optional)')
        self.inVmax = QLineEdit()
        self.labelVmin = QLabel('Set Vmin (Optional)')
        self.inVmin = QLineEdit()

        # Add to layouts
        self.layoutSub.addItem(self.spacer)
        self.layoutSub.addWidget(self.labelGene)
        self.layoutSub.addWidget(self.inGene)
        self.layoutSub.addWidget(self.btnTyped)
        self.layoutSub.addItem(self.spacer)
        self.layoutSub.addWidget(self.labelOR)
        self.layoutSub.addWidget(self.comboGene)
        self.layoutSub.addWidget(self.btnCombo)
        self.layoutSub.addItem(self.spacer)
        self.layoutSub.addWidget(self.labelVmax)
        self.layoutSub.addWidget(self.inVmax)
        self.layoutSub.addWidget(self.labelVmin)
        self.layoutSub.addWidget(self.inVmin)
        self.layoutSub.addItem(self.spacer)

        self.layoutUMAP.addWidget(figUMAP)
        self.layoutUMAP.addWidget(toolbar)
        
        self.layoutMain.addLayout(self.layoutUMAP)
        self.layoutMain.addLayout(self.layoutSub)

        #Initiate Layout
        self.setLayout(self.layoutMain)
    
    def updateUMAP(self, figUMAP):
        try:
            UMAP = FigureCanvasQTAgg(figUMAP)
            toolbar = NavigationToolbar2QT(UMAP, self)
            UMAP.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

            self.clearLayout(self.layoutUMAP)
            self.layoutUMAP.addWidget(UMAP)
            self.layoutUMAP.addWidget(toolbar)

            self.clearLayout(self.layoutMain)
            self.layoutMain.addLayout(self.layoutUMAP)
            self.layoutMain.addLayout(self.layoutSub)

            #Set Layout
            self.setLayout(self.layoutMain)

            pltclose('all')
                
        except Exception as e:
            print('Invalid gene.')
            print(e)
        
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)
        
    