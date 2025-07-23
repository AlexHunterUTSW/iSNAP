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
from scanpy.preprocessing import (
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    normalize_total,
    log1p,
    highly_variable_genes,
    scale,
    neighbors,
    pca
)
from scanpy.tools import (
    umap,
    leiden
)
from scanpy.external.pp import harmony_integrate


print('- Importing Scrublet...')
from scrublet import Scrublet

print('- Importing sklearn')
from sklearn.metrics import silhouette_score

print('- Importing Scipy')
from scipy.sparse import issparse

print('- Importing numpy...')
from numpy import (
    nan,
    float32
)
print('- Importing pandas...')
from pandas import DataFrame

print('- Importing matplotlib...')
from matplotlib.pyplot import (
    close as pltclose,
    figure as pltfigure
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

print('- Importing BBKNN...')
from bbknn import bbknn

print('- Importing PyQt...')
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import ( 
    QWidget, 
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QLineEdit,
    QCheckBox,
    QComboBox,
    QRadioButton,
    QSizePolicy,
    QSpacerItem
)

#############
# Functions #
#############
def calc_qc_metrics(adatas):
    # Mitochondrial genes
    adatas.var['mito'] = adatas.var_names.str.startswith(('MT-','mt-'))
    
    # Rib genes
    adatas.var['RP'] = adatas.var_names.str.startswith(('RPS','RPL'))
    # *** In input, set an option for mice or human: mice is 'MT', human is 'mt'

    # Calculates metrics for qc, No ranking top genes
    calculate_qc_metrics(adatas, percent_top=None, qc_vars = ['mito', 'RP'], inplace=True) 
    adatas.obs['pct_counts_mito'] = adatas.obs['pct_counts_mito'].replace('', nan).fillna(0)
    adatas.obs['pct_counts_RP'] = adatas.obs['pct_counts_RP'].replace('', nan).fillna(0)

    return adatas

def quality_control(adatas, min_cells, min_genes, perMt, perRib, seed, scrublet, output_folder):
    """
    Calculates quality control metrics and filters transcriptomic dataset.

    Args:
        adatas (obs): Object containing raw transcriptomic data.
        min_cells (int): Cell count threshold for filtering genes.
        min_genes (int): Gene count threshold for filtering.
        output_folder (str): Pathway to save .h5ad of adatas.
       

    Returns:
        obs: Quality controlled transcriptomic data.
    """
    # Initial cell ids so that removed cells can be documented in Filtered_Cells.csv
    adatasRaw = adatas.copy()
    initial_cells = adatasRaw.obs_names.copy()

    if scrublet == True:
        # Use raw counts if available (Scrublet needs raw counts)
        
        counts = adatas.raw.X.copy() if adatas.raw is not None else adatas.X.copy()
        
        # Make sure matrix is in correct format
        counts = counts.tocsc().astype(float32)
        if issparse(counts):
            counts = counts.tocsc().astype(float32)
        else:
            counts = counts.astype(float32)
        # Run Scrublet
        scrub = Scrublet(counts, random_state=seed)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()

        # Add results to AnnData object
        adatas.obs['doublet_score'] = doublet_scores
        adatas.obs['predicted_doublet'] = predicted_doublets

        # Filter Doublets
        adatas = adatas[~adatas.obs['predicted_doublet']].copy()

    # Filter w/ metrics
    adatas = adatas[adatas.obs['pct_counts_mito'] <= perMt].copy() # Percentage Mitochondrial Gene
    adatas = adatas[adatas.obs['pct_counts_RP'] <= perRib].copy() # Percentage Ribosomal
    
    filter_cells(adatas, min_genes=min_genes) # Filter cells below minimum gene count threshold
    filter_genes(adatas, min_counts=min_cells) # Filter genes below minimum count threshold

    after_cells = adatas.obs_names.copy()

    removed_cells = [cell for cell in initial_cells if cell not in after_cells]
    
    adatas_filtered = adatasRaw[adatasRaw.obs_names.isin(removed_cells)].copy()

    df_removed = DataFrame()
    df_removed['sample'] = adatas_filtered.obs['sample'].copy()
    df_removed['cell id'] = adatas_filtered.obs['cell_id'].copy()
    df_removed['gene counts'] = adatas_filtered.obs['total_counts'].copy()
    df_removed['pct mito'] = adatas_filtered.obs['pct_counts_mito'].copy()
    df_removed['pct ribo'] = adatas_filtered.obs['pct_counts_RP'].copy()

    outfolder = join(output_folder, 'Filtered Cells Table')

    if not pathexists(outfolder):
        makedirs(outfolder)
    df_removed.to_csv(join(outfolder,f'Filtered Cells.csv'), index=False)
    
    return adatasRaw, adatas


def normalization(adatas, boolHVG, nHVG, seed):
    """
    Normalizes and logs transcriptomic data.

    Args:
        adatas (obs): Object containing quality controlled transcriptomic data.
        norm_target (int): Normalizing target. Tends to be factors of 10. 1 = unit, 100 = percentage. Default = None.

       

    Returns:
        obs: Normalized transcriptomic data.
    """
    adatas.layers["counts_raw"] = adatas.X.copy() # Save raw count data

    normalize_total(adatas,target_sum=1) # Normalize to median total counts
    log1p(adatas) # Logarithmize data
    
    if boolHVG:
        highly_variable_genes(adatas, n_top_genes=nHVG, batch_key="sample")

    scale(adatas, zero_center=False)

    adatas.layers["counts_norm"] = adatas.X.copy() # Save normalized count data

    pca(adatas, random_state=seed) # Principle components analysis

    return adatas

def pcaVarianceRatio(adatas, n_pcs=50, pca_key='pca', dpi=200):
    varRatio = adatas.uns[pca_key]['variance_ratio'][:n_pcs]
    figPCA = pltfigure(constrained_layout=True, dpi=dpi)
    axPCA = figPCA.add_subplot(111)

    for i in range(n_pcs):
        axPCA.text(i, varRatio[i], f'PC{i}', fontsize=6, rotation=90, ha='center', va='center')

    axPCA.set_title('Variance Ratio vs PC Number')
    axPCA.set_ylabel('Variance Ratio')
    axPCA.set_xlabel('PC Number')
    axPCA.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False
    ) # labels along the bottom edge are off
    axPCA.set_xlim(-1, n_pcs+1)
    axPCA.set_ylim(0, max(varRatio)*1.1)
    axPCA.margins(x=2, y=2)

    pltclose('all')

    return figPCA 


def integration(adatas, npcs, nNeigh, intMeth, seed):
    """
    Integrates data with multiple samples. This allows for handling errors that may occur from sample to sample before analysis.

    Args:
        adatas (obs): Transcriptomic data that have been quality controlled and normalized.
       

    Returns:
        obs: Fully processed transcriptomic data that are ready to be analyzed.
    """
    if len(set(adatas.obs['sample'])) > 1:
        if intMeth == 'Harmony':
            print('Integrating using Harmony...')
            harmony_integrate(adatas, key = 'sample', adjusted_basis='X_pca', random_state=seed)
            print('Computing Nearest Neighbors...')
            neighbors(adatas, n_neighbors=nNeigh, n_pcs = npcs, random_state=seed, use_rep='X_pca')
        
        elif intMeth == 'BBKNN':
            print('Integrating using BBKNN...')
            bbknn(adatas, batch_key='sample', pynndescent_n_neighbors=nNeigh, n_pcs=npcs, pynndescent_random_state=seed)
        
        elif intMeth == 'No Integration':
            print('Computing Nearest Neighbors...')
            neighbors(adatas, n_neighbors=nNeigh, n_pcs = npcs, random_state=seed)
        else:
            print('err')
    
    else:
        neighbors(adatas, n_neighbors=nNeigh, n_pcs = npcs, random_state=seed)

    umap(adatas, random_state=seed)

    return adatas


def cluster(adatas, res, seed):
    leiden(adatas, key_added="leiden", resolution=res, random_state=seed)
    leidenParts = [adatas.obs['leiden'], adatas.uns['leiden']]
    return leidenParts

def silhouetteLeiden(adatas):
    return silhouette_score(
        adatas.obsm['X_pca'], 
        adatas.obs['leiden']
        )


##########
# Widget #
##########

class TableNormalize(QWidget):
    def __init__(self):
        super().__init__()

        # Assign Layouts
        self.layoutMain = QVBoxLayout()
        self.layoutParam = QGridLayout()

        # Create Widgets
        self.hspacer =  QSpacerItem(1, 1, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.labelTable = QLabel('Quality Control Results:')
        self.tableQC = QTableWidget()
        self.labelHVG = QLabel('Highly Variable Genes')
        self.checkHVG = QCheckBox()
        self.checkLabelHVG = QLabel('Filter highly variable genes')
        self.inHVG = QLineEdit('3000')
        self.labelInHVG = QLabel('Number of highly variable genes to keep')

        # Set Layout
        self.layoutParam.addWidget(self.checkHVG, 0, 0, alignment = Qt.AlignCenter)
        self.layoutParam.addWidget(self.checkLabelHVG, 0, 1, alignment = Qt.AlignLeft)
        self.layoutParam.addWidget(self.inHVG, 1, 0, alignment = Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelInHVG, 1, 1, alignment = Qt.AlignLeft)
        self.layoutParam.addItem(self.hspacer, 0, 2)
        
        self.layoutMain.addWidget(self.labelTable)
        self.layoutMain.addWidget(self.tableQC)
        self.layoutMain.addWidget(self.labelHVG)
        self.layoutMain.addLayout(self.layoutParam)

        self.setLayout(self.layoutMain)

    # Functions
    def createTable(self, uniqueSamples, numCellsRaw, avgGenesRaw, numCellsQC, avgGenesQC):
        self.tableQC.setColumnCount(5)
        self.tableQC.setRowCount(len(uniqueSamples))

        headers = ['Sample', 'Cell Count Before QC', 'Average Genes Before QC', 'Cell Count After QC', 'Average Genes After QC']
        # Set Column Headers
        self.tableQC.setHorizontalHeaderLabels(headers)

        # Loop for each unique sample to add data for each sample
        for i in range(len(uniqueSamples)):
            self.tableQC.setItem(i, 0, QTableWidgetItem(uniqueSamples[i]))

            # Add raw data
            self.tableQC.setItem(i, 1, QTableWidgetItem(str(numCellsRaw[i])))
            self.tableQC.setItem(i, 2, QTableWidgetItem(f'{avgGenesRaw[i]:.1f}'))

            # Add QC data
            self.tableQC.setItem(i, 3, QTableWidgetItem(str(numCellsQC[i])))
            self.tableQC.setItem(i, 4, QTableWidgetItem(f'{avgGenesQC[i]:.1f}'))

        self.tableQC.resizeColumnsToContents()
        for i in range(self.tableQC.columnCount()):
            self.tableQC.setColumnWidth(i, self.tableQC.columnWidth(i) + 16)
        
        
class PCAInt(QWidget):
    def __init__(self):
        super().__init__()

        # Assign Layouts
        self.layoutMain = QVBoxLayout()
        self.layoutParam = QGridLayout()
        

        # Create Widgets
        self.labelPCA = QLabel('Principal Component Analysis (PCA)')

        self.plotPCA = FigureCanvasQTAgg()
        
        self.hspacer =  QSpacerItem(1, 1, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.inNPC = QLineEdit('50')
        self.labelNPC = QLabel('Number of PCs')
        self.inNeigh = QLineEdit('15')
        self.labelNeigh = QLabel('Number of Neighbors')
        self.labelIntegration = QLabel('Integration')

        self.inInt = QComboBox()
        self.inInt.addItem('Harmony')
        self.inInt.addItem('BBKNN')
        self.inInt.addItem('No Integration')
        self.labelInt = QLabel('Integration Method')

        #Initiate Layout
        self.setLayout(self.layoutMain)

    def setPage(self, figPCA, recPCs, multiSample):
        self.inNPC.setText(str(recPCs))

        plotPCA = FigureCanvasQTAgg(figPCA)
        toolbar = NavigationToolbar2QT(plotPCA, self)
        self.plotPCA.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.clearLayout(self.layoutParam)
        self.layoutParam.addWidget(self.inNPC, 0, 0, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelNPC, 0, 1, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.inNeigh, 1, 0, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelNeigh, 1, 1, alignment=Qt.AlignLeft)
        if multiSample:
            self.layoutParam.addWidget(self.inInt, 2, 0, alignment=Qt.AlignLeft)
            self.layoutParam.addWidget(self.labelInt, 2, 1, alignment=Qt.AlignLeft)
        self.layoutParam.addItem(self.hspacer, 0, 2)

        self.clearLayout(self.layoutMain)
        self.layoutMain.addWidget(self.labelPCA)
        self.layoutMain.addWidget(plotPCA)
        self.layoutMain.addWidget(toolbar)
        self.layoutMain.addLayout(self.layoutParam)
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)


class TestLeiden(QWidget):
    def __init__(self):
        super().__init__()
        self.singleRes = False

        # Assign Layouts
        self.layoutMain = QVBoxLayout()
        self.layoutParam = QGridLayout()

        # Create Widgets
        self.plotUMAP = FigureCanvasQTAgg()
        self.plotUMAP.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        self.labelLeiden = QLabel('Leiden Settings:')
        self.inSingle = QCheckBox()
        self.inSingle.toggled.connect(self.leidenEnable)
        self.labelSingle = QLabel('Test one resolution')
        self.inRes1 = QLineEdit('1')
        self.labelRes1 = QLabel('Resolution for Leiden Clustering')
        self.inRes2 = QLineEdit('0.6')
        self.labelRes2 = QLabel('Resolution for Leiden Clustering')
        self.inRes3 = QLineEdit('0.3')
        self.labelRes3 = QLabel('Resolution for Leiden Clustering')
        self.hspacer =  QSpacerItem(1, 1, QSizePolicy.Expanding, QSizePolicy.Minimum)
        
        # Add widgets to layout
        self.layoutParam.addWidget(self.inSingle, 0, 0, alignment=Qt.AlignCenter)
        self.layoutParam.addWidget(self.labelSingle, 0, 1, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.inRes1, 1, 0, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelRes1, 1, 1, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.inRes2, 2, 0, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelRes2, 2, 1, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.inRes3, 3, 0, alignment=Qt.AlignLeft)
        self.layoutParam.addWidget(self.labelRes3, 3, 1, alignment=Qt.AlignLeft)
        self.layoutParam.addItem(self.hspacer, 0, 2)

        #Initiate Layout
        self.setLayout(self.layoutMain)

    def setPage(self, figUMAP):
        self.plotUMAP = FigureCanvasQTAgg(figUMAP)
        self.toolbar = NavigationToolbar2QT(self.plotUMAP, self)
        self.plotUMAP.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.clearLayout(self.layoutMain)
        self.layoutMain.addWidget(self.plotUMAP)
        self.layoutMain.addWidget(self.toolbar)
        self.layoutMain.addWidget(self.labelLeiden)
        self.layoutMain.addLayout(self.layoutParam)

        pltclose('all')
    
    def leidenEnable(self):
        self.singleRes = self.inSingle.isChecked()

        if self.singleRes:
            self.inRes2.setEnabled(False)
            self.inRes3.setEnabled(False)
        else:
            self.inRes2.setEnabled(True)
            self.inRes3.setEnabled(True)
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)


class SelectLeiden(QWidget):
    def __init__(self):
        super().__init__()

        # Assign Layouts
        self.layoutMain = QVBoxLayout()
        self.layoutUMAP = QHBoxLayout()
        self.layoutParam = QGridLayout()
        
        # Create Widget
        self.vspacer = QSpacerItem(1, 1, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.hspacer =  QSpacerItem(1, 1, QSizePolicy.Expanding, QSizePolicy.Minimum)
        
        self.labelDEGMeth = QLabel('Select Statistical Test for DEG')
        self.comboDEGMeth = QComboBox()
        self.comboDEGMeth.addItem('wilcoxon')
        self.comboDEGMeth.addItem('logreg')
        self.comboDEGMeth.addItem('t-test')
        self.comboDEGMeth.addItem('t-test_overestim_var')

        #Initiate Layout
        self.setLayout(self.layoutMain)

    
    def setPage(self, figList, score, resList, singleRes):
        self.clearLayout(self.layoutUMAP)
        self.clearLayout(self.layoutMain)
        self.clearLayout(self.layoutParam)

        self.plotUMAP = []
        for i in range(len(figList)):
            self.plotUMAP.append(FigureCanvasQTAgg(figList[i]))
            self.plotUMAP[i].setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.layoutUMAP.addWidget(self.plotUMAP[i], stretch=1)
        
        if not singleRes:
            labelLeiden = QLabel('Leiden Plots: Choose One')
            self.rbRes1 = QRadioButton('', self)
            self.labelRes1 = QLabel(f'Resolution {resList[0]} (Silhouette Score: {score[0]:.2f})')
            self.rbRes2 = QRadioButton('', self)
            self.labelRes2 = QLabel(f'Resolution {resList[1]} (Silhouette Score: {score[1]:.2f})')
            self.rbRes3 = QRadioButton('', self)
            self.labelRes3 = QLabel(f'Resolution {resList[2]} (Silhouette Score: {score[2]:.2f})')

            self.layoutParam.addWidget(labelLeiden, 0, 0)
            self.layoutParam.addWidget(self.rbRes1, 1, 0, alignment=Qt.AlignCenter)
            self.layoutParam.addWidget(self.labelRes1, 1, 1, alignment=Qt.AlignLeft)
            self.layoutParam.addWidget(self.rbRes2, 2, 0, alignment=Qt.AlignCenter)
            self.layoutParam.addWidget(self.labelRes2, 2, 1, alignment=Qt.AlignLeft)
            self.layoutParam.addWidget(self.rbRes3, 3, 0, alignment=Qt.AlignCenter)
            self.layoutParam.addWidget(self.labelRes3, 3, 1, alignment=Qt.AlignLeft)
            self.layoutParam.addWidget(self.comboDEGMeth, 4, 0, alignment=Qt.AlignLeft)
            self.layoutParam.addWidget(self.labelDEGMeth, 4, 1, alignment=Qt.AlignLeft)

        else:
            labelScore = QLabel(f'Silhouette Score: {score[0]:.2f}')
            self.layoutParam.addWidget(labelScore, 0, 0)
            self.layoutParam.addWidget(self.comboDEGMeth, 1, 0, alignment=Qt.AlignLeft)
            self.layoutParam.addWidget(self.labelDEGMeth, 1, 1, alignment=Qt.AlignLeft)

        self.layoutParam.addItem(self.hspacer, 0, 2)

        self.layoutMain.addLayout(self.layoutUMAP)
        self.layoutMain.addLayout(self.layoutParam)
        
        pltclose('all')
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)