###################
# Import Packages #
###################
print('- Importing os...')
from os import (
    makedirs
)
from os.path import (
    exists as pathexists,
    basename,
    join
)
print('- Importing Scanpy...')
from scanpy import read_10x_h5, read_10x_mtx, read_h5ad
print('- Importing anndata...')
from anndata import concat as adconcat
print('- Importing scipy...')
from scipy.sparse import csr_matrix

print('- Importing PyQt...')
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import ( 
    QWidget, 
    QVBoxLayout, 
    QHBoxLayout,
    QGridLayout,
    QListWidget, 
    QMessageBox, 
    QLineEdit, 
    QLabel,
    QCheckBox,
    QSpacerItem,
    QSizePolicy
)

print('- Importing matplotlib...')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

####################
# Read in datasets #
####################

def read_data(input_paths, modality):
    """
    Reads all samples.

    Args:
         input_paths (list): Strings of samples' directories. Must include cell_feature_matrix.h5 file. Must include cells.cvs if spatially resolved.
         is_spatial (bool): Whether the data is spatially resolved.
       

    Returns:
        obs: Raw Transcriptomic Data
    """
    adatas = [] # Initialize list of adata

    if modality == 'Xenium':
        for data_path in input_paths:
            dataset_name = basename(data_path) # Get name of sample
            print('Reading sample:', dataset_name)
            adata = read_10x_h5(join(data_path, 'cell_feature_matrix.h5')) # Read Sequencing Data

            adata.obs_names = [f'{cell}_{dataset_name}' for cell in adata.obs_names] # Rename cell name to include sample name
            
            adatas.append(adata) # Append individual adata to list, adatas

        adatas = adconcat(adatas, join='inner', label = 'sample') # Combine sample datasets into one. Only joins in common data. Adds sample label

        adatas.obs['cell_id'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[0]) # Cell ID
        adatas.obs['sample'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[1]).astype('category') # Takes only sample name from cell_sample format for sample labels
    
    elif modality == 'Cell Ranger':
        for data_path in input_paths:
            dataset_name = basename(data_path) # Get name of sample
            print('Reading sample:', dataset_name)  
        
            adata = read_10x_mtx(join(data_path,'outs','filtered_feature_bc_matrix'))

            adata.obs_names = [f'{cell}_{dataset_name}' for cell in adata.obs_names] # Rename cell name to include sample name
            
            adatas.append(adata) # Append individual adata to list, adatas
            
        adatas = adconcat(adatas, join='inner') # Combine sample datasets into one. Only joins in common data. Adds sample label

        adatas.obs['cell_id'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[0]) # Cell ID
        adatas.obs['sample'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[1]).astype('category') # Takes only sample name from cell_sample format for sample labels

    elif modality == 'Single .h5ad file':
        if len(input_paths) == 1:
            data_path = input_paths[0]
            dataset_name = basename(data_path)[:-5]
            print('Reading sample:', dataset_name)  
            adatas = read_h5ad(data_path)
            if not 'cell_id' in adatas.obs.columns or not 'sample' in adatas.obs.columns:
                adatas.obs_names = [f'{cell}_{dataset_name}' for cell in adatas.obs_names]
                adatas.obs['cell_id'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[0]) # Cell ID
                adatas.obs['sample'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[1]).astype('category') # Takes only sample name from cell_sample format for sample labels
            

        else:
            adatas = []
            for data_path in input_paths:
                adata = read_h5ad(data_path)
                if not 'cell_id' in adata.obs.columns or not 'sample' in adata.obs.columns:
                    dataset_name = basename(data_path)[:-5]
                    adata.obs_names = [f'{cell}_{dataset_name}' for cell in adata.obs_names]

                adatas.append(adata)

            adatas = adconcat(adatas, join='inner')

            adatas.obs['cell_id'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[0]) # Cell ID
            adatas.obs['sample'] = adatas.obs_names.to_series().apply(lambda x: x.split("_")[1]).astype('category') # Takes only sample name from cell_sample format for sample labels
            
        

    else:
        raise ValueError('Invalid Modality')
    
    adatas.X = csr_matrix(adatas.X) # Convert count matrix to sparse matrix
    adatas.raw = adatas.copy()

    return adatas

######################
# Save Adatas Object #
######################
def save_adatas(adatas, output_folder, suffix=None):
    """
    Reads all samples.

    Args:
         input_paths (list): Strings of samples' directories. Must include cell_feature_matrix.h5 file. Must include cells.cvs if spatially resolved.
         is_spatial (bool): Whether the data is spatially resolved.

    Returns:
        obs: Raw Transcriptomic Data
    """
    # adatas.X = csr_matrix(adatas.X) # Convert count matrix to sparse matrix
    print('Saving anndata object...')

    outpath = join(output_folder,'Adatas')

    if not pathexists(outpath):
        makedirs(outpath)

    if suffix:
        adatas.write(join(outpath, f'adatas {suffix}.h5ad')) # Save Xenium dataset
        print(f'The anndata object was saved to your output folder as "adatas {suffix}.h5ad".')
    else:
        adatas.write(join(outpath, "adatas.h5ad")) # Save Xenium dataset
        print('The anndata object was saved to your output folder as "adatas.h5ad".')
    

#####################
# Dataset Read Page #
#####################
class ReadQC(QWidget):
    def __init__(self):
        super().__init__()

        # Assign Layouts
        self.layoutMain = QVBoxLayout()
        self.layoutParam = QGridLayout()
        self.layoutScrublet = QHBoxLayout()
        self.layoutViolin = QHBoxLayout()
        

        # Assign Widgets #
        self.hspacer = QSpacerItem(1, 1, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        # Datasets Read
        self.labelRead = QLabel('Datasets Read')
        self.listDatasets = QListWidget()

        self.labelFilter = QLabel('Filtering')
        # Remove Doublets
        self.checkScrublet = QCheckBox()
        self.labelScrublet = QLabel('Remove Doublets')

        # Filter Percentage of mitochondrial genes
        self.inPerMt = QLineEdit('100')
        self.labelPerMt = QLabel('Max % Mitochondrial Genes')
        
        # Filter Percentage of ribosomal genes
        self.inPerRib = QLineEdit('100')
        self.labelPerRib = QLabel('Max % Ribosomal Genes')

        # Filter Minimum Cells
        self.inMinCells = QLineEdit('3')
        self.labelMinCells = QLabel('Minimum Cells')

        # Filter Minimum Genes
        self.inMinGenes = QLineEdit('100')
        self.labelMinGenes = QLabel('Minimum Genes')
        
        # Add Widgets to Layouts
        
        self.setLayout(self.layoutMain)

    def setPage(self, figCells, figGenes, figMT, figRP, inPaths, modality):
        while self.listDatasets.count() != 0:
            self.listDatasets.takeItem(0)

        for i in inPaths:
            self.listDatasets.addItem(basename(i))

        # Create Qt Figure
        figCells = FigureCanvasQTAgg(figCells)
        figGenes = FigureCanvasQTAgg(figGenes)
        figMT = FigureCanvasQTAgg(figMT)
        figRP = FigureCanvasQTAgg(figRP)
        
        self.clearLayout(self.layoutParam)
        self.clearLayout(self.layoutViolin)

        # Add Widgets to Layouts
        self.layoutParam.addWidget(self.checkScrublet, 0, 0, alignment=Qt.AlignmentFlag.AlignCenter)
        self.layoutParam.addWidget(self.labelScrublet, 0, 1, alignment=Qt.AlignmentFlag.AlignLeft)
        self.layoutParam.addWidget(self.inMinCells, 1, 0, alignment=Qt.AlignmentFlag.AlignLeft)
        self.layoutParam.addWidget(self.labelMinCells, 1, 1, alignment=Qt.AlignmentFlag.AlignLeft)
        self.layoutParam.addWidget(self.inMinGenes, 2, 0, alignment=Qt.AlignmentFlag.AlignLeft)
        self.layoutParam.addWidget(self.labelMinGenes, 2, 1, alignment=Qt.AlignmentFlag.AlignLeft)
        if modality != 'Xenium':
            self.layoutParam.addWidget(self.inPerMt, 3, 0, alignment=Qt.AlignmentFlag.AlignLeft)
            self.layoutParam.addWidget(self.labelPerMt, 3, 1, alignment=Qt.AlignmentFlag.AlignLeft)
            self.layoutParam.addWidget(self.inPerRib, 4, 0, alignment=Qt.AlignmentFlag.AlignLeft)
            self.layoutParam.addWidget(self.labelPerRib, 4, 1, alignment=Qt.AlignmentFlag.AlignLeft)
        self.layoutParam.addItem(self.hspacer, 0, 2, alignment=Qt.AlignmentFlag.AlignLeft)

        self.layoutViolin.addWidget(figCells)
        self.layoutViolin.addWidget(figGenes)
        self.layoutViolin.addWidget(figMT)
        self.layoutViolin.addWidget(figRP)
        
        self.clearLayout(self.layoutMain)
        self.layoutMain.addWidget(self.labelRead)
        self.layoutMain.addWidget(self.listDatasets)
        self.layoutMain.addWidget(self.labelFilter)
        self.layoutMain.addLayout(self.layoutViolin)
        self.layoutMain.addLayout(self.layoutParam)

    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)


#################
# Error Message #
#################

class err_NoFiles(QMessageBox):
    def __init__(self):
        super().__init__()

        self.setIcon(QMessageBox.Icon.Critical)
        self.setWindowTitle('Error!')
        self.setText('Going back to folder selection. \nTransciptomic or spatial files not found. \nCheck Folders or Modality.')
        
