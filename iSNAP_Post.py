###################
# Import Packages #
###################
print('- Importing sys')
from sys import float_info

print('- Importing os...')
from os import (
    makedirs,
)

from os.path import (
    exists as pathexists,
    join
)

print('- Importing pandas...')
from pandas import DataFrame, concat

print('- Importing numpy...')
from numpy import log10, where, isinf

print('- Importing matplotlib...')
from matplotlib.pyplot import (
    figure as pltfigure,
    close as pltclose,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.lines import Line2D
from adjustText import adjust_text

print('- Importing PyQt...')
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, 
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QLabel,
    QPushButton,
    QCheckBox,
    QScrollArea,
    QGroupBox,
    QComboBox,
    QLineEdit,
    QSizePolicy,
    QSpacerItem
)
from PyQt5.QtGui import QFont

#############
# Functions #
#############

def dgetoCSV(adatas, mainList, groupAList, groupBList, DGEMethod, output_folder):
    # Extract results
    print('Saving DGE to CSV...')
    result = adatas.uns[f'DGE {DGEMethod}']
    groups = result['names'].dtype.names

    # Create an empty list to collect rows
    records = []

    for group in groups:
        for i in range(len(result['names'][group])):
            records.append({
                'gene': result['names'][group][i],
                'logfoldchange': result['logfoldchanges'][group][i],
                'score': result['scores'][group][i],
                'pval': result['pvals'][group][i],
                'pval_adj': result['pvals_adj'][group][i],
            })
    
    # Convert to dataframe
    deg_table = DataFrame(records)

    addColumn = DataFrame({
        'Included Groups': mainList
    })
    deg_table = concat([deg_table, addColumn], axis=1)

    addColumn = DataFrame({
        'Group A (Focus Group)': groupAList
    })
    deg_table = concat([deg_table, addColumn], axis=1)
    
    addColumn = DataFrame({
        'Group B (Reference Group)': groupBList
    })
    deg_table = concat([deg_table, addColumn], axis=1)
    
    # Save to CSV
    outpath = join(output_folder, 'DGE Comparison Table')

    if not pathexists(outpath):
        makedirs(outpath)

    i=1
    if pathexists(join(outpath, 'DGETable.csv')):
        while pathexists(join(outpath, f'DGETable({i}).csv')):
            i += 1
        deg_table.to_csv(join(outpath, f'DGETable({i}).csv'), index=False)
        
    else:
        deg_table.to_csv(join(outpath, 'DGETable.csv'), index=False)
        
    print(f'.csv File created under: {outpath}')
    

def plotVolcano(adatas, DGEMethod, nTopGenes, labelGenes, pvalThresh, logFCThresh):
        genes = adatas.uns[f'DGE {DGEMethod}']['names']['Group A']
        pval = adatas.uns[f'DGE {DGEMethod}']['pvals_adj'].astype(float)
        logfc = adatas.uns[f'DGE {DGEMethod}']['logfoldchanges'].astype(float)

        logpval = -log10(pval)
        logpval = where(isinf(logpval), -log10(float_info.min), logpval)

        colors = []
        significance = []
        for i in range(len(logpval)):
            if pval[i] < pvalThresh:
                if logfc[i] > logFCThresh:
                    colors.append('red')
                    significance.append(1)
                elif logfc[i] < -logFCThresh:
                    colors.append('blue')
                    significance.append(1)
                else: 
                    colors.append('gray')
                    significance.append(0)
            else:
                colors.append('gray')
                significance.append(0)

        labelGenesIndex = sorted(range(len(logfc)), key=lambda i:abs(logfc[i])*significance[i], reverse=True)[:nTopGenes]

        for gene in labelGenes:
            if gene in genes:
                labelGenesIndex.append(where(genes==gene)[0][0])
            else:
                print(f'{gene} not found in gene pool.')

        topgenes = [genes[i] for i in labelGenesIndex]
        toplogpval = [logpval[i] for i in labelGenesIndex]
        toplogfc = [logfc[i] for i in labelGenesIndex]

        figVol = pltfigure(constrained_layout=True, dpi=200)
        axVol = figVol.add_subplot(111)
        axVol.set_ylabel('-log$_{10}$(Pval$_{\\rm adj}$)')
        axVol.set_xlabel('log$_2$FC')
        axVol.axvline(-logFCThresh, linestyle = '--', linewidth=1, c='lightgray')
        axVol.axvline(logFCThresh, linestyle = '--', linewidth=1, c='lightgray')
        axVol.axhline(-log10(pvalThresh), linestyle = '--', linewidth=1, c='lightgray')
        axVol.scatter(logfc, logpval, c=colors, s=1)
        texts = [axVol.text(toplogfc[i], toplogpval[i], topgenes[i], fontsize=8) for i in range(len(topgenes))]
        adjust_text(texts, arrowprops=dict(arrowstyle='-'), min_arrow_len=0)
        
        downCount = colors.count('blue')
        nsCount = colors.count('gray')
        upCount = colors.count('red')

        legend_elements = [
            Line2D([0], [0], color='w', marker='o', markerfacecolor='red', label=f'Up ({upCount})'),
            Line2D([0], [0], color='w', marker='o', markerfacecolor='blue', label=f'Down ({downCount})'), 
            Line2D([0], [0], color='w', marker='o', markerfacecolor='gray', label=f'NS ({nsCount})')
            ]

        axVol.legend(handles=legend_elements)

        pltclose('all')

        return figVol


###########
# Objects #
###########

class TypeSampleTable(QWidget):
    toDGESignal = pyqtSignal(bool)
    def __init__(self):
        super().__init__()

        self.layoutMain =  QHBoxLayout()
        self.layoutSub = QVBoxLayout()

        self.spacer = QSpacerItem(1, 5, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.labelDGE = QLabel('DGE Analysis:')
        self.btnDGESamples = QPushButton('Compare Between Samples')
        self.btnDGESamples.clicked.connect(self.toDGESamples)
        self.btnDGETypes = QPushButton('Compare Between Cell Types')
        self.btnDGETypes.clicked.connect(self.toDGETypes)

        self.layoutSub.addItem(self.spacer)
        self.layoutSub.addWidget(self.labelDGE)
        self.layoutSub.addWidget(self.btnDGESamples)
        self.layoutSub.addWidget(self.btnDGETypes)
        self.layoutSub.addItem(self.spacer)

        self.setLayout(self.layoutMain)
    
    def toDGESamples(self):
        self.toDGESignal.emit(True)
    def toDGETypes(self):
        self.toDGESignal.emit(False)
    
    def setPage(self, sampleTypeMtx, samples, types):
        self.tableTypeSample = QTableWidget()

        self.tableTypeSample.setColumnCount(len(samples))
        self.tableTypeSample.setRowCount(len(types))

        self.tableTypeSample.setHorizontalHeaderLabels(samples)
        self.tableTypeSample.setVerticalHeaderLabels(types)

        for i in range(len(samples)):
            for j in range(len(types)):
                item = QTableWidgetItem(str(int(sampleTypeMtx[j, i])))
                if j==(len(types)-1) or i == (len(samples)-1):
                    font=QFont()
                    font.setBold(True)
                    item.setFont(font)

                self.tableTypeSample.setItem(j, i, item)
        
        self.tableTypeSample.resizeColumnsToContents()
        for i in range(self.tableTypeSample.columnCount()):
            self.tableTypeSample.setColumnWidth(i, self.tableTypeSample.columnWidth(i) + 16)
        
        self.clearLayout(self.layoutMain)
        self.layoutMain.addWidget(self.tableTypeSample)
        self.layoutMain.addLayout(self.layoutSub)
    
    def clearLayout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)

class DGEWindow(QWidget):
    toVolcanoSignal = pyqtSignal(list, list, list, bool, str, int, list, float, float)
    def __init__(self, bySamples, mainItems, groupItems):
        super().__init__()

        self.compSamples = bySamples

        self.setWindowTitle("iSNAP - DGE Analysis")
        self.setGeometry(100, 100, 1280, 720)
        
        self.layoutMain = QVBoxLayout()
        self.layoutGroups = QHBoxLayout()
        self.layoutGroupA = QVBoxLayout()
        self.layoutGroupB = QVBoxLayout()
        self.layoutMainItems = QVBoxLayout()
        self.layoutGroupAItems = QVBoxLayout()
        self.layoutGroupBItems = QVBoxLayout()
        self.layoutMethod = QHBoxLayout()
        self.layoutTopGenes = QHBoxLayout()
        self.layoutLabelGenes = QHBoxLayout()
        self.layoutPVal = QHBoxLayout()
        self.layoutlogFC = QHBoxLayout()

        self.vspacer = QSpacerItem(1, 1, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.hspacer =  QSpacerItem(1, 1, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.labelGroupA = QLabel('Group A (Experimental Group)')
        self.labelGroupB = QLabel('Group B (Reference Group)')

        self.labelDGEMethod = QLabel('Select Statistical Test:')
        self.comboDGEMethod = QComboBox()
        self.comboDGEMethod.addItem('wilcoxon')
        self.comboDGEMethod.addItem('t-test')
        self.comboDGEMethod.addItem('t-test_overestim_var')

        self.labelTopGenes = QLabel('Number of Genes to Label')
        self.lineTopGenes = QLineEdit('10')

        self.labelLabelGenes = QLabel('Label Genes (ex: gene1, gene2, ...)')
        self.lineLabelGenes = QLineEdit()

        self.labelPVal = QLabel('PVal Adj Threshold')
        self.linePVal = QLineEdit('10e-6')

        self.labelLogFC = QLabel('LogFC Threshold')
        self.lineLogFC = QLineEdit('0.5')

        self.layoutParamLabel = QVBoxLayout()
        self.layoutParamLabel.addWidget(self.labelDGEMethod, alignment=Qt.AlignLeft)
        self.layoutParamLabel.addWidget(self.labelTopGenes, alignment=Qt.AlignLeft)
        self.layoutParamLabel.addWidget(self.labelLabelGenes, alignment=Qt.AlignLeft)
        self.layoutParamLabel.addWidget(self.labelPVal, alignment=Qt.AlignLeft)
        self.layoutParamLabel.addWidget(self.labelLogFC, alignment=Qt.AlignLeft)

        self.layoutParamWidget = QVBoxLayout()
        self.layoutParamWidget.addWidget(self.comboDGEMethod, alignment=Qt.AlignLeft)
        self.layoutParamWidget.addWidget(self.lineTopGenes, alignment=Qt.AlignLeft)
        self.layoutParamWidget.addWidget(self.lineLabelGenes, alignment=Qt.AlignLeft)
        self.layoutParamWidget.addWidget(self.linePVal, alignment=Qt.AlignLeft)
        self.layoutParamWidget.addWidget(self.lineLogFC, alignment=Qt.AlignLeft)

        self.layoutParamAll = QHBoxLayout()
        self.layoutParamAll.addItem(self.hspacer)
        self.layoutParamAll.addLayout(self.layoutParamLabel)
        self.layoutParamAll.addLayout(self.layoutParamWidget)
        self.layoutParamAll.addItem(self.hspacer)


        
        self.btnDGE = QPushButton('Confirm')
        self.btnDGE.clicked.connect(self.toVolcano)

        if bySamples:
            self.labelMain = QLabel('Choose Cell Types to Use')

        else:
            self.labelMain = QLabel('Choose Samples to Use')
            

        self.checksMain = []
        
        for i in range(len(mainItems)):
            self.checksMain.append(QCheckBox(mainItems[i]))
            self.checksMain[i].setChecked(True)
            self.layoutMainItems.addWidget(self.checksMain[i])
        
        self.checksGroupA = [None] * len(groupItems)
        self.checksGroupB = [None] * len(groupItems)

        for i in range(len(groupItems)):
            self.checksGroupA[i] = EntangledCheckBox(groupItems[i], i, True)
            self.checksGroupA[i].entangledSignal.connect(self.disableOtherCheck)
            self.checksGroupA[i].toggled.connect(self.checksGroupA[i].entangleFunc)
            self.layoutGroupAItems.addWidget(self.checksGroupA[i])

            self.checksGroupB[i] = EntangledCheckBox(groupItems[i], i, False)
            self.checksGroupB[i].entangledSignal.connect(self.disableOtherCheck)
            self.checksGroupB[i].toggled.connect(self.checksGroupB[i].entangleFunc)
            self.layoutGroupBItems.addWidget(self.checksGroupB[i])

        self.groupMain = QGroupBox()
        self.groupGroupA = QGroupBox()
        self.groupGroupB = QGroupBox()
        
        self.groupMain.setLayout(self.layoutMainItems)
        self.groupGroupA.setLayout(self.layoutGroupAItems)
        self.groupGroupB.setLayout(self.layoutGroupBItems)

        self.scrollMain = QScrollArea()
        self.scrollGroupA = QScrollArea()
        self.scrollGroupB = QScrollArea()

        self.scrollMain.setWidget(self.groupMain)
        self.scrollMain.setWidgetResizable(True)
        self.scrollMain.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.scrollGroupA.setWidget(self.groupGroupA)
        self.scrollGroupA.setWidgetResizable(True)
        self.scrollGroupA.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.scrollGroupB.setWidget(self.groupGroupB)
        self.scrollGroupB.setWidgetResizable(True)
        self.scrollGroupB.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        # Set Layouts

        self.layoutGroupA.addWidget(self.labelGroupA, alignment=Qt.AlignCenter)
        self.layoutGroupA.addWidget(self.scrollGroupA, alignment=Qt.AlignCenter)

        self.layoutGroupB.addWidget(self.labelGroupB, alignment=Qt.AlignCenter)
        self.layoutGroupB.addWidget(self.scrollGroupB, alignment=Qt.AlignCenter)

        self.layoutGroups.addLayout(self.layoutGroupA)
        self.layoutGroups.addLayout(self.layoutGroupB)

        self.layoutMain.addWidget(self.labelMain, alignment=Qt.AlignCenter)
        self.layoutMain.addWidget(self.scrollMain, alignment=Qt.AlignCenter)
        self.layoutMain.addLayout(self.layoutGroups)
        self.layoutMain.addItem(self.vspacer)
        self.layoutMain.addLayout(self.layoutParamAll)
        self.layoutMain.addWidget(self.btnDGE, alignment=Qt.AlignCenter)

        self.setLayout(self.layoutMain)

        self.show()

    def toVolcano(self):
        try:
            self.btnDGE.setEnabled(False)
            mainList = []
            for i in range(len(self.checksMain)):
                if self.checksMain[i].isChecked():
                    mainList.append(self.checksMain[i].text())
            
            groupAList = []
            groupBList = []
            for i in range(len(self.checksGroupA)):
                if self.checksGroupA[i].isChecked():
                    groupAList.append(self.checksGroupA[i].text())
                if self.checksGroupB[i].isChecked():
                    groupBList.append(self.checksGroupB[i].text())

            DGEMethod = self.comboDGEMethod.currentText()
            
            focusGenes = []
            if self.lineLabelGenes.text().strip() != '':
                focusGenes = [x.strip() for x in self.lineLabelGenes.text().split(',')]

            if self.lineTopGenes.text().strip() == '':
                nTopGenes = 0
            else:
                nTopGenes = int(self.lineTopGenes.text())

            pvalThresh = float(self.linePVal.text())
            logFCThresh = float(self.lineLogFC.text())

            self.toVolcanoSignal.emit(
                mainList, 
                groupAList, 
                groupBList, 
                self.compSamples, 
                DGEMethod, 
                nTopGenes, 
                focusGenes, 
                pvalThresh, 
                logFCThresh
                )

        except Exception as e:
            print(e)

    def disableOtherCheck(self, index, isChecked, isA):
        if isA:
            if isChecked:
                self.checksGroupB[index].setEnabled(False)
            else:
                self.checksGroupB[index].setEnabled(True)
        else:
            if isChecked:
                self.checksGroupA[index].setEnabled(False)
            else:
                self.checksGroupA[index].setEnabled(True)
    

class VolcanoPlot(QWidget):
    def __init__(self, figVol):
        super().__init__()

        self.setWindowTitle("iSNAP - DGE Analysis Volcano Plot")
        self.setGeometry(100, 100, 1280, 720)
        
        self.layoutMain = QVBoxLayout()

        self.volcanoPlot = FigureCanvasQTAgg(figVol)

        self.toolbar = NavigationToolbar2QT(self.volcanoPlot, self)

        self.layoutMain.addWidget(self.volcanoPlot)
        self.layoutMain.addWidget(self.toolbar)

        self.setLayout(self.layoutMain)

        self.show()

        print('Done!')
        

class EntangledCheckBox(QCheckBox):
    entangledSignal = pyqtSignal(int, bool, bool)
    
    def __init__(self, text, index, isA):
        super().__init__()

        self.setText(text)
        self.index = index
        self.isA = isA
        
    def entangleFunc(self):
        self.entangledSignal.emit(self.index, self.isChecked(), self.isA)

            