#!/bin/python3

from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import subprocess as sp
import itertools
import glob, os, re, sys, traceback
solarized = {
    "base03":    "#002b36",
    "base02":    "#073642",
    "base01":    "#586e75",
    "base00":    "#657b83",
    "base0":     "#839496",
    "base1":     "#93a1a1",
    "base2":     "#eee8d5",
    "base3":     "#fdf6e3",
    "yellow":    "#b58900",
    "orange":    "#cb4b16",
    "red":       "#dc322f",
    "magenta":   "#d33682",
    "violet":    "#6c71c4",
    "blue":      "#268bd2",
    "cyan":      "#2aa198",
    "green":     "#859900"
}
colors = itertools.cycle(solarized[c] for c in ["base0", "yellow", "orange", "red",
                                       "magenta", "violet", "blue", "cyan",
                                       "green"])
pg.setConfigOptions(background=solarized['base03'])
pg.setConfigOptions(foreground=solarized['base00'])
## Switch to using white background and black foreground
# pg.setConfigOption('background', 'w')
# pg.setConfigOption('foreground', 'k')
observables = {}

# class ParameterComboBox(QtWidgets.QComboBox):
#     def __init__(self, names):
#         super().__init__()
#         self.addItems(names)

class ChooseDataCombobox(QtWidgets.QComboBox):
    def __init__(self, folderBox):
        super().__init__()
        self.folderBox=folderBox
        self.folderBox.currentIndexChanged.connect(self.updateDataList)
        self.folderBox.activated.connect(self.updateDataList)

    def updateDataList(self):
        print("updating data list...")
        self.blockSignals(True)
        self.clear()
        path = self.folderBox.getCurrentFolder()
        self.files = glob.glob(path+"/*.dat")
        self.names = [os.path.basename(f) for f in self.files]
        self.paths = dict(zip(self.names, self.files));
        self.addItems(self.names)
        self.blockSignals(False)

    def getCurrentPath(self):
        return self.paths[self.currentText()]

class ChooseFolderCombobox(QtWidgets.QComboBox):
    def __init__(self):
        super().__init__()
        self.folders = glob.glob("build/test-results/physics/SoftBody/*")
        self.names = [os.path.basename(f) for f in self.folders]
        self.paths = dict(zip(self.names, self.folders));
        self.addItems(self.names)
    def getCurrentFolder(self):
        return self.paths[self.currentText()]

class CurveControlWidget(QtWidgets.QGroupBox):
    def __init__(self, mw, active=False, title=""):
        super().__init__()
        self.mw = mw
        self.setCheckable(True)
        self.setChecked(active)
        self.setLayout(QtWidgets.QGridLayout())
        self.setTitle(title)
        self.layout().addWidget(QtWidgets.QLabel("X"), 0, 0)
        self.layout().addWidget(QtWidgets.QLabel("Y"), 1, 0)
        self.layout().addWidget(QtWidgets.QLabel("f(X)"), 2, 0)
        self.X = QtWidgets.QComboBox()
        self.Y = QtWidgets.QComboBox()
        self.f = QtWidgets.QLineEdit()
        # self.f.setReadOnly(True)
        self.layout().addWidget(self.X, 0, 1)
        self.layout().addWidget(self.Y, 1, 1)
        self.layout().addWidget(self.f, 2, 1)

    def names(self):
        return list(self.mw.data.dtype.names)


    def updateNames(self):
        self.X.blockSignals(True)
        self.Y.blockSignals(True)
        xItem = self.X.currentText()
        self.X.clear()
        self.X.addItems(self.names())

        yItem = self.Y.currentText()
        if yItem != "f(X)" and yItem not in self.names():
            yItem = self.names()[0]
        self.Y.clear()
        self.Y.addItems(self.names() + ["f(X)"])
        self.X.blockSignals(False)
        self.Y.blockSignals(False)
        if xItem != "" and xItem in self.names():
            self.X.setCurrentText(xItem)
        if yItem != "" and yItem in self.names():
            self.Y.setCurrentText(yItem)

    def xData(self):
        return self.mw.data[self.X.currentText()]

    def yData(self):
        if self.Y.currentText() == "f(X)":
            try:
                print(self.f.text())
                X = self.xData();X
                return eval(self.f.text()) + 0 * self.xData()
            except:
                return 0 * self.xData()
        else:
            return self.mw.data[self.Y.currentText()]
    def xTitle(self):
        return self.names()[self.X.currentIndex()]
    def yTitle(self):
        if self.Y.currentText() == "f(X)":
            return self.f.text()
        return self.names()[self.Y.currentIndex()]


class MainWidget(QtGui.QMainWindow):
    def __init__(self):
        self.app = QtGui.QApplication([])
        super().__init__()
        self.setWindowTitle('Plot Data')
        self.resize(800,800)
        self.cw = QtGui.QWidget()
        self.setCentralWidget(self.cw)
        self.l = QtGui.QVBoxLayout()
        self.cw.setLayout(self.l)
        self.fileChooser = QtWidgets.QWidget()
        self.fileChooser.setLayout(QtWidgets.QFormLayout())
        self.selectFolder = ChooseFolderCombobox()
        self.selectFile = ChooseDataCombobox(self.selectFolder)
        self.fileChooser.layout().addRow("Folder", self.selectFolder)
        self.fileChooser.layout().addRow("File", self.selectFile)
        self.compileButton = QtWidgets.QPushButton("Compile!")
        self.plotConfig = QtWidgets.QWidget()
        self.plotConfig.setLayout(QtWidgets.QHBoxLayout())
        self.curveControlArea = QtWidgets.QScrollArea()
        self.curveControlWidget = QtWidgets.QWidget()
        self.curveControlWidget.setLayout(QtWidgets.QHBoxLayout())
        self.curveControlArea.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum,
            QtGui.QSizePolicy.Minimum,
            ))

        self.curveControls = [CurveControlWidget(self, active=True),
                              CurveControlWidget(self),
                              CurveControlWidget(self),
                              CurveControlWidget(self),
                              CurveControlWidget(self)]
        for w in self.curveControls:
            self.curveControlWidget.layout().addWidget(w)
        self.plotConfig.layout().addWidget(self.curveControlArea)
        self.curveControlArea.setWidget(self.curveControlWidget)
        # self.plotConfig.layout().addItem(QtWidgets.QSpacerItem(10, 1))

        self.chooseCommonObservableX = QtWidgets.QComboBox()
        self.chooseCommonObservableY = QtWidgets.QComboBox()

        self.plotConfig.layout().addWidget(self.fileChooser)
        self.fileChooser.layout().addRow(self.compileButton)
        self.fileChooser.layout().addRow("All X", self.chooseCommonObservableX)
        self.fileChooser.layout().addRow("All Y", self.chooseCommonObservableY)

        self.l.addWidget(self.plotConfig)
        self.pw = pg.PlotWidget(name='Plot1')  ## giving the plots names allows us to link their axes together
        self.l.addWidget(self.pw)
        self.commonsLoaded  = False
        self.connectStuff()
        self.selectFolder.currentIndexChanged.emit(0)
        self.selectFile.currentIndexChanged.emit(0)

    def loadCommonObservablesList(self):
        self.chooseCommonObservableX.clear()
        self.chooseCommonObservableY.clear()
        self.chooseCommonObservableX.addItems(self.observables)
        self.chooseCommonObservableY.addItems(self.observables)

    def connectStuff(self):
        self.selectFile.currentIndexChanged.connect(self.loadData)
        self.compileButton.clicked.connect(self.recompile)
        self.chooseCommonObservableX.currentIndexChanged.connect(self.setSameObservables)
        self.chooseCommonObservableY.currentIndexChanged.connect(self.setSameObservables)

        for w in self.curveControls:
            w.f.editingFinished.connect(self.updatePlots)
            w.X.currentIndexChanged.connect(self.updatePlots)
            w.Y.currentIndexChanged.connect(self.updatePlots)
            w.toggled.connect(self.updatePlots)

    def setSameObservables(self):
        X = self.chooseCommonObservableX.currentText()
        Y = self.chooseCommonObservableY.currentText()
        print("Number of particles: ", self.nParticles)
        for i in range(self.nParticles):
            self.curveControls[i].X.setCurrentText(X + str(i))
            self.curveControls[i].Y.setCurrentText(Y + str(i))

    def recompile(self):
        sp.run("gradle test", shell=True)
        self.loadData()
        self.updatePlots()

    def loadData(self):
        self.data = np.genfromtxt(self.selectFile.getCurrentPath(), names=True)
        self.observables = set(re.sub(r"\d+", "", n) for n in self.data.dtype.names if re.match("\w+\d+", n))
        self.nParticles = max([int(n[-1:]) for n in self.data.dtype.names
                               if re.match("\w+\d+", n)]) + 1
        for w in self.curveControls:
            w.updateNames()
        if not self.commonsLoaded:
            self.loadCommonObservablesList()
            self.commonsLoaded = True

    def updatePlots(self):
        self.pw.clear()
        if not hasattr(self, "data"):
            return
        for i, p in enumerate(self.curveControls):
            if not p.isChecked():
                continue
            # curve.setPen('w')  ## white pen
            pen = pg.mkPen(next(colors))
            self.pw.plot(x=p.xData(), y=p.yData(), symbol='s', symbolPen=pen,
                                 pen=pen, symbolSize=3.5)
            self.pw.setLabel('left', p.yTitle())
            self.pw.setLabel('bottom', p.xTitle())


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    mw = MainWidget()
    mw.show()
    # data = np.genfromtxt("plot.dat", names=True)
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
