#!/bin/python3

from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import subprocess as sp
import itertools
import glob, os, re
## Switch to using white background and black foreground
# pg.setConfigOption('background', 'w')
# pg.setConfigOption('foreground', 'k')
observables = {}

class ParameterComboBox(QtWidgets.QComboBox):
    def __init__(self, names):
        super().__init__()
        self.addItems(names)

class ChooseDataCombobox(QtWidgets.QComboBox):
    def __init__(self, folderBox):
        super().__init__()
        self.folderBox=folderBox
        self.folderBox.currentIndexChanged.connect(self.updateDataList)
        self.updateDataList()
    def updateDataList(self):
        self.clear()
        path = self.folderBox.getCurrentFolder()
        self.files = glob.glob(path+"/*.dat")
        self.names = [os.path.basename(f) for f in self.files]
        self.paths = dict(zip(self.names, self.files));
        self.addItems(self.names)
    def getCurrentPath(self):
        return self.paths.get(self.currentText(), None)

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
        self.names = list(self.mw.data.dtype.names)
        self.X = ParameterComboBox(self.names)
        self.Y = ParameterComboBox(self.names + ["f(X)"])
        self.f = QtWidgets.QLineEdit()
        self.f.editingFinished.connect(self.mw.updatePlots)
        # self.f.setReadOnly(True)
        self.layout().addWidget(self.X, 0, 1)
        self.layout().addWidget(self.Y, 1, 1)
        self.layout().addWidget(self.f, 2, 1)
        self.X.currentIndexChanged.connect(self.mw.updatePlots)
        self.Y.currentIndexChanged.connect(self.mw.updatePlots)
        self.toggled.connect(self.mw.updatePlots)
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
        return self.names[self.X.currentIndex()]
    def yTitle(self):
        if self.Y.currentText() == "f(X)":
            return self.f.text()
        return self.names[self.Y.currentIndex()]


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
        self.selectFile.currentIndexChanged.connect(self.updatePlots)
        self.compileButton = QtWidgets.QPushButton("Compile!")
        self.compileButton.clicked.connect(self.recompile)
        self.plotConfig = QtWidgets.QWidget()
        self.plotConfig.setLayout(QtWidgets.QHBoxLayout())
        self.curveControlArea = QtWidgets.QScrollArea()
        self.curveControlWidget = QtWidgets.QWidget()
        self.curveControlWidget.setLayout(QtWidgets.QHBoxLayout())
        self.curveControlArea.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum,
            QtGui.QSizePolicy.Minimum,
            ))

        self.loadData()
        self.curveControls = [CurveControlWidget(self, active=True),
                              CurveControlWidget(self),
                              CurveControlWidget(self),
                              CurveControlWidget(self),
                              CurveControlWidget(self)]
        for w in self.curveControls:
            self.curveControlWidget.layout().addWidget(w)
        self.chooseCommonObservableX = QtWidgets.QComboBox()
        self.chooseCommonObservableX.addItems(self.observables)
        self.chooseCommonObservableY = QtWidgets.QComboBox()
        self.chooseCommonObservableY.addItems(self.observables)
        self.chooseCommonObservableX.currentIndexChanged.connect(self.setSameObservables)
        self.chooseCommonObservableY.currentIndexChanged.connect(self.setSameObservables)
        assert(len(self.curveControls) >= self.nParticles)
        self.plotConfig.layout().addWidget(self.curveControlArea)
        self.curveControlArea.setWidget(self.curveControlWidget)
        # self.plotConfig.layout().addItem(QtWidgets.QSpacerItem(10, 1))
        self.plotConfig.layout().addWidget(self.fileChooser)
        self.fileChooser.layout().addRow(self.compileButton)
        self.fileChooser.layout().addRow("All X", self.chooseCommonObservableX)
        self.fileChooser.layout().addRow("All Y", self.chooseCommonObservableY)

        self.l.addWidget(self.plotConfig)
        self.pw = pg.PlotWidget(name='Plot1')  ## giving the plots names allows us to link their axes together
        self.l.addWidget(self.pw)
        self.t = QtCore.QTimer()
        # self.t.timeout.connect(a)
        self.t.start(50)

    def setSameObservables(self):
        X = self.chooseCommonObservableX.currentText()
        Y = self.chooseCommonObservableY.currentText()
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

    def updatePlots(self):
        try:
            self.loadData()
        except:
            return
        self.pw.clear()
        flatui = itertools.cycle(["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"])
        for p, c in zip(self.curveControls, flatui):
            if not p.isChecked():
                continue
            # curve.setPen('w')  ## white pen
            pen = pg.mkPen(c)
            self.pw.plot(x=p.xData(), y=p.yData(), symbol='s', symbolPen=pen,
                                 pen=pen, symbolSize=3.5)
            self.pw.setLabel('left', p.yTitle())
            self.pw.setLabel('bottom', p.xTitle())
            # pw.setXRange(0, 2)
            # pw.setYRange(0, 1e-10)
            # curve.setPen(pen, symbol='.')
            # curve.setShadowPen(pen, width=6, cosmetic=True)


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    mw = MainWidget()
    mw.show()
    # data = np.genfromtxt("plot.dat", names=True)
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
