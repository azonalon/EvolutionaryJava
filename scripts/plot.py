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

    def updateDataList(self):
        print("updating data list...")
        self.blockSignals(True)
        self.clear()
        path = self.folderBox.getCurrentFolder()
        self.files = glob.glob(path+"/*.dat")
        self.names = [os.path.basename(f) for f in self.files]
        self.paths = dict(zip(self.names, self.files));
        self.paths[""] = "";
        self.addItems(self.names)
        self.setCurrentText("")
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
        super().__init__(parent=mw)
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
        if xItem != "" and xItem in self.names():
            self.X.setCurrentText(xItem)
        if yItem != "" and yItem in self.names():
            self.Y.setCurrentText(yItem)
        self.X.blockSignals(False)
        self.Y.blockSignals(False)

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

class AnimationControl(QtWidgets.QGroupBox):
    def __init__(self):
        super().__init__()
        self.setCheckable(True)
        self.setChecked(False)
        self.setLayout(QtWidgets.QHBoxLayout())
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.play = QtWidgets.QPushButton("Run")
        self.stepSelect = QtWidgets.QSpinBox()
        self.play.setCheckable(True)
        self.play.toggled.connect(self.togglePlay)
        self.timer = QtCore.QTimer()
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.step)
        self.stepRight = QtWidgets.QPushButton(">")
        self.stepLeft = QtWidgets.QPushButton("<")
        self.stepRight.setAutoRepeat(True)
        self.stepRight.clicked.connect(self.step)
        self.stepLeft.clicked.connect(self.stepBack)

        self.layout().addWidget(self.stepSelect)
        self.layout().addWidget(self.play)
        self.layout().addWidget(self.stepLeft)
        self.layout().addWidget(self.stepRight)
        self.layout().addWidget(self.slider)
    def setTimePoints(self, nPoints):
        self.slider.setMinimum(0)
        self.slider.setMaximum(nPoints-1)
    def active(self):
        return self.isChecked()
    def currentIndex(self):
        return self.slider.value()

    def step(self):
        pos = self.slider.value()
        ds = self.stepSelect.value()
        if pos+ds < self.slider.maximum():
            self.slider.setValue(pos+ds)
        else:
            self.slider.setValue(0)

    def stepBack(self):
        pos = self.slider.value()
        if pos-1 < self.slider.minimum():
            self.slider.setValue(pos-1)
        else:
            self.slider.setValue(self.slider.maximum())
    def togglePlay(self, checked):
        if not checked:
            self.timer.stop()
        else:
            self.timer.start()


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
        self.animate = AnimationControl()
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
            QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum,
            ))

        self.curveControls = []
        self.curveControls = [CurveControlWidget(self, active=True)]
        for w in self.curveControls:
            self.curveControlWidget.layout().addWidget(w)
        self.plotConfig.layout().addWidget(self.curveControlArea)
        # self.plotConfig.layout().addItem(QtWidgets.QSpacerItem(10, 1))

        self.chooseCommonObservableX = QtWidgets.QComboBox()
        self.chooseCommonObservableY = QtWidgets.QComboBox()
        self.toggleAll = QtWidgets.QCheckBox()
        self.toggleAll.stateChanged.connect(self.toggleActive)

        self.plotConfig.layout().addWidget(self.fileChooser)
        self.fileChooser.layout().addRow(self.compileButton)
        self.fileChooser.layout().addRow("Set X", self.chooseCommonObservableX)
        self.fileChooser.layout().addRow("Set Y", self.chooseCommonObservableY)
        self.fileChooser.layout().addRow("Toggle", self.toggleAll)

        self.l.addWidget(self.plotConfig)
        self.pw = pg.PlotWidget(name='Plot1')  ## giving the plots names allows us to link their axes together
        self.l.addWidget(self.animate)
        self.l.addWidget(self.pw)
        self.commonsLoaded  = False
        self.connectStuff()
        self.selectFolder.currentIndexChanged.emit(0)
        self.selectFile.currentIndexChanged.emit(0)
        self.setSameObservables()

    def toggleActive(self, state):
        for w in self.curveControls:
            w.setChecked(state)




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
        self.animate.slider.valueChanged.connect(self.updatePlots)

        # for w in self.curveControls:
        #     w.f.editingFinished.connect(self.updatePlots)
        #     w.X.currentIndexChanged.connect(self.updatePlots)
        #     w.Y.currentIndexChanged.connect(self.updatePlots)
        #     w.toggled.connect(self.updatePlots)

    def setSameObservables(self):
        X = self.chooseCommonObservableX.currentText()
        Y = self.chooseCommonObservableY.currentText()
        print("Number of particles: ", self.nParticles)
        for i in range(self.nParticles):
            if i >= len(self.curveControls):
                w = CurveControlWidget(self, active=False)
                self.curveControls.append(w)
                w.X.setCurrentText(X + str(i))
                w.Y.setCurrentText(Y + str(i))
                self.curveControlWidget.layout().addWidget(w)
                w.f.editingFinished.connect(self.updatePlots)
                w.X.currentIndexChanged.connect(self.updatePlots)
                w.Y.currentIndexChanged.connect(self.updatePlots)
                w.toggled.connect(self.updatePlots)
                w.updateNames()
            else:
                self.curveControls[i].X.setCurrentText(X + str(i))
                self.curveControls[i].Y.setCurrentText(Y + str(i))
                self.curveControls[i].setChecked(True)
        self.curveControlArea.setWidget(self.curveControlWidget)

    def recompile(self):
        sp.run("gradle test", shell=True)
        self.selectFolder.currentIndexChanged.emit(0)
        self.loadData()
        self.updatePlots()

    def loadData(self):
        print("loading data")
        f = self.selectFile.getCurrentPath()
        if f != '':
            self.data = np.genfromtxt(f, names=True)
            self.observables = set(re.sub(r"\d+", "", n) for n in self.data.dtype.names if re.match("\w+\d+", n))
            self.nParticles = max([int(re.sub("\D+", "", n)) for n in self.data.dtype.names
                                   if re.match("\D+\d+", n)]) + 1
            self.animate.setTimePoints(len(self.data["Time"]))
            for w in self.curveControls:
                w.updateNames()
            if not self.commonsLoaded:
                self.loadCommonObservablesList()
                self.commonsLoaded = True
            self.updatePlots()

    def updatePlots(self):
        self.pw.clear()
        t0 = 0
        t1 = -1
        ps=3.5
        colors = itertools.cycle(solarized[c] for c in ["base0", "yellow", "orange", "red",
                                               "magenta", "violet", "blue", "cyan",
                                               "green"])
        if self.animate.active():
            t0 = self.animate.currentIndex()
            t1 = t0 + 1
            ps=10
        if not hasattr(self, "data"):
            return
        for i, p in enumerate(self.curveControls):
            if not p.isChecked():
                continue
            # curve.setPen('w')  ## white pen
            pen = pg.mkPen(next(colors))
            self.pw.plot(x=p.xData()[t0:t1], y=p.yData()[t0:t1], symbol='s', symbolPen=pen,
                                 pen=pen, symbolSize=ps)
            self.pw.setLabel('left', p.yTitle())
            self.pw.setLabel('bottom', p.xTitle())


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    mw = MainWidget()
    mw.show()
    # data = np.genfromtxt("plot.dat", names=True)
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
