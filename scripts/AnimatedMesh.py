#!/bin/python
# -*- coding: utf-8 -*-
"""
Simple examples demonstrating the use of GLMeshItem.

"""
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl

app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: GLMeshItem')
w.setCameraPosition(azimuth=270, elevation=90, pos=(100,100,200))

g = gl.GLGridItem()
g.scale(2,2,1)
w.addItem(g)

import numpy as np

def read2DVertexData(path):
    faces = np.genfromtxt("build/test-results/physics/ElasticModelTests/ball.indices", dtype=int)
    nVertices = int(np.max(faces)) + 1
    data = np.genfromtxt("build/test-results/physics/ElasticModelTests/ball.dat")
    data = np.transpose(data)
    ts, vertexData = data[0], data[1:1+2*nVertices]
    s = vertexData.shape
    s[0]/2
    data.shape[0]/2
    vertices = np.zeros((s[0]+s[0]//2, s[1]))
    vertices[0::3] = vertexData[0::2]
    vertices[1::3] = vertexData[1::2]
    vertices = np.transpose(vertices)
    s = vertices.shape
    newShape = (s[0], s[1]//3, 3)
    vertices = vertices.reshape(newShape)
    return ts, vertices, faces


# Example 4:
# wireframe
# md = gl.MeshData.sphere(rows=4, cols=8)
# m1 = gl.GLMeshItem(meshdata=md,
#                    smooth=False,
#                    drawFaces=False,
#                    drawEdges=True,
#                    edgeColor=(1,1,1,1))
# w.addItem(m1)

ts, vertices, faces = read2DVertexData(
    "build/test-results/physics/ElasticModelTests/ball.dat"
)
md2 = gl.MeshData(vertexes=vertices[0], faces=faces)
m2 = gl.GLMeshItem(meshdata=md2,
                   smooth=False,
                   drawFaces=False,
                   drawEdges=True,
                   edgeColor=(1,1,1,1))
w.addItem(m2)

i = 0
def update():
    global i
    # m4.rotate(1, 1, 1, 0.1)
    if i >= len(vertices):
        i = 0
    # md2.setVertexes(verts=vertices[i])
    md2 = gl.MeshData(vertexes=vertices[i], faces=faces)
    m2.setMeshData(meshdata=md2)
    i+=1
#
timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(30)


## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
