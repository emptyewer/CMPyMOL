import math
import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
from matplotlib.figure import Figure
from PyQt4 import QtGui

class HeatMapFrame(QtGui.QFrame):
    def __init__(self, image, aa, pymol, pmap, parent=None):
        super(HeatMapFrame, self).__init__(parent)
        self.setFrameShape(QtGui.QFrame.NoFrame)
        self.parent = parent
        self.graph_view = GraphView(image, aa, pymol, pmap, self)

    def resizeEvent(self, event):
        self.graph_view.setGeometry(self.rect())

class GraphView(QtGui.QWidget):
    def __init__(self, image, aa, pymol, pmap, parent=None):
        super(GraphView, self).__init__(parent)
        self.aa = aa
        self.pymol = pymol
        self.dpi = 300
        self.pmap = pmap
        self.fig = Figure((4.5, 4.5), dpi=self.dpi)
        self.axes = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.canvas.mpl_connect('button_press_event', self._onpick)
        self.layout = QtGui.QVBoxLayout()
        self.layout.addWidget(self.canvas)
        self.layout.setStretchFactor(self.canvas, 1)
        self.setLayout(self.layout)
        self.dmap = self.cmap_discretize(cm.spectral, 10)
        self.image = image
        self.map = self.axes.imshow(self.image, interpolation='nearest', cmap=self.dmap, aspect='equal', origin='lower')
        self.canvas.show()
        self.set_parameters()

    def cmap_discretize(self, cmap, N):
        """Return a discrete colormap from the continuous colormap cmap.

            cmap: colormap instance, eg. cm.jet.
            N: number of colors.

        Example
            x = resize(arange(100), (5,100))
            djet = cmap_discretize(cm.jet, 5)
            imshow(x, cmap=djet)
        """
        if type(cmap) == str:
            cmap = cm.get_cmap(cmap)
        colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
        colors_rgba = cmap(colors_i)
        indices = np.linspace(0, 1., N + 1)
        cdict = {}
        for ki, key in enumerate(('red', 'green', 'blue')):
            cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki]) for i in xrange(N + 1)]
        # Return colormap object.
        return colors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)

    def _get_clicked_residues(self, event):
        xmin, xmax = self.axes.get_xlim()
        ymin, ymax = self.axes.get_ylim()
        return(int(math.ceil(event.xdata - xmin))-1, int(math.ceil(event.ydata - ymin))-1)

    def _onpick(self, event):
        if self.pmap.use_ca:
            atom = 'CA'
        else:
            atom = 'CB'

        index1, index2 = self._get_clicked_residues(event)
        if self.image[index1][index2] > 0:
            self.pymol.select_aminoacids(self.aa[index1], self.aa[index2], atom, self.pmap.cutoff)

    def set_parameters(self):
        cbar_ax = self.fig.add_axes([0.12, 0.94, 0.75, 0.03])
        cbar_ax.xaxis.labelpad = 3.0
        cbar_ax.tick_params(length=2.0, direction='out', pad=0.0)
        cbar = self.fig.colorbar(self.map, cax=cbar_ax, orientation='horizontal', drawedges=False)
        cbar_ax.xaxis.set_ticks_position('top')
        cbar_ax.xaxis.set_label_position('bottom')
        cbar_ax.xaxis.set_label_text('Contact Counts')
        self.axes.tick_params(axis=u'both', which=u'both', length=0)
        self.axes.xaxis.set_ticks(range(0, len(self.aa), 1))
        self.axes.yaxis.set_ticks(range(0, len(self.aa), 1))
        self.axes.set_xticklabels(self.aa, rotation=90)
        self.axes.set_yticklabels(self.aa)
        self.fig.patch.set_facecolor((0.886, 0.886, 0.886))
        ticks_font = mpl.font_manager.FontProperties(family='times new roman', style='normal', size=12, weight='normal',
                                                     stretch='normal')
        labels = [self.axes.title, self.axes.xaxis.label, self.axes.yaxis.label, cbar.ax.xaxis.label]
        labels += self.axes.get_xticklabels() + self.axes.get_yticklabels() + cbar.ax.get_xticklabels()
        for item in labels:
            item.set_fontproperties(ticks_font)
            item.set_fontsize(4)

    def update_graph(self, image):
        self.axes.clear()
        self.image = image
        self.map = self.axes.imshow(self.image, interpolation='nearest', cmap=self.dmap, aspect='equal', origin='lower')
        self.set_parameters()
        self.canvas.draw()
