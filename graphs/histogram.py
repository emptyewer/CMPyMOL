import math
import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.figure import Figure
from PyQt4 import QtGui

class HistogramFrame(QtGui.QFrame):
    def __init__(self, pymol, pmap, parent=None):
        super(HistogramFrame, self).__init__(parent)
        self.setFrameShape(QtGui.QFrame.NoFrame)
        self.parent = parent
        self.graph_view = GraphView(pymol, pmap, self)

    def resizeEvent(self, event):
        self.graph_view.setGeometry(self.rect())

class GraphView(QtGui.QWidget):
    def __init__(self, pymol, pmap, parent=None):
        super(GraphView, self).__init__(parent)
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
        if self.pmap.use_ca:
            self.xcoor = self.pmap.residue_numbers_ca[self.pmap.parent.current_model]
        else:
            self.xcoor = self.pmap.residue_numbers_cb[self.pmap.parent.current_model]
        self.ycoor = self.pmap.histogram_maps[self.pmap.parent.current_model]
        self.bar = self.axes.bar(self.xcoor, self.ycoor, width=1.0, linewidth=0)
        self.canvas.show()
        self.set_parameters()

    def _get_clicked_residues(self, event):
        xmin, xmax = self.axes.get_xlim()
        return(int(math.ceil(event.xdata - xmin))-1)

    def _onpick(self, event):
        if self.pmap.use_ca:
            atom = 'CA'
            chains = self.pmap.chain_names_ca[self.pmap.parent.current_model]
            residue_numbers = self.pmap.residue_numbers_ca[self.pmap.parent.current_model]
        else:
            atom = 'CB'
            chains = self.pmap.chain_names_cb[self.pmap.parent.current_model]
            residue_numbers = self.pmap.residue_numbers_cb[self.pmap.parent.current_model]

        index = self._get_clicked_residues(event)
        self.pymol.select_density(residue_numbers[index], atom, self.pmap.cutoff, chains[index])

    def set_parameters(self):
        self.axes.tick_params(axis=u'both', which=u'both', length=0)
        self.axes.set_xlim(min(self.xcoor), max(self.xcoor))
        self.axes.set_ylim(0, max(self.ycoor) + 1)
        self.axes.set_xlabel('Residue Number')
        self.axes.set_ylabel('Contact Counts')
        fractions = self.ycoor / max(self.ycoor)
        normalized_colors = colors.Normalize(fractions.min(), fractions.max())
        count = 0
        for rect in self.bar:
            c = cm.jet(normalized_colors(fractions[count]))
            rect.set_facecolor(c)
            count += 1
        self.fig.patch.set_facecolor((0.886, 0.886, 0.886))
        ticks_font = mpl.font_manager.FontProperties(family='times new roman', style='normal', size=12, weight='normal',
                                                     stretch='normal')
        labels = [self.axes.title, self.axes.xaxis.label, self.axes.yaxis.label]
        labels += self.axes.get_xticklabels() + self.axes.get_yticklabels()
        for item in labels:
            item.set_fontproperties(ticks_font)
            item.set_fontsize(4)

    def update_graph(self):
        self.axes.clear()
        if self.pmap.use_ca:
            self.xcoor = self.pmap.residue_numbers_ca[self.pmap.parent.current_model]
        else:
            self.xcoor = self.pmap.residue_numbers_cb[self.pmap.parent.current_model]
        self.ycoor = self.pmap.histogram_maps[self.pmap.parent.current_model]
        self.bar = self.axes.bar(self.xcoor, self.ycoor, width=1.0, linewidth=0)
        self.set_parameters()
        self.canvas.draw()
