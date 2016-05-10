import math
from itertools import groupby
from operator import itemgetter
import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib.patches import Rectangle
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from PyQt4 import QtGui

class ContactMapWidget(QtGui.QWidget):
    def __init__(self, pymol, protein_map, parent=None):
        super(ContactMapWidget, self).__init__(parent)
        self.dpi = 300
        self.pymol = pymol
        self.protein_map = protein_map
        self.parent = parent
        self.gridpsec = gridspec.GridSpec(1, 1)
        self.figure = Figure((8, 8), dpi=self.dpi)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axes = self.figure.add_subplot(self.gridpsec[0])
        self.layoutVertical = QtGui.QVBoxLayout(self)
        self.layoutVertical.setStretchFactor(self.canvas, 1)
        self.layoutVertical.addWidget(self.canvas)
        self.selection_count = 0
        self._ispressed = False
        self.rect = None
        self.glob_ss_rects = []
        self.glob_charged_dots = []
        self.glob_hydrophobic_dots = []
        self.glob_bfactor_dots = []

    def add_rect(self, remove=True):
        if self.rect and remove:
            self.rect.remove()
        self.rect = Rectangle((0, 0), 0, 0, alpha=0.5, facecolor='magenta', edgecolor='none')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.axes.add_patch(self.rect)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def plot(self, image):
        self.map = self.axes.imshow(image, interpolation='nearest', cmap='gray', aspect='equal', origin='lower')
        self.set_parameters()
        self.canvas.draw()
        self.add_rect()

    def on_press(self, event):
        self._ispressed = True
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self, event):
        self._ispressed = False
        self.selection_count += 1
        x0 = min(int(math.floor(self.x0)), int(math.floor(self.x1)))
        y0 = min(int(math.floor(self.y0)), int(math.floor(self.y1)))
        x1 = max(int(math.floor(self.x0)), int(math.floor(self.x1)))
        y1 = max(int(math.floor(self.y0)), int(math.floor(self.y1)))
        if self.parent.show_contact_map:
            map = self.parent.contact_array
        else:
            map = self.parent.variance_array
        xrange = []
        yrange = []
        chainx = []
        chainy = []
        if self.protein_map.use_ca:
            residue_numbers = self.protein_map.residue_numbers_ca
            chain_names = self.protein_map.chain_names_ca
            atom = 'CA'
        else:
            residue_numbers = self.protein_map.residue_numbers_cb
            chain_names = self.protein_map.chain_names_cb
            atom = 'CB'

        for i in range(x0, x1+1):
            for j in range(y0, y1+1):
                if map[j][i] > 0:
                    xrange.append(str(residue_numbers[self.parent.current_model][j]))
                    yrange.append(str(residue_numbers[self.parent.current_model][i]))
                    chainx.append(chain_names[self.parent.current_model][j])
                    chainy.append(chain_names[self.parent.current_model][i])
        self.pymol.select_range(xrange, yrange, set(chainx), set(chainy), atom, self.selection_count)

    def on_motion(self, event):
        if self._ispressed:
            self.x1 = event.xdata
            self.y1 = event.ydata
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.canvas.draw()

    def set_parameters(self):
        self.axes.tick_params(axis=u'both', which=u'both', length=0)
        self.figure.patch.set_facecolor((0.886, 0.886, 0.886))
        ticks_font = mpl.font_manager.FontProperties(family='times new roman', style='normal', size=12, weight='normal',
                                                     stretch='normal')
        labels = [self.axes.title, self.axes.xaxis.label, self.axes.yaxis.label]
        labels += self.axes.get_xticklabels() + self.axes.get_yticklabels()
        for item in labels:
            item.set_fontproperties(ticks_font)
            item.set_fontsize(4)
        self.figure.subplots_adjust(bottom=0.0, left=0.0, top=1.0, right=1.0)

    def _get_ranges(self, data):
        '''
        Returns Starting and Ending Indexes for secondary structure sequence that is passed on (data)
        '''
        ranges = []
        for k, g in groupby(enumerate(data), lambda i_x: i_x[0] - i_x[1]):
            group = map(itemgetter(1), g)
            ranges.append((group[0], group[-1]))
        return ranges

    def _get_alpha_beta(self, ss):
        '''
        Read residue indexes for alpha helix and beta sheets
        '''
        alpha_asg = []
        beta_asg = []
        res_count = 0
        for asg in ss:
            if asg is 'H':
                alpha_asg.append(res_count)
            elif asg is 'E':
                beta_asg.append(res_count)
            res_count += 1
        return (alpha_asg, beta_asg)

    def _get_rectangle_x(self, r, c):
        xs = r[0]
        ys = 0
        xe = r[1]
        ye = len(self.protein_map.residue_numbers_ca[self.parent.current_model])
        rec = Rectangle((xs, ys), width=(xe - xs), height=(ye - ys), alpha=0.25, facecolor=c, edgecolor=c)
        return rec

    def _get_rectangle_y(self, r, c):
        xs = 0
        ys = r[0]
        xe = len(self.protein_map.residue_numbers_ca[self.parent.current_model])
        ye = r[1]
        rec = Rectangle((xs, ys), width=(xe - xs), height=(ye - ys), alpha=0.25, facecolor=c, edgecolor=c)
        return rec

    def _get_dot(self, x, y, c):
        rec = Rectangle((x, y), width=1.0, height=1.0, alpha=1.0, facecolor=c, edgecolor=c)
        return rec

    def overlay_secondary_structure(self):
        if len(self.glob_ss_rects) == 0:
            alpha, beta = self._get_alpha_beta(self.protein_map.secondary_structure)
            alpha_range = self._get_ranges(alpha)
            beta_range = self._get_ranges(beta)
            for a in alpha_range:
                rect = self._get_rectangle_x(a, '#FF1654')
                self.glob_ss_rects.append(rect)
                self.axes.add_patch(rect)
                rect = self._get_rectangle_y(a, '#FF1654')
                self.glob_ss_rects.append(rect)
                self.axes.add_patch(rect)
            for b in beta_range:
                rect = self._get_rectangle_x(b, '#06D6A0')
                self.glob_ss_rects.append(rect)
                self.axes.add_patch(rect)
                rect = self._get_rectangle_y(b, '#06D6A0')
                self.glob_ss_rects.append(rect)
                self.axes.add_patch(rect)
            self.canvas.draw()
            self.add_rect()
        else:
            for r in self.glob_ss_rects:
                r.set_visible(True)
            self.canvas.draw()
            self.add_rect()

    def remove_secondary_structure_overlay(self):
        for r in self.glob_ss_rects:
            r.set_visible(False)
        self.canvas.draw()

    def overlay_charged_interactions(self, array):
        if len(self.glob_charged_dots) == 0:
            positive = ['ARG', 'LYS']
            negative = ['ASP', 'GLU']
            if self.protein_map.use_ca:
                residue_names = self.protein_map.residue_names_ca
            else:
                residue_names = self.protein_map.residue_names_cb
            shape = array.shape
            for i in range(shape[0]):
                for j in range(shape[1]):
                    if array[j][i] > 0:
                        if ((residue_names[self.parent.current_model][i] in positive) and
                                (residue_names[self.parent.current_model][j] in negative)) or \
                                ((residue_names[self.parent.current_model][j] in positive) and
                                (residue_names[self.parent.current_model][i] in negative)):
                            dot = self._get_dot(i, j, '#00A6ED')
                            self.axes.add_patch(dot)
                            self.glob_charged_dots.append(dot)
            self.canvas.draw()
            self.add_rect()
        else:
            for r in self.glob_charged_dots:
                r.set_visible(True)
            self.canvas.draw()
            self.add_rect()

    def remove_charged_interactions(self):
        for r in self.glob_charged_dots:
            r.set_visible(False)
        self.canvas.draw()

    def overlay_hydrophobic_interactions(self, array):
        if len(self.glob_hydrophobic_dots) == 0:
            hydrophobic = ['ALA', 'ILE', 'LEU', 'PHE', 'TRP', 'VAL', 'MET', 'PRO']
            if self.protein_map.use_ca:
                residue_names = self.protein_map.residue_names_ca
            else:
                residue_names = self.protein_map.residue_names_cb
            shape = array.shape
            for i in range(shape[0]):
                for j in range(shape[1]):
                    if array[j][i] > 0:
                        if ((residue_names[self.parent.current_model][i] in hydrophobic) and
                                (residue_names[self.parent.current_model][j] in hydrophobic)):
                            dot = self._get_dot(i, j, '#F0C808')
                            self.axes.add_patch(dot)
                            self.glob_hydrophobic_dots.append(dot)
            self.canvas.draw()
            self.add_rect()
        else:
            for r in self.glob_hydrophobic_dots:
                r.set_visible(True)
            self.canvas.draw()
            self.add_rect()

    def remove_hydrophobic_interactions(self):
        for r in self.glob_hydrophobic_dots:
            r.set_visible(False)
        self.canvas.draw()


    def overlay_bfactors(self, array):
        if len(self.glob_bfactor_dots) == 0:
            if self.protein_map.use_ca:
                bfactors = self.protein_map.bfactors_ca
            else:
                bfactors = self.protein_map.bfactors_cb

            shape = array.shape
            for i in range(shape[0]):
                for j in range(shape[1]):
                    if array[j][i] > 0:
                        bfactor = bfactors[self.parent.current_model][i]
                        if bfactor > self.protein_map.bfactor_cutoff:
                            dot = self._get_dot(i, j, '#FB3640')
                            self.axes.add_patch(dot)
                            self.glob_bfactor_dots.append(dot)
            self.canvas.draw()
            self.add_rect()
        else:
            for r in self.glob_bfactor_dots:
                r.set_visible(True)
            self.canvas.draw()
            self.add_rect()

    def remove_bfactors(self):
        for r in self.glob_bfactor_dots:
            r.set_visible(False)
        self.canvas.draw()

    def reset_overlays(self):
        self.glob_charged_dots = []
        self.glob_ss_rects = []
        self.glob_hydrophobic_dots = []
        self.parent.ss_showing = False
        self.parent.charged_showing = False
        self.parent.hydrophobic_showing = False

    def update_plot(self, image):
        self.axes.clear()
        self.axes.imshow(image, interpolation='nearest', cmap='gray', aspect='equal', origin='lower')
        self.set_parameters()
        self.add_rect(False)
        self.canvas.draw()
        self.reset_overlays()



