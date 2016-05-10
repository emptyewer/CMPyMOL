import os
import sys
import signal
from copy import copy
from PyQt4 import QtCore, QtGui, uic

import functions.connect as connect
import functions.pymolcommands as pymolcommands
import functions.proteinmap as pmap
from functions.about import About_Dialog

import graphs.heatmap as heatmap
import graphs.contactmap as contactmap
import graphs.histogram as histogram
import overlays.overlays as overlays

app = QtGui.QApplication(sys.argv)
form_class, base_class = uic.loadUiType(os.path.join('ui', 'window.ui'))


class CMPyMOL_GUI(QtGui.QMainWindow, form_class):
    def __init__(self, *args):
        super(CMPyMOL_GUI, self).__init__(*args)
        self.setupUi(self)
        self.connection = connect.connect(self)
        # Get PyMOL path if not set already
        self.connection.get_pymol_path()
        # Launch PyMOL
        self.connection.make_pymol_connection()
        # PDB
        self.pymolcommands = pymolcommands.pymolcommands(self)
        self.overlays = overlays.Overlays(self)
        self.protein_map = None
        self.heatmap_fig = None
        self.contactmap_fig = None
        self.histogram_fig = None
        self.current_model = 1
        self.model_limit = 1
        self.ss_showing = False
        self.charged_showing = False
        self.hydrophobic_showing = False
        self.show_contact_map = True
        self.bfactor_showing = False
        self.selected_aa = []


    @QtCore.pyqtSlot()
    def on_abt_btn_clicked(self):
        about_dialog = About_Dialog()
        about_dialog.exec_()

    @QtCore.pyqtSlot()
    def on_quit_btn_clicked(self):
        self.append_status_text(">>> Quitting CMPyMOL. Bye...")
        self.appExit()

    @QtCore.pyqtSlot()
    def on_loadpdb_btn_clicked(self):
        self.pymolcommands.load_pdb()
        self.current_model = 1
        if self.pymolcommands.pdb_path != '':
            self.pymolcommands.reinitialize_pymol()
            self.pymolcommands.send_pdb_to_pymol()
            self.pymolcommands.init_pdb_representation()
            self.protein_map = pmap.pmap(self.pymolcommands.pdb_path, self.distance_cutoff.value(), self.ca.isChecked(),
                                         self)

            for aa in self.protein_map.aminoacids:
                child = self.findChild(QtGui.QCheckBox, aa)
                child.stateChanged.connect(self.on_aa_box_toggled)

            self.model_limit = self.protein_map.model_count
            self.plot_contactmap()
            self.plot_heatmap()
            self.plot_histogram()

            self.frame_label.setText("Frame %d" % self.current_model)
            if self.model_limit == 1:
                self.variance_map_btn.setEnabled(False)
                self.next_btn.setEnabled(False)
                self.prev_btn.setEnabled(False)
            else:
                self.variance_map_btn.setEnabled(True)
                self.next_btn.setEnabled(True)
                self.prev_btn.setEnabled(True)
            self.overlays.calculate_secondary_structure()

    @QtCore.pyqtSlot()
    def on_next_btn_clicked(self):
        if self.protein_map:
            if self.current_model + 1 <= self.model_limit:
                self.current_model += 1
                self.frame_label.setText("Frame %d" % self.current_model)
            self.pymolcommands.change_frame(self.current_model)
            self.update_contactmap()
            self.update_heatmap()


    @QtCore.pyqtSlot()
    def on_prev_btn_clicked(self):
        if self.protein_map:
            if self.current_model - 1 >= 1:
                self.current_model -= 1
                self.frame_label.setText("Frame %d" % self.current_model)
            self.pymolcommands.change_frame(self.current_model)
            self.update_contactmap()
            self.update_heatmap()

    @QtCore.pyqtSlot(float)
    def on_distance_cutoff_valueChanged(self, value):
        if self.protein_map:
            self.protein_map.cutoff = value
            self.update_contactmap()
            self.update_heatmap()

    @QtCore.pyqtSlot(bool)
    def on_ca_toggled(self, val):
        if self.protein_map:
            if self.ca.isChecked():
                self.protein_map.use_ca = True
            else:
                self.protein_map.use_ca = False
            self.update_contactmap()
            self.update_heatmap()

    @QtCore.pyqtSlot()
    def on_aa_box_toggled(self):
        for aa in self.protein_map.aminoacids:
            child = self.findChild(QtGui.QCheckBox, aa)
            if child.isChecked():
                if not aa in self.selected_aa:
                    self.selected_aa.append(aa)
            else:
                if aa in self.selected_aa:
                    self.selected_aa.remove(aa)

        if len(self.selected_aa) > 2:
            for aa in self.selected_aa[:-2]:
                child = self.findChild(QtGui.QCheckBox, aa)
                child.setChecked(False)
                if aa in self.selected_aa:
                    self.selected_aa.remove(aa)

        if len(self.selected_aa) == 2:
            if self.protein_map.use_ca:
                atom = 'CA'
            else:
                atom = 'CB'
            self.pymolcommands.select_aminoacids(self.selected_aa[0], self.selected_aa[1], atom,
                                                 self.distance_cutoff.value())

        if len(self.selected_aa) < 2:
            self.pymolcommands.reset_pdb_view()

    @QtCore.pyqtSlot()
    def on_ss_btn_clicked(self):
        if not self.ss_showing:
            self.contactmap_fig.overlay_secondary_structure()
            self.ss_showing = True
        else:
            self.contactmap_fig.remove_secondary_structure_overlay()
            self.ss_showing = False

    @QtCore.pyqtSlot()
    def on_charged_btn_clicked(self):
        if not self.charged_showing:
            self.contactmap_fig.overlay_charged_interactions(self.contact_array)
            self.charged_showing = True
        else:
            self.contactmap_fig.remove_charged_interactions()
            self.charged_showing = False

    @QtCore.pyqtSlot()
    def on_hydrophobic_btn_clicked(self):
        if not self.hydrophobic_showing:
            self.contactmap_fig.overlay_hydrophobic_interactions(self.contact_array)
            self.hydrophobic_showing = True
        else:
            self.contactmap_fig.remove_hydrophobic_interactions()
            self.hydrophobic_showing = False

    @QtCore.pyqtSlot()
    def on_variance_map_btn_clicked(self):
        if self.show_contact_map:
            self.show_contact_map = False
            self.map_label.setText('Variance Map')
            self.variance_map_btn.setText('Show Contact Map')
        else:
            self.show_contact_map = True
            self.map_label.setText('Contact Map')
            self.variance_map_btn.setText('Show Variance Map')
        self.update_contactmap()

    @QtCore.pyqtSlot(int)
    def on_bfactor_slider_valueChanged(self, value):
        self.bfactor_label.setText('B-factor Cut-Off (%d)' % value)
        if self.protein_map:
            self.protein_map.bfactor_cutoff = value

    @QtCore.pyqtSlot()
    def on_bfactor_slider_sliderReleased(self):
        if self.bfactor_showing:
            self.contactmap_fig.remove_bfactors()
            self.contactmap_fig.glob_bfactor_dots = []
            self.contactmap_fig.overlay_bfactors(self.contact_array)

    @QtCore.pyqtSlot()
    def on_bfactor_btn_clicked(self):
        if not self.bfactor_showing:
            self.contactmap_fig.overlay_bfactors(self.contact_array)
            self.bfactor_showing = True
        else:
            self.contactmap_fig.remove_bfactors()
            self.bfactor_showing = False





    def append_status_text(self, text):
        cursor = self.status_text.textCursor()
        cursor.movePosition(QtGui.QTextCursor.EndOfLine)
        cursor.insertText(text + "\n")

    def rescaled_contact_map(self):
        self.contact_array = copy(self.protein_map.contact_maps_ca[self.current_model]) if self.protein_map.use_ca \
                             else copy(self.protein_map.contact_maps_cb[self.current_model])
        self.variance_array = copy(self.protein_map.variance_maps_ca[self.current_model]) if self.protein_map.use_ca \
                              else copy(self.protein_map.variance_maps_cb[self.current_model])
        self.contact_array[self.contact_array >= self.protein_map.cutoff] = 0.0
        self.contact_array[self.contact_array != 0.0] = 255.0

        self.variance_array *= self.variance_array / self.variance_array.max()
        self.variance_array[self.variance_array > 0.1] = 255.0
        self.variance_array[self.contact_array >= self.protein_map.cutoff] = 0.0

    def plot_contactmap(self):
        if self.contactmap_fig:
            self.contactmap_layout.removeWidget(self.contactmap_fig)
        self.rescaled_contact_map()
        self.contactmap_fig = contactmap.ContactMapWidget(self.pymolcommands, self.protein_map, self)
        if self.show_contact_map:
            self.contactmap_fig.plot(self.contact_array)
        else:
            self.contactmap_fig.plot(self.variance_array)
        self.contactmap_layout.addWidget(self.contactmap_fig)

    def update_contactmap(self):
        self.rescaled_contact_map()
        if self.show_contact_map:
            self.contactmap_fig.update_plot(self.contact_array)
        else:
            self.contactmap_fig.update_plot(self.variance_array)

    def plot_heatmap(self):
        if self.heatmap_fig:
            self.tabWidget.removeTab(1)
        self.protein_map.calculate_heatmap_histogram()
        self.heatmap_fig = heatmap.HeatMapFrame(self.protein_map.heat_maps[self.current_model],
                                                self.protein_map.aminoacids, self.pymolcommands, self.protein_map)
        self.tabWidget.addTab(self.heatmap_fig, "Heat Map")

    def update_heatmap(self):
        self.protein_map.calculate_heatmap_histogram()
        self.heatmap_fig.graph_view.update_graph(self.protein_map.heat_maps[self.current_model])
        self.histogram_fig.graph_view.update_graph()

    def plot_histogram(self):
        if self.histogram_fig:
            self.tabWidget.removeTab(1)
        self.protein_map.calculate_heatmap_histogram()
        self.histogram_fig = histogram.HistogramFrame(self.pymolcommands, self.protein_map)
        self.tabWidget.addTab(self.histogram_fig, "Histogram")

    def appExit(self):
        os.kill(self.connection.pymol_pid, signal.SIGTERM)
        app.quit()
        sys.exit()

if __name__ == '__main__':
    form = CMPyMOL_GUI()
    form.show()
    app.aboutToQuit.connect(form.appExit)
    app.exec_()
