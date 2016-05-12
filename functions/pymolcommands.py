import os
from itertools import groupby
from operator import itemgetter
from PyQt4 import QtGui

class pymolcommands():
    def __init__(self, parent):
        self.parent = parent
        self.connection = self.parent.connection
        self.pdb_path = None

    def load_pdb(self):
        file_filter = "*.pdb"
        self.pdb_path = str(QtGui.QFileDialog.getOpenFileName(QtGui.QFileDialog(),
                                                              caption="Load PDB...",
                                                              filter="PDB files (%s)" % file_filter))
        self.parent.append_status_text(">>> PDB %s loaded" % (os.path.basename(self.pdb_path)))


    def reinitialize_pymol(self):
        self.connection.connected_socket.do('reinitialize')

    def send_pdb_to_pymol(self):
        self.connection.connected_socket.do('load %s' % self.pdb_path)

    def init_pdb_representation(self):
        self.connection.connected_socket.do("as cartoon")
        # self.connection.connected_socket.do("show lines")
        self.connection.connected_socket.do("color gray60")
        self.connection.connected_socket.do("set seq_view, on")
        self.connection.connected_socket.do("set seq_view_format, 1")
        self.connection.connected_socket.do('set seq_view_label_spacing, 2')
        self.connection.connected_socket.do('set dash_color, teal')
        self.connection.connected_socket.do('set label_position, [-2,0,0]')
        self.connection.connected_socket.do('set dash_gap, 0')
        self.connection.connected_socket.do('set dash_width, 0.7')
        self.connection.connected_socket.do('set dash_length, 0.3')
        self.connection.connected_socket.do('set dash_use_shader, 1')
        self.connection.connected_socket.do('set label_color, yellow')
        self.connection.connected_socket.do('set sphere_scale, 0.4')

    def reset_pdb_view(self):
        self.connection.connected_socket.do("as cartoon")
        # self.connection.connected_socket.do("show lines")
        self.connection.connected_socket.do("color gray60")

    def select_aminoacids(self, aa1, aa2, atom, cutoff):
        self.reset_pdb_view()
        self.connection.connected_socket.do('select %s-%s, resname %s and name %s within %f of resname %s and name %s' %
                                            (aa1, atom, aa1, atom, cutoff, aa2, atom))
        self.connection.connected_socket.do('select %s-%s, resname %s and name %s within %f of resname %s and name %s' %
                                            (aa2, atom, aa2, atom, cutoff, aa1, atom))
        self.connection.connected_socket.do('select %s, br. %s-%s' % (aa1, aa1, atom))
        self.connection.connected_socket.do('select %s, br. %s-%s' % (aa2, aa2, atom))
        self.connection.connected_socket.do('show sticks, %s' % aa1)
        self.connection.connected_socket.do('show sticks, %s' % aa2)
        self.connection.connected_socket.do('color tv_red, %s' % aa1)
        self.connection.connected_socket.do('color marine, %s' % aa2)
        self.connection.connected_socket.do('show spheres, %s-%s' % (aa1, atom))
        self.connection.connected_socket.do('show spheres, %s-%s' % (aa2, atom))
        self.connection.connected_socket.do('distance distance-%s-%s-%f, %s-%s, %s-%s, cutoff=%f' %
                                            (aa1, aa2, cutoff, aa1, atom, aa2, atom, cutoff))
        self.connection.connected_socket.do('group %s_%s, %s %s %s-%s %s-%s distance-%s-%s-%f' %
                                            (aa1, aa2, aa1, aa2, aa1, atom, aa2, atom, aa1, aa2, cutoff))
        self.connection.connected_socket.do('zoom animate=1')

    def get_residues_range(self, range):
        text = ''
        int_range = []
        for i in range:
            int_range.append(int(i))

        for k, g in groupby(enumerate(sorted(set(int_range))), lambda (i, x): i - (x)):
            list = map(itemgetter(1), g)
            if len(list) > 1:
                text += "%d-%d," % (list[0], list[-1])
            else:
                text += "%d," % (list[0])
        return text

    def select_range(self, x_range, y_range, chainx, chainy, atom, count):
        self.reset_pdb_view()
        if len(x_range) != 0 and len(y_range) != 0:
            self.connection.connected_socket.do('select Range1_%d, resid %s and chain %s' %
                                                (count, "+".join(x_range), "+".join(chainx)))
            self.connection.connected_socket.do('select Range2_%d, resid %s and chain %s' %
                                                (count, "+".join(y_range), "+".join(chainy)))
            self.connection.connected_socket.do('select Range1_%s_%d, resid %s and name %s and chain %s' %
                                                (atom, count, "+".join(x_range), atom, "+".join(chainx)))
            self.connection.connected_socket.do('select Range2_%s_%d, resid %s and name %s and chain %s' %
                                                (atom, count, "+".join(y_range), atom, "+".join(chainy)))
            self.connection.connected_socket.do('color hotpink, Range1_%d' % count)
            self.connection.connected_socket.do('color limegreen, Range2_%d' % count)
            self.connection.connected_socket.do('show sticks, Range1_%d' % count)
            self.connection.connected_socket.do('show sticks, Range2_%d' % count)
            self.connection.connected_socket.do('show spheres, Range1_%s_%d' % (atom, count))
            self.connection.connected_socket.do('show spheres, Range2_%s_%d' % (atom, count))
            self.connection.connected_socket.do('distance Distance_%d, Range1_%s_%d, Range2_%s_%d, cutoff=%d' %
                                                (count, atom, count, atom, count, self.parent.protein_map.cutoff))
            self.connection.connected_socket.do('zoom Distance_%d, 5, animate=1' % count)
            self.connection.connected_socket.do('group Selection_%d, Range1_%d Range2_%d Range1_%s_%d Range2_%s_%d '
                                                'Distance_%d' % (count, count, count, atom, count, atom, count, count))

            self.parent.sel1_label.setText("Chain: %s\n\nResidues: %s" % (",".join(chainx), self.get_residues_range(
                x_range)))
            self.parent.sel2_label.setText("Chain: %s\n\nResidues: %s" % (",".join(chainy), self.get_residues_range(
                y_range)))

    def select_density(self, index, atom, cutoff, chain):
        self.reset_pdb_view()
        self.connection.connected_socket.do('select index-%d-%s, resid %d and name %s and chain %s' %
                                            (index, atom, index, atom, chain))
        self.connection.connected_socket.do('show spheres, index-%d-%s' % (index, atom))
        # self.connection.connected_socket.do('select index-%d, index %d' % (index, index))
        # self.connection.connected_socket.do('show sticks, index-%d' % index)
        self.connection.connected_socket.do('select Surround-%d, byres((index-%d-%s around %f) and chain %s and name '
                                            '%s)' % (index, index, atom, cutoff, chain, atom))
        self.connection.connected_socket.do('show sticks, Surround-%d' % index)
        self.connection.connected_socket.do('util.cbag Surround-%d' % index)
        self.connection.connected_socket.do('color magenta, index-%d-%s' % (index, atom))
        self.connection.connected_socket.do('center index-%d-%s, animate=1' % (index, atom))
        self.connection.connected_socket.do('zoom Surround-%d, 5, animate=1' % index)


    def change_frame(self, frame):
        self.connection.connected_socket.do('frame %d' % frame)












