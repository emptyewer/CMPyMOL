#!/opt/local/bin/python


'''
CMPyMOL 1.5

https://github.com/emptyewer/CMPyMOL

Author: Venkatramanan Krishnamani (Version 1.5)

'''

#  The MIT License (MIT)
# =======================
#
# CMPyMOL source code in this file is copyrighted, but you are
# free to use and copy it as long as you don't change or remove any of
# the copyright notices.
#
# -----------------------------------------------------------------------------------
# CMPyMOL
# Copyright (C) 2015 Venkatramanan Krishnamani
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be included in all copies
# or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#

import sys
import os
import re
import wx
import sqlite3 as lite
import subprocess
import cPickle as pickle

try:
    import matplotlib as mpl
    mpl.use('WXAgg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.patches import Rectangle
    from matplotlib.widgets import Button
    from matplotlib.widgets import Slider
    from matplotlib.widgets import CheckButtons
except ImportError:
    print '*** Error: This software requires the "matplotlib" module ***'

try:
    if sys.platform == 'darwin':
        import PIL.Image as Image
        import PIL.ImageDraw as ImageDraw
    else:
        import Image
        import ImageDraw
    from PIL import ImageOps
    from operator import itemgetter
    from itertools import groupby
except ImportError:
    print "*** Error: Cannot Import operator or itertools. ***"
    sys.exit()

# Customizing Plot Parameters
mpl.rcParams['figure.figsize'] = [11, 9]
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['grid.color'] = 'white'
mpl.rcParams['axes.labelsize'] = 'medium'
mpl.rcParams['xtick.major.pad'] = 10
mpl.rcParams['ytick.major.pad'] = 10
mpl.rcParams['axes.titlesize'] = 'x-large'
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['axes.labelweight'] = 'bold'

debug = True
progressbar = None
frametext = None
visualtext = None
variance_map_showing = 0
# Secondary Structure
# (Boolean) If secondary Structure calculation is previously performed
ss_calculated = 0
# (Boolean) If Secondary Structure Overlay is Visible
ss_showing = 0
# (matplotlib.widgets.Rectangle) Describes the Secondary Structure Overlay
glob_ss_rects = []
# Charged Interactions
# (matplotlib.image.AxesImage) Charged Interactions Overlay
charge_img_inst = None
# (Boolean) If Charges Interactions Overlay is previously calculated
charge_calculated = 0
# (Boolean) If Charges Interactions Overlay is Visible
charge_showing = 0
# Hydrophobic Interactions
# (matplotlib.image.AxesImage) Hydrophobic Interactions Overlay
hp_img_inst = None
# (Boolean) If Hydrophobic Interactions is previously calculated
hp_calculated = 0
# (Boolean) If Hydrophobic Interactions is Visible
hp_showing = 0
# B-factor
bfac_cutoff = 25                # (Int) Default cutoff for B-factor calculation
bfac_img_inst = None            # (matplotlib.image.AxesImage) B-factor Overlay
# (Boolean) If B-factor Overlay is previously calculted
bfac_calculated = 0
bfac_showing = 0                # (Boolean) If B-factor Overlay is showing
bfac_sliderax = None            # (matplotlib axes) B-factor Slider Axes
# Interacting Aminoacids Overlay
# (List) List of aminoacids chosen from checkedbuttons
buttons_checked = []
# (Boolean) If Aminoacid overlay is previously calculated
aa_calculated = 0
# (matplotlib.image.AxesImage) Aminoacids Overlay
aa_img_inst = None
# Heatmap
# (Boolean) Aminoacid-Aminoacid Contacts previously calculated
heatmap_calculated = 0
# Contact Density Map
condensmap_calculated = 0       # (Boolean) If Contact Density Map Calculated
# Database variables
db_path = ''
db_handle = None
db_cursor = None
# Connection to PyMOL
connectedSocket = ''
# program locations
stridepath = ''
pymolpath = ''
# Contact Map variables
cmap_atom = r'CA'
# (Float) Distance Cutoff for Contact Map Calculation
distance_cutoff = 8
#Image or Calculate
calculate_cmap = 1
image_filepath = ''
# PDB variables
pdb_path = ''
pdb_istrajectory = False
# Trajectory vaiables
model_indexes = []
# PDB properties
residue_index_map = {}
chain_index_map = {}
pdbresiduecount = 0
processed_frames = 0
x_img_length = 0
y_img_length = 0
glob_image = None
glob_ax = None
glob_text = None
progressbar = None
glob_con_ax = None
glob_heat_ax = None
contact_density_closed = True
heat_map_closed = True
pymol_pid = 0

current_model_index = 0
COMAP = None
HEMAP = None
CDENS = None
VARIA = None
FRAME = ''
NAME = ''


class Help:
    '''
    Calculates and displays hydrogen bonding partners network on top of contact map.
    '''

    def showabout(self, event):
        try:
            description = """CMPyMOL is an add-on software to molecular visualization program PyMOL.
It combines the 3D visualization capabilities of PyMOL and 2D protein contact maps with
an interactive interface for scientific analysis."""
            licence = """The MIT License (MIT)
=======================

CMPyMOL source code in this file is copyrighted, but you are
free to use and copy it as long as you don't change or remove any of
the copyright notices.

-----------------------------------------------------------------------------------
CMPyMOL
Copyright (C) 2016 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE."""
            info = wx.AboutDialogInfo()
            info.SetName('CMPyMOL')
            info.SetVersion('1.5')
            info.SetDescription(description)
            info.SetCopyright('(C) 2016 Venkatramanan Krishnamani')
            info.SetWebSite('https://github.com/emptyewer/CMPyMOL')
            info.SetLicence(licence)
            info.AddDeveloper(
                'Venkatramanan Krishnamani\nUniversity of Iowa\nIowa City, IA')
            wx.AboutBox(info)
        finally:
            wx.YieldIfNeeded()

    def showparameters(self, event):
        app = Info(0)
        app.MainLoop()


class InfoFrame(wx.Frame):

    def __init__(self, parent, id, title):
        wx.Frame.__init__(
            self,
            parent,
            id,
            title,
            wx.DefaultPosition,
            wx.Size(
                640,
                490))
        cutoff = "Distance Cut-Off = %.1f angstrom" % (distance_cutoff)
        atom = "Atom for Contact Map Calculation = %s " % (cmap_atom)
        pdb = "Input PDB file = %s" % (pdb_path)
        pymol = "PyMOL Path = %s" % (pymolpath)
        stride = "Stride Path = %s" % (stridepath)
        image_size = "Size of Contact Map = %d pixels x %d pixels" % (
            x_img_length, y_img_length)
        traj_info = "Input is a Trajectory = %r" % (pdb_istrajectory)
        frames_info = "Number of Frames = %d" % (len(model_indexes) + 1)
        current_frame_info = "Current Frame = %d" % (current_model_index + 1)
        database_info = "Database File Name = %s" % (db_path)
        p = subprocess.Popen(['du',
                              '-h',
                              os.path.join(os.environ['HOME'],
                                           '.CMapperDir')],
                             stdout=subprocess.PIPE)
        out, err = p.communicate()
        cache_info = "Cache Size : %s" % (out.split()[0])
        panel = wx.Panel(self, -1)
        wx.StaticBox(panel, label='', pos=(5, 5), size=(630, 455))
        wx.StaticText(panel, -1, pdb, (25, 25), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, atom, (25, 50), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, cutoff, (25, 75), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, pymol, (25, 100), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, stride, (25, 125), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, image_size, (25, 150), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, traj_info, (25, 175), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, frames_info, (25, 200), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, current_frame_info,
                      (25, 225), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, database_info, (25, 250), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, ur'   \u03b1lpha Helix',
                      (65, 280), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, ur'   \u03b2eta Sheets',
                      (65, 305), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, '   Charged Residues',
                      (65, 330), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, '   Hydrophobic Residues',
                      (65, 355), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, '   Residues with B-Factor above slider value',
                      (65, 380), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, '   Aminoacids Interaction Sites (Checkboxes)',
                      (65, 405), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, '   Mouse Selections',
                      (65, 430), style=wx.ALIGN_LEFT)
        wx.StaticText(panel, -1, cache_info, (500, 395), style=wx.ALIGN_LEFT)
        button = wx.Button(
            panel,
            id=wx.ID_ANY,
            label="Clear Cache!",
            pos=(
                505,
                420))
        button.Bind(wx.EVT_BUTTON, self.onButton)
        self.Centre()
        self.Bind(wx.EVT_PAINT, self.OnPaint)

    def onButton(self, evt):
        import glob
        filelist = glob.glob(
            os.path.join(
                os.environ['HOME'],
                '.CMapperDir',
                '*_frame_*'))
        for f in filelist:
            os.remove(f)

    def OnPaint(self, evt):
        import time
        import random
        self.dc = wx.PaintDC(self)
        self.dc.Clear()
        self.dc.BeginDrawing()
        self.dc.SetPen(wx.Pen("BLACK", 1))
        self.dc.SetBrush(wx.Brush((255, 0, 0), wx.SOLID))
        self.dc.DrawRectangle(25, 280, 45, 18)
        self.dc.SetBrush(wx.Brush((0, 255, 0), wx.SOLID))
        self.dc.DrawRectangle(25, 305, 45, 18)
        self.dc.SetBrush(wx.Brush((0, 0, 255), wx.SOLID))
        self.dc.DrawRectangle(25, 330, 45, 18)
        self.dc.SetBrush(wx.Brush((255, 255, 0), wx.SOLID))
        self.dc.DrawRectangle(25, 355, 45, 18)
        self.dc.SetBrush(wx.Brush((0, 255, 255), wx.SOLID))
        self.dc.DrawRectangle(25, 380, 45, 18)
        self.dc.SetBrush(wx.Brush((255, 132, 0), wx.SOLID))
        self.dc.DrawRectangle(25, 405, 45, 18)
        self.dc.SetBrush(wx.Brush((255, 0, 255), wx.SOLID))
        self.dc.DrawRectangle(25, 430, 45, 18)
        self.dc.EndDrawing()
        del self.dc


class Info(wx.App):

    def OnInit(self):
        frame = InfoFrame(None, -1, 'CMPyMOL Parameters and Coloring Legend')
        frame.Show(True)
        self.SetTopWindow(frame)
        return True


class ImageFileCMap(wx.Frame):

    def __init__(self, *args, **kw):
        super(ImageFileCMap, self).__init__(*args, **kw)
        self.InitUI()

    def InitUI(self):
        pnl = wx.Panel(self)
        wx.StaticBox(
            pnl, label='Contact Map Calculation', pos=(
                5, 5), size=(
                300, 60))
        self.rb1 = wx.RadioButton(
            pnl, label=ur' from PDB', pos=(
                20, 30), style=wx.RB_GROUP)
        self.rb2 = wx.RadioButton(pnl, label=ur' precalculated', pos=(150, 30))
        self.rb1.Bind(wx.EVT_RADIOBUTTON, self.SetVal)
        self.rb2.Bind(wx.EVT_RADIOBUTTON, self.SetVal)
        self.rb1.SetValue(True)
        self.btn = wx.Button(
            pnl, label='Continue', pos=(
                112, 80), size=(
                80, -1))
        self.btn.Bind(wx.EVT_BUTTON, self.OnClose)
        self.SetSize((310, 140))
        self.SetTitle('Select CMap Type')
        self.Centre()
        self.Show(True)

    def SetVal(self, e):
        global calculate_cmap
        if self.rb1.GetValue():
            calculate_cmap = 1
        elif self.rb2.GetValue():
            calculate_cmap = 0

    def OnClose(self, e):
        global calculate_cmap
        global image_filepath
        if calculate_cmap == 0:
            wildcard = "Data File (*.dat)|*.dat"
            dlg1 = wx.FileDialog(
                None,
                "Choose the Contact Map Data file...",
                os.getcwd(),
                "",
                wildcard=wildcard)
            try:
                if dlg1.ShowModal() == wx.ID_OK:
                    image_filepath = dlg1.GetPath()
                    print ">>> Contact Map Data File:", image_filepath
                else:
                    print "*** Error: No Contact Map Image file was chosen. Exiting..."
                    sys.exit()
            finally:
                dlg1.Destroy()
                wx.YieldIfNeeded()
        else:
            Options(None)
        self.Destroy()


class Options(wx.Frame):

    def __init__(self, *args, **kw):
        super(Options, self).__init__(*args, **kw)
        self.InitUI()

    def InitUI(self):
        pnl = wx.Panel(self)
        wx.StaticBox(
            pnl, label='Calculate distance matrix based on atom...', pos=(
                5, 5), size=(
                300, 60))
        self.rb1 = wx.RadioButton(
            pnl,
            label=ur' C-\u03b1lpha',
            pos=(
                20,
                30),
            style=wx.RB_GROUP)
        self.rb2 = wx.RadioButton(pnl, label=ur' C-\u03b2eta', pos=(150, 30))
        self.rb1.Bind(wx.EVT_RADIOBUTTON, self.SetVal)
        self.rb2.Bind(wx.EVT_RADIOBUTTON, self.SetVal)
        self.rb1.SetValue(True)
        wx.StaticBox(
            pnl, label='Max. Distance that defines a contact....', pos=(
                5, 80), size=(
                300, 60))
        wx.StaticText(pnl, label='Cut-Off Distance', pos=(20, 105))
        self.spnctl = wx.SpinCtrl(
            pnl, value=str(
                int(distance_cutoff)), pos=(
                135, 100), size=(
                80, -1), min=5, max=15)
        self.spnctl.Bind(wx.EVT_SPINCTRL, self.CutoffVal)
        wx.StaticText(pnl, label=ur' \u00c5', pos=(230, 105))
        self.btn = wx.Button(pnl, label='Done', pos=(122, 155), size=(60, -1))
        self.btn.Bind(wx.EVT_BUTTON, self.OnClose)
        self.SetSize((310, 210))
        self.SetTitle('CMPyMOL Parameters')
        self.Centre()
        self.Show(True)

    def CutoffVal(self, e):
        global distance_cutoff
        distance_cutoff = int(self.spnctl.GetValue())

    def SetVal(self, e):
        global cmap_atom
        if self.rb1.GetValue():
            cmap_atom = r'CA'
        elif self.rb2.GetValue():
            cmap_atom = r'CB'

    def OnClose(self, e):
        self.Destroy()


class Initialize:

    def connectPyMOL(self):
        '''
        Opens a XMLRPC socket and connects to PyMOL so information can be passed to it.
        '''
        import xmlrpclib
        global connectedSocket
        connectedSocket = xmlrpclib.Server('http://localhost:9123')
        print ">>> Established connection with PyMOL"

    def start(self):
        import glob
        # Local methods

        def _cmd_exists(cmd):
            import subprocess
            proc = subprocess.Popen(
                ["which", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            status = False
            command = proc.stdout.read().strip()
            if command != '':
                status = True
            return (status, command)

        # Create working directory
        if not os.path.exists(os.path.join(os.environ['HOME'], '.CMapperDir')):
            os.makedirs(os.path.join(os.environ['HOME'], '.CMapperDir'))

        # Locate and launch PyMOL
        global pymolpath

        pymol_loc = _cmd_exists("pymol")

        if pymol_loc[0]:
            global pymol_pid
            pymolpath = pymol_loc[1]
            pyp = subprocess.Popen(
                [pymolpath, '-R'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            pymol_pid = pyp.pid
            print ">>> Lanching PyMOL from location %s (%d)" % (pymolpath, pymol_pid)
        else:
            print "*** ERROR: Cannot locate PyMOL in the system path."
            if sys.platform.startswith('darwin'):
                print "*** Please install PyMOL using Homebrew: brew install homebrew/science/pymol  ***"
            elif sys.platform.startswith("linux"):
                print "*** Please install PyMOL using: apt-get install pymol. ***"

        # Locate stride
        global stridepath
        stride_loc = _cmd_exists("stride")
        if stride_loc[0]:
            stridepath = stride_loc[1]
            print ">>> Stride located at", stridepath
        else:
            print "*** Warning: Cannot locate STRIDE path. Please install from http://webclu.bio.wzw.tum.de/stride/."
            print "*** Note: Secondary structure calculations will be disabled. ***"

    def getPDB(self):
        global pdb_path
        global pdb_istrajectory
        global distance_cutoff
        app = wx.App()
        pdbwildcard = "PDB (*.pdb)|*.pdb|"
        dlg1 = wx.FileDialog(
            None,
            "Choose the PDB file...",
            os.getcwd(),
            "",
            wildcard=pdbwildcard)
        try:
            if dlg1.ShowModal() == wx.ID_OK:
                pdb_path = dlg1.GetPath()
                print ">>> PDB File:", pdb_path
                f = open(pdb_path, 'r')
                # This for loop checks if the PDB file is a trajectory by simply
                # looking at the number of MODEL entries in the inputfile.
                fcount = 0
                for l in f:
                    if re.match(r'MODEL', l):
                        fcount = fcount + 1
                    if fcount >= 2:
                        print ">>> A PDB trajectory file is provided as input. Models will be loaded as frames in PyMOL."
                        pdb_istrajectory = True
                        break
            else:
                print "*** Error: No PDB file was chosen. Exiting..."
                sys.exit()
        finally:
            dlg1.Destroy()
            wx.YieldIfNeeded()
        # Get Options
        ImageFileCMap(None)
        app.MainLoop()
        # Create database
        global db_handle
        global db_cursor
        global db_path
        db_path = os.path.join(os.path.basename(pdb_path)[
                               :-4] + '_' + str(int(distance_cutoff)) + '_' + cmap_atom + '_CMPyMOL_db.sqlite3')
        if os.path.isfile(db_path):
            pass
        else:
            subprocess.Popen(['touch', db_path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print ">>> Creating a new database..."
        db_handle = lite.connect(db_path, isolation_level=None)
        db_cursor = db_handle.cursor()
        db_cursor.execute(
            "CREATE TABLE IF NOT EXISTS data(id integer primary key, name text, path text, contactmap blob, heatmap blob, contactdensity blob, sumdiffsquares blob, variance blob)")
        db_cursor.execute(
            "CREATE TABLE IF NOT EXISTS average(id integer primary key, sumcontactmap blob, processed integer)")
        db_handle.commit()
        print ">>> Sucessfully connected to the database..."

    def get_model_indexes(self):
        global model_indexes
        if pdb_istrajectory:
            linecount = 1
            f = open(pdb_path, 'r')
            for l in f:
                if re.match(r'^MODEL', l):
                    model_indexes.append(linecount)
                linecount = linecount + 1


class PDBfunctions:
    '''
    Calculates Contact Map and other PDB file operations.
    '''

    def unique(self, seq, idfun=None):
        # order preserving
        if idfun is None:
            def idfun(x): return x
        seen = {}
        result = []
        for item in seq:
            marker = idfun(item)
            if marker in seen:
                continue
            seen[marker] = 1
            result.append(item)
        return sorted(result)

    def map_residues(self, pdb):
        import re
        global residue_index_map
        global chain_index_map
        pFile = open(pdb, "r")
        index = 1
        for line in pFile:
            line_list = list(line.rstrip())
            if re.match(
                r'ATOM', ''.join(
                    line_list[
                        0:6]).strip()) and re.match(
                cmap_atom, ''.join(
                    line_list[
                        12:16]).strip()):
                resid = int(''.join(line_list[22:26]).strip())
                chainid = ''.join(line_list[21:22]).strip()
                chain_index_map[index] = chainid
                residue_index_map[index] = resid
                index = index + 1
        pFile.close()

    def generate_contact_map(self, pdb, index=0):
        import numpy as np
        import re
        import math
        import linecache
        import time

        global pdbresiduecount
        # aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        aa = [
            'GLU',
            'ASP',
            'LYS',
            'ARG',
            'HIS',
            'GLN',
            'PRO',
            'ASN',
            'ALA',
            'THR',
            'SER',
            'VAL',
            'GLY',
            'MET',
            'CYS',
            'ILE',
            'LEU',
            'TYR',
            'PHE',
            'TRP']

        def _calc_dist(x, y):
            return math.sqrt((y[0] - x[0]) ** 2 +
                             (y[1] - x[1]) ** 2 + (y[2] - x[2]) ** 2)

        outpath = ''
        pdbarray = []
        if pdb_istrajectory:
            outfile = os.path.basename(pdb)[
                :-4] + '_frame_' + str(index) + "_%.7f.pdb" % time.time()
            outpath = os.path.join(os.environ['HOME'], '.CMapperDir', outfile)
            w = open(outpath, 'w')
            for i in range(model_indexes[index], model_indexes[index + 1]):
                line = linecache.getline(pdb, i)
                line_list = list(line)
                if re.match(
                    r'^ATOM', line) and ''.join(
                    line_list[
                        17:20]).strip() in aa:
                    w.write(line)
                    if re.match(cmap_atom, ''.join(line_list[12:16]).strip()):
                        pdbarray.append(line)
                        pdbresiduecount = pdbresiduecount + 1
            w.close()
        else:
            pdbhandle = open(pdb, 'r')
            outpath = pdb
            for line in pdbhandle:
                line = line.rstrip()
                line_list = list(line)
                if re.match(
                    r'ATOM', ''.join(
                        line_list[
                            0:6]).strip()) and re.match(
                    cmap_atom, ''.join(
                        line_list[
                        12:16]).strip()) and ''.join(
                            line_list[
                                17:20]).strip() in aa:
                    pdbarray.append(line)
                    pdbresiduecount = pdbresiduecount + 1
            pdbhandle.close()

        # Heatmap Initialize
        heatmap = {}
        heatmap_img = []
        for i in aa:
            for j in aa:
                string = i + '-' + j
                heatmap[string] = []

        # Initialize Condensmap
        condensmap = {}
        for n in range(pdbresiduecount):
            condensmap[n] = 0

        # Initialize Distance Matrix
        matrix = []
        matrix = np.zeros((pdbresiduecount, pdbresiduecount), np.float)

        row = 0
        for line1 in pdbarray:
            line_list1 = list(line1.rstrip())
            col = row + 1
            for line2 in pdbarray[col:]:
                line_list2 = list(line2.rstrip())
                dist = _calc_dist(
                    [
                        float(
                            ''.join(
                                line_list1[
                                    30:38]).strip()), float(
                            ''.join(
                                line_list1[
                                    38:46]).strip()), float(
                            ''.join(
                                line_list1[
                                    46:54]).strip())], [
                        float(
                            ''.join(
                                line_list2[
                                    30:38]).strip()), float(
                            ''.join(
                                line_list2[
                                    38:46]).strip()), float(
                            ''.join(
                                line_list2[
                                    46:54]).strip())])
                if (dist < distance_cutoff):
                    if dist == 0:
                        dist = float(distance_cutoff)
                    matrix[row, col] = 1 / dist
                    matrix[col, row] = 1 / dist
                    condensmap[row] = condensmap[row] + 1
                    key1 = ''.join(line_list1[17:20]).strip(
                    ) + '-' + ''.join(line_list2[17:20]).strip()
                    key2 = ''.join(line_list2[17:20]).strip(
                    ) + '-' + ''.join(line_list1[17:20]).strip()
                    if line_list1[21:22][0] != line_list2[21:22][
                            0]:  # only inter-protein interactions are listed
                        heatmap[key1].append(
                            int(''.join(line_list2[22:26]).strip()))
                        heatmap[key2].append(
                            int(''.join(line_list2[22:26]).strip()))
                col = col + 1
            row = row + 1
        # Get unique interactions
        for key in heatmap.keys():
            sequence = heatmap[key]
            if (len(sequence) != 0):
                heatmap[key] = self.unique(sequence)
        for i in aa:
            for j in aa:
                key = i + '-' + j
                heatmap_img.append(len(heatmap[key]))
        heatmap_img = np.asarray(heatmap_img)
        heatmap_img = heatmap_img.reshape((len(aa), len(aa)))
        rescaled_matrix = np.multiply(
            matrix,
            255 /
            (distance_cutoff)).astype(
            np.uint8)

        sum_matrix_prev = []
        sum_matrix = []
        sum_sq_matrix_prev = []
        sum_sq_matrix = []
        variance_matrix = []

        if index > 0:
            _row = []
            while len(_row) < 1:
                db_cursor.execute(
                    "SELECT * FROM data WHERE id=?", (index - 1,))
                _row = db_cursor.fetchall()
                time.sleep(0.1)
            sum_sq_matrix_prev = pickle.loads(str(_row[0][6]))

            _row = []
            while len(_row) < 1:
                db_cursor.execute("SELECT * FROM average WHERE id=?", (0,))
                _row = db_cursor.fetchall()
                time.sleep(0.1)
            sum_matrix_prev = pickle.loads(str(_row[0][1]))

            sum_matrix = np.add(sum_matrix_prev, rescaled_matrix)
            mean_matrix = np.divide(sum_matrix, index + 1)
            sum_sq_matrix = np.add(
                sum_sq_matrix_prev,
                np.square(
                    np.subtract(
                        rescaled_matrix,
                        mean_matrix)))
            variance_matrix = np.sqrt(
                np.divide(
                    sum_sq_matrix,
                    index)).astype(
                np.uint8)
        else:
            sum_matrix = rescaled_matrix
            sum_sq_matrix = np.zeros(
                (pdbresiduecount, pdbresiduecount), np.uint8)
            variance_matrix = np.zeros(
                (pdbresiduecount, pdbresiduecount), np.uint8)
        db_cursor.execute(
            "INSERT INTO data(id, name, path, contactmap, heatmap, contactdensity, sumdiffsquares, variance) VALUES (?,?,?,?,?,?,?,?)",
            (index,
             os.path.basename(pdb),
             outpath,
             lite.Binary(
                 pickle.dumps(
                     rescaled_matrix,
                     protocol=2)),
                lite.Binary(
                 pickle.dumps(
                     heatmap_img,
                     protocol=2)),
                lite.Binary(
                 pickle.dumps(
                     condensmap,
                     protocol=2)),
                lite.Binary(
                 pickle.dumps(
                     sum_sq_matrix,
                     protocol=2)),
                lite.Binary(
                 pickle.dumps(
                     variance_matrix,
                     protocol=2)),
             ))
        db_handle.commit()
        db_cursor.execute(
            "INSERT OR REPLACE INTO average(id, sumcontactmap, processed) VALUES (?,?,?)",
            (0,
             lite.Binary(
                 pickle.dumps(
                     sum_matrix,
                     protocol=2)),
                index +
                1,
             ))
        db_handle.commit()


class ButtonEvents:
    '''
    Calculates and draws various overlays on top of the contact map.
    '''

    def saveheatmap(self, event):
        global HEMAP
        default = os.path.basename(pdb_path)[:-4] + '_interface_heatmap.dat'
        dlg = wx.FileDialog(
            None,
            "Save data as...",
            os.getcwd(),
            default,
            "Heatmap Data (*.dat)|*.dat",
            wx.SAVE | wx.OVERWRITE_PROMPT)
        result = dlg.ShowModal()
        inFile = dlg.GetPath()
        dlg.Destroy()
        if re.search('.+dat', inFile):
            inFile = inFile[:-4]
        if result == wx.ID_OK:  # Save button was pressed
            out = open(inFile + '.dat', 'w')
            for dslice in HEMAP:
                for e in dslice:
                    out.write(str(e) + '\t')
                out.write('\n')
            out.close()
            print "Heatmap written to..." + inFile + '.dat'
            return True
        elif result == wx.ID_CANCEL:  # Either the cancel button was pressed or the window was closed
            return False

    def ss(self, event):
        '''
        Secondary Structure Calculation and Generate a Overlay image to display on top of the Contact Map. Callback event for Secondary Structure Button.
        '''
        global ss_calculated
        global ss_showing
        global progressbar
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()

        def _show_ss():
            '''
            Toggles Seconday Structure Overlay to Visible
            '''
            for r in glob_ss_rects:
                r.set_visible(True)

        def _hide_ss():
            '''
            Toggles Secondary Structure Overlay to Hidden
            '''
            for r in glob_ss_rects:
                r.set_visible(False)

        def _get_rectangle_x(r, c):
            xs = r[0]
            ys = 0
            xe = r[1]
            ye = x_img_length
            rec = Rectangle(
                (xs,
                 ys),
                width=(
                    xe - xs),
                height=(
                    ye - ys),
                alpha=0.30,
                facecolor=c,
                edgecolor=c)
            return rec

        def _get_rectangle_y(r, c):
            xs = 0
            ys = r[0]
            xe = y_img_length
            ye = r[1]
            rec = Rectangle(
                (xs,
                 ys),
                width=(
                    xe - xs),
                height=(
                    ye - ys),
                alpha=0.30,
                facecolor=c,
                edgecolor=c)
            return rec

        def _get_secondary_structure(pdb):
            import subprocess
            import re

            process = subprocess.Popen(
                [stridepath, pdb], stdout=subprocess.PIPE)
            ss = []
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                if re.match(r'ASG ', line):
                    ss.append(line.split()[5])
            return ss

        def _get_ranges(data):
            '''
            Returns Starting and Ending Indexes for secondary structure sequence that is passed on (data)
            '''
            ranges = []
            for k, g in groupby(enumerate(data), lambda i_x: i_x[0] - i_x[1]):
                group = map(itemgetter(1), g)
                ranges.append((group[0], group[-1]))
            return ranges

        def _get_alpha_beta(ss):
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

        if ss_calculated == 0:
            ss = _get_secondary_structure(FRAME)
            alpha, beta = _get_alpha_beta(ss)
            alpha_range = _get_ranges(alpha)
            beta_range = _get_ranges(beta)
            for a in alpha_range:
                rect = _get_rectangle_x(a, 'red')
                glob_ss_rects.append(rect)
                glob_ax.add_patch(rect)
                rect = _get_rectangle_y(a, 'red')
                glob_ss_rects.append(rect)
                glob_ax.add_patch(rect)
            for b in beta_range:
                rect = _get_rectangle_x(b, 'green')
                glob_ss_rects.append(rect)
                glob_ax.add_patch(rect)
                rect = _get_rectangle_y(b, 'green')
                glob_ss_rects.append(rect)
                glob_ax.add_patch(rect)
            ss_calculated = 1

        if ss_showing == 0:
            _show_ss()
            ss_showing = 1
        else:
            _hide_ss()
            ss_showing = 0

        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()

    def charge(self, event):
        '''
        Calculates and generates an overlay for charge-charge interactions. Callback event for charged interactions button.
        '''
        charged_aa = ['ARG', 'LYS', 'ASP', 'GLU']
        global charge_calculated
        global charge_showing
        global charge_img_inst
        global progressbar
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()

        def _get_charged_sequence(pdb):
            '''
            Retrieves all indexes of all charged aminoacids from the sequence
            '''
            import subprocess
            import re
            process = subprocess.Popen(
                [stridepath, pdb], stdout=subprocess.PIPE)
            charged_seq = []
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                if re.match(r'ASG ', line):
                    if (line.split()[1] in charged_aa):
                        charged_seq.append(1)
                    else:
                        charged_seq.append(0)
            return charged_seq

        def _charge_show():
            '''
            Toggles display of Charged interaction overlay to visible
            '''
            charge_img_inst.set_visible(True)

        def _charge_hide():
            '''
            Toggles display of charged interaction overlay to hidden
            '''
            charge_img_inst.set_visible(False)

        if charge_calculated == 0:

            charges_location = _get_charged_sequence(FRAME)
            charged_img = Image.new(
                'RGBA', (x_img_length, y_img_length), (0, 0, 0, 0))
            charged_img_draw = ImageDraw.Draw(charged_img)
            row = 0
            for r in charges_location:
                col = 0
                for c in charges_location:
                    pix = (row, col)
                    if (r == 1 and c == 1) and row != col and glob_image.getpixel(
                            pix)[0] >= 0.8:
                        #charged_img.putpixel(pix, (0, 0, 255, 255))
                        bbox = (pix[0] - 1, pix[1] - 1, pix[0] + 1, pix[1] + 1)
                        charged_img_draw.ellipse(bbox, fill=(0, 0, 255, 235))
                    col = col + 1
                row = row + 1
            charge_img_inst = glob_ax.imshow(
                charged_img, interpolation='nearest', origin='lower')
            charge_calculated = 1

        if charge_showing == 0:
            _charge_show()
            charge_showing = 1
        else:
            _charge_hide()
            charge_showing = 0

        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()

    def hydrophobic(self, event):
        '''
        Calculated and generates an image overlay for interacting hydrophobic aminoacids. Callback event for Hydrophobic Interactions button.
        '''
        hp_aa = ['VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP']
        global hp_calculated
        global hp_showing
        global hp_img_inst
        global progressbar
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()

        def _get_hp_sequence(pdb):
            '''
            Retrieves indexes of hydrophobic aminoacids from the protein sequence
            '''
            import subprocess
            import re
            process = subprocess.Popen(
                [stridepath, pdb], stdout=subprocess.PIPE)
            hp_seq = []
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                if re.match(r'ASG ', line):
                    if (line.split()[1] in hp_aa):
                        hp_seq.append(1)
                    else:
                        hp_seq.append(0)
            return hp_seq

        def _hp_show():
            '''
            Toggles hydrophobic interaction overly to visible
            '''
            hp_img_inst.set_visible(True)

        def _hp_hide():
            '''
            Toggles hydrophobic interaction overlay to hidden
            '''
            hp_img_inst.set_visible(False)

        if hp_calculated == 0:

            hp_location = _get_hp_sequence(FRAME)
            hp_img = Image.new(
                'RGBA', (x_img_length, y_img_length), (0, 0, 0, 0))
            hp_img_draw = ImageDraw.Draw(hp_img)
            row = 0
            for r in hp_location:
                col = 0
                for c in hp_location:
                    pix = (row, col)
                    if (r == 1 and c == 1) and row != col and glob_image.getpixel(
                            pix)[0] >= 0.8:
                        #hp_img.putpixel(pix, (255, 255, 0, 255))
                        bbox = (pix[0] - 1, pix[1] - 1, pix[0] + 1, pix[1] + 1)
                        hp_img_draw.ellipse(bbox, fill=(255, 255, 0, 235))
                    col = col + 1
                row = row + 1
            hp_img_inst = glob_ax.imshow(hp_img, interpolation='nearest')
            hp_calculated = 1

        if hp_showing == 0:
            _hp_show()
            hp_showing = 1
        else:
            _hp_hide()
            hp_showing = 0

        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()

    def bfactor(self, event):
        '''
        Calculates and generates Bfactor overlay. Callback event for Bfactor button.
        '''
        import re
        global bfac_calculated
        global bfac_showing
        global bfac_img_inst
        global progressbar
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()

        bfacs = []
        if (isinstance(event, int) and bfac_calculated == 1 and bfac_showing == 1) or (
                bfac_calculated == 0 and not isinstance(event, int) and bfac_showing == 0):

            pdb_handle = open(FRAME, "r")
            for line in pdb_handle:
                if (re.search(r'ATOM', line[0:5])):
                    if re.search(cmap_atom, line[12:15]):
                        bfacs.append(int(float(line[60:65])))
            bfac_img = Image.new(
                'RGBA', (x_img_length, y_img_length), (0, 0, 0, 0))
            bfac_img_draw = ImageDraw.Draw(bfac_img)
            row = 0
            for r in bfacs:
                col = 0
                for c in bfacs:
                    pix = (row, col)
                    if (r > bfac_cutoff and c > bfac_cutoff) and row != col and glob_image.getpixel(
                            pix)[0] >= 0.8:
                        bbox = (pix[0] - 1, pix[1] - 1, pix[0] + 1, pix[1] + 1)
                        bfac_img_draw.ellipse(bbox, fill=(0, 255, 255, 235))
                    col = col + 1
                row = row + 1
            if (bfac_calculated == 1):
                bfac_img_inst.set_visible(False)
            bfac_calculated = 1
            bfac_img_inst = glob_ax.imshow(bfac_img, interpolation='nearest')
            bfac_showing = 1
        elif (bfac_calculated == 1 and not isinstance(event, int) and bfac_showing == 1):
            bfac_img_inst.set_visible(False)
            bfac_showing = 0
        elif (bfac_calculated == 1 and not isinstance(event, int) and bfac_showing == 0):
            bfac_img_inst.set_visible(True)
            bfac_showing = 1

        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()

    def bfactor_slider(self, event):
        '''
        Callback event when the value of bfactor cutoff slider changes
        '''
        global bfac_cutoff
        bfac_cutoff = int(event)
        self.bfactor(bfac_cutoff)

    def check(self, event):
        '''
        Aminoacid-Aminoacid interactions overlay callback function.
        '''
        import re
        global buttons_checked
        global aa_img_inst
        global progressbar
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()

        def _get_aa_sequence(pdb, residues):
            '''
            Retrieve amino acid sequence.
            '''
            import subprocess
            process = subprocess.Popen(
                [stridepath, pdb], stdout=subprocess.PIPE)
            aa_seq = []
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                if re.match(r'ASG ', line):
                    aa_seq.append(line.split()[1])
            return aa_seq
        aa = re.sub(r'([a-z])', '', event)
        if (aa in buttons_checked):
            buttons_checked.remove(aa)
        else:
            buttons_checked.append(aa)
        if aa_img_inst is not None:
            aa_img_inst.set_visible(False)

        aa_locations = _get_aa_sequence(FRAME, buttons_checked)
        aa_img = Image.new('RGBA', (x_img_length, y_img_length), (0, 0, 0, 0))
        aa_img_draw = ImageDraw.Draw(aa_img)
        row = 0
        for r in aa_locations:
            col = 0
            for c in aa_locations:
                pix = (row, col)
                if len(buttons_checked) == 2:
                    if ((r == buttons_checked[0] and c == buttons_checked[1]) or (r == buttons_checked[
                            1] and c == buttons_checked[0])) and row != col and glob_image.getpixel(pix)[0] >= 0.8:
                        bbox = (pix[0] - 1, pix[1] - 1, pix[0] + 1, pix[1] + 1)
                        aa_img_draw.ellipse(bbox, fill=(255, 132, 0, 255))
                elif len(buttons_checked) == 1:
                    if (r == buttons_checked[0] and c == buttons_checked[
                            0]) and row != col and glob_image.getpixel(pix)[0] >= 0.8:
                        bbox = (pix[0] - 1, pix[1] - 1, pix[0] + 1, pix[1] + 1)
                        aa_img_draw.ellipse(bbox, fill=(255, 132, 0, 255))
                col = col + 1
            row = row + 1
        aa_img_inst = glob_ax.imshow(aa_img, interpolation='nearest')

        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()


class ContactDensity:
    '''
    Calculates and displays a new plot of contact density along the sequence.
    '''

    def showcdensmap(self, event):

        # aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        aa = [
            'GLU',
            'ASP',
            'LYS',
            'ARG',
            'HIS',
            'GLN',
            'PRO',
            'ASN',
            'ALA',
            'THR',
            'SER',
            'VAL',
            'GLY',
            'MET',
            'CYS',
            'ILE',
            'LEU',
            'TYR',
            'PHE',
            'TRP']
        import numpy as np
        import matplotlib.cm as cm
        import matplotlib.colors as colors

        def _onpick(event):
            if event.inaxes == glob_con_ax:
                connectedSocket.do("as lines")
                connectedSocket.do("color gray40")
                connectedSocket.do('select Selection, resi %s and chain %s' % (
                    residue_index_map[int(event.xdata)], chain_index_map[int(event.xdata)]))
                connectedSocket.do('select Surround, byres((residue %s around %f) and chain %s and name %s)' % (
                    residue_index_map[int(event.xdata)], distance_cutoff, chain_index_map[int(event.xdata)], cmap_atom))
                connectedSocket.do('show sticks, Surround')
                connectedSocket.do('util.cbag Surround')
                connectedSocket.do('color magenta, Selection')
                connectedSocket.do('show sticks, Selection')
                connectedSocket.do('center Selection, animate=1')
                connectedSocket.do('disable Selection')
                connectedSocket.do('disable Surround')

        def _on_close(event):
            global contact_density_closed
            contact_density_closed = True

        def _plot(map):
            global glob_con_ax
            global contact_density_closed
            mpl.rcParams['xtick.labelsize'] = 'medium'
            mpl.rcParams['ytick.labelsize'] = 'medium'
            mpl.rcParams['axes.grid'] = False
            mpl.rcParams['figure.subplot.left'] = 0.06
            mpl.rcParams['figure.subplot.right'] = 0.97
            mpl.rcParams['figure.subplot.top'] = 0.9
            mpl.rcParams['figure.subplot.bottom'] = 0.1
            mpl.rcParams['axes.labelsize'] = 'medium'
            mpl.rcParams['ytick.major.pad'] = 4
            mpl.rcParams['xtick.major.pad'] = 4
            fig = plt.figure(figsize=(16, 9))
            glob_con_ax = plt.gca()
            contact_density_closed = False
            fig.canvas.mpl_connect('close_event', _on_close)
            fig.canvas.set_window_title('Contacts Density Histogram')
            fig.canvas.mpl_connect('button_press_event', _onpick)
            bar_heights = np.array(map.values()).astype(float)
            bar_positions = np.array(map.keys())
            bar_rectangles = plt.bar(bar_positions, bar_heights, width=1.0, linewidth=0)
            plt.xlim(bar_positions.min(), bar_positions.max())
            plt.ylim(0, bar_heights.max() + 1)
            plt.title('Contacts Density Histogram ' + '(cutoff distance = ' + str(distance_cutoff) +
                ur' \u00c5)' +
                '\n\n')
            plt.xlabel('\nResidue Index (continuous)')
            plt.ylabel('Contact Counts\n')
            # Set Colors of the bar
            fractions = bar_heights / bar_heights.max()
            normalized_colors = colors.Normalize(
                fractions.min(), fractions.max())
            count = 0
            for rect in bar_rectangles:
                c = cm.jet(normalized_colors(fractions[count]))
                rect.set_facecolor(c)
                count = count + 1
            ax_text = plt.axes([0.81, -0.02, 0.13, 0.075], frameon=False)
            Button(
                ax_text,
                'Click to select and highlight the residue in PyMol...',
                color=(
                    0.33,
                    0.33,
                    0.33))
            fig.show()
            # This dummy image and drawing it on the canvas is neccessary step
            # to give back CGContextRef to the main window. Otherwise the
            # buttons on main contact map window stops responding.
            progressbar.ax.set_visible(False)
            dummy_img = Image.new(
                'RGBA', (x_img_length, y_img_length), (0, 0, 0, 0))
            glob_ax.imshow(dummy_img)
            glob_ax.figure.canvas.draw()
        if contact_density_closed:
            progressbar.ax.set_visible(True)
            glob_ax.figure.canvas.draw()
            _plot(CDENS)


class InteractionMap:
    '''
    Displays and calculates the interaction map of all residues.
    '''

    def showimap(self, event):
        import numpy as np
        import re
        from matplotlib.widgets import Button

        # aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        aa = [
            'GLU',
            'ASP',
            'LYS',
            'ARG',
            'HIS',
            'GLN',
            'PRO',
            'ASN',
            'ALA',
            'THR',
            'SER',
            'VAL',
            'GLY',
            'MET',
            'CYS',
            'ILE',
            'LEU',
            'TYR',
            'PHE',
            'TRP']

        def _on_close(event):
            global heat_map_closed
            heat_map_closed = True

        def _plot(img, aa):
            global glob_heat_ax
            global heat_map_closed
            mpl.rcParams['figure.subplot.left'] = 0.08
            mpl.rcParams['figure.subplot.right'] = 0.94
            mpl.rcParams['figure.subplot.top'] = 0.91
            mpl.rcParams['figure.subplot.bottom'] = 0.08
            mpl.rcParams['xtick.labelsize'] = 'medium'
            mpl.rcParams['ytick.labelsize'] = 'medium'
            mpl.rcParams['axes.labelsize'] = 'medium'
            mpl.rcParams['ytick.major.pad'] = 8
            mpl.rcParams['xtick.major.pad'] = 8
            mpl.rcParams['axes.grid'] = False
            fig = plt.figure(figsize=(9, 9))
            glob_heat_ax = plt.gca()
            heat_map_closed = False
            fig.canvas.mpl_connect('close_event', _on_close)
            fig.canvas.set_window_title('Aminoacid - Pairwise Heat Map')
            im = plt.imshow(
                img,
                interpolation='nearest',
                cmap='spectral',
                origin='lower')
            plt.title(
                'Aminoacid - Pairwise Heat Map\n' +
                '(cutoff distance = ' +
                str(distance_cutoff) +
                ur' \u00c5)' +
                '\n')
            plt.xticks(range(len(aa)), aa, rotation=90)
            plt.yticks(range(len(aa)), aa)
            mpl.rcParams['ytick.major.pad'] = 2
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(plt.gca())
            cax = divider.append_axes("right", "4%", pad="2%")
            plt.colorbar(im, cax=cax)
            plt.ylabel('Contact Counts')
            callback = ButtonEvents()
            buttaxe = ax = plt.axes([0.79, 0.92, 0.15, 0.04])
            savebutton = Button(
                buttaxe, label='Save Data', hovercolor=(
                    1.0, 1.0, 1.0, 0.2))
            savebutton.active = True
            buttaxe._button = savebutton
            savebutton.on_clicked(callback.saveheatmap)
            fig.show()
            # This dummy image and drawing it on the canvas is neccessary step
            # to give back CGContextRef to the main window. Otherwise the
            # buttons on main contact map window stops responding.
            progressbar.ax.set_visible(False)
            dummy_img = Image.new(
                'RGBA', (x_img_length, y_img_length), (0, 0, 0, 0))
            glob_ax.imshow(dummy_img)
            glob_ax.figure.canvas.draw()
        if heat_map_closed:
            progressbar.ax.set_visible(True)
            glob_ax.figure.canvas.draw()
            _plot(HEMAP, aa)


class Variance:
    '''
    Method for displaying variance map.
    '''

    def show_var_map(self, evt):
        global glob_image
        global variance_map_showing
        if variance_map_showing == 0:
            glob_image = Image.fromarray(VARIA)
            glob_image = glob_image.convert('RGB')
            x_img_length, y_img_length = glob_image.size
            glob_ax.imshow(VARIA, origin='lower', interpolation='nearest')
            visualtext.set_text('Variance Map')
            glob_ax.figure.canvas.draw()
            variance_map_showing = 1
        elif variance_map_showing == 1:
            glob_image = Image.fromarray(COMAP)
            glob_image = glob_image.convert('RGB')
            x_img_length, y_img_length = glob_image.size
            glob_ax.imshow(COMAP, origin='lower', interpolation='nearest')
            visualtext.set_text('Contact Map')
            glob_ax.figure.canvas.draw()
            variance_map_showing = 0


class GUISetup:
    '''
    Setsup GUI and Mouse events to draw rectangular selection and other mouse events.
    '''
    event = None
    image = ''
    selection_count = 1

    def __init__(self):
        self.ax = None
        self.x0 = 0
        self.y0 = 0
        self._x0 = 0
        self._y0 = 0
        self.x1 = 0
        self.y1 = 0
        self.x_ = 0
        self.y_ = 0
        self._x_ = 0
        self._y_ = 0
        self.gray = 0
        self.rect = None
        self._rect = None
        self.bss = None
        self.bcharge = None
        self.bhydrophobic = None
        self.bfactor = None
        self.check = None
        self.heatmap = None
        self.toggle = None
        self.condensmap = None
        self.parameters = None
        self.trackingid = None
        self.new = None
        self.next = None
        self.prev = None
        self.variance = None
        self.about = None
        self.region1 = None
        self.region2 = None

    def _on_close(self, event):
        subprocess.call(["kill", "-9", "%d" % pymol_pid])

    def on_press(self, event):
        if event.inaxes == self.ax:
            self.x0 = event.xdata
            self.y0 = event.ydata
            self._x0 = event.ydata
            self._y0 = event.xdata
            self.trackingid = self.ax.figure.canvas.mpl_connect(
                'motion_notify_event', self.tracking)

    def tracking(self, event):
        if event.inaxes == self.ax:
            self.x_ = event.xdata
            self.y_ = event.ydata
            self._x_ = event.ydata
            self._y_ = event.xdata
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_width(self.x_ - self.x0)
            self.rect.set_height(self.y_ - self.y0)
            self._rect.set_xy((self._x0, self._y0))
            self._rect.set_width(self._x_ - self._x0)
            self._rect.set_height(self._y_ - self._y0)
            self.ax.figure.canvas.draw()

    def on_release(self, event):
        def send_commands(a, b, c, d):
            global distance_cutoff
            connectedSocket.do('set sphere_scale, 0.5')
            region1 = 'select Region_I%s, resi %s-%s and chain %s' % (str(
                self.selection_count),
                residue_index_map[a],
                residue_index_map[b],
                chain_index_map[a])
            region1a = 'select Region_IA, resi %s-%s and chain %s and name %s' % (
                residue_index_map[a], residue_index_map[b], chain_index_map[a], cmap_atom)
            connectedSocket.do(region1)
            connectedSocket.do(region1a)
            connectedSocket.do('show sticks, Region_I%s' %
                               str(self.selection_count))
            connectedSocket.do('color tv_red, Region_IA')
            connectedSocket.do('show sphere, Region_IA')
            # connectedSocket.do('center Region_I, animate=1')
            region2 = 'select Region_II%s, resi %s-%s and chain %s' % (str(
                self.selection_count),
                residue_index_map[c],
                residue_index_map[d],
                chain_index_map[c])
            region2a = 'select Region_IIA, resi %s-%s and chain %s and name %s' % (
                residue_index_map[c], residue_index_map[d], chain_index_map[c], cmap_atom)
            connectedSocket.do(region2)
            connectedSocket.do(region2a)
            connectedSocket.do('show sticks, Region_II%s' %
                               str(self.selection_count))
            connectedSocket.do('color tv_blue, Region_IIA')
            connectedSocket.do('show sphere, Region_IIA')
            # connectedSocket.do('center Region_II, animate=1')
            connectedSocket.do(
                'disable Region_I%s' % str(self.selection_count))
            connectedSocket.do(
                'disable Region_II%s' % str(self.selection_count))
            distcmd = 'distance DIST_' + \
                str(self.selection_count) + ', Region_IA, Region_IIA, cutoff=%d' % (distance_cutoff)
            connectedSocket.do(distcmd)
            #connectedSocket.do('hide labels')
            connectedSocket.do('zoom DIST_' + str(self.selection_count) + ', 5, animate=1')
            connectedSocket.do('delete Region_IA')
            connectedSocket.do('delete Region_IIA')
            groupcmd = 'group SEL_' + str(self.selection_count) + ', Region_I%s Region_II%s DIST_' % (
                str(self.selection_count), str(self.selection_count)) + str(self.selection_count)
            connectedSocket.do(groupcmd)

        if event.inaxes == self.ax:
            self.x1 = event.xdata
            self.y1 = event.ydata
            if (int(self.x0) > int(self.x1)):
                xa1 = int(self.x1)
                xa2 = int(self.x0)
            else:
                xa1 = int(self.x0)
                xa2 = int(self.x1)
            if (int(self.y0) > int(self.y1)):
                xb1 = int(self.y1)
                xb2 = int(self.y0)
            else:
                xb1 = int(self.y0)
                xb2 = int(self.y1)
            # print xa1
            # print xa2
            # print xb1
            # print xb2
            first_selection = []
            second_selection = []
            for a in range(xa1, xa2 + 1):
                for b in range(xb1, xb2 + 1):
                    pix = (a, b)
                    if glob_image.getpixel(pix)[0] > 0:
                        if (a not in first_selection):
                            first_selection.append(a)
                        if (b not in second_selection):
                            second_selection.append(b)

            first_selection.sort()
            second_selection.sort()
            connectedSocket.do("as cartoon")
            connectedSocket.do("show lines")
            connectedSocket.do("color gray40")

            if len(first_selection) == 0 or len(second_selection) == 0:
                self.region1.set_text(
                    'I:%d-%d chain %s' %
                    (residue_index_map[xa1],
                     residue_index_map[xa2],
                     chain_index_map[xa1]))
                self.region2.set_text(
                    'II:%d-%d chain %s' %
                    (residue_index_map[xb1],
                     residue_index_map[xb2],
                     chain_index_map[xb1]))
                send_commands(xa1, xa2, xb1, xb2)
                self.selection_count = self.selection_count + 1
            else:
                self.region1.set_text(
                    'I:%d-%d chain %s' %
                    (residue_index_map[
                        min(first_selection)], residue_index_map[
                        max(first_selection)], chain_index_map[
                        min(first_selection)]))
                self.region2.set_text(
                    'II:%d-%d chain %s' %
                    (residue_index_map[
                        min(second_selection)], residue_index_map[
                        max(second_selection)], chain_index_map[
                        min(second_selection)]))
                send_commands(
                    min(first_selection),
                    max(first_selection),
                    min(second_selection),
                    max(second_selection))
                self.selection_count = self.selection_count + 1
            self.ax.figure.canvas.draw()
            self.ax.figure.canvas.mpl_disconnect(self.trackingid)

    def draw_image(self):
        mpl.rcParams['axes.grid'] = True
        mpl.rcParams['figure.subplot.left'] = 0.0
        mpl.rcParams['figure.subplot.right'] = 1.0
        mpl.rcParams['figure.subplot.top'] = 0.91
        mpl.rcParams['figure.subplot.bottom'] = 0.09
        mpl.rcParams['font.family'] = 'serif'
        mpl.rcParams['font.sans-serif'] = 'Times'

        def _add_buttons(self):
            # Text Placements
            Button(plt.axes([0.855, 0.86, 0.13, 0.075],
                            frameon=False), 'Overlays')
            Button(plt.axes([0.025, 0.86, 0.13, 0.075],
                            frameon=False), 'Plots')
            # Buttons for calculation of overlays and contact density maps and
            # check buttons for aminoacid selections
            self.heatmap = Button(
                ax=plt.axes(
                    [
                        0.025,
                        0.81,
                        0.12,
                        0.065]),
                label='Pairwise\nHeat Map',
                color='#FFFFFF',
                hovercolor=(
                    1.0,
                    1.0,
                    1.0,
                    0.2))
            # Contact Density
            self.condensmap = Button(
                ax=plt.axes(
                    [
                        0.025,
                        0.73,
                        0.12,
                        0.065]),
                label='Contacts\nHistogram',
                color='#FFFFFF',
                hovercolor=(
                    1.0,
                    1.0,
                    1.0,
                    0.2))
            # Variance
            color = '#FFFFFF'
            if pdb_istrajectory:
                pass
            else:
                color = (1.0, 1.0, 1.0, 0.1)
            self.variance = Button(
                ax=plt.axes(
                    [
                        0.025,
                        0.65,
                        0.12,
                        0.065]),
                label='Variance\nContact Map',
                color=color,
                hovercolor=(
                    1.0,
                    1.0,
                    1.0,
                    0.2))
            # Help
            self.parameters = Button(
                ax=plt.axes(
                    [
                        0.025,
                        0.57,
                        0.12,
                        0.065]),
                label='Info',
                color='#FFFFFF',
                hovercolor=(
                    1.0,
                    1.0,
                    1.0,
                    0.2))
            # Legend
            self.about = Button(
                ax=plt.axes(
                    [
                        0.025,
                        0.49,
                        0.12,
                        0.065]),
                label='About',
                color='#FFFFFF',
                hovercolor=(
                    1.0,
                    1.0,
                    1.0,
                    0.2))
            # Next Frame
            color = '#FFFFFF'
            if pdb_istrajectory:
                pass
            else:
                color = (1.0, 1.0, 1.0, 0.1)
            self.next = Button(ax=plt.axes([0.68,
                                            0.035,
                                            0.15,
                                            0.04]),
                               label='Next Frame',
                               color=color,
                               hovercolor=(1.0,
                                           1.0,
                                           1.0,
                                           0.2))
            # Previous Frame
            self.prev = Button(ax=plt.axes([0.17,
                                            0.035,
                                            0.15,
                                            0.04]),
                               label='Previous Frame',
                               color=color,
                               hovercolor=(1.0,
                                           1.0,
                                           1.0,
                                           0.2))
            if pdb_istrajectory:
                pass
            else:
                self.next.active = False
                self.prev.active = False
                self.variance.active = False
            # Secondary Structure button
            color = '#FFFFFF'
            if stridepath == '':
                color = (1.0, 1.0, 1.0, 0.1)
            self.bss = Button(ax=plt.axes([0.86,
                                           0.81,
                                           0.12,
                                           0.065]),
                              label='Secondary\nStructure',
                              color=color,
                              hovercolor=(1.0,
                                          0.0,
                                          0.0,
                                          0.7))
            if stridepath == '':
                self.bss.active = False
            # Charged interactions button
            self.bcharge = Button(plt.axes([0.86,
                                            0.73,
                                            0.12,
                                            0.065]),
                                  'Charged\nInteractions',
                                  color='#FFFFFF',
                                  hovercolor=(0.0,
                                              0.0,
                                              1.0,
                                              0.6))
            # Hydrophobic interactions button
            self.bhydrophobic = Button(
                plt.axes(
                    [
                        0.86,
                        0.65,
                        0.12,
                        0.065]),
                'Hydrophobic\nInteractions',
                color='#FFFFFF',
                hovercolor=(
                    1.0,
                    1.0,
                    0.0,
                    0.6))
            # Bfactor button
            self.bfactor = Button(plt.axes([0.86,
                                            0.57,
                                            0.12,
                                            0.065]),
                                  'B-factor',
                                  color='#FFFFFF',
                                  hovercolor=(0.0,
                                              1.0,
                                              1.0,
                                              0.6))
            # B-factor cutoff slider
            self.bfactor_slider = Slider(
                plt.axes(
                    [
                        0.86,
                        0.525,
                        0.095,
                        0.03],
                    axisbg='lightgoldenrodyellow'),
                '',
                0,
                100,
                valinit=bfac_cutoff,
                valfmt=" %d")
            # Progress Bar Slider
            global progressbar
            progressbar = Button(plt.axes(
                [0.82, 0.053, 0.2, 0.02], axisbg='white', frameon=False), 'Calculating...')
            progressbar.label.set_color((1.0, 0.0, 0.0, 1.0))
            progressbar.ax.set_visible(False)
            # Amino acids checkbuttons
            self.check = CheckButtons(
                plt.axes(
                    [
                        0.855,
                        0.080,
                        0.14,
                        0.43]),
                ('ALAnine',
                 'ARGinine',
                 'ASparagiNe',
                 'ASPartic',
                 'CYSteine',
                 'GLutaminNe',
                 'GLUtamic',
                 'GLYcine',
                 'HIStidine',
                 'IsoLEucine',
                 'LEUcine',
                 'LYSsine',
                 'METhionine',
                 'PHEnylalanine',
                 'PROline',
                 'SERine',
                 'THReonine',
                 'TRyPtophan',
                 'TYRosine',
                 'VALine'),
                (False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False))
            self.check.ax.axes.set_axis_bgcolor((0, 0, 0, 0))
            self.check.ax.axes.set_frame_on(False)
            for rect in self.check.rectangles:
                rect.set_width(0.1)
                rect.set_fill(False)
                rect.set_lw(1.25)
            for text in self.check.labels:
                text.set_ha('left')
            for line in self.check.lines:
                a = (line[0].get_xdata()[0], line[0].get_xdata()[1] + 0.065)
                line[0].set_xdata(a)
                a = (line[1].get_xdata()[0], line[1].get_xdata()[1] + 0.065)
                line[1].set_xdata(a)
        # Initialize Interaction Map calculation class
        imapclass = InteractionMap()
        # Initialize button callback events class
        bcallback = ButtonEvents()
        # Initialize Contact Density Plot class
        condensclass = ContactDensity()
        # Initialize Trajectory Class
        trajClass = Trajectory()
        # Initialize Variance Class
        varClass = Variance()
        # Help Class
        helpclass = Help()
        global glob_image
        global glob_ax
        global x_img_length
        global y_img_length
        global calculate_cmap

        # Main figure window
        fig = plt.figure()
        fig.canvas.set_window_title("CMPyMOL")
        fig.canvas.mpl_connect('close_event', self._on_close)

        if calculate_cmap == 0:
            import numpy as np
            imgarray = np.fromfile(image_filepath, dtype=float, sep=' ')
            imgarray = np.subtract(255, imgarray)
            imgarray = np.subtract(imgarray, imgarray.min())
            imgarray = np.divide(imgarray, imgarray.max())
            imgarray = np.multiply(imgarray, 255)
            imgarray = np.reshape(imgarray, (308, 308))
            img = Image.fromarray(imgarray)
            img = img.transpose(Image.FLIP_TOP_BOTTOM)
            plt.imshow(img, origin='lower', interpolation='none')
            glob_image = img
        else:
            plt.imshow(COMAP, origin='lower', interpolation='none')
            glob_image = Image.fromarray(COMAP)
        glob_image = glob_image.convert('RGB')
        x_img_length, y_img_length = glob_image.size
        plt.xlabel('\nResidue Index')
        plt.ylabel('')
        plt.set_cmap('gray')
        self.ax = plt.gca()
        glob_ax = self.ax
        global frametext, visualtext
        frametext = self.ax.annotate('Frame #1/%s' % str(len(model_indexes) + 1),
                                     xy=(2.0,
                                         1),
                                     xycoords='data',
                                     xytext=(-112,
                                             200),
                                     textcoords='offset points',
                                     bbox=dict(boxstyle="round",
                                               fc="0.8"))
        visualtext = self.ax.annotate('Contact Map',
                                      xy=(2.0,
                                          1),
                                      xycoords='data',
                                      xytext=(-112,
                                              225),
                                      textcoords='offset points',
                                      bbox=dict(boxstyle="round",
                                                fc="0.8"))
        self.region1 = self.ax.annotate('Region I',
                                        xy=(2.0,
                                            1),
                                        xycoords='data',
                                        xytext=(-112,
                                                175),
                                        textcoords='offset points',
                                        bbox=dict(boxstyle="round",
                                                  fc="0.8"))
        self.region2 = self.ax.annotate('Region II',
                                        xy=(2.0,
                                            1),
                                        xycoords='data',
                                        xytext=(-112,
                                                150),
                                        textcoords='offset points',
                                        bbox=dict(boxstyle="round",
                                                  fc="0.8"))
        self.ax.fmt_xdata = self.ax.fmt_ydata = lambda coord: '{:d}'.format(
            coord)
        self.ax.set_autoscaley_on(False)
        self.ax.set_autoscalex_on(False)
        xticks = self.ax.xaxis.get_major_ticks()
        for x in xticks:
            x.label1.set_visible(False)
            x.label2On = True
            x.label2.set_rotation('vertical')
        yticks = self.ax.yaxis.get_major_ticks()
        for y in yticks:
            y.label1.set_visible(False)
        # labels for xaxis
        locs, labels = plt.xticks()
        t = []
        for i in range(len(locs)):
            t.append(str(int(locs[i] + 1)))
        plt.xticks(locs[1:(len(locs) - 1)], t[1:(len(t) - 1)])
        _add_buttons(self)
        # Rectangle for selection added to the top corner of the plot window
        self.rect = Rectangle(
            (0, 0), width=0, height=0, alpha=0.6, fc=(
                0, 0, 0, 0), ec='magenta', aa=True, lw=2)
        self._rect = Rectangle(
            (0, 0), width=0, height=0, alpha=0.6, fc=(
                0, 0, 0, 0), ec='magenta', aa=True, lw=2)
        self.ax.add_patch(self.rect)
        self.ax.add_patch(self._rect)
        # Button callbacks
        self.bss.on_clicked(bcallback.ss)
        self.bcharge.on_clicked(bcallback.charge)
        self.bhydrophobic.on_clicked(bcallback.hydrophobic)
        self.bfactor.on_clicked(bcallback.bfactor)
        self.bfactor_slider.on_changed(bcallback.bfactor_slider)
        self.check.on_clicked(bcallback.check)
        self.heatmap.on_clicked(imapclass.showimap)
        self.condensmap.on_clicked(condensclass.showcdensmap)
        self.parameters.on_clicked(helpclass.showparameters)
        self.about.on_clicked(helpclass.showabout)
        self.next.on_clicked(trajClass.loadNextFrame)
        self.prev.on_clicked(trajClass.loadPrevFrame)
        self.variance.on_clicked(varClass.show_var_map)
        # Mouse events on the plot window registered
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        plt.show()


def loadPDBinitial(model_index):
    import time
    # Send pdb to pymol
    _row = []
    while len(_row) < 1:
        db_cursor.execute("SELECT * FROM data WHERE id=?", (model_index,))
        _row = db_cursor.fetchall()
        time.sleep(0.1)
    global FRAME, COMAP, HEMAP, CDENS, VARIA
    NAME = _row[0][1]
    FRAME = _row[0][2]
    if (len(residue_index_map) == 0 or len(chain_index_map) == 0):
        pdb = PDBfunctions()
        pdb.map_residues(FRAME)
    COMAP = pickle.loads(str(_row[0][3]))
    HEMAP = pickle.loads(str(_row[0][4]))
    CDENS = pickle.loads(str(_row[0][5]))
    VARIA = pickle.loads(str(_row[0][7]))
    connectedSocket.do("reinitialize")
    connectedSocket.do("load " + FRAME + "," +
                       os.path.basename(pdb_path)[:-4] + "," + str(model_index + 1))
    connectedSocket.do("frame " + str(model_index + 1))
    connectedSocket.do("as cartoon")
    connectedSocket.do("show lines")
    connectedSocket.do("color gray40")
    connectedSocket.do("set seq_view, on")
    connectedSocket.do("set seq_view_format, 1")
    connectedSocket.do('set seq_view_label_spacing, 2')
    connectedSocket.do('set dash_color, limegreen')
    connectedSocket.do('set label_position, [-2,0,0]')
    connectedSocket.do('set dash_gap, 0')
    connectedSocket.do('set dash_width, 0.7')
    connectedSocket.do('set dash_length, 0.3')
    connectedSocket.do('set dash_use_shader, 1')
    connectedSocket.do('set label_color, yellow')


def loadPDB(model_index):
    import time
    import numpy as np
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    # Send pdb to pymol
    _row = []
    while len(_row) < 1:
        db_cursor.execute("SELECT * FROM data WHERE id=?", (model_index,))
        _row = db_cursor.fetchall()
        time.sleep(0.1)
    global FRAME, COMAP, HEMAP, CDENS, VARIA
    NAME = _row[0][1]
    FRAME = _row[0][2]
    if (len(residue_index_map) == 0 or len(chain_index_map) == 0):
        pdb = PDBfunctions()
        pdb.map_residues(FRAME)
    COMAP = pickle.loads(str(_row[0][3]))
    HEMAP = pickle.loads(str(_row[0][4]))
    CDENS = pickle.loads(str(_row[0][5]))
    VARIA = pickle.loads(str(_row[0][7]))
    frametext.set_text("Frame #%s/%s" %
                       ((model_index + 1), str(len(model_indexes) + 1)))
    connectedSocket.do("load " + FRAME + "," +
                       os.path.basename(pdb_path)[:-4] + "," + str(model_index + 1))
    connectedSocket.do("frame " + str(model_index + 1))
    global glob_image
    if variance_map_showing == 0:
        glob_image = Image.fromarray(COMAP)
        glob_ax.imshow(COMAP, origin='lower', interpolation='nearest')
        visualtext.set_text('Contact Map')
    elif variance_map_showing == 1:
        glob_image = Image.fromarray(VARIA)
        glob_ax.imshow(VARIA, origin='lower', interpolation='nearest')
        visualtext.set_text('Variance Map')
    glob_image = glob_image.convert('RGB')
    x_img_length, y_img_length = glob_image.size
    glob_ax.figure.canvas.draw()
    if not heat_map_closed:
        glob_heat_ax.imshow(
            HEMAP,
            interpolation='nearest',
            cmap='spectral',
            origin='lower')
        glob_heat_ax.figure.canvas.draw()
    if not contact_density_closed:
        bar_heights = np.array(CDENS.values()).astype(float)
        bar_positions = np.array(CDENS.keys())
        bar_rectangles = glob_con_ax.bar(
            bar_positions, bar_heights, width=1.0, linewidth=0)
        fractions = bar_heights / bar_heights.max()
        normalized_colors = colors.normalize(fractions.min(), fractions.max())
        count = 0
        for rect in bar_rectangles:
            c = cm.jet(normalized_colors(fractions[count]))
            rect.set_facecolor(c)
            count = count + 1
        glob_con_ax.figure.canvas.draw()
    overlays = ButtonEvents()
    global ss_showing
    if ss_showing == 1:
        global ss_calculated
        global glob_ss_rects
        ss_calculated = 0
        for r in glob_ss_rects:
            r.remove()
        glob_ss_rects = []
        overlays.ss('Dummy')
        overlays.ss('Dummy')
    global charge_showing
    if charge_showing == 1:
        global charge_calculated
        charge_calculated = 0
        overlays.charge('Dummy')
        overlays.charge('Dummy')
    global hp_showing
    if hp_showing == 1:
        global hp_calculated
        hp_calculated = 0
        overlays.hydrophobic('Dummy')
        overlays.hydrophobic('Dummy')
    overlays.check('Dummy')


def generate_threadedcontactMap():
    import multiprocessing
    import time
    pdb = PDBfunctions()

    def threading_function():
        global processed_frames
        for index in range(multiprocessing.cpu_count()):
            processed_frames = processed_frames + 1
            p = multiprocessing.Process(
                target=pdb.generate_contact_map, args=(
                    pdb_path, processed_frames - 1,))
            p.start()
            time.sleep(1)

    if current_model_index == 0:
        global processed_frames
        db_cursor.execute("SELECT MAX(id) FROM data")
        _row = db_cursor.fetchall()
        if (_row[0][0] is None):
            processed_frames = 0
        else:
            processed_frames = _row[0][0] + 1

    if pdb_istrajectory:
        if current_model_index + multiprocessing.cpu_count() + 1 >= processed_frames:
            threading_function()
    else:
        if current_model_index > processed_frames - 1:
            pdb.generate_contact_map(pdb_path)


class Trajectory:
    '''
    Handles all calculations related to trajectories
    '''

    def loadNextFrame(self, event):
        global current_model_index
        if current_model_index == len(model_indexes) - 1:
            pass
        else:
            current_model_index = current_model_index + 1
            loadPDB(current_model_index)
            import threading
            t = threading.Thread(target=generate_threadedcontactMap())
            t.setDaemon(1)
            t.start()

    def loadPrevFrame(self, event):
        global current_model_index
        if current_model_index == 0:
            pass
        else:
            current_model_index = current_model_index - 1
            loadPDB(current_model_index)

# Main environment
if __name__ == "__main__":
    if (sys.platform.startswith('darwin') or sys.platform.startswith("linux")):
        init = Initialize()
        init.connectPyMOL()
        init.start()
        init.getPDB()
        init.get_model_indexes()

        import threading
        print ">>> Generating Contact Map... (this will take a while)."
        t = threading.Thread(target=generate_threadedcontactMap())
        t.setDaemon(1)
        t.start()
        print ">>> Trying to load PDB..."
        loadPDBinitial(current_model_index)
        gui = GUISetup()
        gui.draw_image()
    else:
        print "*** ERROR: Windoze operating system is not currently supported. ***"
        sys.exit()
