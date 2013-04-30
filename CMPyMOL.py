#!/opt/local/bin/python


'''
CMapper 1.0

http://pymolwiki.org/index.php/CMapper

Author: Venkatramanan Krishnamani (Version 1.0)

'''

#  The MIT License (MIT)
# =======================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you are
# free to use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# -----------------------------------------------------------------------------------
# CMapper
# Copyright (C) 2013 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
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

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.patches import Rectangle
    from matplotlib.widgets import Button
    from matplotlib.widgets import Slider
    from matplotlib.widgets import CheckButtons
except ImportError:
    print '*** Error: This plugin requires the "matplotlib" module ***'

try:
    import PIL.Image as Image
    from operator import itemgetter
    from itertools import groupby
except ImportError:
    print "*** Error: Cannot Import operator or itertools. ***"
    sys.exit()

#### Customizing Plot Parameters
mpl.rcParams['figure.figsize'] = [11,9]
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['grid.color'] = 'white'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['xtick.major.pad'] = 10
mpl.rcParams['ytick.major.pad'] = 10
mpl.rcParams['axes.titlesize'] = 'x-large'
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['axes.labelweight'] = 'bold'
#### Initializing Global Variables
# Contact Map Image
x_img_length = 0                # (Int) Width of Contact Map 
y_img_length = 0                # (Int) Height of Contact Map
glob_ax = None                  # (matplotlib axes) Contact Map Plot
glob_image = None               # (PIL.Image) Contact Map Image
dist_matrix = []                # (Matrix) Inverse Distance Map
# Secondary Structure
ss_calculated = 0               # (Boolean) If secondary Structure calculation is previously performed
ss_showing = 0                  # (Boolean) If Secondary Structure Overlay is Visible
glob_ss_rects = []              # (matplotlib.widgets.Rectangle) Describes the Secondary Structure Overlay
# Charged Interactions
charge_img_inst = None          # (matplotlib.image.AxesImage) Charged Interactions Overlay
charge_calculated = 0           # (Boolean) If Charges Interactions Overlay is previously calculated
charge_showing = 0              # (Boolean) If Charges Interactions Overlay is Visible
# Hydrophobic Interactions
hp_img_inst = None              # (matplotlib.image.AxesImage) Hydrophobic Interactions Overlay
hp_calculated = 0               # (Boolean) If Hydrophobic Interactions is previously calculated
hp_showing = 0                  # (Boolean) If Hydrophobic Interactions is Visible
# B-factor
bfac_cutoff = 25                # (Int) Default cutoff for B-factor calculation
bfac_img_inst = None            # (matplotlib.image.AxesImage) B-factor Overlay
bfac_calculated = 0             # (Boolean) If B-factor Overlay is previously calculted
bfac_showing = 0                # (Boolean) If B-factor Overlay is showing
bfac_sliderax = None            # (matplotlib axes) B-factor Slider Axes
# Interacting Aminoacids Overlay
buttons_checked = []            # (List) List of aminoacids chosen from checkedbuttons
aa_calculated = 0               # (Boolean) If Aminoacid overlay is previously calculated
aa_img_inst = None              # (matplotlib.image.AxesImage) Aminoacids Overlay       
# Heatmap

heatmap_img = None              # (PIL.Image) Aminoacid-Aminoacid Contacts Image
heatmap_calculated = 0          # (Boolean) Aminoacid-Aminoacid Contacts previously calculated
# Contact Density Map
condensmap = {}                 # (Dictionary) Contact Density Count
condensmap_calculated = 0       # (Boolean) If Contact Density Map Calculated
protein_sequence = ''
progressbar = None
#
residue_index_map = {}
chain_index_map = {}
#PATHS
stridepath = ''
pymolpath = ''
pdbpath = ''
distance_cutoff = 0.0             # (Float) Distance Cutoff for Contact Map Calculation
connectedSocket = ''

class Connect:
    '''

    Opens a XMLRPC socket and connects to PyMOL so information can be passed to it.

    '''
    import xmlrpclib
    global connectedSocket
    connectedSocket = xmlrpclib.Server('http://localhost:9123')

class InitialChecks:
    '''
    Locates the path of programs
    '''
    import sys
    import os
    global connectedSocket
    
    def cmd_exists(cmd):
        import subprocess
        proc = subprocess.Popen(["which", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        status = False
        command = proc.stdout.read().strip()
        if command != '':
            status = True
        return (status, command)


    if (sys.platform == 'darwin' or sys.platform == "linux"):
        import subprocess
        import os
        import wx

        app = wx.App()
        
        global stridepath
        global pymolpath
        global pdbpath
        global distance_cutoff

        stride_loc = cmd_exists("stride")
        if stride_loc[0]:
            stridepath = stride_loc[1]
            print ">>> Stride located at", stridepath
        else:
            print "*** Error: Cannot locate STRIDE path. Please install from http://webclu.bio.wzw.tum.de/stride/."
            print "*** Note: Secondary structure calculations will be disabled."

        pymol_loc = cmd_exists("pymol")
        if pymol_loc[0]:
            pymolpath = pymol_loc[1]
            subprocess.Popen([pymolpath, '-R'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print ">>> Lanching PyMOL from location", pymolpath
        else:
            print "*** Error: Cannot locate PyMOL in the system path."
            if sys.platform == 'darwin':
                print "*** Please install PyMOL using macports http://www.macports.org."
            elif sys.platform == 'linux':
                print "*** Please install PyMOL using apt-get install pymol."

        pdbwildcard = "PDB (*.pdb)|*.pdb|"
        dlg1 = wx.FileDialog(None, "Choose the PDB file...", os.getcwd(), "", wildcard=pdbwildcard)
        try:
            if dlg1.ShowModal() == wx.ID_OK:
                pdbpath = dlg1.GetPath()
                print ">>> PDB File:", pdbpath
                connectedSocket.do("reinitialize")
                connectedSocket.do("load " + pdbpath)
                connectedSocket.do("as cartoon")
                connectedSocket.do("color white")
                connectedSocket.do("set seq_view, on")
                connectedSocket.do("set seq_view_format, 1")
                connectedSocket.do('set seq_view_label_spacing, 2')
            else:
                print "*** Error: No PDB file was chosen. Exiting..."
                sys.exit()
        finally:
            dlg1.Destroy()
            wx.YieldIfNeeded()

        dlg2 = wx.TextEntryDialog(None, "Maximum distance between a pair of C-alpha atoms that defines a contact.","Maximum distance cutoff", "12.0")
        try:
            if dlg2.ShowModal() == wx.ID_OK:
                distance_cutoff = float(dlg2.GetValue())
                print ">>> Maximum distance between C-alpha atoms that define a contact : %s" % distance_cutoff
            else:
                print "*** Warning: No distance cutoff was set. Using default value of 10.0 angstrom."
                print ">>> Maximum distance between C-alpha atoms that define a contact : %s" % distance_cutoff
                distance_cutoff = 12.0
        finally:
            dlg2.Destroy()
            wx.YieldIfNeeded()
    else:
        print "*** ERROR: This operating system is not supported."
        sys.exit()

class ButtonEvents:
    '''
    Calculates and draws various overlays on top of the contact map.
    '''
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

        def _get_rectangle_x(r,c):
            xs = r[0]
            ys = 0
            xe = r[1]
            ye = x_img_length
            rec = Rectangle((xs,ys),width=(xe-xs), height=(ye-ys), alpha=0.30, facecolor=c, edgecolor=c)
            return rec
        def _get_rectangle_y(r,c):
            xs = 0
            ys = r[0]
            xe = y_img_length
            ye = r[1]
            rec = Rectangle((xs,ys),width=(xe-xs), height=(ye-ys), alpha=0.30, facecolor=c, edgecolor=c)
            return rec

        def _get_secondary_structure(pdb):
            import subprocess
            import re

            process = subprocess.Popen([stridepath, pdb], stdout=subprocess.PIPE)
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
            for k, g in groupby(enumerate(data), lambda (i,x):i-x):
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
            ss = _get_secondary_structure(pdbpath)
            alpha,beta = _get_alpha_beta(ss)
            alpha_range = _get_ranges(alpha)
            beta_range = _get_ranges(beta)  
            for a in alpha_range:
                rect = _get_rectangle_x(a,'red')
                glob_ss_rects.append(rect)
                glob_ax.add_patch(rect)
                rect = _get_rectangle_y(a,'red')
                glob_ss_rects.append(rect)
                glob_ax.add_patch(rect)
            for b in beta_range:
                rect = _get_rectangle_x(b,'green')
                glob_ss_rects.append(rect)
                glob_ax.add_patch(rect)
                rect = _get_rectangle_y(b,'green')
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
    
    def charge(self,event):
        '''
        Calculates and generates an overlay for charge-charge interactions. Callback event for charged interactions button.
        '''
        charged_aa = ['ARG','LYS','ASP','GLU']
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
            process = subprocess.Popen([stridepath, pdb], stdout=subprocess.PIPE)
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
            import PIL.ImageDraw as ImageDraw
            charges_location = _get_charged_sequence(pdbpath)
            charged_img = Image.new('RGBA',(x_img_length,y_img_length),(0,0,0,0))
            charged_img_draw = ImageDraw.Draw(charged_img)
            row = 0
            for r in charges_location:
                col = 0
                for c in charges_location:
                    pix = (row, col)
                    if (r == 1 and c == 1) and row != col and glob_image.getpixel(pix)[0] >= 0.8:
                        #charged_img.putpixel(pix, (0, 0, 255, 255))
                        bbox = (pix[0]-1,pix[1]-1,pix[0]+1,pix[1]+1)
                        charged_img_draw.ellipse(bbox,fill=(0, 0, 255, 235))
                    col = col + 1
                row = row + 1
            if mpl.__version__ != '1.2.0':
                charged_img = charged_img.transpose(Image.FLIP_TOP_BOTTOM)
            charge_img_inst = glob_ax.imshow(charged_img, interpolation='nearest', origin='lower')
            charge_calculated = 1
            
        if charge_showing == 0:
            _charge_show()
            charge_showing = 1
        else:
            _charge_hide()
            charge_showing = 0

        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()
        
    def hydrophobic(self,event):
        '''
        Calculated and generates an image overlay for interacting hydrophobic aminoacids. Callback event for Hydrophobic Interactions button.
        '''
        hp_aa = ['VAL','ILE','LEU','MET', 'PHE', 'TRP']
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
            process = subprocess.Popen([stridepath, pdb], stdout=subprocess.PIPE)
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
            import PIL.ImageDraw as ImageDraw
            hp_location = _get_hp_sequence(pdbpath)
            hp_img = Image.new('RGBA',(x_img_length,y_img_length),(0,0,0,0))
            hp_img_draw = ImageDraw.Draw(hp_img)
            row = 0
            for r in hp_location:
                col = 0
                for c in hp_location:
                    pix = (row, col)
                    if (r == 1 and c == 1) and row != col and glob_image.getpixel(pix)[0] >= 0.8:
                        #hp_img.putpixel(pix, (255, 255, 0, 255))
                        bbox = (pix[0]-1,pix[1]-1,pix[0]+1,pix[1]+1)
                        hp_img_draw.ellipse(bbox,fill=(255, 255, 0, 235))
                    col = col + 1
                row = row + 1
            if mpl.__version__ != '1.2.0':
                hp_img = hp_img.transpose(Image.FLIP_TOP_BOTTOM)
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

    def bfactor(self,event):
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
        if (type(event) is int and bfac_calculated == 1 and bfac_showing == 1) or (bfac_calculated == 0 and type(event) is not int and bfac_showing == 0):
            import PIL.ImageDraw as ImageDraw
            pdb_handle = open( pdbpath, "r" )
            for line in pdb_handle:
                if (re.search(r'ATOM', line[0:5])):
                    if re.search(r'CA', line[12:15]):
                        bfacs.append(int(float(line[60:65])))
            bfac_img = Image.new('RGBA',(x_img_length,y_img_length),(0,0,0,0))
            bfac_img_draw = ImageDraw.Draw(bfac_img)
            row = 0
            for r in bfacs:
                col = 0
                for c in bfacs:
                    pix = (row, col)
                    if (r > bfac_cutoff and c > bfac_cutoff) and row != col and glob_image.getpixel(pix)[0] >= 0.8:
                        bbox = (pix[0]-1,pix[1]-1,pix[0]+1,pix[1]+1)
                        bfac_img_draw.ellipse(bbox,fill=(0, 255, 255, 235))
                    col = col + 1
                row = row + 1
            if (bfac_calculated == 1):
                bfac_img_inst.set_visible(False)
            bfac_calculated = 1
            if mpl.__version__ != '1.2.0':
                bfac_img = bfac_img.transpose(Image.FLIP_TOP_BOTTOM)
            bfac_img_inst = glob_ax.imshow(bfac_img, interpolation='nearest')
            bfac_showing = 1
        elif (bfac_calculated == 1 and type(event) is not int and bfac_showing == 1):
            bfac_img_inst.set_visible(False)
            bfac_showing = 0
        elif (bfac_calculated == 1 and type(event) is not int and bfac_showing == 0):
            bfac_img_inst.set_visible(True)
            bfac_showing = 1

        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()
                
    def bfactor_slider(self,event):
        '''
        Callback event when the value of bfactor cutoff slider changes
        '''
        global bfac_cutoff
        bfac_cutoff = int(event)
        self.bfactor(bfac_cutoff)
        
    def check(self,event):
        '''
        Aminoacid-Aminoacid interactions overlay callback function.
        '''
        import re
        global buttons_checked
        global aa_img_inst
        global progressbar
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()

        def _get_aa_sequence(pdb,residues):
            '''
            Retrieve amino acid sequence.
            '''
            import subprocess
            process = subprocess.Popen([stridepath, pdb], stdout=subprocess.PIPE)
            aa_seq = []
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                if re.match(r'ASG ', line):
                    aa_seq.append(line.split()[1])
            return aa_seq
        aa = re.sub(r'([a-z])','',event)
        if (aa in buttons_checked):
            buttons_checked.remove(aa)
        else:
            buttons_checked.append(aa)
        if aa_img_inst != None:
            aa_img_inst.set_visible(False)
        import PIL.ImageDraw as ImageDraw
        aa_locations = _get_aa_sequence(pdbpath, buttons_checked)
        aa_img = Image.new('RGBA',(x_img_length,y_img_length),(0,0,0,0))
        aa_img_draw = ImageDraw.Draw(aa_img)
        row = 0
        for r in aa_locations:
            col = 0
            for c in aa_locations:
                pix = (row, col)
                if len(buttons_checked) == 2:
                    if ((r == buttons_checked[0] and c == buttons_checked[1]) or (r == buttons_checked[1] and c == buttons_checked[0])) and row != col and glob_image.getpixel(pix)[0] >= 0.8:
                        bbox = (pix[0]-1,pix[1]-1,pix[0]+1,pix[1]+1)
                        aa_img_draw.ellipse(bbox,fill=(255, 132, 0, 255))
                elif len(buttons_checked) == 1:
                    if (r == buttons_checked[0] and c == buttons_checked[0]) and row != col and glob_image.getpixel(pix)[0] >= 0.8:
                        bbox = (pix[0]-1,pix[1]-1,pix[0]+1,pix[1]+1)
                        aa_img_draw.ellipse(bbox,fill=(255, 132, 0, 255))
                col = col + 1
            row = row + 1
        if mpl.__version__ != '1.2.0':
            aa_img = aa_img.transpose(Image.FLIP_TOP_BOTTOM)
        aa_img_inst = glob_ax.imshow(aa_img, interpolation='nearest')
        
        progressbar.ax.set_visible(False)
        glob_ax.figure.canvas.draw()

class Help:
    '''
    Calculates and displays hydrogen bonding partners network on top of contact map.
    '''
    def showhelp(self,event):
        import webbrowser
        handle = webbrowser.get()
        handle.open('https://github.com/VenkyKrishnamani/CMapper')

class Load:
    '''
    Loads a new pdb into Pymol
    '''
    def new(self, event):
        global dist_matrix
        global connectedSocket
        global pdbpath
        global distance_cutoff
        import time
        import PIL.Image as Image
        import numpy as np
        import glob
                
        if (sys.platform == 'darwin' or sys.platform == "linux"):
            import subprocess
            import os
            import wx

            app = wx.App()

            pdbwildcard = "PDB (*.pdb)|*.pdb|"
            dlg1 = wx.FileDialog(None, "Choose the PDB file...", os.getcwd(), "", wildcard=pdbwildcard)
            try:
                if dlg1.ShowModal() == wx.ID_OK:
                    pdbpath = dlg1.GetPath()
                    print ">>> PDB File:", pdbpath
                    connectedSocket.do("reinitialize")
                    connectedSocket.do("load " + pdbpath)
                    connectedSocket.do("as cartoon")
                    connectedSocket.do("color white")
                    connectedSocket.do("set seq_view, on")
                    connectedSocket.do("set seq_view_format, 1")
                    connectedSocket.do('set seq_view_label_spacing, 2')
                else:
                    print "*** Error: No PDB file was chosen. Exiting..."
                    sys.exit()
            finally:
                dlg1.Destroy()
                wx.YieldIfNeeded()

            dlg2 = wx.TextEntryDialog(None, "Maximum distance between a pair of C-alpha atoms that defines a contact.","Maximum distance cutoff", "12.0")
            try:
                if dlg2.ShowModal() == wx.ID_OK:
                    distance_cutoff = float(dlg2.GetValue())
                    print ">>> Maximum distance between C-alpha atoms that define a contact : %s" % distance_cutoff
                else:
                    print "*** Warning: No distance cutoff was set. Using default value of 10.0 angstrom."
                    print ">>> Maximum distance between C-alpha atoms that define a contact : %s" % distance_cutoff
                    distance_cutoff = 12.0
            finally:
                dlg2.Destroy()
                wx.YieldIfNeeded()
        else:
            print "*** ERROR: This operating system is not supported."
            sys.exit()

        mouse = MouseMonitor()
        pdb = PDBfunctions()
        print ">>> Generating Contact Map..."

        dist_matrix = pdb.generate_contact_map(pdbpath)
        pdb.map_residues(pdbpath)

        if not os.path.exists(os.path.join(os.environ['HOME'],'.CMapperDir')):
            os.makedirs(os.path.join(os.environ['HOME'],'.CMapperDir'))
        filelist = glob.glob(os.path.join(os.environ['HOME'],'.CMapperDir','_temp_img_*'))
        for f in filelist:
            os.remove(f)
        temp_file = "_temp_img_%.7f.png" % time.time()
        image_file = os.path.join(os.environ['HOME'],'.CMapperDir',temp_file)
        #rescaled_matrix = (dist_matrix/((dist_matrix.max() - dist_matrix.min())/255.0)).astype(np.uint8)
        rescaled_matrix = ((dist_matrix - dist_matrix.min())*(255.0/dist_matrix.max())).astype(np.uint8)
        im = Image.fromarray(rescaled_matrix)
        im.save(image_file)
        mouse.set_data(image_file)

class ContactDensity:
    '''
    Calculates and displays a new plot of contact density along the sequence.
    '''

    global connectedSocket

    def showcdensmap(self,event):
        global condensmap
        global protein_sequence
        global glob_ax
        global condensmap_calculated
        global progressbar
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()
        
        axis = None
        aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

        import numpy as np
        import matplotlib.cm as cm
        import matplotlib.colors as colors
        try:
            import Bio.PDB
            from Bio.PDB.Polypeptide import PPBuilder
        except Exception, err:
            print "Please install a required python package: Biopython"
            sys.stderr.write('>>> IMPORT ERROR: %s\n' % str(err))
            sys.exit()

        def _calc_dist(p1,p2):
            return np.linalg.norm(p1-p2)

        def _onpick(event):
            global axis
            global connectedSocket
            if event.inaxes == axis:
                global connectedSocket
                global residue_index_map
                global chain_index_map
                connectedSocket.do("as cartoon")
                connectedSocket.do("color white")
                connectedSocket.do('select Selection, resi %s and chain %s' % (residue_index_map[int(event.xdata)], chain_index_map[int(event.xdata)]))
                connectedSocket.do('disable Selection')
                connectedSocket.do('show spheres, Selection')
                connectedSocket.do('color tv_green, Selection')
                connectedSocket.do('center Selection, animate=1')
        def _plot(map):
            global axis
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
            axis = plt.gca()
            fig.canvas.set_window_title('Contacts Histogram')
            fig.canvas.mpl_connect('button_press_event', _onpick)
            bar_heights = np.array(map.values()).astype(float)
            bar_positions = np.array(map.keys())
            bar_rectangles = plt.bar(bar_positions, bar_heights, width=1.0, linewidth=0)
            plt.xlim(bar_positions.min(), bar_positions.max())
            plt.ylim(0, bar_heights.max()+1)
            plt.title('Contacts Histogram ' + '(cutoff distance = ' + str(distance_cutoff) + ur' \u00c5)' + '\n\n')
            plt.xlabel('\nResidue Index')
            plt.ylabel('Count\n')
            # Set Colors of the bar
            fractions = bar_heights/bar_heights.max()
            normalized_colors = colors.normalize(fractions.min(), fractions.max())
            count = 0
            for rect in bar_rectangles:
                c = cm.jet(normalized_colors(fractions[count]))
                rect.set_facecolor(c)
                count = count + 1
            ax_text = plt.axes([0.81, -0.02, 0.13, 0.075], frameon=False)
            Button(ax_text, 'Click to select and highlight the residue in PyMol...', color=(0.33,0.33,0.33))
            fig.show()
            # This dummy image and drawing it on the canvas is neccessary step to give back CGContextRef to the main window. Otherwise the buttons on main contact map window stops responding.
            progressbar.ax.set_visible(False)
            dummy_img = Image.new('RGBA',(x_img_length,y_img_length),(0,0,0,0))
            glob_ax.imshow(dummy_img)
            glob_ax.figure.canvas.draw()
        
        if condensmap_calculated == 0:
            parser=Bio.PDB.PDBParser(PERMISSIVE=True, QUIET=True)
            structure=parser.get_structure('Protein', pdbpath)
            ppb=PPBuilder()
            
            for pp in ppb.build_peptides(structure):
                chain_sequence = pp.get_sequence()
                protein_sequence = protein_sequence + chain_sequence
            protein_sequence = list(protein_sequence)

            for n in range(len(protein_sequence)):
                condensmap[n] = 0
            
            print ">>> Calculating Contacts Distribution Map..."
            row = 0
            for r1 in structure.get_residues():
                if (r1.resname in aa):
                    col = 0
                    for r2 in structure.get_residues():
                        if (r2.resname in aa):
                            if (col != row):
                                dist = _calc_dist(r1.child_dict['CA'].get_coord(),r2.child_dict['CA'].get_coord())
                                if dist <= distance_cutoff:
                                    condensmap[row] = condensmap[row] + 1
                            col = col + 1
                    row = row + 1
            condensmap_calculated = 1
        _plot(condensmap)

class InteractionMap:
    '''
    Displays and calculates the interaction map of all residues.
    '''
    
    def showimap(self,event):
        global heatmap_img
        global heatmap_calculated
        global progressbar
        global glob_ax
        progressbar.ax.set_visible(True)
        glob_ax.figure.canvas.draw()

        import numpy as np
        try:
            import Bio.PDB
        except Exception, err:
            print "Please install a required python package: Biopython"
            sys.stderr.write('>>> IMPORT ERROR: %s\n' % str(err))
            sys.exit()
        
        aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

        def _calc_dist(p1,p2):
            return np.linalg.norm(p1-p2)
        
        def _plot(img,aa):
            mpl.rcParams['figure.subplot.left'] = 0.08
            mpl.rcParams['figure.subplot.right'] = 0.94
            mpl.rcParams['figure.subplot.top'] = 0.91
            mpl.rcParams['figure.subplot.bottom'] = 0.08
            mpl.rcParams['xtick.labelsize'] = 'small'
            mpl.rcParams['ytick.labelsize'] = 'small'
            mpl.rcParams['axes.labelsize'] = 'medium'
            mpl.rcParams['ytick.major.pad'] = 8
            mpl.rcParams['xtick.major.pad'] = 8
            mpl.rcParams['axes.grid'] = False
            fig = plt.figure(figsize=(9, 9))
            fig.canvas.set_window_title('Aminoacid - Pairwise Heat Map')
            im = plt.imshow(img, interpolation='nearest', cmap='spectral', origin='lower')
            plt.title('Aminoacid - Pairwise Heat Map\n' + '(cutoff distance = ' + str(distance_cutoff) + ur' \u00c5)' + '\n')
            plt.xticks(range(len(aa)), aa, rotation=90)
            plt.yticks(range(len(aa)), aa)
            mpl.rcParams['ytick.major.pad'] = 2
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(plt.gca())
            cax = divider.append_axes("right", "4%", pad="2%")
            plt.colorbar(im, cax=cax)
            plt.ylabel('Count')
            fig.show()
            # This dummy image and drawing it on the canvas is neccessary step to give back CGContextRef to the main window. Otherwise the buttons on main contact map window stops responding.
            progressbar.ax.set_visible(False)
            dummy_img = Image.new('RGBA',(x_img_length,y_img_length),(0,0,0,0))
            glob_ax.imshow(dummy_img)
            glob_ax.figure.canvas.draw()
            
        if heatmap_calculated == 0:
            parser=Bio.PDB.PDBParser(PERMISSIVE=True, QUIET=True)
            structure=parser.get_structure('Protein', pdbpath)
            heatmap = {}
            #Initializing
            for i in aa:
                hm = {}
                for j in aa:
                    hm[j] = 0.0
                heatmap[i] = hm
            print ">>> Calculating Pairwise Heat Map..."
            heatmap_img = np.zeros(len(aa)*len(aa))
            row = 0
            for r1 in structure.get_residues():
                if r1.resname in aa:
                    col = 0
                    for r2 in structure.get_residues():
                        if r2.resname in aa:
                            if (col > row):
                                dist = _calc_dist(r1.child_dict['CA'].get_coord(),r2.child_dict['CA'].get_coord())
                                if dist <= distance_cutoff:
                                    heatmap[r1.resname][r2.resname] = heatmap[r1.resname][r2.resname] + 1
                                    pix = aa.index(r1.resname)*len(aa) + aa.index(r2.resname)
                                    heatmap_img[pix] = heatmap[r1.resname][r2.resname]
                            col = col + 1
                    row = row + 1
            #heatmap_img = heatmap_img/max(heatmap_img)
            heatmap_img = heatmap_img.reshape((len(aa), len(aa)))
            heatmap_calculated = 1

        _plot(heatmap_img,aa)
                
class MouseMonitor:
    '''
    Records Mouse events and draws rectangular selection and other mouse events.
    '''

    event = None
    image = ''
    xdatalist = []
    ydatalist = []
    read_image = ''
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
        self.width = None
        self.height = None
        self.bss = None
        self.bcharge = None
        self.bhydrophobic = None
        self.bfactor = None
        self.check = None
        self.heatmap = None
        self.toggle = None
        self.condensmap = None
        self.help = None
        self.trackingid = None
        self.new = None
    
    def on_press(self, event):
        if event.inaxes == self.ax:
            self.x0 = event.xdata
            self.y0 = event.ydata
            self._x0 = event.ydata
            self._y0 = event.xdata
            self.trackingid = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.tracking)

    def set_plt_cmap(self, event):
        if self.gray == 1:
            plt.set_cmap('gist_ncar')
            self.gray = 0
        else:
            plt.set_cmap('gray')
            self.gray = 1
        self.ax.figure.canvas.draw()

    def on_release(self, event):
        global connectedSocket
        global residue_index_map
        global chain_index_map
        def send_commands(a, b, c, d):
            region1 =  'select Region_I, resi %s-%s and chain %s' % (residue_index_map[a], residue_index_map[b], chain_index_map[a])
            connectedSocket.do(region1)
            connectedSocket.do('show spheres, Region_I')
            connectedSocket.do('color tv_red, Region_I')
            connectedSocket.do('center Region_I, animate=1')
            region2 =  'select Region_II, resi %s-%s and chain %s' % (residue_index_map[c], residue_index_map[d], chain_index_map[c])
            connectedSocket.do(region2)
            connectedSocket.do('show spheres, Region_II')
            connectedSocket.do('color tv_blue, Region_II')
            connectedSocket.do('center Region_II, animate=1')
            connectedSocket.do('disable Region_I')
            connectedSocket.do('disable Region_II')
            groupcmd = 'group SEL_' + str(self.selection_count) + ', Region_I Region_II'
            connectedSocket.do(groupcmd)

        if event.inaxes == self.ax:
            self.x1 = event.xdata
            self.y1 = event.ydata
            if (int(self.x0) > int(self.x1)):
                xa1 = int(self.x1+1)
                xa2 = int(self.x0+1)
            else:
                xa1 = int(self.x0+1)
                xa2 = int(self.x1+1)
            if (int(self.y0) > int(self.y1)):
                xb1 = int(self.y1+1)
                xb2 = int(self.y0+1)
            else:
                xb1 = int(self.y0+1)
                xb2 = int(self.y1+1)
            first_selection = []
            second_selection = []
            for a in range(xa1,xa2+1):
                for b in range(xb1,xb2+1):
                    pix = (a, b)
                    if glob_image.getpixel(pix)[0] > 0:
                        if (a not in first_selection):
                            first_selection.append(a)
                        if (b not in second_selection):
                            second_selection.append(b)
            
            first_selection.sort()
            second_selection.sort()
            connectedSocket.do("as cartoon")
            connectedSocket.do("color white")

            if len(first_selection) == 0 or len(first_selection) == 0:
                send_commands(xa1, xa2, xb1, xb2)
                self.selection_count = self.selection_count + 1
            else:
                send_commands(min(first_selection), max(first_selection), min(second_selection), max(second_selection))
                self.selection_count = self.selection_count + 1
            self.ax.figure.canvas.mpl_disconnect(self.trackingid)
            

    def tracking(self,event):
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

    def set_data(self,im):
        global x_img_length
        global y_img_length
        global glob_image
        
        self.image = im
        self.read_image = mpimg.imread(self.image)
        glob_image = Image.open(self.image)
        glob_image = glob_image.convert('RGB')
        self.xdatalist = range(self.read_image[1,:].size)
        self.ydatalist = range(self.read_image[:,1].size)
        self.width = self.read_image[1,:].size
        self.height = self.read_image[:,1].size
        x_img_length = self.width
        y_img_length = self.height
        self.draw_image()

    def draw_image(self):
        mpl.rcParams['axes.grid'] = True
        global glob_ax
        global bfac_sliderax
        def _add_buttons(self):
            # Text Placements
            ax_text = plt.axes([0.855, 0.86, 0.13, 0.075], frameon=False)
            Button(ax_text, 'Overlays')
            ax_text = plt.axes([0.73, -0.02, 0.13, 0.075], frameon=False)
            Button(ax_text, 'Click and drag to select and highlight a region in PyMol...', color=(0.33,0.33,0.33))
            ax_text = plt.axes([0.025, 0.86, 0.13, 0.075], frameon=False)
            Button(ax_text, 'Plots')
            ax_text = plt.axes([0.025, 0.51, 0.13, 0.075], frameon=False)
            Button(ax_text, '--------------------')
            ax_text = plt.axes([0.86, 0.47, 0.12, 0.065], frameon=False)
            Button(ax_text, '-----------------')
            # Buttons for Overlay Legends
            b = Button(ax=plt.axes([0.23, 0.047, 0.073, 0.03]), label=ur'\u03b1'+'-Helix', color=(1.0,0.0,0.0,0.7))
            b.label.set_fontsize(11.0)
            b = Button(ax=plt.axes([0.31, 0.047, 0.073, 0.03]), label=ur'\u03b2'+'-Sheet', color=(0.0,1.0,0.0,0.7))
            b.label.set_fontsize(11.0)
            b = Button(ax=plt.axes([0.39, 0.047, 0.073, 0.03]), label='Charges', color=(0.0,0.0,1.0,0.6))
            b.label.set_fontsize(11.0)
            b = Button(ax=plt.axes([0.47, 0.047, 0.073, 0.03]), label='Non-Polar', color=(1.0,1.0,0.0,0.6))
            b.label.set_fontsize(11.0)
            b = Button(ax=plt.axes([0.55, 0.047, 0.073, 0.03]), label='B-factor', color=(0.0,1.0,1.0,0.6))
            b.label.set_fontsize(11.0)
            b = Button(ax=plt.axes([0.63, 0.047, 0.073, 0.03]), label='Custom', color=(1.0,0.52,0.0,1.0))
            b.label.set_fontsize(11.0)
            b = Button(ax=plt.axes([0.71, 0.047, 0.073, 0.03]), label='Selection', color=(1.0,0.0,1.0,0.6))
            b.label.set_fontsize(11.0)
            # Buttons for calculation of overlays and contact density maps and check buttons for aminoacid selections
            ax_heatmap = plt.axes([0.025, 0.81, 0.12, 0.065])
            self.heatmap = Button(ax=ax_heatmap, label='Pairwise\nHeat Map', color='#FFFFFF', hovercolor=(1.0,1.0,1.0,0.2))
            # Contact Density
            ax_condensmap = plt.axes([0.025, 0.73, 0.12, 0.065])
            self.condensmap = Button(ax=ax_condensmap, label='Contacts\nHistogram', color='#FFFFFF', hovercolor=(1.0,1.0,1.0,0.2))
            # Toggle of colormap
            ax_toggle = plt.axes([0.025, 0.65, 0.12, 0.065])
            self.toggle = Button(ax=ax_toggle, label='Toggle\nMap Coloring', color='#FFFFFF', hovercolor=(1.0,1.0,1.0,0.2))
            # Help
            ax_help = plt.axes([0.025, 0.57, 0.12, 0.065])
            self.help = Button(ax=ax_help, label='Help', color='#FFFFFF', hovercolor=(1.0,1.0,1.0,0.2))
            # Load Another
            ax_new = plt.axes([0.025, 0.09, 0.12, 0.065])
            self.new = Button(ax=ax_new, label='Load New', color='#FFFFFF', hovercolor=(1.0,1.0,1.0,0.2))
            #Secondary Structure button
            ax_ss = plt.axes([0.86, 0.81, 0.12, 0.065])
            self.bss = Button(ax=ax_ss, label='Secondary\nStructure', color='#FFFFFF', hovercolor=(1.0,0.0,0.0,0.7))
            # Charged interactions button
            ax_charge = plt.axes([0.86, 0.73, 0.12, 0.065])
            self.bcharge = Button(ax_charge, 'Charged\nInteractions', color='#FFFFFF', hovercolor=(0.0,0.0,1.0,0.6))
            # Hydrophobic interactions button
            ax_hydrophobic = plt.axes([0.86, 0.65, 0.12, 0.065])
            self.bhydrophobic = Button(ax_hydrophobic, 'Hydrophobic\nInteractions', color='#FFFFFF', hovercolor=(1.0,1.0,0.0,0.6))
            # Bfactor button
            ax_bfactor = plt.axes([0.86, 0.57, 0.12, 0.065])
            self.bfactor = Button(ax_bfactor, 'B-factor', color='#FFFFFF', hovercolor=(0.0,1.0,1.0,0.6))
            # B-factor cutoff slider
            ax_bfactor_slider = plt.axes([0.86, 0.525, 0.095, 0.03], axisbg='white')
            self.bfactor_slider = Slider(ax_bfactor_slider, '', 0, 100, valinit=bfac_cutoff, valfmt=" %d")
            # Progress Bar Slider
            global progressbar
            ax_pbar = plt.axes([0.82, 0.053, 0.2, 0.02], axisbg='white', frameon = False)
            progressbar = Button(ax_pbar, 'Calculating...')
            progressbar.label.set_color((1.0,0.0,0.0,1.0))
            progressbar.ax.set_visible(False)
            # Amino acids checkbuttons
            ax_aa = plt.axes([0.855, 0.080, 0.14, 0.43])
            self.check = CheckButtons(ax_aa,('ALAnine','ARGinine','ASparagiNe','ASPartic','CYSteine','GLutaminNe','GLUtamic','GLYcine','HIStidine','IsoLEucine','LEUcine','LYSsine','METhionine','PHEnylalanine','PROline','SERine','THReonine','TRyPtophan','TYRosine','VALine'), (False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False))
            self.check.ax.axes.set_axis_bgcolor((0,0,0,0))
            self.check.ax.axes.set_frame_on(False)
            for rect in self.check.rectangles:
                rect.set_width(0.1)
                rect.set_fill(False)
                rect.set_lw(1.25)
            for text in self.check.labels:
                text.set_ha('left')
            for line in self.check.lines:
                a = (line[0].get_xdata()[0],line[0].get_xdata()[1]+0.065)
                line[0].set_xdata(a)
                a = (line[1].get_xdata()[0],line[1].get_xdata()[1]+0.065)
                line[1].set_xdata(a)
        #Initialize Interaction Map calculation class
        imapclass = InteractionMap()
        # Initialize button callback events class
        bcallback = ButtonEvents()
        # Initialize Contact Density Plot class
        condensclass = ContactDensity()
        # Help Class
        helpclass = Help()
        load = Load()
        # Main figure window
        mpl.rcParams['figure.subplot.left'] = 0.0
        mpl.rcParams['figure.subplot.right'] = 1.0
        mpl.rcParams['figure.subplot.top'] = 0.91
        mpl.rcParams['figure.subplot.bottom'] = 0.09
        plt.close("all")
        fig = plt.figure()
        fig.canvas.set_window_title(self.image)
        plt.imshow(self.read_image, origin='lower', interpolation='nearest')    
        plt.title('Contact Map ' + '(cutoff distance = ' + str(distance_cutoff) + ur' \u00c5)' + '\n')
        plt.xlabel('\nResidue Number')
        plt.ylabel('\nResidue Number')
        self.ax = plt.gca()
        self.ax.fmt_xdata = self.ax.fmt_ydata = lambda coord: '{:d}'.format(coord)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.set_plt_cmap(self)
        glob_ax = self.ax
        self.ax.set_autoscaley_on(False)
        self.ax.set_autoscalex_on(False)
        # Add buttons to the plot window
        _add_buttons(self)
        # Button callbacks
        self.bss.on_clicked(bcallback.ss)
        self.bcharge.on_clicked(bcallback.charge)
        self.bhydrophobic.on_clicked(bcallback.hydrophobic)
        self.bfactor.on_clicked(bcallback.bfactor)
        self.bfactor_slider.on_changed(bcallback.bfactor_slider)
        self.check.on_clicked(bcallback.check)
        self.heatmap.on_clicked(imapclass.showimap)
        self.toggle.on_clicked(self.set_plt_cmap)
        self.condensmap.on_clicked(condensclass.showcdensmap)
        self.help.on_clicked(helpclass.showhelp)
        self.new.on_clicked(load.new)
        # Rectangle for selection added to the top corner of the plot window
        self.rect = Rectangle((0,0), width=0, height=0, alpha=0.6, fc=(0,0,0,0), ec='magenta', aa=True, lw=2)
        self._rect = Rectangle((0,0), width=0, height=0, alpha=0.6, fc=(0,0,0,0), ec='magenta', aa=True, lw=2)
        self.ax.add_patch(self.rect)
        self.ax.add_patch(self._rect)
        # Mouse events on the plot window registered
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        #self.trackingid = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.tracking)
        plt.show()

class PDBfunctions:
    '''
    Calculates Contact Map and other PDB file operations.
    '''
    def unique(seq, idfun=None): 
       # order preserving
       if idfun is None:
           def idfun(x): return x
       seen = {}
       result = []
       for item in seq:
           marker = idfun(item)
           if marker in seen: continue
           seen[marker] = 1
           result.append(item)
       return sort(result)
        
    def generate_contact_map(self,pdb):
        global distance_cutoff
        import numpy as np
        try:
            import Bio.PDB
        except Exception, err:
            print "Please install a required python package: Biopython"
            sys.stderr.write('>>> IMPORT ERROR: %s\n' % str(err))
            sys.exit()
        
        def _calc_dist(p1,p2):
            return np.linalg.norm(p1-p2)

        parser=Bio.PDB.PDBParser(PERMISSIVE=True, QUIET=True)
        structure=parser.get_structure('Protein', pdb)
        count = 0
        for atom in structure.get_atoms():
            if 'CA' in atom.name:
                count = count + 1
        matrix = np.zeros((count+1, count+1), np.float)
        aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        row = 0
        for r1 in structure.get_residues():
            if r1.resname in aa:
                col = 0
                for r2 in structure.get_residues():
                    if r2.resname in aa:
                        dist = _calc_dist(r1.child_dict['CA'].get_coord(),r2.child_dict['CA'].get_coord())
                        if (dist <= distance_cutoff):
                            #matrix[row, col] = dist
                            if dist == 0:
                                dist = distance_cutoff
                            matrix[row, col] = 1/dist
                        col = col + 1
                    else:
                        continue
                row = row + 1
            else:
                continue
        return matrix
    
    def map_residues(self,pdb):
        global residue_index_map
        global chain_index_map
        pFile = open( pdb, "r" )
        index = 1
        for line in pFile:
            line = line.rstrip()
            resid = line[22:26].strip()
            chainid = line[21:22].strip()
            if (line[12:16].strip() == 'CA'):
                chain_index_map[index] = chainid
                residue_index_map[index] = resid
                index = index  + 1

### Main function that initiates the calculation.
def _contact_map_visualizer():
    global dist_matrix
    import time
    import PIL.Image as Image
    import numpy as np
    import glob

    Connect()
    InitialChecks()

    mouse = MouseMonitor()
    pdb = PDBfunctions()
    print ">>> Generating Contact Map..."

    global pdbpath
    dist_matrix = pdb.generate_contact_map(pdbpath)
    pdb.map_residues(pdbpath)

    if not os.path.exists(os.path.join(os.environ['HOME'],'.CMapperDir')):
        os.makedirs(os.path.join(os.environ['HOME'],'.CMapperDir'))
    filelist = glob.glob(os.path.join(os.environ['HOME'],'.CMapperDir','_temp_img_*'))
    for f in filelist:
        os.remove(f)
    temp_file = "_temp_img_%.7f.png" % time.time()
    image_file = os.path.join(os.environ['HOME'],'.CMapperDir',temp_file)
    #rescaled_matrix = (dist_matrix/((dist_matrix.max() - dist_matrix.min())/255.0)).astype(np.uint8)
    rescaled_matrix = ((dist_matrix - dist_matrix.min())*(255.0/dist_matrix.max())).astype(np.uint8)
    im = Image.fromarray(rescaled_matrix)
    im.save(image_file)
    mouse.set_data(image_file)

### Main environment
if __name__ == "__main__":  
    _contact_map_visualizer()

