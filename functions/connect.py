import os
import re
import sys
import xmlrpclib
import subprocess
from PyQt4 import QtGui

class connect():
    """Functions to find location of and make connection to PyMOL
    and handle errors related to the connection"""

    def __init__(self, parent):
        # Create CMPyMOL directory
        home = os.path.expanduser('~')
        print home
        self.working_directory = os.path.join(home, '.CMPyMOLDir')
        if not os.path.exists(os.path.join(home, '.CMPyMOLDir')):
            os.makedirs(os.path.join(home, '.CMPyMOLDir'))
        self.pymol_path = ''
        self.pymol_pid = 0
        self.connected_socket = xmlrpclib.Server('http://localhost:9123')
        self.parent = parent

    def get_pymol_path(self):
        try:
            input = open(os.path.join(self.working_directory, 'pymol_path'))
            self.pymol_path = input.readline().rstrip()
            if len(self.pymol_path) > 1 and os.path.exists(self.pymol_path):
                return
            else:
                pass
        except IOError:
            pass

        file_filter = '*.app'
        location = '/Applications'
        if sys.platform == 'win32':
            file_filter = '*.exe'
            location = os.path.join('C:\Program Files')
        elif sys.platform == 'linux':
            file_filter = '*.*'
            location = os.path.join('~')

        temp_path = str(QtGui.QFileDialog.getOpenFileName(QtGui.QFileDialog(),
                                                          caption="Locate PyMOL or MacPyMOL executable",
                                                          directory=os.path.expanduser(location),
                                                          filter="Executables (%s)" % file_filter))

        if sys.platform == 'darwin' and os.path.basename(temp_path) == 'MacPyMOL.app':
            self.pymol_path = os.path.join(temp_path, 'Contents', 'MacOS', 'MacPyMOL')
        elif sys.platform == 'linux':
            self.pymol_path = temp_path
        elif sys.platform == 'win32':
            self.pymol_path = re.sub(r'/', r'\\', temp_path)
        self.write_pymol_path()

    def write_pymol_path(self):
        output = open(os.path.join(self.working_directory, 'pymol_path'), 'w')
        output.write('%s' % self.pymol_path)
        output.close()

    def make_pymol_connection(self):
        pyp = subprocess.Popen([self.pymol_path, '-R'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.pymol_pid = pyp.pid
        self.parent.append_status_text(">>> Launching PyMOL from location %s (pid: %d)" % (self.pymol_path,
                                                                                           self.pymol_pid))
