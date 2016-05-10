import os
import sys
import re
import subprocess

class Overlays():
    def __init__(self, parent):
        self.parent = parent

    def calculate_secondary_structure(self):
        if sys.platform == 'darwin':
            stride_path = os.path.join(os.curdir, 'tools', 'stride_osx')
        elif sys.platform == 'win32':
            stride_path = os.path.join(os.curdir, 'tools', 'stride_win32.exe')
        else:
            self.parent.ss_btn.setEnabled(False)
            return
        process = subprocess.Popen([stride_path, self.parent.protein_map.pdb_path], stdout=subprocess.PIPE)
        ss = []
        while True:
            line = process.stdout.readline()
            if not line:
                break
            if re.match(r'ASG ', line):
                ss.append(line.split()[5])
        self.parent.protein_map.secondary_structure = ss

