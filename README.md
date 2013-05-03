CMPyMOL
=======
CMPyMOL is an add-on software to molecular visualization program PyMOL. It combines the 3D visualization capabilities of PyMOL and 2D protein contact maps with an interactive interface for scientific analysis. This software is freely distributed under the MIT license for Linux and Mac OS X platforms.

Author
======
Venkatramanan Krishnamani (venks@andrew.cmu.edu)

Website
======
http://emptyewer.github.io/CMPyMOL/

Prerequisites
=============
1. Python 2.7
2. Python library dependencies (using pip or easy_install)
 
 a) matplotlib
 ```Shell
 pip install matplotlib
 ```
 b) python imaging library (PIL)
  ```Shell
 pip install pil
 ```
 c) numpy (numeric python)
  ```Shell
 pip install numpy
 ```
3. PyMOL. (It is recommended that the user make sure the directory of installed executable is listed in the $PATH variable.)
4. Stride secondary structure assignment tool. This program can be downloaded from http://webclu.bio.wzw.tum.de/stride/ and easily compiled into an executable. It is recommemded that the executable or its directory is listed in the system variable $PATH. NOTE: If this executable is not installed properly, the secondary structure calculation will be disabled in CMPyMOL.

Mac OS X
--------
Users can install the python libraries using "easy_install" or "pip". It is recommened that the user downloads Enthough Canopy python distribution and management package from https://www.enthought.com/products/canopy/. This package includes a robust python library management software and a python IDE.

PyMOL 1.5.x can be installed using MacPorts http://www.macports.org. NOTE: This automatically adds the executable into $PATH.

Linux
-----
The python dependencies and PyMOL can be installed using apt-get (aptitude) or similar package management system.


Installation
============
There is no need for installation of the script. Optionally, a standalone executable can be complied using "pyinstaller" or "py2exe" or "py2app" packages depending on the users operating system.

Usage
=====
```Shell
python CMPyMOL.py
```

This command will automatically invoke the PyMOL executable and the user is led through the rest of the program with a series of pop-up windows.

Software
========
Clicking (left) and draging a selection of points on the displayed contact map will highlight the corresponding residues involved in that contact in the PyMOL window (as red and blue colored spheres). In addition, several structural/biochemical properties can be overlayed on top of the contact map (listed below). The contact-map data can also be plotted in other representations.

Overlays
--------
1. Secondary structure of the protein as translucent strips over the contact map. This button won't be active if secondary structure calculation program stride is not installed. (Button: Secondary Structure)
2. Contact points where a Charge-Charge interaction occurs. (Button: Charged Interactions)
3. Residues that interact via hydrophobic interaction. (Button: Hydrophobic Interactions)
4. Contact regions that have a b-factor that is higher than a certain cutoff (Button: B-factor). The b-factor cutoff can be varied using a slider (Slider).
5. Highlights a contact point/region where the selected pair of residues are in contact (selected by checking the checkbox). If only one aminoacid is selected from the list, interaction site of the selected aminoacid with another one of the same type is highlighted. (List of checkboxes for each aminoacid)

Plots
-----
1. Pairwise Heat Map - Plots a 20x20 matrix of count of pairwise aminoacid interactions.
2. Contacts Histogram - Plots the number of contacts formed by a given residue as a bar graph. Selecting a particular bar highlights the corresponding residue in the PyMOL window.

Requests and Disclaimer
=======================
Users are welcome to send me an email to request the addition of a particular feature. I am actively developing this software.

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/496c5edce682fd47dca759c644857cea "githalytics.com")](http://githalytics.com/emptyewer/CMPyMOL)
