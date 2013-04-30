CMPyMOL
=======
CMPyMOL combines the 3D visualization capabilities of PyMol and protein contact map into a rich and intuitive tool for scientific analysis. This software is freely distributed under the MIT license for Linux and Mac OS X platforms.

Website
=======
http://pymolwiki.org/index.php/CMPyMOL

Author
======
Venkatramanan Krishnamani (venks@andrew.cmu.edu)

Prerequisites
=============
1. Python 2.7
2. Python framework dependency
	a. matplotlib
	b. python imaging library (PIL)
	c. numpy
	d. biopython
3. PyMOL executable. It is recommended that the user make sure the directory of installed executable is listed in the $PATH variable.
4. Stride executable. This can be downloaded from http://webclu.bio.wzw.tum.de/stride/ and easily compiled into an executable. It is recommemded that the executable or its directory is listed in the system variable $PATH. NOTE: If this executable is not installed properly, the secondary structure calculation will be automatically disabled.

Mac OS X
--------
Users can install the python dependencies easily using "easy_install" or "pip". It is recommened that the user downloads Enthough Canopy python distribution and management package from https://www.enthought.com/products/canopy/. It comes with a robust package management software and a python IDE.

PyMOL 1.5.x can be installed using MacPorts http://www.macports.org. This automatically adds the executable into $PATH.

Linux
-----
The python dependencies and PyMOL can be installed using apt-get (aptitude) or similar package management system.


Installation
============
There is no installation of the script is required required. Optionally a standalone executable can be complied using "pyinstaller" or "py2exe" or "py2app" packages depending on the users operating system.

Usage
=====
"python CMPyMOL.py"

The above command will automatically invoke PyMOL. The user is led through the rest of the program with a series of pop-up window.

Software
========
Left-clicking and draging anywhere on the contact map displayed will highlight the corresponding residues in PyMOL (spheres). In addition, there are several different structural/biochemical properties overlays that can be overlayed on top of the contact map. The user can plot the contact-map data in other representations.

Overlays
--------
1. Secondary Structure of the protein as strips of color (Secondary Structure). This button won't be active if secondary structure calculation program stride is not installed.
2. Contact regions where a Charge-Charge interaction occurs (Charged Interactions)
3. Residues that interact via a hydrophobic interaction (Hydrophobic Interactions)
4. Contact regions that have a b-factor that is higher than a certain cutoff (B-factor). The cutoff can be varied using a slider (slider below the button).
5. Highlights a contact point/region where a pair of residues are in contact. If only one aminoacid is selected from the list, interaction site of the selected aminoacid with another on of the same type is highlighted.

Plots
-----
1. Pairwise Heat Map - Plots a 20x20 matrix of number of pairwise aminoacid interactions as a colored heatmap.
2. Contacts Histogram - Plots the number of residue wise contacts of a given protein.

Requests and Disclaimer
=======================
Users are welcome to send me an email (above) to request the addition of a particular feature. I am actively developing this software. Please do understand that I have other committments, so the requested feature may not be able to implemented promptly.


