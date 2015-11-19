CMPyMOL
=======
CMPyMOL is an add-on software to molecular visualization program PyMOL. It combines the protein 3D visualization capabilities of PyMOL and the protein's 2D contact map with an interactive interface for scientific analysis. This software is freely distributed under the MIT license for Linux and Mac OS X platforms.

Author
======
Venkatramanan Krishnamani (venky.krishna@me.com)

Website
======
http://emptyewer.github.io/CMPyMOL/

Version History
===============

1.5
---
1. Added support for reading multi-model PDB files. (Multi-Model PDB file for NMR structure or Trajectory from MD simulations.)
2. Supports displaying variance of contact points from a series of contact maps generated from multi-model PDB files.
3. CMPyMOL stores the calculated contact maps, heat maps and contact density information in a local SQLite database for fast and easy subsequent access.
4. Cleaner GUI.
5. Parallelized the code for contact map calculation for multi-frame PDB files.


Install Instructions for Linux
==============================

The python dependencies and PyMOL can be installed using ```apt-get``` (aptitude) or ```pip``` or a similar package management system.


0. PyMOL
 ```Shell
 apt-get install pymol
 ```
1. Python 2.7 (provided by default in the distribution)
2. Python external library dependencies (install using pip or easy_install)

 a) wxpython
 ```Shell
 apt-get install python-wxgtk2.8
 ```
 b) matplotlib
 ```Shell
 pip install matplotlib
 ```
 c) python imaging library (PIL)
  ```Shell
 pip install pil
 ```
 d) numpy (numeric python)
  ```Shell
 pip install numpy
 ```
3. Stride secondary structure assignment tool. This program can be downloaded from http://webclu.bio.wzw.tum.de/stride/ and compiled into a stand-alone executable. It is recommended that the Stride executable or its installation directory is added to the $PATH environment variable. NOTE: If this executable is not detected in the $PATH variable, the secondary structure calculation will be disabled in CMPyMOL.

Install Instructions for Mac OS X
=================================

Users can install the python libraries using ```easy_install``` or ```pip``` or ```homebrew```. It is recommended that the user use Enthought Canopy python distribution and management package downloaded from https://www.enthought.com/products/canopy/. This package includes a robust python library management software and a python IDE.

0. PyMOL
 ```Shell
 brew install homebrew/science/pymol
 ```
1. Python 2.7 (provided by default in the distribution)
2. Python external library dependencies (install using ```pip``` or ```easy_install``` or ```brew```)

 a) wxpython
 ```Shell
 brew install wxmac
 ```
 ```Shell
 brew install wxpython
 ```
 b) matplotlib
 ```Shell
 pip install matplotlib
 ```
 c) python imaging library (PIL)
  ```Shell
 pip install pil
 ```
 d) numpy (numeric python)
  ```Shell
 pip install numpy
 ```
3. Stride secondary structure assignment tool. This program can be downloaded from http://webclu.bio.wzw.tum.de/stride/ and compiled into a stand-alone executable. It is recommended that the Stride executable or its installation directory is added to the $PATH environment variable. NOTE: If this executable is not detected in the $PATH variable, the secondary structure calculation will be disabled in CMPyMOL.


Usage
=====
```Shell
python /<path to CMPyMOL directory>/CMPyMOL.py
```

This command will automatically invoke the PyMOL executable and the user is led through the rest of the program with a series of pop-up windows.

Software
========
Clicking (left) and dragging a selection of contact points on the displayed contact map will highlight the corresponding residues in the PyMOL window (as red and blue colored atoms in spheres representation). In addition, several structural/biochemical properties can be overlaid on top of the contact map. The contact-map data can also be plotted in other representations. The calculated contact-map, heat-map and contact density information is stored in a local SQL database. Any subsequent access of the same PDB with matching parameters will be read from the database for fast access. The code for calculating contact map for trajectory files is parallelized for efficiency.

Input
------
1. Single-frame PDB files (local)
2. Multi-frame PDB trajectory files (local)
3. Multi-frame trajectory files should have the following format.

```Shell
MODEL X
.
.
.
ATOM ...
ATOM ...
.
.
.
ENDMDL
```
NOTE: The PDB can include REMARKS, CRYST and other standard PDB information entries. The MODEL line is essential for the software to work properly (ENDMDL is optional).

Overall Interface
-----------------
![alt tag](https://raw.githubusercontent.com/emptyewer/CMPyMOL/master/images/cmpymol.png)

Overlays
--------

![alt tag](https://raw.githubusercontent.com/emptyewer/CMPyMOL/master/images/main.png)

1. Secondary structure of the protein is overlaid as translucent strips over the contact map. This button won't be active if secondary structure calculation program stride is not found in the system path ($PATH). (Button: Secondary Structure)

2. Contact points where a Charge-Charge interaction occurs are highlighted. (Button: Charged Interactions)

3. Residues that interact via hydrophobic interaction are highlighted. (Button: Hydrophobic Interactions)

4. Contact regions that have a B-factor that is higher than a certain cutoff are highlighted (Button: B-factor). The b-factor cutoff can be varied using a slider (Slider).

5. Highlights a contact point/region where the pair of selected residues are in contact (selected by checking the checkboxes). If only one aminoacid is selected from the list, interaction site of the selected aminoacid with another one of the same type is highlighted. (List of checkboxes for each aminoacid)

Plots
-----
1. Pairwise Heat Map - Plots a 20x20 matrix of pairwise aminoacid interaction count.
![alt tag](https://raw.githubusercontent.com/emptyewer/CMPyMOL/master/images/pairwise.png)
2. Contacts Histogram - Plots the number of contacts around a given residue. Selecting a particular bar highlights the corresponding residue in the PyMOL window.
![alt tag](https://raw.githubusercontent.com/emptyewer/CMPyMOL/master/images/contact-density.png)
3. Variance Contact Map - For Multi-frame PDB files (trajectory), this button toggles the displays the variance contact map starting from the initial frame until the current frame. This view can be used to identifying the dynamic regions in a protein.

Word of Caution
===============
When using a multi-frame PDB file, the contact-map for the next frame(s) are being pre-calculated in the background (depending on the number of free CPU cores available). Clicking on "Next Frame" in rapid succession may lead to undesired results and/or a crash.

In the event of a crash, delete the database that is created in the working directory and relaunch the program.

Requests
========
Users are welcome to send me an email to request the addition of a specific feature or to report a bug.

