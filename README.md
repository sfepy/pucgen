PUCGen - Periodic Unit Cell Generator
=====================================

A python script for genereting periodic unit cells.

Requirements
------------

* [Gmsh](http://gmsh.info/) - three-dimensional finite element mesh generator
* [FreeCAD](https://www.freecadweb.org) - 3D CAD/CAE parametric modeling application
* [PyQt5](https://riverbankcomputing.com/software/pyqt/intro) - bindings for the Qt application framework (including OpenGL module)
* [Python bindings for VTK](https://vtk.org/download)

On Ubuntu based Linux distributions users can use the following command to install the required packages:

    apt install gmsh freecad python-pyqt5 python-pyqt5.qtopengl python-vtk6

Installation
------------

* Download the code from the git repository:

      git clone git://github.com/vlukes/mumpspy

or

* Use [pip](https://pypi.org/project/pip/):

      pip install git+git://github.com/vlukes/mumpspy

Usage
-----

* Command line:

      python pucgen.py <input_file.txt>

* GUI:

      python pucgen_gui.py

Examples:
---------

Sample input `example1.puc`:
```
BaseCell;size=(1, 1, 1);el_size=0.1;mat_id=1
SphericalInclusion;radius=0.3;central_point=(0, 0, 0);es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=2
CylindricalChannel;radius=0.1;central_point=(0, 0, 0);direction=x;es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=2
CylindricalChannel;radius=0.15;central_point=(0, 0, 0);direction=y;es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=2
CylindricalChannel;radius=0.2;central_point=(0, 0, 0);direction=z;es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=2
```

Sample input `example2.puc`:
```
BaseCell;size=(1, 1, 1);el_size=0.1;mat_id=1
SandwichLayer;thickness=0.1;central_point=(0, 0, 0);direction=x;es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=2
SandwichLayer;thickness=0.1;central_point=(0, 0, 0);direction=y;es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=2
SandwichLayer;thickness=0.1;central_point=(0, 0, 0);direction=z;es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=2
EllipsoidalInclusion;radius=(0.4, 0.4, 0.4);central_point=(0, 0, 0);direction=(1, 0, 0);es_dmin=1.1;es_dmax=1.3;es_in=0.5;mat_id=1
```
