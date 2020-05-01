#!/bin/bash

git clone https://github.com/schrodinger/pymol-open-source.git
cd pymol-open-source
prefix=$HOME/pymol-open-source-build
# sudo apt-get install build-essential python-dev python-pmw libglew-dev freeglut3-dev libpng-dev libfreetype6-dev libxml2-dev libmsgpack-dev python-pyqt5.qtopengl libglm-dev libnetcdf-dev
python setup.py build install --home=$prefix
