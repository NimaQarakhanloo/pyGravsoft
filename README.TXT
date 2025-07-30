# $Id: README.TXT 291 2009-10-14 07:39:42Z tjansson $
##########################
## 1.1 ABOUT pyGravsoft
##########################
PyGravsoft provides the graphical user interface to the program in Gravsoft.
pyGravsoft is developed by:
Thomas R. N. Jansson (info@tjansson.dk) - www.tjansson.dk
C. C. Tscherning (cct@gfy.ku.dk) - www.gfy.ku.dk/~cct/

The src/ folder contains the source code for the fortran programs.
The bin/ folder contains the binary fortran programs.
The data/ folder contains example data files for use in the programs.
The doc/ folder contains the documentation of the programs. 

##########################
## 1.2 INSTALLATION
##########################
PyGravsoft depends on the Python binaries which can be downloaded here:
http://www.python.org/download/ 
Pygravsoft should work out-of-the-box on all modern Linux systems with python
and the python-tk packages installed. The program is also tested on Windows 
XP and Solaris 8 but PyGravsoft should work on all the machines where python 
is available. The program should work on other *NIX's as Mac OS X and Solaris. 

## 1.2.1 Linux/Unix
Most modern Linux distributions already have python install and there is no
installation. However in some cases the graphical toolkit is not installed and
on fedora for instance. Under fedora this would be resolved by running "yum
install tkinter".

Running the program is done by writing the following:
username@machine:~/pyGravsoft$./launcher.py
or 
username@machine:~/pyGravsoft$python launcher.py

The program uses a absolute path to the binary fortran programs. This path 
is given in the first line of geomodule.py and should be changed upon a new 
installation.

## 1.2.2 Windows
Windows does not have python installed by default and this should be downloaded
from the webpage:
http://www.python.org/download/
After the installation .py files will be associated with python programs and a
double-click on a the launcher.py file will start the program.

##########################
## 1.3 LICENSE
##########################

Copyright 2007 Thomas Jansson IT and university of Copenhagen.

This file is part of pyGravsoft. pyGravsoft constitutes all the python files
related to the graphical user interface. pyGravsoft is free software: you can
redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

pyGravsoft is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
pyGravsoft. If not, see <http://www.gnu.org/licenses/>.
