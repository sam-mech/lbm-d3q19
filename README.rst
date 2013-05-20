================
lbm-d3q19 readme
================
lbm-d3q19 is a simple, three-dimensional lattice-Boltzmann code to approximate
the Navier-Stokes equations of incompressible fluid flow.

lbm-d3q19 is licensed under the GNU General Public License, v.3.
See license.txt for more information.

Requirements
------------
The build requirements are:
  * A C compiler, from e.g. the `GNU Compiler Collection 
    <http://gcc.gnu.org/>`_ (GCC)

Obtaining lbm-d3q19
-------------------
The best way to keep up to date with subsequent updates, bugfixes and 
development, is to use the Git version control system. To obtain a local 
copy, execute::

 git clone https://github.com/anders-dc/lbm-d3q19.git

Build and run instructions
--------------------------
lbm-d3q19 is built using `make`, the platform-specific C compiler.
Execute the following commands from the root directory::

 make

