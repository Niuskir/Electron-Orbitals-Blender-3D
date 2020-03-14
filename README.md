Electron Orbital Visualization (probability clouds) generated in Blender 3D

Blender is the free and open source 3D creation suite. It supports the entirety of the 3D pipeline—modeling, rigging, animation, simulation, rendering, compositing and motion tracking, even video editing and game creation:https://www.blender.org/

The blend file above, which as-is can be directly loaded into Blender 2.82a, includes a Python script which calculates a defined number of probability cloud iso surfaces and generates the corresponding meshes. Diffent configurations of orbitals + corresonding probability levels can be generated based combinations of n, l & m + number of contours + mesh resolution.

The Python script uses Hydrogen atom wavefunctions to provide visualizations of the quantum equations of the atomic orbitals and is based on the great work done by Damon Allen Ph.D., Nick Polfer Ph.D., Corey Stedwell Ph.D and Nathan Roehr Ph.D. which can be found here: 

https://github.com/damontallen/Orbitals/blob/master/Hydrogen%20Orbitals%20(Feb%2018,%202014)%20(dynamic%20entry).ipynb 

and for the Marching Cubes solution to generate the iso meashes on Blender i used the great work done by Robert Forsman, Tom Sapiens and Paul Bourke which can be found here:

https://github.com/mutantbob/blender-marching-cubes

This Python script is not an add-on and must be run from within the Blender text editor.

Values for n, l, m, isolevel, number of contours can be changed in the script for various electron probabilty isosurfaces.

If the variable what = 'single' only one Blender object will be created based on the n, l & m values in lines 738. 740 & 743.
If variable what is unequal to 'single' multiple Blender objects will be created based on the all combinations of the values of n, l & m as defined in the "for" loops in lines 749, 753 and 756.

The number of contours (isolevel surfaces) is defined in line 726 in the variable n_o_c. The progressive "see-trough" capability of a Blender object with multiple isosurfaces (to show the probability electron "cloud") can be created using a Blender material which defines transparency. There are many videos on Youtube showing how to create such a transparent material in Blender. 

Each Blender object generated will be placed in the correct location in 3D space as defined by the n, l & m grid contained in the blend file.   

The script requires the Python modules sympy and mpmath to be installed in the following Blender folder (they do not come standard in Blender):
C:\Program Files\Blender Foundation\Blender 2.82\2.82\python\lib\site-packages

Sympy can be downloaded here: https://www.sympy.org/en/index.html
Mpmath can be downloaded here: http://mpmath.org/

The sw versions i ran this script with are: Blender 2.82a, sympy 1.5.1 and mpmath 1.1.0. This script ran without problems in the Blender 2.83 Alpha version (as of 03/14/20) as well (don't forget to add sympy and mpmath).   

The script calculates the orbital cloud(s) using the proper scientific formulas and then uses the Marching Cubes computer graphics algoritm to visualize the electron orbitals at various isosurfaces (probability levels).

I am not an experienced Python developer and not a scientist so am sure there are better and proper and faster ways to get things done. The script is complete in the sense that it generates single or multiple orbital(s) in a mesh object in Blender and you can play around changing n.l.m, # of contours & isolevel. Rendering can be done in EEVEE or Cycles..

This is a hobby project and was done just for the fun of it.

If you are interested is this kind of stuff using Blender, this can help you get started and if you know Blender, or want to learn Blender, the full capability of this software is available. 

Hope this is of value to someone.

Youtube video of some generated orbitals: 
https://youtu.be/ir7imCX-Ulw
