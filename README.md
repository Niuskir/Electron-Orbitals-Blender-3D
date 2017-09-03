# Electron-Orbitals-Blender-3D
Electron Orbital Visualization (probability clouds) generated in Blender 3D

Included a Blender Python script which calculates probabilty cloud isosurfaces and generates the corresponding meshes at various probabilty levels (isolevels). Various orbitals can be generated based on in the script settable combinations of n, l & m + number of contours + mesh resolution.

I have created this pyhton script by combining great work which can be found here: https://github.com/damontallen/Orbitals/blob/master/Hydrogen%20Orbitals%20(Feb%2018,%202014)%20(dynamic%20entry).ipynb and here: https://github.com/mutantbob/blender-marching-cubes

This Python script is not an add-on and must be run from the Blender text editor.

Values for n, l, m, isolevel, number of contours can be changed in the script for various electron probabilty isosurfaces.

You will need NUMPY, which comes standard with Blender, and SYMPY which you must install separately in the following Blender folder:
C:\Program Files\Blender Foundation\Blender\2.78\python\lib\site-packages
SYMPY is dependent on one or two additional Python extensions but you will get a message on which ones in the Blender console when you run the script. The easiest way i know about to get SYMPY and it's dependents is to install the Anaconda3 environment (https://docs.continuum.io/) and copy the SYMPY + it's dependents from C:\Anaconda3\Lib\site-packages to the Blender site packages folder specifies above.

The sw versions i ran this script with are: Blender 2.78c, NUMPY 14.1 SYMPY 1.0

The script calculates the orbital cloud(s) using the proper scientific formulas and then uses the Marching Cubes computer graphics algoritm to visualize the electron orbitals of the Hydrogen Atom in various isosurfaces.
This is my first Python script so am sure it there are many proper and faster ways to get things done for the parts i have written. The script is completed in the sense that it generates the orbital in a mesh (used Blender 2.78) and you can play around changing n.l.m, # of countours & isolevel. Rendering shows you the various probabilty isosurfaces going inwards.
This is a hobby project and done just for the fun of it.
If you are interested is this kind of stuff with Blender this can help you started and with Blender the sky is the limit.
Hope this is of value to someone.
08/30/17 I have added another script BlenderHydrogenOrbitalVisualization.py which creates a coud of verticies instead of a mesh for you to play with.
