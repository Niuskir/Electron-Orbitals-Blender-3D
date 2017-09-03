i# Electron-Orbitals-Blender-3D
Electron Orbital Visualization (probability clouds) generated in Blender 3D

Included a Python script (Electron_Orbitals_v4) in Blender (Electron_Orbitals_v4.blend) which calculates a defined number of probability cloud iso surfaces and generates the corresponding meshes. Diffent configurations of orbitals + corresonding probability levels can be generated based combinations of n, l & m + number of contours + mesh resolution.

I have created this pyhton script by combining great work which can be found here: 

https://github.com/damontallen/Orbitals/blob/master/Hydrogen%20Orbitals%20(Feb%2018,%202014)%20(dynamic%20entry).ipynb 

and here: 

https://github.com/mutantbob/blender-marching-cubes

This Python script is not an add-on and must be run from the Blender text editor.

Values for n, l, m, isolevel, number of contours can be changed in the script for various electron probabilty isosurfaces.

You will need NUMPY, which comes standard with Blender, and SYMPY which you must install separately in the following Blender folder:
C:\Program Files\Blender Foundation\Blender\2.78\python\lib\site-packages
SYMPY is dependent on one or two additional Python extensions but you will get a message on which ones are missing in the Blender console when you run the script. The easiest way i found to get SYMPY and it's dependents is to install the Anaconda3 environment (https://docs.continuum.io/) and copy the SYMPY + it's dependents from C:\Anaconda3\Lib\site-packages to the Blender site packages folder specified above.

The sw versions i ran this script with are: Blender 2.78c, NUMPY 14.1 SYMPY 1.0

The script calculates the orbital cloud(s) using the proper scientific formulas and then uses the Marching Cubes computer graphics algoritm to visualize the electron orbitals at various isosurfaces (probability levels).
I am not an experienced Python developer and not a scientist so am sure it there are many proper and faster ways to get things done. The script is complete in the sense that it generates the orbital in a mesh object in Blender 3D and you can play around changing n.l.m, # of contours & isolevel. Rendering (Cycles) shows you the various probability isosurfaces going into the cloud(s) due to a simple Transparency material.

This is a hobby project and was done just for the fun of it.

If you are interested is this kind of stuff with Blender this can help you started and if you know Blender well the sky is the limit in tems of what you can do withthe mesh/object.

Hope this is of value to someone.

I have added another script Electron_Orbitals_v1.py which creates a coud of verticies instead of a mesh. This requires still some work to get something viewable.
