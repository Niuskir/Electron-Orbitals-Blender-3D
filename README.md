Electron Orbital Visualization (probability clouds) generated in Blender 3D

Blender is the free and open source 3D creation suite. It supports the entirety of the 3D pipeline—modeling, rigging, animation, simulation, rendering, compositing and motion tracking, even video editing and game creation:https://www.blender.org/

The blend file above, which as-is can be directly loaded into Blender 2.82a, includes a Python script which calculates a defined number of probability cloud iso surfaces and generates the corresponding meshes. Diffent configurations of orbitals + corresonding probability levels can be generated based combinations of n, l & m + number of contours + mesh resolution.

The Python script uses Hydrogen atom wavefunctions to provide visualizations of the quantum equations of the atomic orbitals and is based on the great work done by Damon Allen Ph.D., Nick Polfer Ph.D., Corey Stedwell Ph.D and Nathan Roehr Ph.D. which can be found here: 

https://github.com/damontallen/Orbitals/blob/master/Hydrogen%20Orbitals%20(Feb%2018,%202014)%20(dynamic%20entry).ipynb 

and for the Marching Cubes solution to generate the iso meashes on Blender i used the great work done by Robert Forsman, Tom Sapiens and Paul Bourke which can be found here:

https://github.com/mutantbob/blender-marching-cubes

This Python script is not an add-on and must be run from within the Blender text editor.

Values for n, l, m, isolevel, number of contours can be changed in the script for various electron probabilty isosurfaces.

If the variable what in line 746 is 'single' only one Blender object will be created based on the n, l & m values in lines 752, 755 & 757.
If variable "what" is unequal to 'single', multiple Blender objects will be created based on the all combinations of the values of n, l & m as defined in the "for" loops in lines 763, 767 and 770. Please note that if you run for example "for n in range(1,8)" it will take hours (depending on your PC computing power) to generate all the combinations of probability clouds in Blender objects. 

The number of contours (isolevel surfaces) is defined in line 726 in the variable n_o_c. The progressive "see-trough" capability of a Blender object with multiple isosurfaces (to show the probability electron "cloud") can be created using a Blender material which defines transparency. There are many videos on Youtube showing how to create such a transparent material in Blender. 

Each Blender object generated will be placed in the correct location in 3D space as defined by the n, l & m grid contained in the blend file.   

The script requires the Python modules sympy and mpmath to be added to the blender python installation folder. This can be done by doing the following in Windows:

1) Navigate using the File Explorer to the Blender python install folder, in my case with Blender 3.0 this was C:\Program Files\Blender Foundation\Blender 3.0\3.0\python\bin.  
2) Copy the address in the File Explorerer address bar using Ctrl C 
3) Launch the Windows Command Prompt as Administrator
4) Enter cd and then Ctrl V to paste the address you copied above so you end up with "cd C:\Program Files\Blender Foundation\Blender 3.0\3.0\python\bin" and press enter. 
5) Enter "python -m pip install --upgrade pip" to update to the latest pip version
6) Enter "python -m pip install mpmath --upgrade --force" to install mpmath
7) Enter "python -m pip install sympy --upgrade --force" to install sympy. If needed this will also replace mpmath to the correct matching version. 

You now should have all the python modules needed to run the blender script. You need to repeat these steps with each new install of Blender. Unfortunately i don't know the correct steps on how to add sympy and mpmath in Blender using Linux.

The sw versions i ran this script with are: Blender 3.0.1, sympy 1.9 and mpmath 1.2.1. 

The script calculates the orbital cloud(s) using the proper scientific formulas and then uses the Marching Cubes computer graphics algoritm to visualize the electron orbitals at various isosurfaces (probability levels).

I am not an experienced Python developer and not a scientist so am sure there are better and proper and faster ways to get things done. The script is complete in the sense that it generates single or multiple orbital(s) in a mesh object in Blender and you can play around changing n.l.m, # of contours & isolevel. Rendering can be done in EEVEE or Cycles.

Please note that if you don't know Blender and also have limited computer knowledge you most likely will not be able to run the script succesfully.   

This is a hobby project and was done just for the fun of it.

If you are interested is this kind of stuff using Blender, this can help you get started and if you know Blender, or want to learn Blender, the full capability of this software is available free of charge and there is and increadible amount of free eduction and support available in many places on the interbet.  

Hope this is of value to someone.

Youtube video of some generated orbitals: 
https://youtu.be/ir7imCX-Ulw
