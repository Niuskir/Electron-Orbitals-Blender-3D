import bpy, bmesh
import numpy as np
import sympy

#"I" is sympy's imaginary number
from sympy import symbols,I,latex,pi,diff 
from sympy.utilities.lambdify import lambdastr
from sympy import factorial as fac
from sympy.functions import Abs,sqrt,exp,cos,sin
from sympy import re, im, simplify
import warnings # in order to suppress divide_by_zero warnings...

#display the latex representation of a symbolic variable by default.
from sympy import init_printing 
init_printing(use_unicode=True)

a_0,z,r=symbols("a_0,z,r")
n,m,l=symbols("n,m,l",integer=True)
int_m=symbols("int_m",integer=True)
theta,phi = symbols("\\theta,\\phi",real=True)

#The variables will used with lambdify...
angle_theta, angle_phi, radius = symbols("angle_theta,angle_phi,radius",real=True)


print("numpy version:   %s"%np.__version__)
print("sympy version:   %s"%sympy.__version__)

def createMesh(name, origin, verts, edges, faces):
    # Create mesh and object
    me = bpy.data.meshes.new(name+'Mesh')
    ob = bpy.data.objects.new(name, me)
    ob.location = origin
    ob.show_name = True
    # Link object to scene
    bpy.context.scene.objects.link(ob)

    # Create mesh from given verts, edges, faces. Either edges or
    # faces should be [], or you ask for problems
    me.from_pydata(verts, edges, faces)

    # Update mesh with new data
    me.update(calc_edges=True)
    return ob

def P_l(l,theta): #valid for l greater than equal to zero
    """Legendre polynomial"""
    if l>=0:
        eq=diff((cos(theta)**2-1)**l,cos(theta),l)
    else:
        print("l must be an integer equal to 0 or greater")
        raise ValueError
    return 1/(2**l*fac(l))*eq

def P_l_m(m,l,theta):
    """Legendre polynomial"""
    eq = diff(P_l(l,theta),cos(theta),Abs(m))
    result = sin(theta)**Abs(m)*eq #note 1-cos^2(theta) = sin^2(theta)
    return result

def Y_l_m(l,m,phi,theta):
    """Spherical harmonics"""
    eq = P_l_m(m,l,theta)
    if m>0:
        pe=re(exp(I*m*phi))*sqrt(2)
    elif m<0:
        pe=im(exp(I*m*phi))*sqrt(2)
    elif m==0:
        pe=1
    return abs(sqrt(((2*l+1)*fac(l-Abs(m)))/(4*pi*fac(l+Abs(m))))*pe*eq)

def L(l,n,rho):
    """Laguerre polynomial"""
    _L = 0.
    for i in range((n-l-1)+1): #using a loop to do the summation 
        _L += ((-i)**i*fac(n+l)**2.*rho**i)/(fac(i)*fac(n-l-1.-i)*\
                                          fac(2.*l+1.+i))
    return _L

def R(r,n,l,z=1.,a_0=1.):
    """Radial function"""
    rho = 2.*z*r/(n*a_0)
    _L = L(l,n,rho)
    _R = (2.*z/(n*a_0))**(3./2.)*sqrt(fac(n-l-1.)/\
         (2.*n*fac(n+l)**3.))*exp(-z/(n*a_0)*r)*rho**l*_L
    return _R

def Psi(r,n,l,m,phi,theta,z=1,a_0=1):
    """Wavefunction"""
    _Y = Y_l_m(l,m,phi,theta)
    _R = R(r,n,l)
    return _R*_Y

def P(r,n,l,m,phi,theta):
    """Returns the symbolic equation probability of the location 
    of an electron"""
    return Psi(r,n,l,m,phi,theta)**2*r**2

def display_orbital(n,l,m_,no_of_contours = 16,Opaque=0.5):
    """Diplays a 3D view of electron orbitals"""
    #The plot density settings (don't mess with unless you are sure)
    rng = 12*n*1.5 #This determines the size of the box 
    _steps = 100j#           (it needs to be bigger with n).
    steps = _steps.imag

    _x,_y,_z = np.ogrid[-rng:rng:_steps,-rng:rng:_steps,-rng:rng:_steps]

    P_tex = "" #initialize the LaTex string of the probabilities
    
    #Validate the quantum numbers
    assert(n>=1), "n must be greater or equal to 1"       #validate the value of n
    assert(0<=l<=n-1), "l must be between 0 and n-1"      #validate the value of l
    assert(-l<=max(m_)<=l), "p must be between -l and l"  #validate the value of p
    assert(-l<=min(m_)<=l), "p must be between -l and l"  #validate the value of p
    
        
    for m in m_:
        #Determine the probability equation symbolically and convert
        #it to a string
        prob = lambdastr((radius,angle_phi,angle_theta), P(radius,n,l,m,
                                                           angle_phi,
                                                           angle_theta))
        
        #record the probability equation as a LaTex string
        P_eq = simplify(P(r,n,l,m,phi,theta))
        P_tex+="$$P ="+latex(P_eq)+"$$ \n\n " 
        
        if '(nan)' in prob: #Check for errors in the equation
            print("There is a problem with the probability function.")
            raise ValueError
        
        #Convert the finctions in the probability equation from the sympy  
        #library to the numpy library to allow for the use of matrix 
        #calculations
        prob = prob.replace('sin','np.sin') #convert to numpy
        prob = prob.replace('cos','np.cos') #convert to numpy
        prob = prob.replace('Abs','np.abs') #convert to numpy
        prob = prob.replace('pi','np.pi')   #convert to numpy
        prob = prob.replace('exp','np.exp') #convert to numpy
        
        #convert the converted string to a callable function
        Prob = eval(prob)

        #generate a set of data to plot the contours of.
        w = Prob(r_fun(_x,_y,_z),phi_fun(_x,_y,_z),theta_fun(_x,_y,_z))
                    
        #Remove nan's in grid w               
        w[np.isnan(w)] = 0                
        
        #Determine minimum and maximum value (probability density) in the grid w 
        minw = np.nanmin(w)
        maxw = np.nanmax(w)
#        print(minw, maxw)
        #Determine value slices 
        colorSlice = (maxw - minw) / (no_of_contours)
        
        #Create contour lookup array 
        pr = np.linspace(minw, maxw, no_of_contours)
#        print(pr)
            
        buildGridAndCreate(int(steps), pr, w)

        #########
        #select object as active and change toggle to editmode
        bpy.context.scene.objects.active = bpy.data.objects['Orbital']
        bpy.ops.object.mode_set(mode='EDIT')        
        # Get the active mesh
        ob = bpy.context.edit_object
        me = ob.data


        # Get a BMesh representation
        bm = bmesh.from_edit_mesh(me)

        # Modify the BMesh, can do anything here...
        for v in bm.verts:
            v.co.x += 1.0

        # Show the updates in the viewport
        # and recalculate n-gon tessellation.
        bmesh.update_edit_mesh(me, True)

        bpy.ops.object.mode_set(mode='OBJECT')
        ########
               
        
        #Information used for the 2D slices below
        limits = []
        lengths = []
        for cor in (_x,_y,_z):
            limit = (np.min(cor),np.max(cor))
            limits.append(limit)
            #print(np.size(cor))
            lengths.append(np.size(cor))
            #print(limit)
        return (limits, lengths, _x, _y, _z, P_tex)
        

def buildGridAndCreate(steps, pr, w):

    # mesh arrays
    verts = []
    edges = []
    faces = []
    origin = (0,0,0) 
    # mesh variables    
    scale = 0.1
     
    #fill verts array
    for i in range (0,(steps)):
        for j in range(0,(steps)):
            for k in range(0,(steps)):
     
                colorindex = np.searchsorted(pr, w[i,j,k])
                if colorindex == 2:
                    vert = (i*scale,j*scale,k*scale) 
                    verts.append(vert)
    
#    numX = 50
#    for i in range (0, numX):
#        A = i
#        B = i+1
#        edge = (A,B)
#        edges.append(edge)

 
    #fill faces array

#    numX = 2
#    for i in range (0, numX):
#        A = i
#        B = i+1
#        C = (i+2)
# 
#        face = (A,B,C)
#        faces.append(face)

     
    ob1 = createMesh('Orbital', origin, verts, [], faces)
    
    # Move object to center
    bpy.data.objects['Orbital'].select = True
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')
    bpy.ops.object.location_clear()




#Start the calculation

r_fun = lambda _x,_y,_z: (np.sqrt(_x**2+_y**2+_z**2))
theta_fun = lambda _x,_y,_z: (np.arccos(_z/r_fun(_x,_y,_z)))
phi_fun = lambda _x,_y,_z: (np.arctan(_y/_x)*(1+_z-_z))

#Delete generated Blender object
bpy.ops.object.select_all(action='DESELECT')
if bpy.data.objects.get("Orbital") is not None:
    bpy.data.objects['Orbital'].select = True
    bpy.ops.object.delete(use_global=False)

n = 4
l = 3
m_ = [1]
n_o_c = 64
opacity = 0.5
P_tex = "" #initialize the LaTex string of the probabilities

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    limits, lengths, _x, _y, _z, P_tex = display_orbital(n,l,m_,n_o_c,opacity)

txt = r"$$\textbf{The symbolic expression for the resulting probability equation is:}$$ "
print(txt+P_tex)