print("Proyecto final ASEC JPI v1.0")

#INSTALACIÓN DE LIBRERÍAS
import numpy as np
import openseespy.opensees as ops
import pandas as pd
from openseespy.opensees import *
import openseespy.postprocessing.ops_vis as opsv
import openseespy.postprocessing.Get_Rendering as opsgr 
import matplotlib.pyplot as plt

#SISTEMA DE UNIDADES
# Unidades Base
m = 1
kg = 1
s = 1
# Otras Unidades
cm = 0.01*m
N = kg*m/s**2
kgf = 9.81*N
Pa = N/m**2
MPa = 10**6*Pa
inch = 2.54*cm
ft = 12*inch
ksi = 6894757.2932*Pa
kip = ksi*inch**2
psi = 6894.76*Pa
# Constantes Físicas
g = 9.81*m/s**2

#PROPIEDADES DE MATERIALES Y SECCIONES
fc = 210*kg/cm**2
E = 151*fc**0.5*kgf/cm**2
G = 0.5*E/(1+0.2)
# Viga
b,h = 25*cm, 50*cm
Av = b*h
Izv = b*h**3/12
Iyv = b**3*h/12
aa, bb = max(b,h),min(b,h)
β= 1/3-0.21*bb/aa*(1-(bb/aa)**4/12)
Jxxv = β*bb**3*aa
# Columna
a = 60*cm
Ac = a**2
Izc = a**4/12
Iyc = a**4/12
β= 1/3-0.21*1.*(1-(1.)**4/12)
Jxxc = β*a**4
# Densidad del concreto
ρ = 2400*kg/m**3



#CREACIÓN DE MODELO
def GeoModel(dx, dy, h, nx, ny, nz, Lm):
    from numpy import zeros
    Lx, Ly, Lz = dx*nx, dy*ny, h*nz
    NN = (nx+1)*(ny+1)*(nz+1)+8*(nz+1)
    Nodes = zeros((NN,5))
    # Creando los nodos y asignando coordenadas
    c = 0
    for i in range(nz+1):
        for j in range(ny+1+2):
            for k in range(nx+1+2):
                #Esquinas
                if k == 0 and j == 0:
                    Nodes[c] = [c,0,0,i*h,0.25]
                    c = c + 1
                elif k == 0 and j == ny+2:
                    Nodes[c] = [c,0,ny*dy,i*h,0.25]
                    c = c + 1
                elif k == nx+2 and j == 0:
                    Nodes[c] = [c,nx*dx,0,i*h,0.25]
                    c = c + 1
                elif k == nx+2 and j == ny+2:
                    Nodes[c] = [c,nx*dx,ny*dy,i*h,0.25]
                    c = c + 1
                #Lm en x
                elif k == 1 and j == 0:
                    Nodes[c] = [c,Lm,0,i*h,0.25]
                    c = c + 1
                elif k == 1 and j == ny+2:
                    Nodes[c] = [c,Lm,ny*dy,i*h,0.25]
                    c = c + 1
                elif k == nx+1 and j == 0:
                    Nodes[c] = [c,nx*dx-Lm,0,i*h,0.25]
                    c = c + 1
                elif k == nx+1 and j == ny+2:
                    Nodes[c] = [c,nx*dx-Lm,ny*dy,i*h,0.25]
                    c = c + 1
                #Lm en y
                elif k == 0 and j == 1:
                    Nodes[c] = [c,0,Lm,i*h,0.25]
                    c = c + 1
                elif k == nx+2 and j == 1:
                    Nodes[c] = [c,nx*dx,Lm,i*h,0.25]
                    c = c + 1
                elif k == 0 and j == ny+1:
                    Nodes[c] = [c,0,ny*dy-Lm,i*h,0.25]
                    c = c + 1
                elif k == nx+2 and j == ny+1:
                    Nodes[c] = [c,nx*dx,ny*dy-Lm,i*h,0.25]
                    c = c + 1
                #Lm en x que no son puntos
                elif k == 1 and j != 0 and j != ny+2:
                    c = c
                elif k == nx+1 and j != 0 and j != ny+2:
                    c = c
                #Lm en y que no son puntos
                elif j == 1 and k != 0 and k != nx+2:
                    c = c
                elif j == ny+1 and k != 0 and k != nx+2:
                    c = c
                #Medios
                elif k == 0 and j != 0 and j != 1 and j != ny+1 and j != ny+2:
                    Nodes[c] = [c,0,(j-1)*dy,i*h,0.50]
                    c = c + 1
                elif k == nx+2 and j != 0 and j != 1 and j != ny+1 and j != ny+2:
                    Nodes[c] = [c,nx*dx,(j-1)*dy,i*h,0.50]
                    c = c + 1
                elif j == 0 and k != 0 and k != 1 and k != nx+1 and k != nx+2:
                    Nodes[c] = [c,(k-1)*dx,0,i*h,0.50]
                    c = c + 1
                elif j == ny+2 and k != 0 and k != 1 and k != nx+1 and k != nx+2:
                    Nodes[c] = [c,(k-1)*dx,ny*dy,i*h,0.50]
                    c = c + 1
                #Centros
                else:
                    Nodes[c] = [c,(k-1)*dx,(j-1)*dy,i*h,1.00]
                    c = c+1
    Nodes[:((nx+1)*(ny+1)+8),4]=0


    NE = (nx*(ny+1)+ny*(nx+1)+(nx+1)*(ny+1)-4)*nz+8*nz
    Elems = zeros((NE,6))
    # Creando las conexiones de los elementos verticales
    c = 0
    nn = (nx+1)*(ny+1)+8
    #Muros
    for i in range(nz):
        Elems[c] = [c,0+nn*i,1+nn*i,nn*(i+1),1+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,nx+1+nn*i,nx+2+nn*i,nx+1+nn*(i+1),nx+2+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,0+nn*i,nx+3+nn*i,nn*(i+1),nx+3+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,nx+2+nn*i,nx+4+nn*i,nx+2+nn*(i+1),nx+4+nn*(i+1),3]
        c = c + 1

        Elems[c] = [c,nx+5+(nx+1)*(ny-1)+nn*i,nx+5+(nx+1)*(ny-1)+2+nn*i,nx+5+(nx+1)*(ny-1)+nn*(i+1),nx+5+(nx+1)*(ny-1)+2+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,nx+5+(nx+1)*(ny-1)+1+nn*i,nx+5+(nx+1)*(ny-1)+1+nx+3+nn*i,nx+5+(nx+1)*(ny-1)+1+nn*(i+1),nx+5+(nx+1)*(ny-1)+1+nx+3+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,nn-1-(nx+2)+nn*i,nn-(nx+2)+nn*i,nn-1-(nx+2)+nn*(i+1),nn-(nx+2)+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,nn-2+nn*i,nn-1+nn*i,nn-2+nn*(i+1),nn-1+nn*(i+1),3]
        c = c + 1


    #Columnas
    for i in range(nz):
        for j in range(ny+1):
            for k in range(nx+1):
                if j == 0 and k != 0 and k != nx:
                    Elems[c] = [c,k+1+nn*i,k+1+nn*(i+1),0,0,1]
                    c = c + 1
                elif j != 0 and j != ny:
                    Elems[c] = [c,k+nx+5+(nx+1)*(j-1)+nn*i,k+nx+5+(nx+1)*(j-1)+nn*(i+1),0,0,1]
                    c = c + 1
                elif j == ny and k != 0 and k != nx:
                    Elems[c] = [c,k+nx+5+(nx+1)*(j-1)+3+nn*i,(k+nx+5+(nx+1)*(j-1)+3)+nn*(i+1),0,0,1]
                    c = c + 1
    # Creando las conexiones de los elementos horizontales en x
    for i in range(nz):
        for j in range(ny+1):
            for k in range(nx):
                if j == 0:
                    Elems[c] = [c,k+1+nn*(i+1),k+2+nn*(i+1),0,0,2]
                    c = c + 1
                elif j != 0 and j != ny:
                    Elems[c] = [c,k+nx+5+(nx+1)*(j-1)+nn*(i+1),k+nx+5+1+(nx+1)*(j-1)+nn*(i+1),0,0,2]
                    c = c + 1
                elif j == ny:
                    Elems[c] = [c,k+nx+5+(nx+1)*(j-1)+3+nn*(i+1),k+nx+5+(nx+1)*(j-1)+4+nn*(i+1),0,0,2]
                    c = c + 1

    # Creando las conexiones de los elementos horizontales en y
    for i in range(nz):
        for j in range(nx+1):
            for k in range(ny):
                if j == 0 and k==0:
                    Elems[c] = [c,nx+3+nn*(i+1),nx+5+nn*(i+1),0,0,2]
                    c = c + 1
                elif j == 0 and k!=0:
                    Elems[c] = [c,nx+5+(nx+1)*(k-1)+nn*(i+1),nx+5+(nx+1)*k+nn*(i+1),0,0,2]
                    c = c + 1
                elif j != 0 and j != nx and k==0:
                    Elems[c] = [c,j+1+nn*(i+1),j+1+nx+4+nn*(i+1),0,0,2]
                    c = c + 1
                elif j != 0 and j != nx and k !=0 and k != ny-1:
                    Elems[c] = [c,j+nx+5+(nx+1)*(k-1)+nn*(i+1),j+nx+5+(nx+1)*k+nn*(i+1),0,0,2]
                    c = c + 1
                elif j != 0 and j != nx and k == ny-1:
                    Elems[c] = [c,j+nx+5+(nx+1)*(k-1)+nn*(i+1),j+nx+5+(nx+1)*k+3+nn*(i+1),0,0,2]
                    c = c + 1
                elif j == nx and k==0:
                    Elems[c] = [c,nx+4+nn*(i+1),nx+4+(nx+1)+nn*(i+1),0,0,2]
                    c = c + 1
                elif j == nx and k != 0 and k != ny-1:
                    Elems[c] = [c,nx+4+(nx+1)+(nx+1)*(k-1)+nn*(i+1),nx+4+(nx+1)+(nx+1)*k+nn*(i+1),0,0,2]
                    c = c + 1
                elif j == nx and k == ny-1:
                    Elems[c] = [c,nx+4+(nx+1)+(nx+1)*(k-1)+nn*(i+1),nx+4+(nx+1)+(nx+1)*(k-1)+2+nn*(i+1),0,0,2]
                    c = c + 1

    # Creando centro de diafragmas
    Diap = zeros((nz,4))
    for i in range(nz):
        Diap[i] = [i+1000,Lx/2.0,Ly/2.0,h*(i+1)]
    #
    return Nodes, Elems, Diap



ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

# Generamos la malla
RigidDiaphragm = 'ON'
dx, dy, dz = 5*m, 5*m, 3*m
nx, ny, nz = 4, 4, 10

#Propiedades de los muros
#Muros
Lm = 1*m
t = 0.25*m

ops.uniaxialMaterial('Elastic', 1, E) #Concreto
ops.uniaxialMaterial('Elastic', 2, 2*10**6*kgf*cm**2) #Acero
#ops.uniaxialMaterial('Concrete01', 1, 210, 0.003, 280, 0.005) #Concreto

a = 20*cm
muf = int(round(Lm/a))
tt = np.zeros(muf)
tt[:] = t
ww = np.zeros(muf)
ww[:] = Lm/(muf)
ρρ = np.zeros(muf)
ρρ[:] = 0
conc = np.zeros(muf)
conc[:] = int(1)
acero = np.zeros(muf)
acero[:] = int(2)

# Nodos del Modelo
Nodes, Elems, Diap = GeoModel(dx,dy,dz,nx,ny,nz,Lm)

#CREAMOS NODOS DEL MODELO
# Creamos los nodos
for Ni in Nodes:
  ops.node(int(Ni[0]), *Ni[1:4])

# Definimos diafragmas rígidos
if RigidDiaphragm == 'ON':
  dirDia = 3 # perpendicular al plano del diafragma
  for Nd in Diap:
    ops.node(int(Nd[0]), *Nd[1:4])
    ops.fix(int(Nd[0]),*[0,0,1,1,1,0])
    NodesDi = []
    for Ni in Nodes:
      if Ni[3]==Nd[3]:
        NodesDi.append(int(Ni[0]))
    ops.rigidDiaphragm(dirDia,int(Nd[0]),*NodesDi)

#ASIGNAMOS RESTRICCIONES EN LA BASE
# Restricciones
ops.fixZ(0.0, *[1,1,1,1,1,1], '-tol', 1e-6)

#EJES LOCALES
#Establecemos transformación geométrica
ops.geomTransf('PDelta', int(1), *[1, 0, 0])
ops.geomTransf('Linear', int(2), *[1,-1, 0])

#DEFINIMOS ELEMENTOS CON SUS PROPIEDADES
# Creamos los elementos
for Ele in Elems:
  if int(Ele[5]) == 1: # 1 Columna
    ops.element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Ac, E, G, Jxxc, Iyc, Izc, int(Ele[5]),'-mass', ρ*Ac)
  elif int(Ele[5]) == 2: # 2 Viga
    ops.element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Av, E, G, Jxxv, Iyv, Izv, int(Ele[5]),'-mass', ρ*Av)#*(dx-a)/dx)
  elif int(Ele[5]) == 3: # 3 Muro
    ops.element('MVLEM_3D', int(Ele[0]), int(Ele[1]), int(Ele[2]), int(Ele[3]), int(Ele[4]), muf, '-thick', *tt[:], '-width', *ww[:], '-rho', *ρρ[:], '-matConcrete', *conc[:], '-matSteel', *acero[:], '-matShear', int(1),'-Density', ρ)


#PLOTEO DEL MODELO
#opsgr.plot_model("nodes", "elements")
opsgr.plot_model()
plt.show()