print("Proyecto final ASEC JPI v1.0")

#INSTALACIÓN DE LIBRERÍAS
import numpy as np
import openseespy.opensees as ops
import pandas as pd
from openseespy.opensees import *
import openseespy.postprocessing.ops_vis as opsv
import openseespy.postprocessing.Get_Rendering as opsplt
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
fc = 210*kgf/cm**2
E = 15100*(fc/(kgf/cm**2))**0.5*kgf/cm**2
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
def GeoModel(dx, dy, h, nx, ny, nz, Lmx, Lmy):
    from numpy import zeros
    Lx, Ly, Lz = dx*nx, dy*ny, h*nz
    NN = (nx+1)*(ny+1)*(nz+1)+8*(nz+1)
    Nodes = zeros((NN,5))
    #Calculando fracciones de masas en los nodos
    mnodes = zeros(5)
    if dx > dy and Lmx<(dx-dy)/2 and Lmy<dy/2:
        mnodes[0] = (Lmx**2+Lmy**2)/(8*dx*dy) #nudo esquina
        mnodes[1] = 0.25 - (dy-Lmx)**2/(8*dx*dy) #nudo extremo de muro en x
        mnodes[2] = 0.50 - Lmx/(4*dx) #nudo lateral en x
        mnodes[3] = ((dy+Lmy)**2-3*Lmy**2)/(8*dx*dy) #nudo extremo de muro en y
        mnodes[4] = 0.50 - (2*dy-Lmy)*Lmy/(8*dx*dy) #nudo lateral en y
    elif dx < dy and Lmx<dx/2 and Lmy<(dy-dx)/2:
        mnodes[0] = (Lmx**2+Lmy**2)/(8*dx*dy) #nudo esquina
        mnodes[1] = ((dx+Lmx)**2-3*Lmx**2)/(8*dx*dy) #nudo extremo de muro en x
        mnodes[2] = 0.50 - (2*dx-Lmx)*Lmx/(8*dx*dy) #nudo lateral en x
        mnodes[3] = 0.25 - (dx-Lmy)**2/(8*dx*dy) #nudo extremo de muro en y
        mnodes[4] = 0.50 - Lmy/(4*dy) #nudo lateral en y
    elif dx == dy and Lmx<dx/2 and Lmy<dy/2:
        mnodes[0] = (Lmx**2+Lmy**2)/(8*dx*dy) #nudo esquina
        mnodes[1] = 0.25 - (Lmx**2+dy**2-(dx+dy-Lmx)*Lmx)/(8*dx*dy) #nudo extremo de muro en x
        mnodes[2] = 0.50 - (dx+dy-Lmx)*Lmx/(8*dx*dy) #nudo lateral en x
        mnodes[3] = 0.25 - (Lmy**2+dx**2-(dx+dy-Lmy)*Lmy)/(8*dx*dy) #nudo extremo de muro en y
        mnodes[4] = 0.50 - (dx+dy-Lmy)*Lmy/(8*dx*dy) #nudo lateral en y

    # Creando los nodos y asignando coordenadas
    c = 0
    for i in range(nz+1):
        for j in range(ny+1+2):
            for k in range(nx+1+2):
                #Esquinas
                if k == 0 and j == 0:
                    Nodes[c] = [c,0,0,i*h,mnodes[0]]
                    c = c + 1
                elif k == 0 and j == ny+2:
                    Nodes[c] = [c,0,ny*dy,i*h,mnodes[0]]
                    c = c + 1
                elif k == nx+2 and j == 0:
                    Nodes[c] = [c,nx*dx,0,i*h,mnodes[0]]
                    c = c + 1
                elif k == nx+2 and j == ny+2:
                    Nodes[c] = [c,nx*dx,ny*dy,i*h,mnodes[0]]
                    c = c + 1
                #Lm en x
                elif k == 1 and j == 0:
                    Nodes[c] = [c,Lmx,0,i*h,mnodes[1]]
                    c = c + 1
                elif k == 1 and j == ny+2:
                    Nodes[c] = [c,Lmx,ny*dy,i*h,mnodes[1]]
                    c = c + 1
                elif k == nx+1 and j == 0:
                    Nodes[c] = [c,nx*dx-Lmx,0,i*h,mnodes[1]]
                    c = c + 1
                elif k == nx+1 and j == ny+2:
                    Nodes[c] = [c,nx*dx-Lmx,ny*dy,i*h,mnodes[1]]
                    c = c + 1
                #Lm en y
                elif k == 0 and j == 1:
                    Nodes[c] = [c,0,Lmy,i*h,mnodes[3]]
                    c = c + 1
                elif k == nx+2 and j == 1:
                    Nodes[c] = [c,nx*dx,Lmy,i*h,mnodes[3]]
                    c = c + 1
                elif k == 0 and j == ny+1:
                    Nodes[c] = [c,0,ny*dy-Lmy,i*h,mnodes[3]]
                    c = c + 1
                elif k == nx+2 and j == ny+1:
                    Nodes[c] = [c,nx*dx,ny*dy-Lmy,i*h,mnodes[3]]
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
                #Medios en zonas de esquina
                elif k == 0 and (j == 2 or j == ny):
                    Nodes[c] = [c,0,(j-1)*dy,i*h,mnodes[4]]
                    c = c + 1
                elif k == nx+2 and (j == 2 or j == ny):
                    Nodes[c] = [c,nx*dx,(j-1)*dy,i*h,mnodes[4]]
                    c = c + 1
                elif j == 0 and (k == 2 or k == nx):
                    Nodes[c] = [c,(k-1)*dx,0,i*h,mnodes[2]]
                    c = c + 1
                elif j == ny+2 and (k == 2 or k == nx):
                    Nodes[c] = [c,(k-1)*dx,ny*dy,i*h,mnodes[2]]
                    c = c + 1
                #Medios centrales
                elif k == 0 and j != 0 and j != 1 and j != 2 and j != ny and j != ny+1 and j != ny+2:
                    Nodes[c] = [c,0,(j-1)*dy,i*h,0.50]
                    c = c + 1
                elif k == nx+2 and j != 0 and j != 1 and j != 2 and j != ny and j != ny+1 and j != ny+2:
                    Nodes[c] = [c,nx*dx,(j-1)*dy,i*h,0.50]
                    c = c + 1
                elif j == 0 and k != 0 and k != 1 and k != 2 and k != nx and k != nx+1 and k != nx+2:
                    Nodes[c] = [c,(k-1)*dx,0,i*h,0.50]
                    c = c + 1
                elif j == ny+2 and k != 0 and k != 1 and k != 2 and k != nx and k != nx+1 and k != nx+2:
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
        Elems[c] = [c,0+nn*i,1+nn*i,1+nn*(i+1),nn*(i+1),3] #3 muro en x #Asignar los nodos de forma antihoraria
        c = c + 1
        Elems[c] = [c,nx+1+nn*i,nx+2+nn*i,nx+2+nn*(i+1),nx+1+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,0+nn*i,nx+3+nn*i,nx+3+nn*(i+1),nn*(i+1),4] #4 muro en y
        c = c + 1
        Elems[c] = [c,nx+2+nn*i,nx+4+nn*i,nx+4+nn*(i+1),nx+2+nn*(i+1),4]
        c = c + 1

        Elems[c] = [c,nx+5+(nx+1)*(ny-1)+nn*i,nx+5+(nx+1)*(ny-1)+2+nn*i,nx+5+(nx+1)*(ny-1)+2+nn*(i+1),nx+5+(nx+1)*(ny-1)+nn*(i+1),4]
        c = c + 1
        Elems[c] = [c,nx+5+(nx+1)*(ny-1)+1+nn*i,nx+5+(nx+1)*(ny-1)+1+nx+3+nn*i,nx+5+(nx+1)*(ny-1)+1+nx+3+nn*(i+1),nx+5+(nx+1)*(ny-1)+1+nn*(i+1),4]
        c = c + 1
        Elems[c] = [c,nn-1-(nx+2)+nn*i,nn-(nx+2)+nn*i,nn-(nx+2)+nn*(i+1),nn-1-(nx+2)+nn*(i+1),3]
        c = c + 1
        Elems[c] = [c,nn-2+nn*i,nn-1+nn*i,nn-1+nn*(i+1),nn-2+nn*(i+1),3]
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
Lmx, Lmy = 1*m, 1*m
t = 0.25*m

ops.uniaxialMaterial('Elastic', 1, E) #Concreto Axial
#ops.nDMaterial('ElasticIsotropic', int(1), E, 0.2) #Concreto Axial
ops.uniaxialMaterial('Elastic', 2, 2*10**6*kgf/cm**2) #Acero
ops.uniaxialMaterial('Elastic', 3, G) #Concreto Cortante

#Muros en x
a = 20*cm
mufx = int(round(Lmx/a))
ttx = np.zeros(mufx)
ttx[:] = t
wwx = np.zeros(mufx)
wwx[:] = Lmx/(mufx)
ρρx = np.zeros(mufx)
ρρx[:] = 0 #Cuantía vertical en muros
concx = np.zeros(mufx)
concx[:] = int(1)
acerox = np.zeros(mufx)
acerox[:] = int(2)

#Muros en y
a = 20*cm
mufy = int(round(Lmy/a))
tty = np.zeros(mufy)
tty[:] = t
wwy = np.zeros(mufy)
wwy[:] = Lmy/(mufy)
ρρy = np.zeros(mufy)
ρρy[:] = 0 #Cuantía vertical en muros
concy = np.zeros(mufy)
concy[:] = int(1)
aceroy = np.zeros(mufy)
aceroy[:] = int(2)

# Nodos del Modelo
Nodes, Elems, Diap = GeoModel(dx,dy,dz,nx,ny,nz,Lmx,Lmy)
#print(Elems[:8,:])
#text= np.savetxt('nodos.txt',Nodes)
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
  elif int(Ele[5]) == 3: # 3 Muro en x
    ops.element('MVLEM_3D', int(Ele[0]), int(Ele[1]), int(Ele[2]), int(Ele[3]), int(Ele[4]), mufx, '-thick', *ttx[:], '-width', *wwx[:], '-rho', *ρρx[:], '-matConcrete', *concx[:], '-matSteel', *acerox[:], '-matShear', int(3), '-Poisson', 0.2, '-Density', ρ)
    #ops.element('quad',int(Ele[0]), int(Ele[1]), int(Ele[2]), int(Ele[3]), int(Ele[4]), t, 'PlaneStrain', int(1),0,ρ)
  elif int(Ele[5]) == 4: # 4 Muro en y
    ops.element('MVLEM_3D', int(Ele[0]), int(Ele[1]), int(Ele[2]), int(Ele[3]), int(Ele[4]), mufy, '-thick', *tty[:], '-width', *wwy[:], '-rho', *ρρy[:], '-matConcrete', *concy[:], '-matSteel', *aceroy[:], '-matShear', int(3), '-Poisson', 0.2, '-Density', ρ)
    #ops.element('quad',int(Ele[0]), int(Ele[1]), int(Ele[2]), int(Ele[3]), int(Ele[4]), t, 'PlaneStrain', int(1),0,ρ)

#PLOTEO DEL MODELO
#opsplt.plot_model("nodes", "elements")
#opsplt.plot_model()
#plt.show()

#ASIGNACIÓN DE MASAS Y MODOS DE VIBRACIÓN
# Aplicando Cargas vivas y muertas
wLive = 250*kg/m**2
wLosa = 300*kg/m**2
wAcab = 100*kg/m**2
wTabi = 150*kg/m**2
wTotal = 1.0*(wLosa+wAcab+wTabi)+0.25*wLive
#print(wTotal)
#
Carga = wTotal*dx*dy*m**2
for Ni in Nodes:
    ops.mass(int(Ni[0]),Ni[4]*Carga,Ni[4]*Carga,0.0)

# Obtenemos los modos
Nmodes = 12

# ploteamos el modo 1
#opsplt.plot_modeshape(1, 2000)
#plt.savefig('./modo_1.png',dpi=300)
#opsplt.plot_modeshape(2, 2000)
#opsplt.plot_modeshape(3, 2000)
#plt.show()
vals = ops.eigen(Nmodes)
Tmodes = np.zeros(len(vals))
for i in range(Nmodes):
    Tmodes[i] = 2*np.pi/vals[i]**0.5
    print("T[%i]: %.5f"%(i+1,Tmodes[i]))

#ANÁLISIS PARA OBTENER LA MATRIZ DE MASAS
# Realizamos un análisis para obtener la matriz de Masas
ops.wipeAnalysis()
ops.system('FullGeneral')
ops.numberer("Plain")
ops.constraints('Transformation') 
ops.algorithm('Linear')
ops.analysis('Transient')
ops.integrator('GimmeMCK',1.0,0.0,0.0)
ops.analyze(1,0.0) 

# Obtenemos la matriz de Masas
N = ops.systemSize()         # Número de Grados de Libertad
Mmatrix = ops.printA('-ret')
Mmatrix = np.array(Mmatrix).reshape((N,N))
MF = Mmatrix[-3*nz:,-3*nz:]

#MASAS EFECTIVAS
Tags = ops.getNodeTags()
# print(Tags)
modo = np.zeros((Nmodes,3*nz))
for j in range(1,Nmodes+1):
    ind = 0
    for i in Tags[-nz:]:
        temp = ops.nodeEigenvector(i,j)
        modo[j-1,[ind,ind+1,ind+2]] = temp[0],temp[1],temp[-1]
        ind = ind + 3

# Definimos valores iniciales
Ux,Uy,Rz = np.zeros(3*nz),np.zeros(3*nz),np.zeros(3*nz)
Ux[0::3]=1
Uy[1::3]=1
Rz[2::3]=1
SUMx, SUMy, SUMr = 0., 0., 0.
ni = 0

Mx = sum(sum(MF[0::3,0::3]))
My = sum(sum(MF[1::3,1::3]))
Mr = sum(sum(MF[2::3,2::3]))

df3 = pd.DataFrame(columns=['Modo','T(s)','FPRx','FPRy','FPRr','SumUx','SumUy','SumRz'])
for j in range(1,Nmodes+1):
    FPx = modo[j-1].T@MF@Ux
    FPy = modo[j-1].T@MF@Uy
    FPr = modo[j-1].T@MF@Rz
    FPRx = FPx**2/Mx
    FPRy = FPy**2/My
    FPRr = FPr**2/Mr
    SUMx = SUMx + FPRx
    SUMy = SUMy + FPRy
    SUMr = SUMr + FPRr
    #
    if min(SUMx,SUMy,SUMr)>=0.90 and ni==0:
        ni = j
    df3 = df3.append({'Modo':j, 'T(s)':Tmodes[j-1],'FPRx':FPRx,'FPRy':FPRy,'FPRr':FPRr,'SumUx':SUMx,
                    'SumUy':SUMy,'SumRz':SUMr}, ignore_index=True)
print(df3.round(5))
print('N° mínimo de Modos a considerar:',ni)



#plt.show()



#ANÁLISIS ESTÁTICO
#@title
def espectro_E030(T,Z=0.45,U=1.5,S=1.0,Tp=0.4,Tl=2.5,R=1):
    from numpy import zeros
    n = len(T)
    E030 = zeros(n)
    for i in range(n):
        if T[i]>=0 and T[i]<0.2*Tp:
            E030[i]=2.5#1+7.5*T[i]/Tp
        elif T[i]>=0.2*Tp and T[i]<Tp:
            E030[i]=2.5
        elif T[i]>=Tp and T[i]<Tl:
            E030[i] = 2.5*(Tp/T[i])
        elif T[i]>=Tl:
            E030[i] = 2.5*(Tp*Tl/T[i]**2)
        else:
            print("El periodo no puede ser negativo!")
    return E030*Z*U*S/R

def get_static_loads(coef,p,h,T):
    from numpy import zeros
    n = len(h)
    V = coef*sum(p)
    F = zeros(n)
    #
    if T > 0.0 and T <= 0.5:
        k=1.0
    elif T>0.5:
        k = 0.75+0.5*T
    else:
        print('El periodo es negativo!')
    #
    div = 0.
    for i in range(n):
        div = div + p[i]*h[i]**k
    #
    for i in range(n):
        F[i] = p[i]*h[i]**k/div*V
    return F,k

#Análisis Estático en X
np.set_printoptions(precision=3,linewidth=300,suppress=True)
H = np.arange(1,nz+1)*dz
P = sum(MF[0::3,0::3])*9.80665 # Peso por nivel
#print(H,P)
Ro = 7.0
E030 = espectro_E030(Tmodes,Z=0.45,U=1.0,S=1.0,Tp=0.4,Tl=2.5,R=Ro)
F, k = get_static_loads(E030[0],P,H,Tmodes[0])
CR = E030[0]/(0.45*1.*1.)
print('C/R=',CR)
print(E030[0],k)

ops.wipeAnalysis()
opsplt.createODB('Dual 10-1','Sismo_X',Nmodes=ni)

##Aplicamos fuerzas nodales
ops.timeSeries('Linear',1)
ops.pattern('Plain',1,1)
Le = ny*dy*0.05
for i in range(nz):
    #print(int(Diap[i][0]))
    ops.load(int(Diap[i][0]),F[i],0.,0.,0.,0.,F[i]*Le)


#Realizamos el análisis

ops.constraints('Transformation')
ops.numberer('Plain')
ops.system('FullGeneral')
ops.algorithm('Linear')
ops.integrator('LoadControl',1)
ops.analysis('Static')
ops.analyze(1)
ops.wipe()
#Calculando cortantes
#VS = np.cumsum(F[::-1])[::-1]

#Resultados del análisis estático en X
# Desplazamientos
#df1 = pd.DataFrame(columns=['Nivel','Vx(kN)','UxMax(cm)','UyMax(cm)','DriftX(‰)','DriftY(‰)'])
#tempX, tempY = 0., 0.
#for i in range(nz):
#    desX = ops.nodeDisp(int(Diap[i][0]),1)
#    desY = ops.nodeDisp(int(Diap[i][0]),2)
#    rotZ = ops.nodeDisp(int(Diap[i][0]),6)
#    desX = desX + abs(rotZ*ny*dy/2)
#    desY = desY + abs(rotZ*nx*dx/2)
#    desX, desY = desX*0.75*Ro, desY*0.75*Ro
#    driftX = 1000.*(desX-tempX)/dz
#    driftY = 1000.*(desY-tempY)/dz 
#    tempX, tempY = desX, desY
#    df1 = df1.append({'Nivel':i+1,'Vx(kN)':VS[i]/1000,'UxMax(cm)':desX*100,'UyMax(cm)':desY*100,
#                    'DriftX(‰)':driftX,'DriftY(‰)':driftY}, ignore_index=True)
#print('\nANÁLISIS ESTÁTICO EN X')
#print(df1.round(4))

#Se plotea la deformación obtenida

opsplt.plot_deformedshape(Model="Dual 10-1", LoadCase="Sismo_X")
plt.show