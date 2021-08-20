import numpy as np
#=================================================================
#SISTEMA DE UNIDADES
# Unidades Base
m = 1
kg = 1
s = 1
# Otras Unidades
cm = 0.01*m
N = kg*m/s**2
kN = 1000*N
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
#=================================================================
#GEOMETRÍA DE ESTRUCTURA
dxy = 5*m #separación de las grillas en X e Y
dx, dy, dz = dxy, dxy, 3*m

nxy = 4 #número de grillas en X e Y
nx, ny = nxy, nxy

nzmin=5 #cantidad de pisos mínima
nzmax=10 #cantidad de pisos máxima
dnz=1 #intervalo de pisos
nnz=int((nzmax-nzmin)/dnz+1) #número de intervalos para nz (cantidad de pisos)
Nz=np.zeros(nnz) #lista de nz (cantidad de pisos)
nzi=nzmin
for i in range(nnz):
    Nz[i]=nzi
    nzi=nzi+dnz

#GEOMETRÍA DE ELEMENTOS ESTRUCTURALES
lstf=np.empty((0, 5))

acmin=0.25 #ancho mínimo de columna
acmax=0.80 #ancho máximo de columna
dac=0.05 #intervalo de ancho de columna
nac=int((acmax-acmin)/dac+1) #número de intervalos para ac (ancho de columna)
ac=np.zeros(nac) #lista de ac (ancho de columna)
aci=acmin
for i in range(nac):
    ac[i]=aci
    aci=aci+dac
# print(ac)

bvmin=0.25 #ancho mínimo de viga
bvmax=0.40 #ancho máximo de viga
dbv=0.05 #intervalo de nacho de viga
nbv=int((bvmax-bvmin)/dbv+1) #número de intervalos para bv (andho de viga)
bv=np.zeros(nbv) #lista de bv (ancho de viga)
bvi=bvmin
for i in range(nbv):
    bv[i]=bvi
    bvi=bvi+dbv
# print(bv)

tmmin=bvmin #espesor mínimo de muro
tmmax=bvmax #espesor máximo de muro
dtm=0.05 #intervalo de espesor de muro
ntm=int((tmmax-tmmin)/dtm+1) #número de intervalos para tm (espesor de muro)
tm=np.zeros(ntm) #lista de tm (ancho de muro)
tmi=tmmin
for i in range(ntm):
    tm[i]=tmi
    tmi=tmi+dtm
# print(tm)

Lmmin=4*tmmin #Longitud mínima de muro
Lmmax=dxy/2 #longitud máxima de muro
dLm=0.25 #intervalo de longitud de muro
nLm=int((Lmmax-Lmmin)/dLm+1) #número de intervalos para Lm (longitud de muro)
Lm=np.zeros(nLm) #lista de Lm (longitud de muro)
Lmi=Lmmin
for i in range(nLm):
    Lm[i]=Lmi
    Lmi=Lmi+dLm
# print(Lm)

hvmin=0.25 #peralte mínimo de viga
hvmax=0.8 #peralte máximo de viga
dhv=0.05 #intervalo de lperalte de viga
nhv=int((hvmax-hvmin)/dhv+1) #número de intervalos para hv (peralte de viga)
hv=np.zeros(nhv) #lista de hv (peralte de viga)
hvi=hvmin
for i in range(nhv):
    hv[i]=hvi
    hvi=hvi+dhv
# print(hv)
     

    
lstf=np.empty((0, 6))

total=nac*nbv*ntm*nLm*nhv
# print(total)
total1=int(total/nLm)
#print(total1)
lst1=np.zeros((total1,6))

lst2=np.zeros((total,6))

cont=0
for j in range(nLm):
    for k in range(nac):
        for ii in range(nbv):
            for jj in range(ntm):
                for kk in range(nhv):
                    a = ac[k]
                    b = bv[ii]
                    t = tm[jj]
                    Lmx, Lmy = Lm[j] ,Lm[j]
                    h = hv[kk]
                    #Para nz = 1
                    nz = 1

                    VolVx = (b*h*dx*(nx*(ny+1)-4))*nz
                    VolVy = (b*h*dy*(ny*(nx+1)-4))*nz
                    VolCol = (a*a*dz*(nx+1)*(ny+1)+8-12)*nz
                    VolMurosx = t*Lmx*dz*4*nz
                    VolMurosy = t*Lmy*dz*4*nz
                    Volumen = round(VolVx+VolVy+VolCol+VolMurosx+VolMurosy,2)

                    lst1[cont]=[a,b,t,Lmx,h,Volumen]
                    cont=cont+1
    lst1[lst1[:, 5].argsort()]
    #print(lst1)
    #print(cont)
    cont=0
    for i in range(total1):
        lst2[i+j*total1]=lst1[i]


#print(total)
                
for i in range(total):
    #CONDICIONALES
    #tm tiene como minimo bv
    if lst2[i,2] >= lst2[i,1]:
        #Lm tiene como minimo 4tm
        if lst2[i,3] >= 4*lst2[i,2]:
            #hv tiene como minimo 2bv
            if lst2[i,4] >= 2*lst2[i,1]:
                
                b=lst2[i,1]
                h=lst2[i,4]
                Izv = b*h**3/12
                a=lst2[i,0]
                Izc = a**4/12
                
                if Izv <= Izc:
                    
                    lstf=np.append(lstf,[lst2[i]],axis=0)

nlstf=len(lstf)

#print (lstf)

print('Número de modelos a analizar por cantidad de pisos (de 5 a 10 pisos) = %.0f'%(nlstf))
np.savetxt('lista_filtro.txt',lstf)
                    