import numpy as np

dxy=5
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
     

    
lstf=np.empty((0, 5))
total=nac*nbv*ntm*nLm*nhv
# print(total)
lst1=np.zeros((total,5))

for i in range(0,int(total/nac)):
    for j in range(0,nac):
        lst1[j+i*nac,0]=ac[j]
        
for i in range(0,int(total/nbv)):
    for j in range(0,nbv):
        lst1[j+i*nbv,1]=bv[j]

for i in range(0,int(total/ntm)):
    for j in range(0,ntm):
        lst1[j+i*ntm,2]=tm[j]
        
for i in range(0,int(total/nLm)):
    for j in range(0,nLm):
        lst1[j+i*nLm,3]=Lm[j]
        
for i in range(0,int(total/nhv)):
    for j in range(0,nhv):
        lst1[j+i*nhv,4]=hv[j]
        
    
#CONDICIONALES
for i in range(0,total):
    #tm tiene como minimo bv
    if lst1[i,2] >= lst1[i,1]:
        #Lm tiene como minimo 4tm
        if lst1[i,3] >= 4*lst1[i,2]:
            #hv tiene como minimo 2bv
            if lst1[i,4] >= 2*lst1[i,1]:
                
                b=lst1[i,1]
                h=lst1[i,4]
                Izv = b*h**3/12
                a=lst1[i,0]
                Izc = a**4/12
                
                if Izv <= Izc:
                    
                    lstf=np.append(lstf,[lst1[i]],axis=0)
                    
                

        
        






