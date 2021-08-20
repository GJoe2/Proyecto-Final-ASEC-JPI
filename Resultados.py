import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Leemos los resultados
df6=pd.read_csv('Modelos_analizados.csv',usecols= ['Nmodel','a(m)','b(m)','h(m)','t(m)','Lm(m)','V(kN)','Δmax(‰)','VolConc(m3)'])
#print(df6)

df6.loc[:,'Lm(m)'] = round(df6.loc[:,'Lm(m)'],2)
#print(df6)
#Buscando el modelo con las secciones óptimas
#Este será el que tiene mayor deriva (Δmax(‰)) y menor volumen de concreto (VolConc(m3))

Ropt = df6.loc[:,'Δmax(‰)']/df6.loc[:,'VolConc(m3)'] #índice
#print(Ropt)
df6=df6.assign(Ropt=Ropt)
#print(df6)

#Separamos los Lm y ploteamos
Lmmin=4*0.25 #Longitud mínima de muro
Lmmax=5/2 #longitud máxima de muro
dLm=0.25 #intervalo de longitud de muro
nLm=int((Lmmax-Lmmin)/dLm+1) #número de intervalos para Lm (longitud de muro)
LM=np.zeros(nLm) #lista de Lm (longitud de muro)
Lmi=Lmmin
for i in range(nLm):
    LM[i]=Lmi
    Lmi=Lmi+dLm
#print(Lm)


#Ploteamos
#Δmax(‰) VS VolConc(m3)
vecX = np.array(df6.loc[:,'Δmax(‰)'])
vecY = np.array(df6.loc[:,'VolConc(m3)'])
limx = 0.9*vecX.min()
liMx = 1.1*vecX.max()
limy = 0.9*vecY.min()
liMy = 1.1*vecY.max()


for i in range(nLm) :
    dflmi = df6['Lm(m)'] == LM[i]
    dflm = df6[dflmi]
    #print(dflm)

    X = np.array(dflm.loc[:,'Δmax(‰)'])
    Y = np.array(dflm.loc[:,'VolConc(m3)'])

    #xopt[i] = dflm['Ropt'].argmax()
    #yopt[i] = dflm['Ropt'].max()


    plt.plot(X,Y,'.',label='Modelo Lm=%.2fm'%(LM[i]),lw = 0.8)



#x = df6['Ropt'].argmax()
#y = df6['Ropt'].max()
#print(y)


#
#plt.plot(vecX,vecY,'r.',label='Modelo analizado',lw = 0.8)
plt.legend()
plt.xlabel('Δmax(‰)')
plt.ylabel('VolConc(m3)')
plt.axis([limx,liMx,limy,liMy])
plt.yticks(np.arange(limy, liMy, (liMy-limy)/10))
plt.show()
#plt.savefig('./drift_volu.png')

