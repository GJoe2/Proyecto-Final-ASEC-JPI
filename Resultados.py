import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Leemos los resultados
df6=pd.read_csv('Modelos_analizados.csv',usecols= ['Nmodel','a(m)','b(m)','h(m)','t(m)','Lm(m)','V(kN)','Δmax(‰)','VolConc(m3)'])
#print(df6)

#Buscando el modelo con las secciones óptimas
#Este será el que tiene mayor deriva (Δmax(‰)) y menor volumen de concreto (VolConc(m3))

Ropt = df6.loc[:,'Δmax(‰)']/df6.loc[:,'VolConc(m3)']
#print(Ropt)
df6=df6.assign(Ropt=Ropt)
print(df6)
#df6=pd.DataFrame(columns=['Nmodel','a(m)','b(m)','h(m)','t(m)','Lm(m)','V(kN)','Δmax(‰)','VolConc(m3)'])
#print(df6)
#x = df6['Ropt'].argmax()
#y = df6['Ropt'].max()
#print(y)
#Ploteamos
#Δmax(‰) VS VolConc(m3)
#vecX = np.array(df6.loc[:,'Δmax(‰)'])
#vecY = np.array(df6.loc[:,'VolConc(m3)'])
#limx = 0.9*vecX.min()
#liMx = 1.1*vecX.max()
#limy = 0.9*vecY.min()
#liMy = 1.1*vecY.max()
#
#plt.plot(vecX,vecY,'r.',label='Modelo analizado',lw = 0.8)
#plt.legend()
#plt.xlabel('Δmax(‰)')
#plt.ylabel('VolConc(m3)')
#plt.axis([limx,liMx,limy,liMy])
#plt.yticks(np.arange(limy, liMy, (liMy-limy)/10))
#plt.savefig('./drift_volu.png')

