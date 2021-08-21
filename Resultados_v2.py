import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

#Leemos los resultados
df6=pd.read_csv('Modelos_analizados_v2.csv',usecols= ['#Niveles','Nmodel','a(m)','b(m)','h(m)','t(m)','Lm(m)','V(kN)','Δmax(‰)','VolConc(m3)'])
#print(df6)

df6.loc[:,'Lm(m)'] = round(df6.loc[:,'Lm(m)'],2)
#print(df6)
#Buscando el modelo con las secciones óptimas
#Este será el que tiene mayor deriva (Δmax(‰)) y menor volumen de concreto (VolConc(m3))

Ropt = df6.loc[:,'Δmax(‰)']/df6.loc[:,'VolConc(m3)'] #índice
#print(Ropt)
df6=df6.assign(Ropt=Ropt)
#print(df6)

#Separamos los #Lm y ploteamos
nzmin=5 #cantidad de pisos mínima
nzmax=10 #cantidad de pisos máxima
dnz=1 #intervalo de pisos
nnz=int((nzmax-nzmin)/dnz+1) #número de intervalos para nz (cantidad de pisos)
Nz=np.zeros(nnz) #lista de nz (cantidad de pisos)
nzi=nzmin
for i in range(nnz):
    Nz[i]=nzi
    nzi=nzi+dnz

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
vecX = np.array(df6.loc[:,'Ropt'])
vecY = np.array(df6.loc[:,'#Niveles'])
limx = 0.9*vecX.min()
liMx = 1.1*vecX.max()
limy = 0.9*vecY.min()
liMy = 1.1*vecY.max()

df7 = pd.DataFrame(columns=['#Niveles','Nmodel_opt','a(m)','b(m)','h(m)','t(m)','Lm(m)','V(kN)','Δmax(‰)','VolConc(m3)','Ropt']) #cuadro de modelos óptimos
for i in range(nLm) :
    dflmi = df6['Lm(m)'] == LM[i]
    dflm = df6[dflmi]
    #print(dflm)

    X = np.array(dflm.loc[:,'Ropt'])
    Y = np.array(dflm.loc[:,'#Niveles'])

    annotations = np.array(dflm.loc[:,'Nmodel'])

    plt.plot(X,Y,'.-',label='Modelo Lm=%.2fm'%(LM[i]),lw = 0.8)
    for j, label in enumerate(annotations):
        plt.annotate('%.0f'%(label),(X[j],Y[j]))

for i in range(nnz) :
    dfnz = df6['#Niveles'] == Nz[i]
    dfnz = df6[dfnz]
    opt = dfnz['Ropt'].idxmax()

    df7 = df7.append({'#Niveles':'%.0f'%Nz[i],'Nmodel_opt':'%.0f'%dfnz.loc[opt,'Nmodel'],'a(m)':dfnz.loc[opt,'a(m)'],'b(m)':dfnz.loc[opt,'b(m)'],'h(m)':dfnz.loc[opt,'h(m)'],'t(m)':dfnz.loc[opt,'t(m)'],'Lm(m)':dfnz.loc[opt,'Lm(m)'],'V(kN)':dfnz.loc[opt,'V(kN)'],'Δmax(‰)':dfnz.loc[opt,'Δmax(‰)'],'VolConc(m3)':dfnz.loc[opt,'VolConc(m3)'],'Ropt':dfnz.loc[opt,'Ropt']}, ignore_index=True)




#
plt.legend()
plt.xlabel('Ropt = Δmax(‰) / VolConc(m3)')
plt.ylabel('#Niveles')
plt.axis([limx,liMx,limy,liMy])
xx=(vecX.max()-vecX.min())/5
plt.xticks(np.arange(vecX.min(), vecX.max()+xx, xx))
plt.yticks(np.arange(vecY.min(), vecY.max()+1, 1))

adjust_text(texts, only_move={'points':'y', 'texts':'y'}, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
plt.show()
plt.savefig('./Niveles_Ropt.png')

print(df7.round(4))
