import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

#función adjust_text recuperado de: https://github.com/Phlya/adjustText

#Leemos los resultados
df6=pd.read_csv('ENTREGABLE\Modelos_analizados.csv',usecols= ['#Niveles','Nmodel','a(m)','b(m)','h(m)','t(m)','Lm(m)','V(kN)','Δmax(‰)','VolConc(m3)'])
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
dLm=0.10 #intervalo de longitud de muro
nLm=int((Lmmax-Lmmin)/dLm+1) #número de intervalos para Lm (longitud de muro)
LM=np.zeros(nLm) #lista de Lm (longitud de muro)
Lmi=Lmmin
for i in range(nLm):
    LM[i]=round(Lmi,2)
    Lmi=Lmi+dLm
#print(Lm)

#Ploteamos
#Δmax(‰) VS VolConc(m3)
vecX = np.array(df6.loc[:,'Ropt'])
vecY = np.array(df6.loc[:,'#Niveles'])
limx = 0.8*vecX.min()
liMx = 1.05*vecX.max()
limy = 0.9*vecY.min()
liMy = 1.025*vecY.max()

texts = []
df7 = pd.DataFrame(columns=['#Niveles','Nmodel_opt','a(m)','b(m)','h(m)','t(m)','Lm(m)','V(kN)','Δmax(‰)','VolConc(m3)','Ropt']) #cuadro de modelos óptimos

fig, axs = plt.subplots(figsize=(17.5, 10), dpi=100)


for i in range(nLm) :
    dflmi = df6['Lm(m)'] == LM[i]
    dflm = df6[dflmi]
    #print(dflm)

    X = np.array(dflm.loc[:,'Ropt'])
    Y = np.array(dflm.loc[:,'#Niveles'])


    annotations = np.array(dflm.loc[:,'Nmodel'])

    plt.plot(X,Y,'.-',label='Lm=%.2fm'%(LM[i]),lw = 0.8)
    for j, label in enumerate(annotations):
        texts.append(plt.text(X[j],Y[j],'%.0f'%(label),fontsize=8))
 
for i in range(nnz) :
    dfnz = df6['#Niveles'] == Nz[i]
    dfnz = df6[dfnz]
    opt = dfnz['Ropt'].idxmax()

    df7 = df7.append({'#Niveles':'%.0f'%Nz[i],'Nmodel_opt':'%.0f'%dfnz.loc[opt,'Nmodel'],'a(m)':dfnz.loc[opt,'a(m)'],'b(m)':dfnz.loc[opt,'b(m)'],'h(m)':dfnz.loc[opt,'h(m)'],'t(m)':dfnz.loc[opt,'t(m)'],'Lm(m)':dfnz.loc[opt,'Lm(m)'],'V(kN)':dfnz.loc[opt,'V(kN)'],'Δmax(‰)':dfnz.loc[opt,'Δmax(‰)'],'VolConc(m3)':dfnz.loc[opt,'VolConc(m3)'],'Ropt':dfnz.loc[opt,'Ropt']}, ignore_index=True)
    



#
plt.title('Ropt vs #Niveles para cada Lm', fontsize=16)
plt.xlabel('Ropt = Δmax(‰) / VolConc(m3)', fontsize = 12)
plt.ylabel('#Niveles', fontsize = 12)
plt.axis([limx,liMx,limy,liMy])
xx=(vecX.max()-vecX.min())/5
plt.xticks(np.arange(vecX.min(), vecX.max()+xx, xx), fontsize = 10)
plt.yticks(np.arange(vecY.min(), vecY.max()+1, 1), fontsize = 10)

plt.legend(ncol=2,title='N° de modelo para cada Lm:' ,loc=1, fontsize = 10,fancybox=True)
adjust_text(texts, arrowprops=dict(arrowstyle='-'),expand_text=(2, 2),expand_points=(1.5, 1.5),expand_objects=(1.5, 1.5))

#plt.show()
plt.savefig('ENTREGABLE./Niveles_Ropt.png')

print(df7.round(4))
df7.to_csv('ENTREGABLE\Modelos_óptimos.csv')
