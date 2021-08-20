acmin=0.25
acmax=0.80
dac=0.05
nac=int((acmax-acmin)/dac+1)
bvmin=0.25
bvmax=0.40
dbv=0.05
nbv=int((bvmax-bvmin)/dbv+1)

ac=np.zeros((1,nac))
for i in range(0,nac):
    ac[0,i]=acmin
    acmin=acmin+dac

bv=np.zeros((1,nbv))
for i in range(0,nbv):
    bv[0,i]=bvmin
    bvmin=bvmin+dbv

dtm=0.05
Lmmax=2.5
dLm=0.10
hvmax=0.8
dhv=0.10