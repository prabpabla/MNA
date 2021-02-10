#--------import required modules-----------
import locale
import matplotlib.pyplot as plt
import numpy as np
import math as mt
import scipy as sc

from matrixFinal import *

frequency = 10e3
w = 2.0*mt.pi*frequency
Is = 1e-15
VT = 0.026

#bascially doing this x = ([10; 8; 0.8; 10; -0.009; -0.009]); from matlab
x = [10, 8, 0.8, 10, -0.009, -0.009]
x = np.reshape(x,(len(x),1))

error1=1
error2=1
index = 0

while (True):
    val = (Is*(mt.exp(x[2][0]/VT)-1))

    print("error",error1,"   ",error2)

    fx = [0,0,val,0,0,0]
    fx= np.reshape(fx,(len(fx),1))
    phi = np.add((np.dot(G,x)),(np.subtract(fx,b)))
    
    dfx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((Is/VT)*mt.exp(x[2][0]/VT)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    dfx = np.reshape(dfx,(6,6))
    
    M = np.add(G,dfx)
    
    P,L,U = sc.linalg.lu(M)
    Z = np.linalg.solve(L,np.dot(P,phi))
    y = np.linalg.solve(U,Z) 

    
    #delx = np.subtract(x,y) 
    delx = np.dot(y,-1)
    x = np.add(delx,x)
    error1 = sc.linalg.norm(delx)
    error2 = sc.linalg.norm(phi)
    
    
    
    if (error1<1e-12 and error2<1e-12):
        break
    
    index = index+1
    print(index)



print("Initial Solution of X\n",x)


h = 0.3e-6;
T = 0.0001;
time = np.arange(0,(4*T),h)

An = np.add(G,np.dot(C,(2/h)))
xArray = []
xValueArray1 = []
xValueArray2 = []

useFirst = False
if(useFirst == True):
    for i in range(0,len(time)-1):
        val = (Is*(mt.exp(x[2]/VT)-1))
        K = np.dot(np.subtract(np.dot(C,(2/h)),G),x)
        
        f = [0,0,val,0,0,0]
        f = np.reshape(f,(len(f),1))
        
        bn = [0, 0, 0, 0, (25*mt.sin(w*time[i])), 10]
        bn = np.reshape(bn,(len(bn),1))
        
        bn1 = [0, 0, 0, 0,(25*mt.sin(w*time[i+1])), 10]    
        bn1 = np.reshape(bn1,(len(bn1),1))
    
    
        while (True):
            val = (Is*(mt.exp(x[2][0]/VT)-1))
            fx = [0,0,val,0,0,0]  
            fx= np.reshape(fx,(len(fx),1))
            dfx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((Is/VT)*mt.exp(x[2][0]/VT)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];                   
            dfx = np.reshape(dfx,(6,6))
            phi = np.subtract(np.subtract(np.add(np.subtract(np.add(np.dot(An,x),fx),K),f),bn),bn1)
            
            M= np.add(An,dfx)
    
            P,L,U = sc.linalg.lu(M)
            Z = sc.linalg.solve(L,np.dot(P,phi))
            y = sc.linalg.solve(U,Z) 
            
            
            #delx = np.subtract(x,y)
            #error1 = sc.linalg.norm(np.subtract(delx,x))
            #error2 = sc.linalg.norm(phi)
            #x = delx
            
            delx = np.dot(y,-1)
            x = np.add(delx,x)
            error1 = sc.linalg.norm(delx)
            error2 = sc.linalg.norm(phi)
            
            
            if (error1<1e-12 and error2<1e-12):
                break
            
            xArray.append(x)
            xValueArray1.append(x[0])
            xValueArray2.append(x[2])
            
            index = index+1
            
else:
    counter = 0
    for i in range(0,len(time)-1):
        k = np.dot(np.subtract(np.dot(C,(2/h)),G),x)
        
        fxn = [0,0,(Is*(np.exp(x[2][0]/VT)-1)),0,0,0]
        fxn = np.reshape(fxn,(len(fxn),1))
        
        bn = [0, 0, 0, 0, (25*np.sin(w*time[i])), 10]
        bn = np.reshape(bn,(len(bn),1))
        
        bnplusone = [0, 0, 0, 0,(25*np.sin(w*time[i+1])), 10]    
        bnplusone = np.reshape(bnplusone,(len(bnplusone),1))

        while(True):
            A = np.dot(np.add(np.dot(C,(2/h)),G),x)
            fxnplusone = [0,0,(Is*(np.exp(x[2][0]/VT)-1)),0,0,0]
            fxnplusone = np.reshape(fxnplusone,(len(fxnplusone),1))
            
            dfxnplusone = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((Is/VT)*np.exp(x[2][0]/VT)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];                   
            dfxnplusone = np.reshape(dfxnplusone,(6,6))
            
            flux = np.subtract(np.subtract(np.add(np.subtract(np.add(A,fxnplusone),k),fxn),bn),bnplusone)
            J = np.add(np.add(G,np.dot(C,(2/h))),dfxnplusone)    
            
            P,L,U = sc.linalg.lu(J)
            Z = sc.linalg.solve(L,(np.dot(P,flux)))
            y = sc.linalg.solve(U,Z) 
            
            deltax = y
            deltax = np.dot(-1,deltax)
            
            x = np.add(x,deltax)
            counter = counter+1
            
           
            
            if(sc.linalg.norm(deltax)<1e-12 and sc.linalg.norm(flux)<1e-12):
                break;
                
            xArray.append(x)
            xValueArray1.append(x[0])
            xValueArray2.append(x[2])


print("Final Solution of X\n",x)  
plt.plot(xValueArray1)
plt.plot(xValueArray2)
plt.show()    
































'''
h = 0.01e-9 #time step
time = 0.0
itteration = 3000


value = 0.0
f = []
for i in range(0,itteration):  
    if time<=2e-9:
        value=0.0
    elif time>2e-9 and time<=3e-9:
        value=(5/1e-9) * time - 10
    elif time>3e-9 and  time <=13e-9:
        value=5.0
    elif time > 13e-9 and time <= 14e-9:
        value=((-5/1e-9)) * time + 70
    else: 
        value=0.0
        
    f.append(value)     
    time = time + h





C_div_H = np.divide(C,h)
G_plus_C_div_H = np.add(C_div_H,G) # A
'''






'''
Vout=[]
#PLOT:
fmin = 0#Hz
fmax = 10000#Hz
Nrpt = 1000#Number of frequency points
#g_tuple = len(G)
#print (g_tuple)
numberofNodes = len(G)

F = np.linspace(fmin, fmax, Nrpt)

for n in range(0, Nrpt):
    w = 2 * mt.pi * F[n]
    s = (1j * w)
      


    A = (G + s * C)

    X = sc.linalg.solve(A,b)
    #X=np.mat(A.T)*np.mat(b)
    #print(Vout)
    Vout.append(X[7])
#print((F))
#print(X)
#print(np.absolute(Vout))
F=list(map(int,F))
plt.plot(F, np.absolute(Vout),'o')
plt.show()
'''
