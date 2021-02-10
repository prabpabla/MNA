#--------import required modules-----------
import locale
import numpy as np
from opentest import *
n=0

value=list(map(float,value))
node1=list(map(int,node1))
node2=list(map(int,node2))
refnode1=list(map(int,refnode1))
refnode2=list(map(int,refnode2))
refn=0

n1 = max(node1)
n2 = max(node2)

#---------setting n to the highest node value---------
if (n1 > n2):
        n = n1
elif(n1 < n2):
        n = n2
elif(n1 == n2):
        n = n1
nmin=n

#------incrementing n to acomodate for inductors-----
for index in name:
        if ((index[0:1] == "L") or (index[0:1] == "E") or (index[0:4] == "VCVS")):
                n+=1
nmax=n
#--------------create G, C, and b matrix-------------
G = np.zeros((n,n))
C = np.zeros((n,n))
b = np.zeros((n,1))

print ("Min Node: ",nmin)
print("Max Node: ",nmax)
#-------------Put values in matrix --------------------
for index in range(len(name)):


        #--------------resistor stamp-----------------
        if (name[index][0:1] == "R"):

                mn1 = node1[index]
                mn2 = node2[index]
                g = 1/(value[index])

                if(mn1 == 0):
                        mn2-=1
                        G[mn2,mn2] = G[mn2,mn2] + g
                elif(mn2==0):
                        mn1-=1
                        G[mn1,mn1] = G[mn1,mn1] + g
                else:
                        mn1-=1
                        mn2-=1
                        G[mn1,mn2] = G[mn1,mn2] - g
                        G[mn1,mn1] = G[mn1,mn1] + g
                        G[mn2,mn2] = G[mn2,mn2] + g
                        G[mn2,mn1] = G[mn2,mn1] - g

        #--------------capacitor stamp--------------
        if((name[index][0:1] == "C") and (name[index][0:3] != "CUR")):
                mn1 = node1[index]
                mn2 = node2[index]
                c=value[index]

                if(mn1 == 0):
                        mn2-=1
                        C[mn2,mn2] = C[mn2,mn2] + c
                elif(mn2==0):
                        mn1-=1
                        C[mn1,mn1] = C[mn1,mn1] + c
                else:
                        mn1-=1
                        mn2-=1
                        C[mn1,mn2] = C[mn1,mn2] - c
                        C[mn1,mn1] = C[mn1,mn1] + c
                        C[mn2,mn2] = C[mn2,mn2] + c
                        C[mn2,mn1] = C[mn2,mn1] - c

        #----inductor currrnt flow from node 1 to node 2----
        if(name[index][0:1] == "L"):
                mn1 = node1[index]
                mn2 = node2[index]
                l=value[index]

                if(mn1 == 0):
                        mn2-=1
                        G[mn2,nmin] = G[mn2,nmin]-1
                        G[nmin,mn2] = G[nmin,mn2]-1
                        C[nmin,nmin] = C[nmin,nmin]-l
                elif(mn2==0):
                        mn1-=1
                        G[mn1,nmin] = G[mn1,nmin]+1
                        G[nmin,mn1] = G[nmin,mn1]+1
                        C[nmin,nmin] = C[nmin,nmin]-l
                else:
                        mn1-=1
                        mn2-=1
                        G[mn2,nmin] = G[mn2,nmin]-1
                        G[nmin,mn2] = G[nmin,mn2]-1
                        G[mn1,nmin] = G[mn1,nmin]+1
                        G[nmin,mn1] = G[nmin,mn1]+1
                        C[nmin,nmin] = C[nmin,nmin]-l
                nmin+=1

        #--- current source current flow from node 2 to node 1---
        if(name[index][0:3] == "CUR"):
                mn1 = node1[index]
                mn2 = node2[index]
                cur=value[index]

                if(mn1 == 0):
                        mn2-=1
                        b[mn2,0] = b[mn2,0]-cur
                elif(mn2==0):
                        mn1-=1
                        b[mn1,0] = b[mn1,0]+cur
                else:
                        mn1-=1
                        mn2-=1
                        b[mn2,0] = b[mn2,0]-cur
                        b[mn1,0] = b[mn1,0]+cur

        # --- current flow from node 1 to node 2----
        if(name[index][0:1] == "E"):
                #print("AAAAAAAAAAAAAAAAAAA")
                mn1 = node1[index]
                mn2 = node2[index]
                E=value[index]

                if(mn1 == 0):
                        mn2-=1
                        G[mn2,nmin] = G[mn2,nmin]-1
                        G[nmin,mn2] = G[nmin,mn2]-1
                        b[nmin,0] = b[nmin,0]+E
                elif(mn2==0):
                        mn1-=1
                        G[mn1,nmin] = G[mn1,nmin]+1
                        G[nmin,mn1] = G[nmin,mn1]+1
                        b[nmin,0] = b[nmin,0]+E
                else:
                        mn1-=1
                        mn2-=1
                        G[mn2,nmin] = G[mn2,nmin]-1
                        G[nmin,mn2] = G[nmin,mn2]-1
                        G[mn1,nmin] = G[mn1,nmin]+1
                        G[nmin,mn1] = G[nmin,mn1]+1
                        b[nmin,0] = b[nmin,0]+E
                nmin+=1

        # --- current flow from +(n1) to -(n2) ----
        if(name[index][0:4] == "VCVS"):
                mn1 = node1[index]
                mn2 = node2[index]
                mn3 = refnode1[refn]
                mn4 = refnode2[refn]
                Ea=value[index]
                refn+=1
                if(mn1 == 0):
                        mn2-=1
                        mn3-=1
                        mn4-=1
                        G[mn2,nmin] = -1
                        G[nmin,mn2] = -1
                        G[nmin,mn3] = G[nmin,mn3]-Ea
                        G[nmin,mn4] = G[nmin,mn4]+Ea

                elif(mn2==0):
                        mn1-=1
                        mn3-=1
                        mn4-=1 
                        G[mn1,nmin] = 1
                        G[nmin,mn1] = 1
                        G[nmin,mn3] = G[nmin,mn3]-Ea
                        G[nmin,mn4] = G[nmin,mn4]+Ea

                elif(mn3==0):
                        mn1-=1
                        mn2-=1
                        mn4-=1
                        G[mn1,nmin] = 1
                        G[nmin,mn1] = 1
                        G[mn2,nmin] = -1
                        G[nmin,mn2] = -1
                        G[nmin,mn4] = G[nmin,mn4]+Ea

                elif(mn4==0):
                        mn1-=1
                        mn2-=1
                        mn3-=1
                        G[mn1,nmin] = 1
                        G[nmin,mn1] = 1
                        G[mn2,nmin] = -1
                        G[nmin,mn2] = -1
                        G[nmin,mn3] = G[nmin,mn3]-Ea

                else:
                        mn1-=1
                        mn2-=1
                        mn3-=1
                        mn4-=1
                        G[mn1,nmin] = 1
                        G[nmin,mn1] = 1
                        G[mn2,nmin] = -1
                        G[nmin,mn2] = -1
                        G[nmin,mn3] = G[nmin,mn3]-Ea
                        G[nmin,mn4] = G[nmin,mn4]+Ea
                nmin+=1

        # --- current flow from n1 to n2 or n3 to n4----
        if(name[index][0:4] == "VCCS"):
                mn1 = node1[index]
                mn2 = node2[index]
                mn3 = refnode1[refn]
                mn4 = refnode2[refn]
                Ivccs = value[index]
                refn+=1
                if(mn1 == 0):
                        mn2-=1
                        mn3-=1
                        mn4-=1
                        G[mn2,mn3] = G[mn2,mn3]+Ivccs
                        G[mn2,mn4] = G[mn2,mn4]-Ivccs
                elif(mn2==0):
                        mn1-=1
                        mn3-=1
                        mn4-=1
                        G[mn1,mn3] = G[mn1,mn3]+Ivccs
                        G[mn1,mn4] = G[mn1,mn4]-Ivccs

#----------------- Testing ----------------------

print("printing G \n", G)
print("printing C \n", C)
print("printing b \n", b)
#A = np.mat(G)+np.mat(C)
#x=np.mat(A.T)*np.mat(b)
#print(A)
#print(x)
