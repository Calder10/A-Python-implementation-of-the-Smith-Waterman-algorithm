"""
    Implementazione dell'algoritmo di Waterman Smith
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import importlib

def blos_score(s1,s2,ms):
    pkg = importlib.import_module('Bio.SubsMat.MatrixInfo')
    X=getattr(pkg,ms)
    A=np.zeros((len(s2),len(s1)))
    for j in range (len(s2)):
        for i in range (len(s1)):
            if (s2[j],s1[i]) in X :
                A[j][i]=X[(s2[j],s1[i])]
            else:
                A[j][i]=X[(s1[i],s2[j])]
    return A


def my_min(M):
    min=A[0][0]
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if (M[i][j]<min):
                min=M[i][j]
    return min

def my_max(M):
    max=0
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if (M[i][j]>=max):
                max=M[i][j]
                i_max=i
                j_max=j
    return (max,i_max,j_max)

def norm_matrix(A):
    m=my_min(A)
    for i in range (0,len(s2)):
        for j in range(0,len(s1)):
            A[i][j]=A[i][j]-m
    return A

def ws(A,s1,s2,d):
    F=np.zeros((len(s2)+1,len(s1)+1))
    P = np.empty((len(s2)+1,len(s1)+1),dtype=np.dtype('U25'))
    F[0][0]=0
    P[0][0]="End"

    for j in range(1,len(s1)+1):
        F[0][j]=0
        P[0][j]="End"

    for i in range(1,len(s2)+1):
        F[i][0]=0
        P[i][0]='End'

    for i in range(1,len(s2)+1):
        for j in range(1,len(s1)+1):
            #print(str(i)+" "+str(j))
            v=[(A[i-1][j-1]+F[i-1][j-1]),(F[i][j-1]-d),(F[i-1][j]-d),0]
            #print(v)
            F[i][j]=max(v)
            aus=v.index(max(v))
            if aus==0:
                P[i][j]="Diag"
            if aus==1:
                P[i][j]="Left"
            if aus==2:
                P[i][j]="Up"
            if aus==3:
                P[i][j]='End'

    s1_a=[]
    s2_a=[]
    score,i,j=my_max(F)
    X_P=[]
    Y_P=[]
    while True:
        if P[i][j]=="End":
            break
        if P[i][j]=='Diag':
            X_P.append(i)
            Y_P.append(j)
            s1_a.append(s1[j-1])
            s2_a.append(s2[i-1])
            i=i-1
            j=j-1
        elif (P[i][j]=='Up'):
            X_P.append(i)
            Y_P.append(j)
            s1_a.append("-")
            s2_a.append(s2[i-1])
            i=i-1

        elif (P[i][j]=='Left'):
            X_P.append(i)
            Y_P.append(j)
            s2_a.append("-")
            s1_a.append(s1[j-1])
            j=j-1


    s1_a.reverse()
    s2_a.reverse()
    s1_a=' '.join(map(str,s1_a))
    s2_a=' '.join(map(str,s2_a))
    F=F[1:len(s2)+1,1:len(s1)+1]
    return s1_a,s2_a,score,X_P,Y_P,F


def plot(S,s1,s2,x,y,d,ms):
    """
        La funzione prende in input:
        S->Matrice con i relativi score
        s1,s2->Le 2 sequenze che si stanno allineando
        x,y-> acisse e ordinate dei punti che formano il percorso dell'allineamento

        Fornisce in output un grafico che ci mostra il percorso migliore ottenuto
    """
    t=s1+"_"+s2+" Delta="+str(d)+" MS="+ms
    plt.figure(figsize=(10,10))
    plt.title(t)
    plt.plot(y,x,'r')
    f=sns.heatmap(S,annot=True,cmap="Blues_r",linewidths=.5,xticklabels=s1,yticklabels=s2)
    plt.show()

print("Inserisci il nome della matrice di sostituzione da utilizzare-->")
ms=input()
print("Inserisci la sequenza 1--->")
s1=input()
print("Inserisci la sequenza 2--->")
s2=input()
print("Inserisci il valore di delta--->")
d=int(input())
A=blos_score(s1,s2,ms)
s1_a,s2_a,score,X_P,Y_P,F=ws(A,s1,s2,d)
print("Allineamento:"+"\n"+s1_a +" \n"+ s2_a +"\n"+ "Score: "+str(score))
plot(F,s1,s2,X_P,Y_P,d,ms)

"""
TFDERILGVQTYWAECLA
QTFWECIKGDNATY

QERTY
QERS

AFGIVHKLIVS
AFGIHKIVS

"""
