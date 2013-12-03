'''
Created on Jun 29, 2009

@author: lee
'''
#from numpy.ma.core import zeros

from numpy import *

#gaussian elimination based solver
def solver(m, sol, eps=1.0 / (10 ** 10)):
      (h, w) = (len(m), len(m[0]))
      for y in range(0, h):
            m[y].append(sol[y])
            
      w = w + 1
      for y in range(0, h):
            maxrow = y
            for y2 in range(y + 1, h):# find pivot
                  if abs(m[y2][y]) > abs(m[maxrow][y]):
                        maxrow = y2
                        
            temp = m[y]
            m[y] = m[maxrow]
            m[maxrow]=temp

            if abs(m[y][y]) <= eps:#check singularity
                  return []
            for y2 in range(y + 1, h):#eliminate column y
                  c = m[y2][y] / m[y][y]
                  for x in range(y, w):
                        m[y2][x] -= m[y][x] * c
      for y in range(h - 1, 0 - 1, -1):#back substitution
            c = m[y][y]
            for y2 in range(0, y):
                  for x in range(w - 1, y - 1, -1):
                        m[y2][x] -= m[y][x] * m[y2][y] / c
            m[y][y] /= c
            for x in range(h, w):# normalize row y
                  m[y][x] /= c
                  
      rtr = []
      for t in range(0,h):
            rtr.append(m[t][w-1])
            
      return rtr


from random import *
def divide(X,Y,split):
    n = len(X)
    mix=range(0,n)
    shuffle(mix)
      
    size = int(split*n)
      
    trainX = [0]*size
    trainY = [0]*size
    testX =  [0]*(n-size)
    testY =  [0]*(n-size)
      
    for i in range(0,size):
        trainX[i] = X[mix[i]]
        trainY[i] = Y[mix[i]]
      
    for i in range(size,n):
        testX[i-size] = X[mix[i]]
        testY[i-size] = Y[mix[i]]
            

            
    return [array(trainX),array(trainY),array(testX),array(testY)]

def nipals(X,Y,a, tol=10 ** -10):
      n = X.shape[0]
      m = X.shape[1]
      p = Y.shape[1]
      
      
      P = array(zeros((m,a)))
      Q=array(zeros((p,a)))
      W=array(zeros((m,a)))
      B=array(zeros((a,a)))
      T=array(zeros((n,a)))
      
      Yres = array(Y)
      #Determine principal components (SVD?)
      for h in range(0,a):
            tOld = zeros((n,1))
            u = Yres[:,0]
            while True:
                  wOld = dot(X.T,u)
                  w = wOld/linalg.norm(wOld)
                  t = dot(X,w)#/(w'*w);
                  
                  qOld = dot(Yres.T,t)
                  q = qOld/linalg.norm(qOld)
                  u = dot(Yres,q)
                  # Test for convergence
                  err = linalg.norm(tOld - t)/linalg.norm(t)
                  
                  if err < tol:
                  # Calculate the final values for this component
                        pOld = dot(X.T,t)
                        p = pOld/linalg.norm(pOld)
                        # Store component vectors into the arrays
                        P[:,h] = p
                        Q[:,h] = q
                        W[:,h] = w
                        T[:,h] = t
                        # One difference from Geladi's algorithm is that
                        # we never calculate the X-residual, so we need
                        # to calculate the regression coefficients using
                        # least squares method (See Eq. 14 in Geladi.)
                        
                        top = dot(T[:,0:h+1].T,T[:,0:h+1]) 
                        bottom =dot( T[:,0:h+1].T,u)
                        B[0:h+1,h]=solver(top.tolist(),bottom.tolist())

                        #B[0:h+1,h]= linalg.lstsq(top,bottom)[0]
                        
                        Ypred = dot(dot(T,B),Q.T)
                        Yres = Y - Ypred
                        break
                  else:
                        tOld = t
      return P,Q,W,B


def readDataMats():
    f = file("grades",'r')
    l = f.readline()[:-1].split('\t')
    Y = []
    data = [float(p) for p in l]
    Y.append(data)
    Y = array(Y).T

    f = file("data",'r')
    l = f.readline()[:-1].split('\t')
    X = []
    while len(l)>1:
        data = [float(p) for p in l]
        X.append(data)
        l =  f.readline()[:-1].split('\t')
    X = array(X).T



    genelist = []
    f = file("genes",'r')
    l = f.readline()[:-1].split('\t')
    while len(l)>1:
        l =  f.readline()[:-1]
        genelist.append(l)
    return X,Y,genelist





X,Y,genelist = readDataMats()
[trainX,trainY,testX,testY]=divide(X,Y,.50)


#standard usage
#P,Q,W,B = nipals(X,Y,15,10**-10)
#test = X
#Ypred = dot(dot(dot(test,W),B),Q.T)

P,Q,W,B = nipals(trainX,trainY,15,10**-32)

Ypred = dot(dot(dot(testX,W),B),Q.T)


import pylab
zippedYY = sorted([(testY[i][0],Ypred[i][0])  for i in xrange(len(testY)) ])

pylab.plot(zippedYY)

pylab.show()



#gotta remember how to scale by loading vectors
#noreps = set()
#print B
#for j in range(0,15):
#    for i in xrange(len(W[:,j])):
#        if abs(W[i,j]*B[j])>.06: 
#            noreps.add(genelist[i])
#for gen in noreps:print gen.split("\t")[0]



