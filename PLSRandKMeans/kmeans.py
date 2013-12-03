from random import shuffle,randint
from copy import deepcopy


'''
    Euclidean for now
'''
def dist(X,Y):
    d = 0.0
    for i in xrange(len(X)):
        d=d+(X[i]-Y[i])*(X[i]-Y[i])
    return d**.5
    
def least(q,D,fnc=dist):
    l = 1000000.0
    for i in xrange(len(D)):
        if fnc(q,D[i])<l:
            lp = i
            l=fnc(q,D[i])
    return lp
    

'''
    assign to clusters
'''
def assignClusters(A,means):
    clusters = [[] for i in xrange(len(means))]
    for i in xrange(len(A)):
        clusters[least(A[i],means)].append(A[i])
    return clusters

'''
    update means
'''
def kmeansUpdate(clusters,means,dim):
    for i in xrange(len(clusters)):
        means[i]=[0.0 for k in xrange(dim)]
        for j in xrange(len(clusters[i])):
            l = len(clusters[i][j])
            for k in xrange(l):
                #print means[i]
                means[i][k] = means[i][k]+clusters[i][j][k]
            for k in xrange(l):
                means[i][k]= means[i][k] / float(l)
    return means

def kmeans(A,k,dim,eps = .000001):
    R = range(len(A))
    shuffle(R)
    means = [A[r] for r in R[:k]]
    print means
    prev = deepcopy(means)
    clusters = assignClusters(A,means)
    means = kmeansUpdate(clusters,means,dim)
    dif = sum(map(dist,prev, means))
    print "error = ",dif
    while dif>eps:
        prev = deepcopy(means)
        clusters = assignClusters(A,means)
        means = kmeansUpdate(clusters,means,dim)
        dif = sum(map(dist,prev, means))
        print "error = ",dif
    return means, dif
    


def readDataMats():
    from numpy import array
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
    
def classifier(centroids,data):
    ret = [0]*len(data)
    retvals = [0.0]*len(data)
    for i in xrange(len(data)):
        mindist = dist(data[i],centroids[0])
        argmindist = 0
        for j in xrange(1,len(centroids)):
            if dist(data[i],centroids[j]) < mindist:
                mindist = dist(data[i],centroids[j])
                argmindist = j
        ret[i] = argmindist
        retvals[i] = mindist
    return ret,retvals




    
X,Y,genelist = readDataMats()
dim = len(X[0])
k = 3
means,dif = kmeans(X,k,dim)
ret,retvals = classifier(means,X)
counts = [0,0,0]
for i in range(len(Y)):
    counts[ret[i]]+=Y[i][0]
    print ret[i]
    
print counts



#random set of data in 2d
#k=3
#dim = 100
#npts = 20000
#A = [[randint(0,100)  for d in xrange(dim)] for i in xrange(npts)]
#means,dif = kmeans(A,k,dim)
#print means


    

