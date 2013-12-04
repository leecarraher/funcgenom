from random import shuffle,randint
from copy import deepcopy



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
    return [trainX,trainY,testX,testY]
    
def dist(X,Y):
    '''
        Euclidean for now
    '''
    d = 0.0
    for i in xrange(len(X)):
        d=d+(X[i]-Y[i])*(X[i]-Y[i])
    return d
    
def least(q,D,fnc=dist):
    l = fnc(q,D[0])
    lp = 0

    for i in xrange(1,len(D)):
        tmp=fnc(q,D[i])
        if tmp<l:
            lp = i
            l = tmp
    return lp
    


def assignClusters(A,means,clusters):
    '''
        assign to clusters
    '''
    swaps = 0
    newclusters = [list() for i in xrange(k)]
    for i in xrange(len(clusters)):
        for j in xrange(len(clusters[i])):
            arglst = least(A[clusters[i][j]],means)
            newclusters[arglst].append(clusters[i][j])
            swaps += int(arglst != i)
    return newclusters,swaps


def kmeansUpdate(A,clusters,dim):
    '''
        update means
    '''
    means  = []
    for cluster in clusters:
        mean=[0.0 for k in xrange(dim)]
        l = len(cluster)
        for point in cluster:
            for d in xrange(dim):
                mean[d] = mean[d]+(A[point][d] / float(l))
        means.append(mean)
    return means

def kmeans(A,k,dim,maxiters = 1000):
    #some data storage structures
    R = range(len(A))
    shuffle(R)
    clusters = []
    part = len(R)/k
    for i in xrange(k):
        clusters.append( R[i*part:(i+1)*(part)])
    means = kmeansUpdate(A,clusters,dim)
    clusters,swaps = assignClusters(A,means,clusters)
    while swaps>2 and maxiters>0:
        maxiters-=1
        means = kmeansUpdate(A,clusters,dim)
        clusters,swaps = assignClusters(A,means,clusters)
        print "swaps = ",swaps
    return means, swaps
    
    
    
def readCDTFile(filename):
    '''
        read a cdt file as output by genomics portals
        assums genelists are the same, but not neccessarily
        the same order
    '''
    f = file(filename,'r')
    gene = []
    expression = []
    raw = f.readline()[:-1].split("\t")
    
    #first 4 idx are labels
    samples =[k for k in raw[4:]]
    #eweight garbage
    raw = f.readline()
    
    #first gene expression sample
    raw = f.readline()[:-1].split("\t")

    while len(raw)>1:
        gene.append(raw[1]+raw[2])
        expression.append([float(k) for k in raw[4:]])
        #centroids.append([float(k) for k in raw[-2:]])
        raw = f.readline()[:-1].split("\t")
        
    #sort expressions by gene name so lists are aligned
    gene,expression = zip(*sorted(zip(gene,expression)))
        
    return samples,gene,expression

    
def classifier(centroids,data):
    '''
        Apply nearest centroid classifier to data
        returns nearest centroid idx and distance
    '''
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

def findlabels(Y,tildaY,k,labels):
    '''
        find which cluster has the most of a particular label
        return a mapping from cluster id->label
    '''
    counts = [ [0]*labels for i in range(k)]
    ret = {}
    for i in xrange(len(Y)):
        counts[tildaY[i]][Y[i]]+=1
    for i in xrange(k):
        clu = counts[i]
        mxlbl = clu[0]
        argmx = 0
        for ct in xrange(1,labels):
            if clu[ct]>mxlbl:
                mxlbl=clu[ct]
                argmx = ct
        ret[i] = argmx
    return ret
            

#from numpy import array

#read data
Y,genelist,X = readCDTFile("../train.cdt")

#need to transpose for kmeans
X = zip(*X) #transpose X, numpy .T prevents us from using pypy

Y = [int(y[0])-1 for y in Y]
[trainX,trainY,testX,testY]=divide(X,Y,.50)

dim = len(trainX[0])

labels = 3 # number of labels, comes from dataset
k = labels*2 # k parameter for k-means
means,dif = kmeans(trainX,k,dim)

#assign centroid labels based on max labels from the training set
ret,retvals = classifier(means,trainX)
labelMap = findlabels(trainY,ret,k,labels)

#assign cluster ids to unseen dataset
ret,retvals = classifier(means,testX)
ctequal = 0

#accumulate correct matches
for i in range(len(testY)):
    ctequal += (testY[i]==labelMap[ret[i]])
    
print "Classification Accuracy:",
print float(ctequal)/float(len(testY))
    


    

