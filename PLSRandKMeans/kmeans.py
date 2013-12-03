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
    return [array(trainX),array(trainY),array(testX),array(testY)]


'''
    Euclidean for now
'''
def dist(X,Y):
    d = 0.0
    if len(X)!=len(Y):
        print "error incompatible dimensions"
    for i in xrange(len(X)):
        d=d+(X[i]-Y[i])*(X[i]-Y[i])
    return d**.5
    
def least(q,D,fnc=dist):
    l = fnc(q,D[0])
    lp = 0
    for i in xrange(1,len(D)):
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
            

from numpy import array

#read data
Y,genelist,X = readCDTFile("train.cdt")

#need to transpose for kmeans
Y = array([int(y[0])-1 for y in Y]).T
X = array(X).T

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
    


    

