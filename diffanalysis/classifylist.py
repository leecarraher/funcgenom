#euclidean distance
def dist(a,b):return sum((a[i]-b[i])**2 for i in xrange(len(a)))**.5


def classifier(centroids,data):
    '''
        Apply nearest centroid classifier to data
        returns nearest centroid idx and distance
    '''
    ret = [0]*len(data)
    retvals = [0.0]*len(data)
    #loop over vectors
    for i in xrange(len(data)):
        #Argmin over centroids
        mindist = dist(data[i],centroids[0])
        argmindist = 0
        for j in xrange(1,len(centroids)):
            if dist(data[i],centroids[j]) < mindist:
                mindist = dist(data[i],centroids[j])
                argmindist = j
        ret[i] = argmindist
        retvals[i] = mindist
    return ret,retvals

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


'''
    Use expert identified gene centroid labels for grade 1 and 3 primary breast 
    cancer tumors to classify grade 2 identified gene samples as either grade 1
    or 3
'''

#read training set data    
samples,gene,expression = readCDTFile("train.cdt")

#compute class centroids
from pylab import array
expression =  array(expression).T

#only two centroids, split them by samplename 1GSMXXX or 3GSMXXX
expression1 = [0.0]*len(gene)
expression3 = [0.0]*len(gene)
for e in xrange(len(expression)):
    if samples[e][0] == '1':
        for j in xrange(len(expression[e])):
            expression1[j]+=expression[e][j]/float(len(samples))
    else:
        for j in xrange(len(expression[e])):
            expression3[j]+=expression[e][j]/float(len(samples))

#list of centroids
centroids = [expression1,expression3]

#read test set data
samples,gene,expression = readCDTFile("test.cdt")
#use centroids to classify
classes,vals = classifier(array(centroids),array(expression).T)
for i in range(len(classes)):
    #groups are 1&3 so apply affine xfrm
    print str(samples[i]) + '->' + str(classes[i]*2+1) 


