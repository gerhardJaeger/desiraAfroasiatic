from numpy import *
import pandas as pd
from Bio import pairwise2
from scipy import stats
from multiprocessing import Process,Manager
ncores=20 # change according to your hardware

data = pd.read_csv('afroasiaticMatrix.csv',index_col=0)
concepts = array(data.columns[-40:])


# exclude reconstructed languages
protoLanguages = [l for l in data.index if 'PROTO' in l]

data = data[~data.index.isin(protoLanguages)]

for c in concepts:
    data[c] = [str(x) for x in data[c].values]


################################################
# Infrastructure for weighted sequence alignment
################################################

f = open('sounds41.txt')
sounds = array([x.strip() for x in f.readlines()])
f.close()


f = open('pmi-world.txt','r')
l = f.readlines()
f.close()
logOdds = array([x.strip().split() for x in l],double)



lodict = dict()
for i in xrange(len(sounds)):
    for j in xrange(len(sounds)):
        lodict[sounds[i],sounds[j]] = logOdds[i,j]

f = open('gapPenalties.txt')
gp1,gp2 = array([x.strip() for x in f.readlines()],double)
f.close()


def sscore(a,b,lodict,gp1,gp2):
    """a,b: ASJP strings
    lodict: logodds dictionary
    gp1,gp2: gap penalties
    return PMI score of a/b
    """
    out = pairwise2.align.globalds(a,b,lodict,gp1,gp2)
    if len(out)==0: return nan
    return out[0][2]



def scoreNW(x,y,lodict,gp1,gp2):
    """x,y: sequences of ASJP strings, separated by '-'
    lodict: logodds dictionary
    gp1,g2: gap penalties
    returns maximal PMI score for the Cartesian product of x and y"""
    if '0' in [x,y]: return nan
    x1=x.split('-')
    y1=y.split('-')
    return max([sscore(xx,yy,lodict,gp1,gp2) for xx in x1 for yy in y1])



def pmiSims40(l1,l2,data=data):
    """data: dataframe with languages as rows
             and concepts as columns
       l1,l2: two languages from data.index
       returns a 40-items vector of calibrated PMI similarities
       for the words from l1/l2 for the 40 concepts"""
    l1List = data[concepts].ix[l1].values
    l2List = data[concepts].ix[l2].values
    simMtr = array([[scoreNW(x,y,lodict,gp1,gp2) for x in l1List]
                          for y in l2List])
    dg = diag(simMtr).copy()
    fill_diagonal(simMtr,nan)
    cmpr = simMtr[isnan(simMtr)==False]
    ranks = array([stats.gmean(1.+arange(sum(cmpr>x),1.+sum(cmpr>=x)))
                   for x in dg],double)
    stc = -log(ranks/(1+len(cmpr)))
    stc[isnan(dg)] = -1
    return stc


taxa = array(data.index)


lpairs = array([(l1,l2)
                for i,l1 in enumerate(taxa)
                for j,l2 in enumerate(taxa)
                if i<j])

packages = array_split(lpairs,ncores)


manager = Manager()
return_dict = manager.dict()

def doWork(i,pck):
    return_dict[i] = array([pmiSims40(l1,l2,data) for l1,l2 in pck])

jobs = []
for i,pck in enumerate(packages):
    p = Process(target=doWork,args=(i,pck))
    p.start()
    jobs.append(p)

for p in jobs:
    p.join()

results = concatenate([return_dict[i] for i in xrange(ncores)])

output = pd.DataFrame(lpairs)

output.columns = ['language1','language2']

for i,c in enumerate(concepts):
    output[c] = results[:,i]

output.to_csv('conceptwiseSimilarities.csv')
