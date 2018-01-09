import os,subprocess
from numpy import *
from multiprocessing import Pool,cpu_count,Process,Manager
from ete2 import Tree
import tempfile
import pandas as pd

cores = 40 # change according to your hardware

def phylip_output(m,taxa=[],file='data.phy'):
    "writes a distance matrix into a phylip file"
    if len(taxa) != len(m):
        taxa = xrange(len(m))
        mx = 10
    else:
        mx = max([len(x) for x in taxa])
    f = open(file,'w')
    f.write(str(len(taxa))+'\n')
    for nm, row in zip(taxa,m):
        f.write(str(nm).ljust(mx))
        for cell in row:
            f.write(' '+str(cell))
        f.write('\n')
    f.close()


data = pd.read_csv('conceptwiseSimilarities.csv',index_col=0)

cronbach = pd.read_csv('cronbachValues.csv',header=None,squeeze=True,index_col=0)

taxa = array(cronbach[cronbach>=0.6].index)

concepts = array(data.columns[-40:])

data = data[data.language1.isin(taxa)]
data = data[data.language2.isin(taxa)]

nrMts = []
for c in concepts:
    cMtx = zeros((len(taxa),len(taxa)))
    ix1 = list(pd.match(data.language1,taxa))
    ix2 = list(pd.match(data.language2,taxa))
    cMtx[ix1,ix2] = data[c].values
    cMtx[ix2,ix1] = data[c].values
    nrMts.append(cMtx)

matrices = zeros((len(taxa),len(taxa),len(nrMts)))

for c in xrange(40):
    matrices[:,:,c] = nrMts[c]

minSim = -sqrt(40)-.001
maxSim = (log(40*(39)+1)-1)*sqrt(40)+.001


def aggregate(v):
    v = v[v!=-1]
    if len(v)==0:
        sim = 0.
    else:
        sim = (mean(v)-1)*sqrt(len(v))
    return (maxSim-sim)/(maxSim-minSim)


def weightedAggregate(v,weights):
    w = weights[v!=-1]
    w /= sum(w)
    v = v[v!=-1]
    if len(v)==0:
        sim = 0.
    else:
        sim = (dot(v,w)-1)*sqrt(len(v))
    return (maxSim-sim)/(maxSim-minSim)

def simplex(n):
    x = concatenate([[0],sort(random.rand(n-1)),[1]])
    return x[1:] - x[:n]

bbs = array([simplex(40) for i in xrange(1000)])

bbs = r_[[array([1./40]*40)],bbs]

def ldistNWPV_bbs(l1,l2,bbs=bbs):
    return array([weightedAggregate(matrices[l1,l2],rsp) for rsp in bbs])



def fillPackage(pck):
    out = pck
    for p in out:
        p[2:] = ldistNWPV_bbs(int(p[0]),int(p[1]))
    return out

print "start parallel computing"


# initializing cores cores
pool = Pool(cores)



n = len(matrices[0])
nPairs = n*(n+1)/2 # number of non-identical pairs of languages

# longMatrix is a type double numpy array
# rows correspond to language pairs
# the first two columns are initialized with language numbers
# the last column is initialized with 0.0
longMatrix = zeros((nPairs,len(bbs)+2))
k = 0
for i in range(n):
    for j in range(i+1):
        longMatrix[k,0]=i
        longMatrix[k,1]=j
        k += 1
# split longMatrix into cores sub-matrices, each to be processed
# by a different core
packages = array_split(longMatrix,cores)

# for each k: longMatrix[k,2] is filled with ldistNWPV40(longMatrix[k,0],longMatrix[k,1])
longMatrix = vstack(pool.map(fillPackage,packages))

pool.close()
pool.join()

dists = zeros((n,n,len(bbs)))
for p in longMatrix:
    i = int(p[0])
    j = int(p[1])
    dists[i,j] = dists[j,i] = p[2:]

def fastme(dists,taxa,suffix=''):
    outfile = tempfile.mkstemp(suffix=suffix)[1]
    distsfile = tempfile.mkstemp(suffix=suffix)[1]
    phylip_output(dists,taxa=taxa,file=distsfile)
    p = subprocess.Popen('fastme -I /dev/null -i '+ distsfile +' -o '+ outfile +' -T 1 >/dev/null',shell=True)
    os.waitpid(p.pid, 0)
    f = open(outfile)
    nwk = f.readlines()
    f.close()
    p = subprocess.Popen('rm -f '+ outfile,shell=True)
    p = subprocess.Popen('rm -f '+ distsfile,shell=True)
    nwk = concatenate([array(x.strip(),'c') for x in nwk]).tostring()
    return nwk



packages = array_split(arange(len(bbs)),cores)

manager = Manager()

return_dict = manager.dict()

def njP(idx,pck):
    return_dict[idx] = [fastme(-log(1-dists[:,:,i]),taxa,str(idx)) for i in pck]

jobs = []
for i in xrange(cores):
    p = Process(target=njP,args=(i,packages[i]))
    jobs.append(p)
    p.start()

for p in jobs:
    p.join()

trees = concatenate([return_dict[i] for i in xrange(cores)])

treesRerooted = []
for tr in trees:
    t = Tree(tr)
    t.set_outgroup('LATE_EGYPTIAN')
    treesRerooted.append(t.write(format=1))

treesRerooted[0].write(outfile='afroasiatic.tre')

f = open('afroasiatic-bbs.tree','w')
for t in treesRerooted[1:]:
    f.write(t+'\n')
f.close()

