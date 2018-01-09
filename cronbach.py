from numpy import *
import pandas as pd

# this script computes Cronbach's alpha for all languages in the sample

data = pd.read_csv('conceptwiseSimilarities.csv',index_col=0)

concepts = array(data.columns[-40:])

taxa = unique(data[['language1','language2']].values.flatten())

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


def cronbach(x):
    itemwise = sum(apply_along_axis(var,0,x))
    total = var(apply_along_axis(sum,1,x))
    return 1.*len(x)/(len(x)-1)*(1-itemwise/total)



matricesFilled = matrices.copy()

for i in xrange(len(matricesFilled)):
    for j in xrange(len(matricesFilled)):
        x = matricesFilled[i,j]
        fill = mean(x[x>-1])
        matricesFilled[i,j,x==-1] = fill


crList = []
for i in xrange(len(taxa)):
    x = matricesFilled[i]
    x = delete(x,i,0)
    cr = cronbach(x)
    crList += [cr]
crList = pd.Series(crList,index=taxa)

crList.sort_values().to_csv('cronbachValues.csv')
