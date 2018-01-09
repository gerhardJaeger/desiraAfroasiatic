# reads list-afroasiatic.txt into a matrix
# which is saved as afroasiaticMatric.csv.
# each row corresponds to one language
# the first four columns hold language name, genus, family, and numberOfSpeakers/ExtinctionStatus
# the remaining columns hold the words for the 40 concepts
# diacritics and preceding segments are removed
# if several translations are given, they are combined into one string, separated by '-'
# missing values are indicated as '0'
# entries marked as loans are kept



from numpy import *
from collections import defaultdict
import re


# ASJP is based on a 100-item Swadesh list, but for most languages
# only the 40 most stable concepts are covered. "concepts" holds
# the indices of these 40 concepts within the 100-item list
concepts = array([1,2,3,11,12,18,19,21,22,23,25,28,30,31,34,39,40,41,
                  43,44,47,48,51,53,54,57,58,61,66,72,74,75,77,82,85,
                  86,92,95,96,100])


# mspli splits string s into sub-strings, where each symbol in the
# string seps can serve as separator
def msplit(s,seps):
    r = [s]
    for sp in seps:
        r = concatenate([ri.split(sp) for ri in r])
    return list(r)


f = open('list-afroasiatic.txt')
rawlist = [x.strip() for x in f.readlines() if x.strip()!='']
f.close()


conceptDict = dict(array([x.split() for x in rawlist[2:42]]))

conceptNames = array([conceptDict[str(i)] for i in concepts])

rawlist = [re.sub(r"(, |,)","-",x) for x in rawlist]



lIndices = [n for n in xrange(len(rawlist)) if '{' in rawlist[n]]

lnames = defaultdict(lambda:0)
for n in lIndices:
    language = msplit(rawlist[n],['.','{','|'])[:3]
    d = dict()
    d['family'] = language[1]
    d['genus'] = language[2]
    metadata = rawlist[n+1]
    d['pop'] = re.sub(r"[a-zA-Z]",'',metadata).split()[-1]
    lnames[language[0]] = d


swadeshDict = defaultdict(lambda:'0')
for n in lIndices:
    language = msplit(rawlist[n],['.','{','|'])[0]
    remaining = [i for i in lIndices if i > n]
    if len(remaining)>0:
        nextLIndex = min(remaining)
    else: nextLIndex = -1
    for x in rawlist[n+2:nextLIndex]:
        c = int(double(x.split()[0]))
        if c in concepts:
            x = re.sub(r"\([^\)]*\)","",x)
            word = x.split()[2].strip(',')
            if not word in ['xxx','XXX']:
                if not '???' in word:
                    swadeshDict[language,c] = word#.split(',')[0]


languages = [l.split('{')[0] for l in array(rawlist)[lIndices]]

f = file('afroasiaticMatrix.csv',"w")
f.write('language,family,genus,pop,'+','.join(conceptNames)+'\n')
for l in languages:
    lEdited = re.sub(r"\(","_",l)
    lEdited = re.sub(r"\)","",lEdited)
    s = lEdited+','+lnames[l]['family']+','+lnames[l]['genus']+','+lnames[l]['pop']+','
    for c in concepts:
        word = swadeshDict[l,c]
        word = re.sub(r"~~","~",word)
        word = re.sub(r"\%","",word)
        word = re.sub(r"\*","",word)
        word = re.sub(r"\"","",word)
        word = re.sub(r".~","",word)
        word = re.sub(r"(.)(.)(.)\$",r"\2",word)
        s += word+','
    f.write(s.strip(',')+'\n')
f.close()
