#Author: Martin Basterrechea

from Bio import AlignIO
from Bio.Align import AlignInfo
from math import ceil
from Bio.SubsMat import MatrixInfo

#Notes: use diversity(genetics) as a value?

filename = "KOG0019.fasta"

#Parameters
params = {
    'conserved':  0.5, # Default is 0.5
    'highly conserved': 0.85, # Default is 0.85
    'cont-non-conserved': 8, #Default is 8 All longer stretches of non-conserved are discarded.
    'final-length-1': 15, #Default is 15
    'gaps': 'none', #Default is none
    'final-length-2' : 10
    }

def allPairs(L):
    L = list(L)
    while L:
        i = L.pop()
        for j in L: yield i, j

def scoreMatch(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def filterBlocks(filename, params):
    #Read alignment file
    alignment = AlignIO.read(open("testfiles/" + filename ), "fasta")

    length = alignment.get_alignment_length()

    #Step 1 in algorithm
    conserved = (len(alignment) * params['conserved']) + 1
    hConserved = len(alignment) * params['highly conserved']
    blosum = MatrixInfo.blosum62
    info = []

    for n in range(length):
        s = alignment[:, n]
        chars = {}
        chars['-'] = 0 #For gap counting
        for i in s:
            if i in chars:
                chars[i] += 1
            else:
                chars[i] = 1

        mostCommon = ('-', 1)
        for key in chars.keys():
            if chars[key] > mostCommon[1]:
                mostCommon = (key, chars[key])
        status = 'X' # Non conserved
        if mostCommon[1] > conserved:
            status = 'C' # Conserved
        if mostCommon[1] > hConserved:
            status = 'H' # Highly conserved


        score = 0
        count = 0
        for p in allPairs(s):
            if not p[0] in '-*' and not p[1] in '-*':
                score += scoreMatch( p, blosum )
                count += 1
        if count == 0:
            avgScore = 0
        else:
            avgScore = score / count
        #Count gaps
        gaps = chars['-']

        info.append({
        'TS': score,
        'AS': avgScore,
        'S': {1 : status },
        'MC': mostCommon,
        'G': gaps
        })

    # Step 2

    count = 0
    for i in range(len(info)):
        n = info[i]
        if n['S'][1] == 'X':
            count += 1
            if count == params['cont-non-conserved']:
                for j in range(count):
                    info[i-j]['S'][2] = 'X'
            elif count > params['cont-non-conserved']:
                n['S'][2] = 'X' # Invalid
            else:
                n['S'][2] = 'V' # Valid

        else:
            n['S'][2] = 'V' # Valid
            count = 0

    #Step 3 First from the beginning to the end, then backwards
    block = False
    for i in range(len(info)):
        n = info[i]
        if n['S'][2] == 'V' and not block:
            if n['S'][1] == 'H':
                n['S'][3] = 'V'
                block = True
            else:
                n['S'][3] = 'X'
        elif n['S'][2] == 'X':
            block = False
            n['S'][3] = 'X'
        else:
            n['S'][3] = 'V'
    block = False
    for i in range(len(info),0, -1):
        n = info[i-1]
        if n['S'][2] == 'V' and not block:
            if n['S'][1] == 'H':
                n['S'][3] = 'V'
                block = True
            else:
                n['S'][3] = 'X'
        elif n['S'][2] == 'X':
            block = False
            n['S'][3] = 'X'
        elif n['S'][2] == 'V' and block:
            n['S'][3] = 'V'

    #Step 4

    count = 0
    for i in range(len(info)):
        n = info[i]
        if n['S'][3] == 'V':
            count += 1
            if count == params['final-length-1']: #Block is longer, save previous positions
                for j in range(count):
                    info[i-j]['S'][4] = 'V'
            elif count > params['final-length-1']:
                n['S'][4] = 'V' # Valid
            else: #Count is less than FL1
                n['S'][4] = 'X' # Invalid

        else:
            n['S'][4] = 'X' # Invalid
            count = 0

    #Step 5
    if params['gaps'] == 'none':
        cutting = True
        for i in range(len(info)):
            n = info[i]
            if n['S'][4] == 'V' and n['G'] > 0 and not cutting:
                n['S'][5] = 'X' #No gaps allowed
                cutting = True
                for j in range(i):
                    p = info[i-j]
                    if p['S'][5] == 'V' and p['S'][1] == 'X': #If previous is valid and non conserved, take it out
                        p['S'][5] = 'X'
                    elif p['S'][5] == 'V' and p['S'][1] != 'X':
                        break #We found the last non-conserved block, break loop
            elif n['S'][4] == 'V' and n['G'] == 0 and not cutting:
                n['S'][5] = 'V'
            elif n['S'][4] == 'V' and cutting:
                if n['S'][1] == 'X':
                    n['S'][5] = 'X'
                elif n['G'] > 0:
                    n['S'][5] = 'X'
                else:
                    n['S'][5] = 'V'
                    cutting = False
            elif n['S'][4] == 'X':
                n['S'][5] = 'X'
    else:
        for n in info:
            n['S'][5] = n['S'][4]


    #Step 6

    count = 0
    for i in range(len(info)):
        n = info[i]
        if n['S'][5] == 'V':
            count += 1
            if count == params['final-length-2']: #Block is longer, save previous positions
                for j in range(count):
                    info[i-j]['S'][6] = 'V'
            elif count > params['final-length-2']:
                n['S'][6] = 'V' # Valid
            else: #Count is less than FL2
                n['S'][6] = 'X' # Invalid

        else:
            n['S'][6] = 'X' # Invalid
            count = 0
    return(alignment, info)


def printAlign(alignment, info):
    #Print valid columns
    length = alignment.get_alignment_length()
    for n in range(length):
        s = alignment[:, n]
        #print('{s} |  MC: {MC[1]} | TS: {TS} | AS: {AS} | {S} '.format(s = s,**info[n]))
        if info[n]['S'][6] == 'V':
            print("{} | {}".format(s,n))
        else:
            print("{} | {} - Deleted".format(s,n))

def calculateValidBlocks(alignment, info): #using index 0
    validBlocks = []
    firstValue = 0
    valid = False
    for n in range(alignment.get_alignment_length()):
        if info[n]['S'][6] == 'V':
            if not valid:
                valid = True
                firstValue = n
        if info[n]['S'][6] == 'X':
            if valid:
                valid = False
                validBlocks.append( ( firstValue, n) )
    return validBlocks



def compareGblocks(alignment, filename):
    gblocks = AlignIO.read(open("gblocksfiles/" + filename + '-gb'), "fasta")

    gLength = gblocks.get_alignment_length()
    lLength = alignment.get_alignment_length()
    if gLength != lLength:
        print('Sequence lenght is not the same for filename {}, LA is {} and GB is {}'.format(filename, lLength, gLength))

    for n in range(gLength):
        g = gblocks[:, n]
        l = alignment[:, n]
        if g != l:
            print('{}: Column {} is not the same: G: {} and L: {}'.format(filename, n, g, l))

def writeAlign(alignment, validBlocks, filename):
    finalAlignment = False
    for block in validBlocks:
        if not finalAlignment:
            finalAlignment = alignment[:, block[0]:block[1]]
        else:
            finalAlignment += alignment[:, block[0]:block[1]]
    #AlignIO.write([finalAlignment], "my_example.phy", "fasta")
    return finalAlignment



alignment, info = filterBlocks(filename, params)
#printAlign(alignment, info)
validBlocks = calculateValidBlocks(alignment, info)

finalAln = writeAlign(alignment, validBlocks, filename)

compareGblocks(finalAln, filename)