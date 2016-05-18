#Author: Martin Basterrechea

from Bio import AlignIO
from Bio.Align import AlignInfo
from math import ceil
from Bio.SubsMat import MatrixInfo
from sys import argv
import argparse

version = '0.0.3'

def parseArguments(version): #Parse arguments
    defaultParams = {
        'conserved':  0.5, # Default is 0.5
        'highly-conserved': 0.85, # Default is 0.85
        'cont-non-conserved': 8, #Default is 8 All longer stretches of non-conserved are discarded.
        'final-length-1': 15, #Default is 15
        'gaps': 'none', #Default is none
        'final-length-2' : 10,
        'infmt': 'fasta',
        'outfmt': 'fasta'
        }

    desc = "Squeeze out relevant blocks from alignments"

    parser = argparse.ArgumentParser(description = desc )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s ' + version
        )

    parser.add_argument(
        'infile',
        type=argparse.FileType('r')
        )

    parser.add_argument(
        '-o', '--outfile',
        default = False,
        help = 'Name of file containing the output alignment.',
        metavar = 'OUT'
        )

    parser.add_argument(
        '-c', '--conserved',
        type = float,
        default = defaultParams['conserved'],
        help = 'A position will be considered conserved if this fraction of sequences plus one have the same value. Default: 0.5',
        metavar = 'CONS'
        )

    parser.add_argument(
        '-d', '--hconserved',
        type = float,
        default = defaultParams['highly-conserved'],
        help = 'A position will be considered highly conserved if this fraction of sequences plus one have the same value. Default: 0.85',
        metavar = 'HCONS'
        )

    parser.add_argument(
        '-e', '--contnoncon',
        type = int,
        default = defaultParams['cont-non-conserved'],
        help = 'Blocks of non-conserved positions longer than this will be removed. Default: 8',
        metavar = 'CONTNON'
        )

    parser.add_argument(
        '-f', '--finallength1',
        type = int,
        default = defaultParams['final-length-1'],
        help = 'Blocks of valid positions should be at least this long before removing gaps. Default: 15',
        metavar = 'FL1'
        )

    parser.add_argument(
        '-g', '--gaps',
        choices = ['all', 'none'],
        default = defaultParams['gaps'],
        help = "Allow gaps. Default: none."
        )

    parser.add_argument(
        '-i', '--finallength2',
        type = int,
        default = defaultParams['final-length-2'],
        help = 'Blocks of valid positions should be at least this long. Default: 10',
        metavar = 'FL2'
        )

    parser.add_argument(
        '--infmt', '-j',
        choices = ['fasta', 'clustal', 'phylip', 'emboss', 'stockholm'],
        default = defaultParams['infmt'],
        help = 'Input format. Default: fasta'
        )

    parser.add_argument(
        '--outfmt', '-p',
        choices = ['fasta', 'clustal', 'phylip', 'stockholm'],
        default = defaultParams['outfmt'],
        help = 'Output format. Default: fasta'
        )

    parser.add_argument(
        '--debug',
        action = 'store_const',
        const = True,
        default = False,
        help = 'Print some information regarding sequence filtering.'
        )

    parser.add_argument(
        '-D',
        action = 'store_const',
        const = True,
        default = False,
        help = "Print accepted blocks positions on output's sequence description. Available for fasta and stockholm formats."
        )

    parser.add_argument(
        '-V',
        default = False,
        help = 'Print valid block ranges to VALIDFILE',
        metavar = 'VALIDFILE'
    )

    parser.add_argument(
        '-I',
        default = False,
        help = 'Print discarded block ranges to INVALIDFILE',
        metavar = 'INVALIDFILE'
    )
    args = parser.parse_args()
    filename = args.infile.name


    params = {
        'conserved':  args.conserved, # Default is 0.5
        'highly-conserved': args.hconserved, # Default is 0.85
        'cont-non-conserved': args.contnoncon, #Default is 8 All longer stretches of non-conserved are discarded.
        'final-length-1': args.finallength1, #Default is 15
        'gaps': args.gaps, #Default is none
        'final-length-2' : args.finallength2,
        'filename' : filename
        }
    return(params, args)

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

def filterBlocks( params ):
    #Read alignment file
    alignment = AlignIO.read(open( params['filename'] ), "fasta")

    length = alignment.get_alignment_length()

    #Step 1 in algorithm
    conserved = (len(alignment) * params['conserved']) + 1
    hConserved = len(alignment) * params['highly-conserved']
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

        status2 = 'X' # Non conserved
        if mostCommon[1] >= conserved:
            status2 = 'C' # Conserved
        if mostCommon[1] >= hConserved:
            status2 = 'H' # Highly conserved

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
        'S': {1 : status, 1.5 : status2 },
        'MC': mostCommon,
        'G': gaps
        })

    # Step 2

    count = 0
    for i in range(len(info)):
        n = info[i]
        if n['S'][1] == 'X':
            count += 1
            if count == params['cont-non-conserved'] + 1:
                for j in range(count):
                    info[i-j]['S'][2] = 'X'
            elif count > params['cont-non-conserved'] + 1:
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
                        p['S'][5] = 'X'   #Possible bug in GBlocks, it doesnt go back
                        a = 0
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
    short = True #Short description, long is for debug
    for n in range(length):
        s = alignment[:, n]
        if short: print('{s} |  MC: {MC[1]} | TS: {TS} | AS: {AS} | {S} '.format(s = s,**info[n]))
        else:
            if info[n]['S'][6] == 'V':
                print("{} | {}".format(s,n))
            else:
                print("{} | {} - Deleted".format(s,n))

def calculateValidBlocks(alignment, info): #using index 1
    index = 1 #change to 0 to use index 0.
    validBlocks = []
    invalidBlocks = []
    firstValue = 0
    firstValueInvalid = 0
    valid = False
    length = alignment.get_alignment_length()
    for n in range(length):
        if info[n]['S'][6] == 'V':
            if not valid:
                valid = True
                firstValue = n
                invalidBlocks.append( ( firstValueInvalid + index, n - 1 + index) )
            if valid and n == length - 1:
                validBlocks.append( ( firstValue + index , n + index ) )
        if info[n]['S'][6] == 'X':
            if valid:
                valid = False
                validBlocks.append( ( firstValue + index , n - 1 + index ) )
                firstValueInvalid = n
            if not valid and n == length - 1:
                invalidBlocks.append( ( firstValueInvalid + index , n + index ) )

    return ( validBlocks, invalidBlocks )



def compareGblocks(alignment): #This function compares the alignment to GBlocks output
    gblocks = AlignIO.read(open("gblocksfiles/" + params['filename'] + '-gb'), "fasta")

    gLength = gblocks.get_alignment_length()
    lLength = alignment.get_alignment_length()
    if gLength != lLength:
        print('Sequence length is not the same for filename {}, LA is {} and GB is {}'.format(params['filename'], lLength, gLength))

    for n in range(max(gLength, lLength)):
        try:
            g = gblocks[:, n]
        except IndexError:
            g = 'NA'
        try:
            l = alignment[:, n]
        except IndexError:
            l = 'NA'
        if g != l:
            print('{}: Column {} differs: G: {} and L: {}'.format(params['filename'], n, g, l))

def writeAlign(alignment, blocks, filename):

    detailed = args.D

    outfmt = args.outfmt

    validBlocks = blocks[0]
    invalidBlocks = blocks[1]

    validBlockString = ''

    details = []
    for block in validBlocks:
        details.append(str(block))
    validBlockString = ', '.join(details)


    if args.V:
        with open(args.V, 'w') as o:
            for block in validBlocks:
                o.write( '{}\t{}\n'.format( block[0] , block[1] ) )

    if args.I:
        with open(args.I, 'w') as o:
            for block in invalidBlocks:
                o.write( '{}\t{}\n'.format( block[0] , block[1] ) )


    finalAlignment = False
    for block in validBlocks:
        if not finalAlignment:
            finalAlignment = alignment[:, block[0]:block[1]]
        else:
            finalAlignment += alignment[:, block[0]:block[1]]
    if detailed:
        for record in finalAlignment:
            record.description = validBlockString
    AlignIO.write([finalAlignment], filename, outfmt)
    return finalAlignment

params, args = parseArguments(version)

alignment, info = filterBlocks(params)

if (args.debug): printAlign(alignment, info)
blocks = calculateValidBlocks(alignment, info)

if not args.outfile:
    outfile = params['filename'].split('/')[-1].split('.')[:-1]
    outfile = '.'.join(outfile) + '.' + args.outfmt

writeAlign(alignment, blocks, outfile)
