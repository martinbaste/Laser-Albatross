#Author: Martin Basterrechea

from Bio import AlignIO
from Bio.Align import AlignInfo
from math import ceil, floor
from Bio.SubsMat import MatrixInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import sys
import argparse
import json

version = '0.1.1'

chmGroup = { #Taken from http://www.ebi.ac.uk/Tools/msa/clustalw2/help/faq.html#24
    'small' : 'AVFPMILW',
    'acidic' : 'DE',
    'basic' : 'RK',
    'rest' : 'STYHCNGQ' # Hydroxyl + sulfhydryl + amine + G
}

def parseArguments(version): #Parse arguments
    defaultParams = {
        'conserved':  0.5, # Default is 0.5
        'highly-conserved': 0.85, # Default is 0.85
        'cont-non-conserved': 8, #Default is 8 All longer stretches of non-conserved are discarded.
        'final-length-1': 15, #Default is 15
        'gaps': 'none', #Default is none
        'final-length-2' : 10,
        'infmt': 'fasta',
        'outfmt': 'fasta',
        'windowSize' : 5,
        'windowScore' : 7,
        }

    desc = "Squeeze out relevant blocks from alignments."

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
        '-m', '--mode',
        choices = ['gblocks', 'window'],
        default = 'gblocks',
        help = "Use GBlocks or Window mode."
        )

    parser.add_argument(
        '-o', '--outfile',
        default = False,
        help = 'Name of file containing the output alignment.',
        metavar = 'FILE'
        )

    parser.add_argument(
        '-c', '--conserved',
        type = float,
        default = defaultParams['conserved'],
        help = 'A position will be considered conserved if this fraction of sequences plus one have the same value. Default: 0.5',
        metavar = 'FLOAT'
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
        choices = ['all', 'partial', 'none'],
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
        '-w', '--windowsize',
        type = int,
        default = defaultParams['windowSize'],
        help = 'Window size (for window mode). Default: {}'.format(defaultParams['windowSize']),
        metavar = 'VALUE'
        )

    parser.add_argument(
        '-x', '--windowscore',
        type = float,
        default = defaultParams['windowScore'],
        help = 'Window score above this threshold will be taken out (for window mode). Default: {}'.format(defaultParams['windowScore']),
        metavar = 'VALUE'
        )

    parser.add_argument(
        '--infmt', '-j',
        choices = ['fasta', 'clustal', 'phylip', 'emboss', 'stockholm'],
        default = defaultParams['infmt'],
        help = 'Input format. Accepted: fasta, clustal, phylip, emboss, stockholm. Default: fasta',
        metavar = 'INFMT'
        )

    parser.add_argument(
        '--outfmt', '-p',
        choices = ['fasta', 'clustal', 'phylip', 'stockholm'],
        default = defaultParams['outfmt'],
        help = 'Output format. Accepted: fasta, clustal, phylip, stockholm. Default: fasta',
        metavar = 'OUTFMT'
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

    parser.add_argument(
        '-X',
        default = False,
        help = 'Print discarded blocks to FILE',
        metavar = 'FILE'
        )
    parser.add_argument(
        '-H',
        action = 'store_const',
        default = False,
        const = True,
        help = 'Create HTML output file.'
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
        'filename' : filename,
        'windowScore' : args.windowscore,
        'windowSize' : args.windowsize,
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

def filterBlocks(params):
    #Read alignment file

    alignment = AlignIO.read(open( params['filename'] ), args.infmt )
    #except ValueError as err:
    #    print('Error while opening the file: {}'.format(err), file = sys.stderr)
    #    exit()
    length = alignment.get_alignment_length()

    #Step 1 in algorithm
    conserved = (len(alignment) * params['conserved']) + 1
    hConserved = len(alignment) * params['highly-conserved']
    blosum = MatrixInfo.blosum62
    info = []
    allGaps = params['gaps'] == 'all' #Some things have to be calculated
                                        #differently if there are all gaps
    for n in range(length):
        s = alignment[:, n]
        chars = {}
        chars['-'] = 0 #For gap counting
        charGr = { 's': 0, 'a' : 0, 'b' : 0, 'r' : 0, 'o' : 0}
        for i in s:
            if i in chars:
                chars[i] += 1
            else:
                chars[i] = 1
            if i in chmGroup['small']:
                charGr['s'] += 1
            elif i in chmGroup['acidic']:
                charGr['a'] += 1
            elif i in chmGroup['basic']:
                charGr['b'] += 1
            elif i in chmGroup['rest']:
                charGr['r'] += 1
            else:
                if not (allGaps and i == '-'):
                    charGr['o'] += 1

        mostCommon = ('-', 0)
        for key in chars.keys():
            if chars[key] > mostCommon[1]:
                if not (allGaps and key == '-'):
                    mostCommon = (key, chars[key])

        # In this case gaps don't count towards the "conserved" requirement, so
        # if there are gaps, the required frequency for an aminoacid to be
        # "conserved" is less.
        numSeq = len(alignment)
        if allGaps:
            numSeq -= chars['-']
        if allGaps:
            conserved = (numSeq * params['conserved']) + 1
            hConserved = numSeq * params['highly-conserved']
        status = 'X' # Non conserved
        if mostCommon[1] >= conserved:
            status = 'C' # Conserved
        if mostCommon[1] >= hConserved:
            status = 'H' # Highly conserved

        heterozygosity = 1
        for key in chars.keys():
            if key != '-':
                heterozygosity -= (chars[key]/numSeq)**2

        heterozygosityGroup = 1
        for key in charGr.keys():
            if key != '-':
                heterozygosityGroup -= (charGr[key]/numSeq)**2

        #Count gaps
        gaps = chars['-']

        #Only "none" for gaps makes the position count as non-conserved.
        if params['gaps'] == 'none':
            if gaps:
                status = 'X'


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


        info.append({
        'H' : (heterozygosity + heterozygosityGroup)/2,
        'TS': score,
        'AS': float(avgScore),
        'S': {1 : status },
        'MC': mostCommon,
        'G': gaps
        })
    if args.mode == 'gblocks':
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

        #Step 5 #Filter out gaps if none are allowed
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
    elif args.mode == 'window':
        wSize = params['windowSize']
        for n in range(length - wSize):
            scoreSum = 0
            for i in range(wSize):
                scoreSum += info[n+i]['H']
            avgScore = scoreSum / wSize
            if avgScore * 10 >= params['windowScore']:
                for i in range(wSize):
                    info[n+i]['S'][6] = 'X'
            if not 6 in info[n]['S']:
                if params['gaps'] == 'none' and info[n]['G'] != 0:
                    info[n]['S'][6] = 'X'
                else:
                    info[n]['S'][6] = 'V'
            if n == length - wSize - 1:
                for i in range(wSize + 1):
                    if not 6 in info[n+i]['S']:
                        if params['gaps'] == 'none' and info[n+i]['G'] != 0:
                            info[n+i]['S'][6] = 'X'
                        else:
                            info[n+i]['S'][6] = 'V'
            info[n]['H'] = avgScore




    return(alignment, info)

def getInitialJson(alignment):
     seqArray = []
     for record in alignment:
         seqArray.append({'name': record.id ,'seq': str(record.seq) })
     return json.dumps(seqArray, sort_keys=True, indent=4)

def printAlign(alignment, info):
    #Print valid columns
    length = alignment.get_alignment_length()
    short = True #Short description, long is for debug
    count = 0
    freq = [0]*10
    f = open('numbers','w')
    for n in range(length):
        freq[floor(info[n]['H']*10)] += 1
        f.write(str(info[n]['H']) + '\n')
        s = alignment[:, n]
        if info[n]['S'][6] == 'V':
            count += 1
        if short: print('{s} |  H: {H} | MC: {MC[1]} | TS: {TS} | AS: {AS:.2} | {S} '.format(s = s,**info[n]), file = sys.stderr)
        else:
            if info[n]['S'][6] == 'V':
                print("{} | {}".format(s,n), file = sys.stderr)
            else:
                print("{} | {} - Deleted".format(s,n), file = sys.stderr)
    print(count, file = sys.stderr)
    for n in range(len(freq)):
        freq[n] = freq[n]/length
    print(freq)
    f.close()

def calculateMetadata(alignment, info):
    metadata = {}
    #Calculate valid blocks using index 1 and valid string
    index = 1 #change to 0 to use index 0.
    validBlocks = []
    invalidBlocks = []
    firstValue = 0
    firstValueInvalid = -1
    valid = False

    validString = ''
    scoreString = ''

    length = alignment.get_alignment_length()
    for n in range(length):
        scoreString += str(int(floor(info[n]['H']*10)))
        if info[n]['S'][6] == 'V':
            validString += 'X'
            if not valid:
                valid = True
                firstValue = n
                if n != 0:
                    invalidBlocks.append( ( firstValueInvalid + index, n - 1 + index) )
            if valid and n == length - 1:
                validBlocks.append( ( firstValue + index , n + index ) )
        if info[n]['S'][6] == 'X':
            validString += '.'
            if valid:
                valid = False
                validBlocks.append( ( firstValue + index , n - 1 + index ) )
                firstValueInvalid = n
            if not valid and n == length - 1:
                invalidBlocks.append( ( firstValueInvalid + index , n + index ) )
    metadata['blocks'] = ( validBlocks, invalidBlocks )
    metadata['validString'] = validString
    metadata['scoreString'] = scoreString


    return metadata

def writeAlign(alignment, metadata, filename, initialJson):

    detailed = args.D

    outfmt = args.outfmt
    blocks = metadata['blocks']
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
            finalAlignment = alignment[:, block[0]-1:block[1]]
        else:
            finalAlignment += alignment[:, block[0]-1:block[1]]
    if detailed:
        for record in finalAlignment:
            record.description = validBlockString

    if args.H: #JSON



        seqArray = []
        for record in finalAlignment:
            seqArray.append({'name': record.id ,'seq': str(record.seq) })

        jsonSeq = json.dumps(seqArray, sort_keys=True, indent=4)
        with open(filename + '.html', 'w') as h:
            h.write(htmlOutput(initialJson, jsonSeq))


    outFilename = filename + '.' + args.outfmt
    if finalAlignment:
        AlignIO.write([finalAlignment], outFilename, outfmt)
    else:
        print("All positions were removed, try lowering the thresholds.", file = sys.stderr)

    if args.X:
        invalidBlockString = ''

        details = []
        for block in invalidBlocks:
            details.append(str(block))
        invalidBlockString = ', '.join(details)
        invalidAlignment = False
        for block in invalidBlocks:
            if not invalidAlignment:
                invalidAlignment = alignment[:, block[0]:block[1]]
            else:
                invalidAlignment += alignment[:, block[0]:block[1]]
        if detailed:
            for record in invalidAlignment:
                record.description = invalidBlockString
        AlignIO.write([invalidAlignment], args.X, outfmt)


    return finalAlignment

def htmlOutput(initialJson, seqJson):
    html = """<!DOCTYPE html>
    <html>
      <head>
        <meta charset="UTF-8">
        <link type="text/css" rel="stylesheet" href="https://cdn.rawgit.com/wilzbach/msa/master/css/msa.css">
      </head>
      <body>

        <script src="https://cdn.bio.sh/msa/0.4/msa.min.gz.js"></script>


        <div id='initialDiv'></div>
        <div id='finalDiv'></div>

        <script>
          var gffParser = msa.io.gff;
          var xhr = msa.io.xhr;

          var seqs = sequences@
          var initialSeqs = sequences2@
          var rootDiv = document.getElementById('initialDiv');
          /* global rootDiv */
          // set your custom properties

          var opts = {
            seqs : initialSeqs,
            el : rootDiv,
            vis: {
                labelId: false
            }
          };


          // init msa
          var i = new msa.msa(opts);

          renderMSA(i);
          function renderMSA(m) {

              // the menu is independent to the MSA container
              var menuOpts = {vis : {labelId : false} };
              menuOpts.el = document.getElementById('div');
              menuOpts.msa = m;
              menuOpts.menu = "small";
              var defMenu = new msa.menu.defaultmenu(menuOpts);
              m.addView("menu", defMenu);

              // call render at the end to display the whole MSA
              m.render();
          }

          var finalDiv = document.getElementById('finalDiv');

          var opts = {
            seqs : seqs,
            el : finalDiv,
            vis: {
                labelId: false
            }
          };
          var f = new msa.msa(opts);

          renderMSA(f);

        </script>
      </body>
    </html>"""
    return html.replace("sequences2@", initialJson).replace("sequences@", seqJson)


def main():
    global args
    params, args = parseArguments(version)

    alignment, info = filterBlocks(params)
    metadata = calculateMetadata(alignment, info)

    initialJson = ''
    if args.H:
        # Generate JSON for HTML output, add the "valid blocks" sequence to the initial one.
        validSeq = Seq(metadata['validString'])
        validSeqRecord = SeqRecord(seq = validSeq, id = 'Vs', name = 'Valid Blocks')
        scoreSeq = Seq(metadata['scoreString'])
        scoreSeqRecord = SeqRecord(seq = scoreSeq, id = 'Sc', name = 'Heterozygosity Score')
        jsonAlignment = MultipleSeqAlignment([validSeqRecord, scoreSeqRecord])
        for record in alignment:
            jsonAlignment.append(record)
        initialJson = getInitialJson(jsonAlignment)

    if (args.debug): printAlign(alignment, info)

    blocks = metadata['blocks']

    if not args.outfile:
        outfile = params['filename'].split('/')[-1].split('.')[:-1]
        outfile = '.'.join(outfile)

    writeAlign(alignment, metadata, outfile, initialJson)

main()
