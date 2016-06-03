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
import datetime

version = '0.1.1'

chmGroup = { #Taken from http://www.ebi.ac.uk/Tools/msa/clustalw2/help/faq.html#24
    'small' : 'AVFPMILW',
    'acidic' : 'DE',
    'basic' : 'RK',
    'rest' : 'STYHCNGQ' # Hydroxyl + sulfhydryl + amine + G
}

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
    'windowScore' : 6.5,
    'mode' : 'window',
    'H' : False,
    'debug' : False,
    'detailed' : False,
    'nucleotide' : False
    }

def parseArguments(version): #Parse arguments


    desc = "Squeeze out relevant blocks from alignments."

    parser = argparse.ArgumentParser(description = desc )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s ' + version
        )

    parser.add_argument(
        'infile',
        nargs = '+',
        type=argparse.FileType('r')
        )

    parser.add_argument(
        '-m', '--mode',
        choices = ['gblocks', 'window'],
        default = defaultParams['mode'],
        help = "Use GBlocks or Window mode."
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
        help = 'Blocks of valid positions should be at least this long. Default: {}'.format(defaultParams['final-length-2']),
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
        help = 'Input format. Accepted: fasta, clustal, phylip, emboss, stockholm. Default: {}'.format(defaultParams['infmt']),
        metavar = 'INFMT'
        )

    parser.add_argument(
        '--outfmt', '-p',
        choices = ['fasta', 'clustal', 'phylip', 'stockholm'],
        default = defaultParams['outfmt'],
        help = 'Output format. Accepted: fasta, clustal, phylip, stockholm. Default: {}'.format(defaultParams['outfmt']),
        metavar = 'OUTFMT'
        )

    parser.add_argument(
        '--debug',
        action = 'store_const',
        const = True,
        default = defaultParams['debug'],
        help = 'Print some information regarding sequence filtering.'
        )

    parser.add_argument(
        '-D',
        action = 'store_const',
        const = True,
        default = defaultParams['detailed'],
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
        default = defaultParams['H'],
        const = True,
        help = 'Create HTML output file.'
        )
    parser.add_argument(
        '-N',
        action = 'store_const',
        default = defaultParams['nucleotide'],
        const = True,
        help = 'Analyse nucleotide sequence. Default is aminoacid.'
        )
    args = parser.parse_args()
    filenames = []
    for files in args.infile:
        filenames.append(files.name)

    params = {
        'conserved':  args.conserved, # Default is 0.5
        'highly-conserved': args.hconserved, # Default is 0.85
        'cont-non-conserved': args.contnoncon, #Default is 8 All longer stretches of non-conserved are discarded.
        'final-length-1': args.finallength1, #Default is 15
        'gaps': args.gaps, #Default is none
        'final-length-2' : args.finallength2,
        'filenames' : filenames,
        'windowScore' : args.windowscore,
        'windowSize' : args.windowsize,
        'mode' : args.mode,
        'infmt': args.infmt,
        'outfmt': args.outfmt,
        'H' : args.H,
        'debug' : args.debug,
        'detailed' : args.D,
        'validRange' : args.V,
        'invalidRange' : args.I,
        'discardedBlocks' : args.X,
        'nucleotide' : args.N
        }
    return params

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

def filterBlocks(filename, params = {}):
    #Read alignment file
    if params == {}:
        params = defaultParams
    else:
        for param in defaultParams.keys():
            if not (param in params):
                params[param] = defaultParams[param]

    alignment = AlignIO.read( open(filename, 'rU'), params['infmt'] )
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
        if not params['nucleotide']:
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



        if not params['nucleotide']:
            heterozygosityGroup = 1
            for key in charGr.keys():
                if key != '-':
                    heterozygosityGroup -= (charGr[key]/numSeq)**2
            heterozygosity = (heterozygosity + heterozygosityGroup)/2

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
        'H' : heterozygosity,
        'TS': score,
        'AS': float(avgScore),
        'S': {1 : status },
        'MC': mostCommon,
        'G': gaps
        })
    if params['mode'] == 'gblocks':
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
    elif params['mode'] == 'window':
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
    for n in range(length):
        freq[floor(info[n]['H']*10)] += 1
        s = alignment[:, n]
        if info[n]['S'][6] == 'V':
            count += 1
        if short: print('{s} |  H: {H} | MC: {MC[1]} | TS: {TS} | AS: {AS:.2} | {S} '.format(s = s,**info[n]), file = sys.stderr)
        else:
            if info[n]['S'][6] == 'V':
                print("{} | {}".format(s,n), file = sys.stderr)
            else:
                print("{} | {} - Deleted".format(s,n), file = sys.stderr)


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

def writeAlign(alignment, metadata, filename, initialJson, info, params):

    detailed = params['detailed']

    outfmt = params['outfmt']
    blocks = metadata['blocks']
    validBlocks = blocks[0]
    invalidBlocks = blocks[1]

    validBlockString = ''

    details = []
    for block in validBlocks:
        details.append(str(block))
    validBlockString = ', '.join(details)


    if params['validRange']:
        with open(params['validRange'], 'w') as o:
            for block in validBlocks:
                o.write( '{}\t{}\n'.format( block[0] , block[1] ) )

    if params['invalidRange']:
        with open(params['invalidRange'], 'w') as o:
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

    if params['H']: #JSON



        seqArray = []
        if finalAlignment:
            for record in finalAlignment:
                seqArray.append({'name': record.id ,'seq': str(record.seq) })

        jsonSeq = json.dumps(seqArray, sort_keys=True, indent=4)
        with open(filename + '.html', 'w') as h:
            h.write(htmlOutput(initialJson, jsonSeq, info, params, filename, False))


    outFilename = filename + '.' + outfmt
    if finalAlignment:
        AlignIO.write([finalAlignment], outFilename, outfmt)
    else:
        print("All positions were removed, try lowering the thresholds.", file = sys.stderr)

    if params['discardedBlocks']:
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
        AlignIO.write([invalidAlignment], params['discardedBlocks'], outfmt)


    return finalAlignment

def htmlOutput(initialJson, seqJson, info, params, infile, CGI = True):
    x = []
    y = []
    p = [] # Present
    for n in range(len(info)):
        x.append(n+1)
        y.append(int(floor(info[n]['H']*10)))
        if info[n]['S'][6] == 'V':
            p.append(0)
        else:
            p.append(params['windowScore'])
    script = """
    <script>
        SCOREPLOT = document.getElementById('scoreplot');
        trace1 = {{
        x: {0} ,
        y: {1},
        type : 'scatter',
        name : 'Score',
        marker: {{
            color : '#2980b9'
        }}
        }};

        trace2 = {{
        x: {0} ,
        y: {2},
        type: 'bar',
        name : 'Discarded',
        marker: {{
            color : '#bdc3c7'
        }}
        }};

        data = [trace1, trace2];
        Plotly.plot( SCOREPLOT, data , {{margin : {{t : 0}} }});
    </script>""".format(x, y, p)

    head = """
<!DOCTYPE html>
<head>
    <meta charset="utf-8">
    <title>Laser Albatross - Martín Basterrechea</title>
    <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">
    <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/grids-responsive.css">

    <link rel="stylesheet" href="http://purecss.io/combo/1.18.13?/css/layouts/blog.css">
    <link rel="stylesheet" href="https://cdn.rawgit.com/wilzbach/msa/master/css/msa.css">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <style>
      .post {
        padding-bottom: 0;
      }
      #scoreplot {
        width: 100%;
        height: 200px;
      }
      .button-green {
        color: white;
        background: #3498db;
        border-radius: 4px;
      }
      .sidebar {
        background: #2980b9;
      }
      .smenubar_alink {
        background: #3498db !important;
      }
      label {
        color: #777;
      }
      .table-center {
        font-size: 11px;
        margin: 40px 0 0 auto;
      }
    </style>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>"""
    html = """
<body>
    <div id="layout" class="pure-g">
        <div class="sidebar pure-u-1 pure-u-md-1-4">
          <div class="header">
            <h1 class="brand-title">
              Laser Albatross
            </h1>
            <nav class="nav">
              <ul class="nav-list">
                <li class="nav-item">
                  <a class="pure-button" href="/">Home</a>
                </li>
                <li class="nav-item">
                  <a class="pure-button" href="/about.html">About</a>
                </li>
              </ul>
            </nav>
            table@
          </div>
        </div>
        <div class="content pure-u-1 pure-u-md-3-4">
          <div class="posts">
            <section class="post">
              <header class="post-header">
                <h2 class="content-subhead">Output alignment <a class="pure-button button-green" href="@download">Download</a> </h2>
              </header>
              <div class="post-description">
                <div id='finalDiv'></div>
              </div>
            </section>
            <section class="post">
              <header class="post-header">
                <h2 class="content-subhead">Score plot</h2>
              </header>
              <div class="post-description">
                <div id="scoreplot"></div>
              </div>
            </section>
            <section class="post">
              <header class="post-header">
                <h2 class="content-subhead">Original alignment</h2>
              </header>
              <div class="post-description">
                <div id='initialDiv'></div>
              </div>
            </section>

          </div>



        </div>
    </div>
    <script src="https://cdn.bio.sh/msa/0.4/msa.min.gz.js"></script>

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
        },
        colorscheme : {'scheme' : 'clustal2'}
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
        },
        colorscheme : {'scheme' : 'clustal2'}
      };
      var f = new msa.msa(opts);

      renderMSA(f);

    </script>
    plotscript
  </body>
</html>"""

    if params['mode'] == 'window':
        paramkeys = ['time', 'input file', 'infmt', 'mode', 'gaps', 'windowScore', 'windowSize']
    else:
        paramkeys = ['time', 'input file', 'infmt', 'mode', 'gaps', 'conserved', 'highly-conserved']

    table = '<div><table class="pure-table pure-table-horizontal table-center">\n<tr><th>Parameter</th><th>Value</th></tr>'
    now = datetime.datetime.now()
    for param in paramkeys:
        if param == 'time':
            table+= '<tr><td>time:</td><td>' + now.strftime("%Y-%m-%d %H:%M") + '</td>\n'
        elif param == 'input file':
            #table+= '<tr><td>input file:</td><td>' + infile.split('/')[-1].split('-')[-2] + '</td>\n'
            table+= '<tr><td>input file:</td><td>' + infile + '</td>\n'
        else:
            table+= '<tr><td>' + param + '</td><td>' + str(params[param]) + '</td>\n'
    table += '</table></div>'
    html = head + html.replace("sequences2@", initialJson).replace("sequences@", seqJson).replace("plotscript", script).replace("table@", table)
    if not CGI:
        html = html.replace("@download", infile.split('/')[-1] + '.' + params['outfmt'])
    return head + html.replace("sequences2@", initialJson).replace("sequences@", seqJson).replace("plotscript", script).replace("table@", table)


def getCGI(filename, params = {}): #Return HTML output for CGI script
    alignment, info = filterBlocks(filename, params)
    #Delete file
    metadata = calculateMetadata(alignment, info)

    initialJson = ''
    validSeq = Seq(metadata['validString'])
    validSeqRecord = SeqRecord(seq = validSeq, id = 'Vs', name = 'Valid Blocks')
    scoreSeq = Seq(metadata['scoreString'])
    scoreSeqRecord = SeqRecord(seq = scoreSeq, id = 'Sc', name = 'Heterozygosity Score')
    jsonAlignment = MultipleSeqAlignment([validSeqRecord, scoreSeqRecord])
    for record in alignment:
        jsonAlignment.append(record)
    initialJson = getInitialJson(jsonAlignment)

    validBlocks = metadata['blocks'][0]
    invalidBlocks = metadata['blocks'][1]

    printAlign(alignment, info)

    finalAlignment = False
    for block in validBlocks:
        if not finalAlignment:
            finalAlignment = alignment[:, block[0]-1:block[1]]
        else:
            finalAlignment += alignment[:, block[0]-1:block[1]]

    seqArray = []
    if finalAlignment:
        for record in finalAlignment:
            seqArray.append({'name': record.id ,'seq': str(record.seq) })

    jsonSeq = json.dumps(seqArray, sort_keys=True, indent=4)

    outfmt = params['outfmt']
    outfile = filename.split('/')[-1].split('.')[:-1]
    outfile = 'tmp/' + '.'.join(outfile) + '-out' + '.' + outfmt


    AlignIO.write([finalAlignment], outfile, outfmt)

    return(htmlOutput(initialJson, jsonSeq, info, params, filename), outfile)


def main():

    params = parseArguments(version)

    for filename in params['filenames']:
        alignment, info = filterBlocks(filename, params)
        metadata = calculateMetadata(alignment, info)

        initialJson = ''
        if params['H']:
            # Generate JSON for HTML output, add the "valid blocks" sequence to the initial one.
            validSeq = Seq(metadata['validString'])
            validSeqRecord = SeqRecord(seq = validSeq, id = 'Valid blocks', name = 'Valid Blocks')
            scoreSeq = Seq(metadata['scoreString'])
            scoreSeqRecord = SeqRecord(seq = scoreSeq, id = 'Score', name = 'Heterozygosity Score')
            jsonAlignment = MultipleSeqAlignment([validSeqRecord, scoreSeqRecord])
            for record in alignment:
                jsonAlignment.append(record)
            initialJson = getInitialJson(jsonAlignment)

        if params['debug']: printAlign(alignment, info)

        blocks = metadata['blocks']

        outfile = filename.split('/')[-1].split('.')[:-1]
        outfile = '.'.join(outfile) + '-out'

        writeAlign(alignment, metadata, outfile, initialJson, info, params)

if __name__ == "__main__":
    main()
