from multiprocessing import Pool, shared_memory
import gzip
from optparse import OptionParser
import sys, gzip, os
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-p","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
parser.add_option("-t","--threadN",action = 'store',type = 'int',dest = 'threadN',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None or opt.threadN == None:
    print('Basic usage')
    print('')
    print('     python find_teloRead.py -t 24 -p Test.HiFi.01C')
    print('')
    sys.exit()

batchN = opt.threadN
batchIDX_LIST = range(0, batchN)
fileName = opt.INPUT
tmpDir = 'tmp'
MINIMUM_COUNT = 50

if not os.path.exists(tmpDir):
    os.makedirs(tmpDir)

def find_teloMotif(sequence):
    motif = 'AACCCT'
    motif_rev = 'AGGGTT'
    motifLen = len(motif)

    result_LIST = []
    seqLen = len(sequence)
    for idx in range(seqLen - (motifLen - 1)):
        subSequence = sequence[idx:idx + motifLen]
        if subSequence == motif:
            result_LIST += [(idx + 1, '+')]
        elif subSequence == motif_rev:
            result_LIST += [(idx + 1, '-')]

    return result_LIST

def read_fastq(batchIDX):
    fin  = gzip.open(f'{fileName}.fastq.gz', 'rt')
    fout = gzip.open(f'{tmpDir}/{fileName}.tmp_{batchIDX:03d}.fastq.gz' , 'wt')

    data_LIST = []
    qualityCount_LIST = [0] * 200
    for lineIDX, line in enumerate(fin):
        if lineIDX > 10000: break
        if (lineIDX>>2)%batchN != batchIDX: continue
        
        if lineIDX%4 == 0:
            readName = line.rstrip('\n')
        elif lineIDX%4 == 1:
            sequence = line.rstrip('\n')
        elif lineIDX%4 == 2:
            pass
        elif lineIDX%4 == 3:
            qualities = line.rstrip('\n')

            telo_LIST = find_teloMotif(sequence)

            if len(telo_LIST) > MINIMUM_COUNT:
                fout.write(readName + '\n')
                fout.write(sequence + '\n')
                fout.write('+' + '\n')
                fout.write(qualities + '\n')
    fin.close()
    fout.close()

with Pool(processes=batchN) as pool:
    result_LIST = pool.map(read_fastq, batchIDX_LIST)

fout  = gzip.open(f'{fileName}.telo.fastq.gz', 'wt')
for batchIDX in batchIDX_LIST:
    fin = gzip.open(f'{tmpDir}/{fileName}.tmp_{batchIDX:03d}.fastq.gz','rt')
    for line in fin:
        fout.write(line)
    fin.close()
fout.close()