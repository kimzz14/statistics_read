################################################################################################
last_modified = '2024-11-21'
version = '0.0.1'
print('Last Modified:', last_modified, 'Version:', version)
################################################################################################
from optparse import OptionParser
from joblib import Parallel, delayed
import gzip, sys
import numpy as np

#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'input',help = "")
parser.add_option("-o","--output",action = 'store',type = 'string',dest = 'output',help = "")
parser.add_option("-s","--sampleN",action = 'store',type = 'int',dest = 'sampleN',help = "")
parser.add_option("-t","--totalN",action = 'store',type = 'int',dest = 'totalN',help = "")
(opt, args) = parser.parse_args()
if opt.input == None or opt.output == None or opt.sampleN == None:
    print('Basic usage')
    print('')
    print('     python fastqSampler.py -s 10000 -t 1000000 (optional) -i test -o test.sub')
    print('')
    sys.exit()

totalN = 0
if opt.totalN == None:
    print('Since the total number of reads was not entered, it will be calculated.')
    fin = gzip.open(opt.input + '_1.fastq.gz', 'rt')
    for line in fin:
        totalN += 1
    fin.close()
    totalN = int(totalN / 4)
else:
    totalN = opt.totalN

sampleN = opt.sampleN

print('totalN:', totalN, 'sampleN', sampleN)

random_LIST = np.random.choice(np.arange(totalN + 1), sampleN, replace=False)
random_LIST = np.sort(random_LIST)

fout = open(opt.output + '.readIDX', 'w')
for random in random_LIST:
    fout.write(str(random) + '\n')
fout.close()


random_LIST = np.append(random_LIST, random_LIST[-1])
random_LIST = random_LIST * 4

def run_single(infile, outfile):
    fin = gzip.open(infile, 'rt')
    fout = gzip.open(outfile, 'wt')
    lineN = 0
    writeN = 0
    randomIDX = 0
    for lineIDX, line in enumerate(fin):
        if lineIDX == random_LIST[randomIDX]:
            if randomIDX == sampleN:
                break
            randomIDX += 1
            writeN = 4
        
        if writeN > 0:
            fout.write(line)
            lineN += 1
            writeN -= 1
    fout.close()
    fin.close()

Parallel(n_jobs=2)(delayed(run_single)(infile, outfile) for infile, outfile in [(opt.input + '_1.fastq.gz', opt.output + '_1.fastq.gz'), (opt.input + '_2.fastq.gz', opt.output + '_2.fastq.gz')])

print('done!')