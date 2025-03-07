from multiprocessing import Pool, shared_memory
import gzip
import numpy as np

from optparse import OptionParser
import sys, gzip
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
parser.add_option("-t","--threadN",action = 'store',type = 'int',dest = 'threadN',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None or opt.threadN == None:
    print('Basic usage')
    print('')
    print('     python statistics_read.py -t 24 -i test.fastq(.gz)')
    print('')
    sys.exit()

threadN = opt.threadN
batchN = 256
batchIDX_LIST = range(0, batchN)
fileName = opt.INPUT


#######################################################################################################
def batch_runner(func_name):
    result_LIST = None
    qualityCount_LIST = [0] * 200
    with Pool(processes=threadN) as pool:
        result_LIST = pool.map(func_name, batchIDX_LIST)

    data_LIST = []
    for _data_LIST, _qualityCount_LIST in result_LIST:

        for idx, qualityCount in enumerate(_qualityCount_LIST):
            qualityCount_LIST[idx] += qualityCount

        data_LIST += _data_LIST

    data_LIST = np.array(data_LIST)
    qualityCount_LIST = np.array(qualityCount_LIST)
    return data_LIST, qualityCount_LIST

#######################################################################################################
def read_fastq(batchIDX):
    if fileName.endswith('gz') == True:
        fin = gzip.open(fileName, 'rt')
    else:
        fin = open(fileName, 'rt')

    data_LIST = []
    qualityCount_LIST = [0] * 200
    for lineIDX, line in enumerate(fin):
        if (lineIDX>>2)%batchN != batchIDX: continue

        if lineIDX%4 == 1:
            sequence = line.rstrip('\n')
            sequenceLen = len(sequence)
            continue

        if lineIDX%4 == 3:
            quality_LIST = line.rstrip('\n')
            quality_LIST = list(map(ord, quality_LIST))

            for quality in quality_LIST:
                qualityCount_LIST[quality] += 1
            if len(quality_LIST) == 0:
                qualityMean = 0
            else:
                qualityMean = int(sum(quality_LIST) / len(quality_LIST)) - 33
            data_LIST += [[sequenceLen, qualityMean]]
            continue
    fin.close()
    return data_LIST, qualityCount_LIST



#######################################################################################################
dataSet_LIST, qualityCount_LIST = batch_runner(read_fastq)
sequenceLen_LIST = dataSet_LIST[:, 0]
qualityMean_LIST = dataSet_LIST[:, 1]

#
sequenceN = dataSet_LIST.shape[0]
meanLength = sequenceLen_LIST.mean()
minLength = sequenceLen_LIST.min()
maxLength = sequenceLen_LIST.max()

meanQuality = qualityMean_LIST.mean()
minQuality = qualityMean_LIST.min()
maxQuality = qualityMean_LIST.max()

#
totalLength_Q00 = sequenceLen_LIST.sum()
totalLength_Q10 = sequenceLen_LIST[(qualityMean_LIST >= 10)].sum()
totalLength_Q20 = sequenceLen_LIST[(qualityMean_LIST >= 20)].sum()
totalLength_Q30 = sequenceLen_LIST[(qualityMean_LIST >= 30)].sum()
totalLength_Q35 = sequenceLen_LIST[(qualityMean_LIST >= 35)].sum()
totalLength_Q40 = sequenceLen_LIST[(qualityMean_LIST >= 40)].sum()

#
sequenceLen_LIST.sort()

N_LIST = []
N = 0
subLength = 0
for idx, length in enumerate(sequenceLen_LIST[::-1]):
    subLength += length
    while float(subLength / totalLength_Q00) >= (float(N)/100):
        N_LIST += [(length, idx + 1)]
        N += 1

#qualityCount
fout = open(fileName + '.qualityCount', 'w')
for idx in range(33, 133):
    fout.write(str(idx - 33) + '\t' + str(qualityCount_LIST[idx]) + '\n')
fout.close()

#log
fout = open(fileName + '.log', 'w')
fout.write(f'Number of sequences:'  + '\t' + str(sequenceN) + '\n')
fout.write(f'Total bases(Q00):\t{totalLength_Q00}\n')
fout.write(f'Total bases(Q10):\t{totalLength_Q10} ({totalLength_Q10/totalLength_Q00 * 100:.2f}%)\n')
fout.write(f'Total bases(Q20):\t{totalLength_Q20} ({totalLength_Q20/totalLength_Q00 * 100:.2f}%)\n')
fout.write(f'Total bases(Q30):\t{totalLength_Q30} ({totalLength_Q30/totalLength_Q00 * 100:.2f}%)\n')
fout.write(f'Total bases(Q35):\t{totalLength_Q35} ({totalLength_Q35/totalLength_Q00 * 100:.2f}%)\n')
fout.write(f'Total bases(Q40):\t{totalLength_Q40} ({totalLength_Q40/totalLength_Q00 * 100:.2f}%)\n')
fout.write(f'N50(L50):\t{N_LIST[50][0]}({N_LIST[50][1]})\n')
fout.write(f'N90(L90):\t{N_LIST[90][0]}({N_LIST[90][1]})\n')
fout.write(f'mean length:\t{meanLength:0.4f}\n')
fout.write(f'Min length:\t{minLength}\n')
fout.write(f'Max length:\t{maxLength}\n')
fout.write(f'mean quality:\t{meanQuality:0.4f}\n')
fout.write(f'min quality:\t{minQuality:0.4f}\n')
fout.write(f'max quality:\t{maxQuality:0.4f}\n')
fout.close()


#####################################################################################################################s
sequenceLen_LIST.sort()
sequenceLen_LIST_shm = shared_memory.SharedMemory(create=True, size=sequenceLen_LIST.nbytes)
sequenceLen_LIST_shared = np.ndarray(sequenceLen_LIST.shape, dtype=sequenceLen_LIST.dtype, buffer=sequenceLen_LIST_shm.buf)
np.copyto(sequenceLen_LIST_shared, sequenceLen_LIST)



def parallel_histogram(batchIDX, bin_func):
    histo_LIST = np.zeros((bin_func(maxLength) + 1, 2), dtype=int)
    for idx, length in enumerate(sequenceLen_LIST_shared):
        if idx%batchN != batchIDX: continue
        bin = bin_func(length)
        histo_LIST[bin][0] += 1
        histo_LIST[bin][1] += length
    return histo_LIST

#1kb count
def bin_func(length):
    return int(length/1000)

histogram_LIST = None
with Pool(processes=threadN) as pool:
    histogram_LIST = pool.starmap(parallel_histogram, [(batchIDX, bin_func) for batchIDX in batchIDX_LIST])
histogram = np.sum(histogram_LIST, axis=0)

fout = open(fileName + '.count.1kb', 'w')
for countIDX, (count, sum) in enumerate(histogram):
    fout.write('\t'.join(map(str,[countIDX, count, sum])) + '\n')
fout.close()

#log count
import math
def bin_func(length):
    return int(math.log2(length))

histogram_LIST = None
with Pool(processes=threadN) as pool:
    histogram_LIST = pool.starmap(parallel_histogram, [(batchIDX, bin_func) for batchIDX in batchIDX_LIST])
histogram = np.sum(histogram_LIST, axis=0)

fout = open(fileName + '.count.log2', 'w')
for countIDX, (count, sum) in enumerate(histogram):
    fout.write('\t'.join(map(str,[countIDX, count, sum])) + '\n')
fout.close()

sequenceLen_LIST_shm.close()
sequenceLen_LIST_shm.unlink()