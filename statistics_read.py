#####################################################################################
class Statistics_Length:
    def __init__(self, length_LIST):
        self.length_LIST = length_LIST

        self.sortedLength_LIST = sorted(self.length_LIST, reverse=True)
        self.totalLength = sum(self.sortedLength_LIST)

        self.maxLength = self.sortedLength_LIST[0]
        self.minLength = self.sortedLength_LIST[-1]

    def get_total(self):
        return self.totalLength
    
    def get_lengthN(self):
        return len(self.sortedLength_LIST)

    def get_max(self):
        return self.maxLength
    
    def get_min(self):
        return self.minLength

    def get_N(self, N):
        sumLength = 0
        for idx, length in enumerate(self.sortedLength_LIST):
            sumLength += length
            if float(sumLength) / self.totalLength >= (float(N)/100):
                return length, idx+1

    def get_count(self, func):
        count_LIST = [0]*(int(func(self.maxLength)) + 1)
        sum_LIST   = [0]*(int(func(self.maxLength)) + 1)
        for length in self.sortedLength_LIST:
            countIDX = int(func(length))
            count_LIST[countIDX] += 1
            sum_LIST[countIDX]   += length
        return count_LIST, sum_LIST
#####################################################################################
from optparse import OptionParser
import sys, gzip
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None:
    print('Basic usage')
    print('')
    print('     python statistics_read.py -i test.fastq.gz')
    print('')
    sys.exit()

infile = opt.INPUT

#read file
length_LIST = []

fin = gzip.open(infile, 'rt')
for lineIDX, line in enumerate(fin):
    if lineIDX%4 == 1:
        sequence = line.rstrip('\n')
        length_LIST += [len(sequence)]
fin.close()

#Statistics_Length
statistics_Length = Statistics_Length(length_LIST)

#log
fout = open(opt.INPUT + '.log', 'w')
fout.write('Number of sequences:'  + '\t' + str(statistics_Length.get_lengthN()) + '\n')
fout.write('Total bases:'  + '\t' + str(statistics_Length.get_total()) + '\n')
fout.write('N50(L50):'  + '\t' + '\t'.join(map(str,(statistics_Length.get_N(50)))) + '\n')
fout.write('N90(L90):'  + '\t' + '\t'.join(map(str,(statistics_Length.get_N(90)))) + '\n')
fout.write('Min length:'  + '\t' + str(statistics_Length.get_min()) + '\n')
fout.write('Max length:'  + '\t' + str(statistics_Length.get_max()) + '\n')
fout.close()

#N1 ~ N100, L1 ~ N100
fout = open(infile + '.' + 'NL', 'w')
for N in range(1, 101):
    fout.write('N' + str(N) + '(L' + str(N) + ')' + '\t' + '\t'.join(map(str,(statistics_Length.get_N(N)))) + '\n')
fout.close()

#1kb count
fout = open(opt.INPUT + '.count.1kb', 'w')
count_LIST, sum_LIST = statistics_Length.get_count(lambda x: x/1000)
for countIDX, (count, sum) in enumerate(zip(count_LIST, sum_LIST)):
    fout.write('\t'.join(map(str,[countIDX, count, sum])) + '\n')
fout.close()

import math
#log count
fout = open(opt.INPUT + '.count.log2', 'w')
count_LIST, sum_LIST = statistics_Length.get_count(lambda x: math.log2(x))
for countIDX, (count, sum) in enumerate(zip(count_LIST, sum_LIST)):
    fout.write('\t'.join(map(str,[countIDX, count, sum])) + '\n')
fout.close()