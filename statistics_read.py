from lib.JobTimer.JobTimer import JobTimer
from optparse import OptionParser
import sys, gzip, pragzip

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

jobTimer = JobTimer()


#count line from gzip
lineN = 0

fin = pragzip.open(infile)
while chunk := fin.read(1024*1024):
    lineN += chunk.count(b'\n')
fin.close()

#read file
qualCount_LIST  = [0]*93
readLength_LIST = []

jobTimer.reset()
fin = gzip.open(infile, 'rt')
for lineIDX, line in enumerate(fin):
    if (lineIDX)%int(lineN / 20) == 0:
        jobTimer.check()
        percentage = float(lineIDX+1)/lineN
        print("Read File... [{0:6.2f}%] remainTime: {1}".format(percentage*100, jobTimer.remainTime(percentage)))

    if lineIDX%4 == 0:
        pass
    elif lineIDX%4 == 1:
        read = line.rstrip('\n')
        readLength = len(read)
        readLength_LIST += [readLength]
        pass
    elif lineIDX%4 == 2:
        pass
    elif lineIDX%4 == 3:
        pass
    else:
        qual_LIST = line.rstrip('\n')
        for qual in qual_LIST:
            qualCount_LIST[ord(qual) - 33] += 1

jobTimer.check()
percentage = 1.0
print("Read File... [{0:6.2f}%] remainTime: {1}".format(percentage*100, jobTimer.remainTime(percentage)))
fin.close()

#write qual
fout = open(opt.INPUT + '.qual', 'w')
for qualIDX, qualCount in enumerate(qualCount_LIST):
    fout.write('\t'.join([qualIDX, qualCount]) + '\n')
fout.close()


#
readLength_LIST = sorted(readLength_LIST, reverse=True)

readN = len(readLength_LIST)

totalReadLength = sum(readLength_LIST)
sumReadLength = 0

maxReadLength = max(readLength_LIST)
minReadLength = min(readLength_LIST)

for readLength in readLength_LIST:
    sumReadLength += readLength
    if float(sumReadLength) / totalReadLength > 0.5:
        N50 = readLength
        break

print('Number of reads:'  + '\t' + str(readN))
print('Total bases:'  + '\t' + str(totalReadLength))
print('Read length N50:'  + '\t' + str(N50))
print('Min read length:'  + '\t' + str(minReadLength))
print('Max read length:'  + '\t' + str(maxReadLength))


fout = open(opt.INPUT + '.logCount', 'w')

import math

logCount_LIST = [0]*(int(math.log2(maxReadLength)) + 1)
logSum_LIST = [0]*(int(math.log2(maxReadLength)) + 1)

for readLength in readLength_LIST:
    logCountIDX = int(math.log2(readLength))
    logCount_LIST[logCountIDX] += 1
    logSum_LIST[logCountIDX] += readLength

for logCountIDX, (logCount, logSum) in enumerate(zip(logCount_LIST, logSum_LIST)):
    fout.write('\t'.join(map(str,[logCountIDX, logCount, logSum])) + '\n')

fout.close()

fout = open(opt.INPUT + '.count', 'w')

count_LIST = [0]*(int((maxReadLength)/1000) + 1)
sum_LIST   = [0]*(int((maxReadLength)/1000) + 1)

for readLength in readLength_LIST:
    countIDX = int((readLength)/1000)
    count_LIST[countIDX] += 1
    sum_LIST[countIDX] += readLength

for countIDX, (count, sum) in enumerate(zip(count_LIST, sum_LIST)):
    fout.write('\t'.join(map(str,[countIDX, count, sum])) + '\n')

fout.close()