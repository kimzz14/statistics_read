###############################################################################
import time, math
class JobTimer:
    def __init__(self):
        self.tic = None

    def reset(self):
        self.tic = time.time()

    def check(self):
        self.toc = time.time()

    def remainTime(self, percentage):
        remainTime = (self.toc - self.tic)/percentage*(1.0-percentage)

        day = math.floor(remainTime/60/60/24)
        hour = math.floor((remainTime - day*60*60*24)/60/60)
        minute = math.floor((remainTime - day*60*60*24 - hour*60*60)/60)
        second = math.floor((remainTime - day*60*60*24 - hour*60*60 - minute*60))

        timeText = str(day).rjust(3, ' ') + ' days, ' + str(hour).rjust(2,'0') + ':' + str(minute).rjust(2,'0') + ':' + str(second).rjust(2,'0')

        return timeText
###############################################################################

from optparse import OptionParser
import sys, gzip

#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
parser.add_option("-o","--output",action = 'store',type = 'string',dest = 'OUTPUT',help = "")
parser.add_option("-m","--min",action = 'store',type = 'int',dest = 'MIN',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None or opt.OUTPUT == None or opt.MIN == None:
    print('Basic usage')
    print('')
    print('     python filter_read.py -m 9000 -i test.fastq.gz -o pass.fastq.gz')
    print('')
    sys.exit()

infile = opt.INPUT
outfile = opt.OUTPUT
minLength = opt.MIN


jobTimer = JobTimer()


#count line from gzip
lineN = 0
fin = gzip.open(infile, 'rt')
for line in fin:
    lineN += 1
fin.close()

#read file
qualCount_LIST  = [0]*100
readLength_LIST = []

jobTimer.reset()
fin = gzip.open(infile, 'rt')
fout = gzip.open(outfile, 'wt')
for lineIDX, line in enumerate(fin):
    if (lineIDX)%int(lineN / 20) == 0:
        jobTimer.check()
        percentage = float(lineIDX+1)/lineN
        print("Read File... [{0:6.2f}%] remainTime: {1}".format(percentage*100, jobTimer.remainTime(percentage)))
        sys.stdout.flush()

    if lineIDX%4 == 0:
        readID = line
    elif lineIDX%4 == 1:
        readSeq = line
    elif lineIDX%4 == 2:
        readSpace = line
    elif lineIDX%4 == 3:
        readQual = line

        if len(readSeq.rstrip('\n')) >= minLength:
            fout.write(readID)
            fout.write(readSeq)
            fout.write(readSpace)
            fout.write(readQual)


    else:
        print('bug')        


jobTimer.check()
percentage = 1.0
print("Read File... [{0:6.2f}%] remainTime: {1}".format(percentage*100, jobTimer.remainTime(percentage)))
sys.stdout.flush()
fin.close()
fout.close()