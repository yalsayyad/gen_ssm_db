
dateSet = set()
import random


def makeNewFile(date):
    with  open('shell_%s.cfg'%(date),'w') as f:
        headerStr = """Notification = never
getenv = true
Executable = /astro/net/pogo3/yusra/ssm/gordon_copyback/ssm_src/gen_ephems_monthly.csh
Initialdir = /astro/net/pogo3/yusra/ssm/gordon_copyback/workingDir/%s
Universe = vanilla
requirements = ((Arch == "X86_64")) && Machine!= "ullr.astro.washington.edu"

Output = /astro/net/pogo3/yusra/ssm/gordon_copyback/workingDir/%s/output.out
Error = /astro/net/pogo3/yusra/ssm/gordon_copyback/workingDir/%s/error.out
Log = /astro/net/pogo3/yusra/ssm/gordon_copyback/workingDir/%s/log.out

"""%(date, date, date,  date)
        f.write(headerStr)

def printFile(fileroot, date):
    with open('shell_%s.cfg'%(date),'a') as f:
        f.write('Arguments = %s %s %i\n'%(fileroot, date, random.random()*30))
        f.write('Queue\n')

dateSet = set()
with open('redo_Jan14.dat', 'rb') as f:
    for line in f:
        ex, fileroot, yr = line.split()
        if yr not in dateSet:
            dateSet.add(yr)
            makeNewFile(yr)
            printFile(fileroot.replace('.des', ''), yr)
        else:
            printFile(fileroot.replace('.des', ''), yr)
