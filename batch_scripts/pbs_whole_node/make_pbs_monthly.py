"""python script to generate pdb files, adapted from file form ljones to run on athena/minerva
"""

import sys
import os
from optparse import OptionParser
import numpy as np

# the root location for local (node) directories
scratchpartition = "/oasis/scratch/yusra/"
NPROC = 15

def parse_options(argv):
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-u", "--username", action="store", type="string", dest="username", help="send email to username@astro.", default="ljones")
    parser.add_option("-w", "--walltime", action="store", type="string", dest="walltime", help="estimated time for job to run", default="12:00:00")
    parser.add_option("-n", "--nodedir", action="store", type="string", dest="nodedir", help="copy files to /state/partition1/nodedir", default=None)
    parser.add_option("-f", "--commandfile", action="store", type="string", dest="commandfile", help="use tcsh script in 'commandfile' as executable portion of pbs script", default=None)
    parser.add_option("-d", "--desDirectory", action="store", type="string", dest="desDirectory", help="directory with des files", default=None)
    (options, argv) = parser.parse_args()
    return options

def write_headervars(pbsout, username, jobname, walltime):
    """ Write some typical PBS header file information """
    print >>pbsout, "#!/bin/bash"
    print >>pbsout, "#PBS -q normal"
    print >>pbsout, "#PBS -A nsa101"
    print >>pbsout, "#PBS -m a"
#    print >>pbsout, "#PBS -S /bin/tcsh"
    print >>pbsout, "#PBS -N %s"  %(jobname)
    # set an email address for job notification
    print >>pbsout, "#PBS -M %s@astro.washington.edu" %(username)
    print >>pbsout, "#PBS -V"
    # consolidate error and output logs into one file
    print >>pbsout, "#PBS -j oe"
    # set the walltime
    print >>pbsout, "#PBS -l walltime=%s" %(walltime)
    # and you could change this if desired/required
    print >>pbsout, "#PBS -l nodes=1:ppn=16:native"
    print >>pbsout, "#PBS -V"
    print >>pbsout, "cd /oasis/scratch/yusra/$PBS_JOBID"
    print >>pbsout, "### ---------------------------------------"
    print >>pbsout, "### Beginning of Executable Sections "
    print >>pbsout, "### ---------------------------------------"
    return

def write_logginglines(pbsout, nodedir):
    """ Write some lines that will record information useful for logging """
    print >>pbsout, " "
    print >>pbsout, "### ---------------------------------------"
    print >>pbsout, "### Create logging information"
    print >>pbsout, "### ---------------------------------------"
    print >>pbsout, " "
    # create some logging info    
    print >>pbsout, "set pbs_job_id = `echo $PBS_JOBID | awk -F . '{print $1}'`"
    print >>pbsout, "set num_procs = `wc -l < $PBS_NODEFILE`"
    print >>pbsout, "set master_node_id = `hostname`"
    # Print useful diagnostic information to log files
    print >>pbsout, "echo This is job $pbs_job_id"
    print >>pbsout, "echo The master node of this job is $master_node_id"
    print >>pbsout, "echo The directory from which this job was submitted is $PBS_O_WORKDIR"
    if nodedir != None:
        print >>pbsout, "echo The directory in which this job will run is %s" %(scratchpartition + nodedir)
    else:
        print >>pbsout, "echo No local node directory was indicated - job will run in `echo $PBS_O_WORKDIR`"
    print >>pbsout, "echo This job is running on  $num_procs processors"
    print >>pbsout, "echo This job is starting at `date`"
    print >>pbsout, "echo ---"
    return


def write_node_cleandirlines(pbsout, nodedir):
    """Add the tcsh job commands to remove the directories from the node"""
    """ (for normal script operation - nonfailure mode)"""
    # simply delete the directories we created here (no copy back)
    # I haven't implemented copy back to original dir because was having
    #  issues with NFS deleting files. Consider scp instead.
    print >>pbsout, "### ---------------------------------------"
    print >>pbsout, "### DELETE the local (node) directories and all files"
    print >>pbsout, "###   (does not delete parent directories if created"
    print >>pbsout, "### ---------------------------------------"
    print >>pbsout, "### echo Now deleting files in %s" % (scratchpartition + nodedir)
    print >>pbsout, "### /bin/rm -rf %s" % (scratchpartition + nodedir)
    print >>pbsout, "echo ---"
    print >>pbsout, "echo PBS job finished at `date`"
    print >>pbsout, " "
    print >>pbsout, "###"
    return

if __name__ == "__main__":
    # jobname, pbsfilename, walltime and job command can be input from command line
    options = parse_options(sys.argv)
    username = options.username
    walltime = options.walltime
    commandfile = options.commandfile
    nodedir = options.nodedir
    desDirectory = options.desDirectory
    for month in np.arange(49343, 53000, 30):
        for i in range(13):
            desDirectory = '/oasis/scratch/yusra/temp_project/orbits/49353/S1_%02d' % (i)
            famSubset = os.path.basename(os.path.normpath(desDirectory))
            pbsfilename = 'S1_%02d_%i.pbs' % (i, month)
            try:
                pbsout = open(pbsfilename, 'w')
            except IOError:
                print "Could not open %s for writing PBS script" % (pbsfilename)
                sys.exit()
            # write the header PBS stuff
            jobname = pbsfilename.replace('.pbs', '')
            write_headervars(pbsout, username, jobname, walltime)
            # write the lines to give PBS some useful logging info
            write_logginglines(pbsout, nodedir)
            print >>pbsout, """ls %s/*_49353.des | parallel -j+0 '/oasis/scratch/yusra/temp_project/ssm_src/gen_ephems_monthly.csh {/.} %i' """%(desDirectory, month)
            print >>pbsout, 'rm fort.10'
            # write the remove local directories stuff
            if nodedir is not None:
                write_node_cleandirlines(pbsout, nodedir)

            pbsout.close()
            print "Created PBS file %s" % (pbsfilename)