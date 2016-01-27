#!/bin/csh

echo working directory is
pwd

set fileroot = $1
set STARTEPOCH = $2

#### CONFIGURATION ###
set TIMELENGTH = 30
set ORBIT_DIR = "/oasis/scratch/yusra/temp_project/orbits/"
set COEFF_DIR = "/oasis/scratch/yusra/temp_project/coeffs/"
set GEN_COEF_DIR = "/oasis/scratch/yusra/temp_project/ssm_src/"
set oorb_dir = "/oasis/scratch/yusra/temp_project/oorb/main/"
set FAILURES_FILE = "/oasis/scratch/yusra/temp_project/FAILURES.lis"
set NANS_FILE = "/oasis/scratch/yusra/temp_project/failed_nans.txt"

# User specific aliases and functions
# Not necessary to load lsst stack
setenv LSST_DEVEL /home/yusra/lsst_devel
source /oasis/scratch/yusra/temp_project/lsstStack/loadLSST.bash


#SET NECESSARY ENVIRONMENT VARIABLES
setenv PATH /oasis/scratch/yusra/temp_project/oorb/main:$PATH
setenv OORB_CONF /oasis/scratch/yusra/temp_project/oorb/main/krughoff.conf
setenv OORB_DATA /oasis/scratch/yusra/temp_project/oorb/data
setenv LD_LIBRARY_PATH /oasis/scratch/yusra/temp_project/oorblib:$LD_LIBRARY_PATH
setenv PYTHONPATH /oasis/scratch/yusra/temp_project/oorblib:$PYTHONPATH

echo $PYTHONPATH

# Lets go!
set foldername =  `echo $fileroot | awk '{split($0,a,"."); print a[1]}'`
set middle =  `echo $fileroot | awk '{split($0,a,"_"); print a[2]}'`
set famtype =  `echo $fileroot | awk '{split($0,a,"_"); print a[1]}'`
set endname = `echo $fileroot | awk '{split($0,a,"_"); print a[3]}'`
set epoch =  `echo $endname | awk '{split($0,a,"."); print a[1]}'`

echo epoch is $epoch
echo foldername is $foldername
echo middle is $middle
echo famtype is $famtype
echo endname is $endname

set plainfileroot = $famtype"_"$middle

echo startepoch is $STARTEPOCH
echo timelength is $TIMELENGTH
echo $ORBIT_DIR/$epoch/$foldername/$fileroot".des" $plainfileroot".des"

#Copy orbit file over to node
echo rsync -pt $ORBIT_DIR/$epoch/$foldername/$fileroot".des" $plainfileroot".des"
rsync -pt $ORBIT_DIR/$epoch/$foldername/$fileroot".des" $plainfileroot".des"

# Generate coefficients
echo python $GEN_COEF_DIR/gen_coeff_flexible_recursive.py $plainfileroot".des" $STARTEPOCH 30 14 $TIMELENGTH
python $GEN_COEF_DIR/gen_coeff_flexible_recursive.py $plainfileroot".des" $STARTEPOCH 30 14 $TIMELENGTH

echo 'python done'

if(-e $plainfileroot".des.done.txt") then
    #gzip and rsync back to shared storage
    echo $plainfileroot".des.coef_vartime_14.dat"
    echo $plainfileroot".des.resid_sum_vartime_14.dat"
    gzip $plainfileroot".des.coef_vartime_14.dat"
    gzip $plainfileroot".des.resid_sum_vartime_14.dat"

    rsync -pt $plainfileroot".des.resid_sum_vartime_14.dat.gz" $COEFF_DIR/$STARTEPOCH/$foldername/$plainfileroot".des.resid_sum_vartime_14.dat.gz" 
    rsync -pt $plainfileroot".des.coef_vartime_14.dat.gz" $COEFF_DIR/$STARTEPOCH/$foldername/$plainfileroot".des.coef_vartime_14.dat.gz"
    cat $plainfileroot".des.failed_14.dat" >> $NANS_FILE
else
    echo $fileroot, $STARTEPOCH, $TIMELENGTH, `date` >> $FAILURES_FILE 
endif

ls -lht .

rm $plainfileroot".des"
rm $plainfileroot".des.done.txt"
rm $plainfileroot".des.coef_vartime_14.dat.gz"
rm $plainfileroot".des.resid_sum_vartime_14.dat.gz"
rm $plainfileroot".des.failed_14.dat"
#rm fort.10
#ls -lht .

