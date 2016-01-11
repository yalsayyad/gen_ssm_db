# gen_ssm_db

## Install Dependencies

```
git clone https://github.com/EUPSForge/oorb.git
cd oorb
./configure gfortran opt
cd python
make pyoorb
cd ..
```

### OpenOrb Libraries

This produces two libraries:
 * python/pyoorb.so
 * lib/liboorb.so or ../lib/liboorb.dylib

Note: do not be tricked by the file called python/liboorb.f90. This is not the liboorb you are looking for.

Make pyoorb is visible to your $PYTHONPATH and the liboorb.* is visible to your $LD_LIBRARY_PATH or $DYLD_LIBRARY_PATH:

export LD_LIBRARY_PATH=$PWD/lib/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$PWD/lib/:$DYLD_LIBRARY_PATH
export PYTHONPATH=$PWD/python/:$PYTHONPATH


### OpenOrb configuration

You also need to control the integration configuration paramters. You can use the default configuration
in main/oorb.conf or point to any file you choose. The configuration used in CATSIM's 2011 ephemeris generation is in gen_ssm_db/config/oorb.conf.

csh:
```
setenv OORB_DATA $PWD/data
setenv OORB_CONF $PWD/main/oorb.conf
```
bash:
```
export OORB_DATA=$PWD/data
export OORB_CONF=$PWD/main/oorb.conf
```

You will also need the planetary ephemeris file: de405.dat.
Put this in data/

## Run unit tests
git clone https://github.com/yalsayyad/gen_ssm_db.git
cd gen_ssm_db
export PYTHONPATH=$PWD:$PYTHONPATH

python 



