#!/bin/sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/boost-current/lib
/usr/local/epd-7.0-2/bin/python /homes/peterhm/gbd/dmco_fit_posterior.py $1 $2 $3 $4 $5