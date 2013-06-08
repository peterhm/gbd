#!/bin/bash

#$ -S /bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/boost-current/lib
/usr/local/epd-7.2-2/bin/python dmco_fit_posterior.py "$1" "$2" "$3" "$4" "$5" "$6" "$7"