#!/bin/bash

#$ -S /bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/boost-current/lib
/usr/local/epd-7.3-2/bin/python /homes/peterhm/gbd/dmco_plot_fits_pdf.py $1 $2 $3 $4 