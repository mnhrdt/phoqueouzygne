#!/bin/csh
#
#  don't include the .-files
#
setenv COPYFILE_DISABLE true
#
#   create the tar file
#
rm -r */bin/*
rm -r */lib/*
tar -cvf /topex/ftp/pub/gmtsar/GMTSAR_V9.6.tar .
#tar -cvf ../GMTSAR_V9.6.tar .
#
