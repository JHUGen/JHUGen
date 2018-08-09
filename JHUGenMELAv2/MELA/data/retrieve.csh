#!/bin/tcsh -f
# Script to retrieve the libmcfm.so library from the link specified in download.url

cd `dirname $0`/$1

if (! -e libmcfm_701.so) then
  cat download.url | xargs wget --no-check-certificate -q
endif
