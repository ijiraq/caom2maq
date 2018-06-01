#!/bin/bash
source activate astroconda
# for camera in HSC FCS CIA
for camera in SUP
do
   for year in 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017
   do 
      python subaru2caom2.py --repo tardis  --verbose ${camera} ${year} >& ${camera}_${year}_err.txt  & 
   done 
done
