import numpy as np
from astroquery.eso import Eso
import sys


target_names = sys.argv[1].split(',')
print ''
print 'Check the observations in ESO Archive...'  
for star in target_names:
    print '\n'+star+':\n',
    found=False
    # Search in all the instruments in ESO Archive	
    for inst in Eso.list_instruments():
    	# 
        x=Eso.query_instrument(inst,column_filters={'target':star,'dp_cat':'SCIENCE','box':'00 05 00'})
        try:
            if len(x)>0:
                print inst.upper()+'('+str(len(x))+')',
                found=True
        except:
            continue
    if found:
    	print ''
    else:
    	print 'No observations found'
	    	
