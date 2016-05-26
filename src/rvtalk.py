#!/usr/bin/python
'''
Brief Description

Detailed Description

@package rvtalk
@author ndelee
@version \e \$Revision$
@date \e \$Date$

Usage: rvtalk.py
'''

# -----------------------------
# Standard library dependencies
# -----------------------------

# -------------------
# Third-party imports
# -------------------
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
# -----------------
# Class Definitions
# -----------------

# --------------------
# Function Definitions
# --------------------

# -------------
# Main Function
# -------------
def rvtalk_main():
    rg = np.array([115,56,9])
    ms = np.array([71,41,45])
    rc = np.array([18,5,0])
    sg = np.array([9,10,3])    
    
    pp = PdfPages('rvtalk.pdf')
    ind = np.arange(3)    # the x locations for the groups
    width = 0.5
    
    pl.bar(ind,rg,width,color='r')
    pl.bar(ind,ms,width,color='purple',bottom=rg)
    pl.bar(ind,rc,width,color='g',bottom=rg+ms)
    pl.bar(ind,sg,width,color='b',bottom=rg+ms+rc)
    
    
    pl.title("APOGEE Gold Sample Companion Types")
    pl.ylabel("Number of Companion Candidates")
    pl.xticks(ind + width/2., ('Binaries', 'Brown Dwarfs','Planets'))
    pl.xlim(-.25,2.75)
    pl.legend(('Red Giant', 'Main Seq.','Red Clump','Subgiant'))
    pl.savefig("rvtalk.png")
    pp.savefig()
    pp.close()

if __name__ == '__main__':
    rvtalk_main()
    
##
#@mainpage
 #@copydetails  rvtalk
    