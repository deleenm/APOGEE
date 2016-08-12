#!/usr/bin/python
'''
Brief Description

Detailed Description

@package schedule
@author deleenm
@version \e \$Revision$
@date \e \$Date$

Usage: schedule.py
'''

# -----------------------------
# Standard library dependencies
# -----------------------------
from __future__ import print_function, division
import sys
import time
import datetime
# -------------------
# Third-party imports
# -------------------
from astropy.table import Table
from astropy.time import Time
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import astropy as astpy
# -----------------
# Class Definitions
# -----------------

# --------------------
# Function Definitions
# --------------------
def read_master(masterfile):
    #Create Table
    mytable = Table.read(masterfile,format='ascii')
    print ("Read MasterTable Read Done")
    return(mytable)

def apolst(intimes):
    
    t = Time(intimes,format='jd')
    #Deals with there being a slight difference betwee utc and ut1 which are too small for me to care
    t.delta_ut1_utc = 0.
    #Location of APO 32.789278, -105.820278)
    lst = t.sidereal_time('apparent',longitude=-105.820278)
    return(lst.deg / 15.0)

def master_proj_lst(masterfile,startmjd=0.0,endmjd=None,withmeng=False,pergood=0.45):
    
    mjd_list = list()
    lst_list = list()
    len_list = list()
    
    master_tab = read_master(masterfile)
    
    #Use specified date range
    if(endmjd != None):
        master_tab = master_tab[(master_tab['MJD']  <= (endmjd + 2400000.0))]
    master_tab = master_tab[(master_tab['MJD'] >= (startmjd + 2400000.0))]
    
    # Define APOGEE-II observing parameters
    par = {'exposure': 67, 'overhead': 20, 'ncarts': 9, 'maxz': 3, 'moon_threshold': 15, 
               'sn_target': 3136}
    
    hours=0
    nslots=0
    for row in range(len(master_tab)):
        mjd_list.append(np.floor(master_tab['MJD'][row]) - 2400000.0)
        #See if apogee time
        if (withmeng == True):
            engint = 2
        else:
            engint = 1
        
        if(master_tab['Eng'][row] <= engint):
            #Check for consistency
            if (master_tab['MJD_start_bright'][row] > master_tab['MJD_end_bright'][row]):
                print("{} Bright time start and stop wrong!".format(master_tab['MJD'][row]))
            
            # Define APOGEE-II blocks for tonight
            bright_time = (master_tab['MJD_end_bright'][row] - master_tab['MJD_start_bright'][row]) * 24
            hours = hours + bright_time
            nslots = min([int(round(bright_time / ((par['exposure'] + par['overhead']) / 60))), par['ncarts']])
            if nslots != 0:
                times = [master_tab['MJD_start_bright'][row] + (par['exposure'] + par['overhead']) 
                         / 60 / 24 * x for x in range(nslots)]
                lengths = [(par['exposure'] + par['overhead']) / 60 for x in range(nslots)]
        
                # If APOGEE-II starts the split night
                if (master_tab['MJD_start_bright'][row] < master_tab['MJD_start_dark'][row] 
                                or master_tab['MJD_start_dark'][row] == 0):
                    # Determine whether we should add another exposure (leftover time > 15min)
                    if len(times) < par['ncarts'] and bright_time - sum(lengths) > 0:
                        times.append(master_tab['MJD_start_bright'][row] + len(times) 
                                     * (par['exposure'] + par['overhead']) / 60 / 24)
                        lengths.append(bright_time - sum(lengths))
                    # Because APOGEE-II is first, the last exposure will not have overhead
                    else: lengths[-1] = par['exposure'] / 60
        
                # APOGEE-II ends the night
                else:
                    # Determine whether we can add an exposure (leftover time > 15min)
                    if len(times) < par['ncarts'] and bright_time - sum(lengths) > par['overhead'] / 60:
                        times.append(master_tab['MJD_start_bright'][row] + len(times) * (par['exposure'] + par['overhead']) / 60 / 24)
                        lengths.append(bright_time - sum(lengths))

            else:
                times = []
                lengths = []
            
        #Now determine the LSTs for each visit
        if nslots != 0:
            times = np.array(times)
            lengths = np.array(lengths)
    
            midlst = apolst(times + (lengths/24.0/2.0))
            lst_list.extend(midlst)
            len_list.extend(lengths)
            
 
    return(np.array(lst_list),np.array(len_list),hours)

# -------------
# Main Function
# -------------
def schedule_main():
    (oldlst,oldlen,oldhours) = master_proj_lst('../schedule/Sch_base.jan16.v2h.txt',
                                               startmjd=57618,endmjd=57945)
    (newlst,newlen,newhours) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=57618,endmjd=57945)
    
    pp = PdfPages('../schedule/schedule_diff.pdf')
    
    #Plot both visit distributions
    pl.hist(oldlst,bins=24,range=(0,24))
    pl.hist(newlst,bins=24,color='r')
    pl.xlabel('LST (Hours)')
    pl.ylabel('# Visits')
    pl.legend(("Jan16",'Aug16'))
    pp.savefig()
    pl.clf()

    #Plot Visit Difference
    (oldarr, oldbin) = np.histogram(oldlst, 24)
    (newarr, newbin) = np.histogram(newlst, 24)
    lstdiffarr = oldarr-newarr
    center = (oldbin[:-1] + oldbin[1:]) / 2
    pl.bar(center,lstdiffarr,width=1)
    pl.xlabel('LST (Hours)')
    pl.ylabel('Difference in # Visits')
    pp.savefig()
    pl.clf()

    #Plot Hourly difference
    #digitize gives the left right bin I want left bin edge
    oldbin_indice = (np.digitize(oldlst,bins=range(24)) -1)
    newbin_indice = (np.digitize(newlst,bins=range(24)) -1)
    
    #Sum up the hours appropriately
    oldhourhist = list()
    newhourhist = list()
    
    for i in range(24):
        oldhourhist.append(np.sum(oldlen[(oldbin_indice == i)]))
        newhourhist.append(np.sum(newlen[(newbin_indice == i)]))
    
    hrdiffarr = np.array(oldhourhist)-np.array(newhourhist)
    center = (oldbin[:-1] + oldbin[1:]) / 2
    pl.bar(center,hrdiffarr,width=1)
    pl.xlabel('LST (Hours)')
    pl.ylabel('Difference in # Hours')
    pp.savefig()
    pl.clf()
    pp.close()

    #Print out some useful info
    print("Old Hours: {} New Hours: {}".format(oldhours,newhours))
    print("Diff Hours: {} Diff Visits:{} ".format(oldhours-newhours,np.sum(lstdiffarr)))

    print('LST Diff')
    for i in range(24):
        print('{} {:.1f}'.format(i,-1.0*hrdiffarr[i]))

    print("Finished")

if __name__ == '__main__':
    schedule_main()
    
##
#@mainpage
 #@copydetails  schedule
    