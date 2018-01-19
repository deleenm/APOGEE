#!/usr/bin/env python
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
from data_plots import read_mjd
from current_visits import read_planned
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

def read_lst():
     mytable = Table.read('../lstdist.txt',format='ascii')
     return mytable

def apolst(intimes):
    
    t = Time(intimes,format='jd')
    #Deals with there being a slight difference betwee utc and ut1 which are too small for me to care
    t.delta_ut1_utc = 0.
    #Location of APO 32.789278, -105.820278)
    lst = t.sidereal_time('apparent',longitude=-105.820278)
    return(lst.deg / 15.0)

def master_proj_lst(masterfile,startmjd=0.0,endmjd=None,withmeng=False):
    
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
    for row in range(len(master_tab)):
        mjd_list.append(np.floor(master_tab['MJD'][row]) - 2400000.0)
        nslots=0
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
            nslots = min([int(bright_time / ((par['exposure'] + par['overhead']) / 60)), par['ncarts']])
            if nslots != 0:
                times = [master_tab['MJD_start_bright'][row] + (par['exposure'] + par['overhead']) 
                         / 60 / 24 * x for x in range(nslots)]
                lengths = [(par['exposure'] + par['overhead']) / 60 for x in range(nslots)]
        
                # If APOGEE-II starts the split night
                if (master_tab['MJD_start_bright'][row] < master_tab['MJD_start_dark'][row] 
                                or master_tab['MJD_start_dark'][row] == 0):
                    # Determine whether we should add another exposure (leftover time > 15min)
                    if len(times) < par['ncarts'] and bright_time - sum(lengths) > par['overhead'] / 60:
                        times.append(master_tab['MJD_start_bright'][row] + len(times) 
                                     * (par['exposure'] + par['overhead']) / 60 / 24)
                        lengths.append(bright_time - sum(lengths))
                    # Because APOGEE-II is first, the last exposure will not have overhead
                    else: 
                        pass
                        #lengths[-1] = par['exposure'] / 60
        
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

def sched_diff():
    #Year 3 differences
    (begdate,enddate)= (57618,57945)
    #(begdate,enddate)= (57645,57648)
    
    (oldlst,oldlen,oldhours) = master_proj_lst('../schedule/Sch_base.jan16.v2h.txt',
                                               startmjd=begdate,endmjd=enddate)
    (newlst,newlen,newhours) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=begdate,endmjd=enddate)
    
    
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
    print("Old Hours: {:.2f} New Hours: {:.2f}".format(oldhours,newhours))
    print("Old Visit Hours: {:.2f} New Visit Hours: {:.2f}".format(np.sum(oldlen),np.sum(newlen)))
    print("Diff Hours: {:.2f} Diff Visit Hours: {:.2f} Diff Visits:{} ".format(oldhours-newhours,
                                                                     np.sum(oldlen)-np.sum(newlen),
                                                                     np.sum(lstdiffarr)))

#    print("Old Lengths: {}".format(oldlen))
#    print("New Lengths: {}".format(newlen))
#    print("Old Times: {}".format(oldhours))
#    print("New Times: {}".format(newhours))
#    print("Old LST: {}".format(oldlst))
#    print("New LST: {}".format(newlst))
    
    print('LST Diff')
    for i in range(24):
        print('{} {:.2f}'.format(i,-1.0*hrdiffarr[i]))

def predict_mjd(end1,good_weather,weng=False):
    beg1  = 0
    
    (lst,length,hours) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg1,endmjd=end1,withmeng=weng)
    
    edge_efficiency = 0.985 #Deals with slack at end of night
    efficiency = good_weather*edge_efficiency
    #Only take visits with lengths greater than 33+20 minutes
    print('Projected Visits {:.0f} through {}'.format(len(lst[(length > 0.88)])*good_weather,end1))
    print('Projected Visits (Blanton) {:.0f} through {}'.format(hours*60/87*efficiency,end1))
    

def predict_visits():
    #Year 1 
    (beg1,end1)= (56840,57210)
    #Year 2
    (beg2,end2)= (57250,57581)
    #Year 3
    (beg3,end3)= (57618,57945)
    #Year 4
    (beg4,end4)= (57983,58302)
    #Year 5
    (beg5,end5)= (58340,58666)
    #Year 6
    (beg6,end6)= (58710,59032)
    
    pp = PdfPages('../proj_lst.pdf')

    weng = False   
    (lsty1,leny1,hoursy1) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg1,endmjd=end1,withmeng=weng)
    (lsty2,leny2,hoursy2) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg2,endmjd=end2,withmeng=weng)
    (lsty3,leny3,hoursy3) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg3,endmjd=end3,withmeng=weng)
    (lsty4,leny4,hoursy4) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg4,endmjd=end4,withmeng=weng)
    (lsty5,leny5,hoursy5) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg5,endmjd=end5,withmeng=weng)
    (lsty6,leny6,hoursy6) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg6,endmjd=end6,withmeng=weng)
    (lstfull,lenfull,hoursfull) = master_proj_lst('../schedule/Sch_base.Aug16.RM_ELG.txt',
                                               startmjd=beg1,endmjd=end6,withmeng=weng)
    

    #Only take visits with lengths greater than 33+20 minutes
    (histy1,binsy1) = np.histogram(lsty1[(leny1 > 0.88)],bins=24)
    (histy2,binsy2) = np.histogram(lsty2[(leny2 > 0.88)],bins=24)
    (histy3,binsy3) = np.histogram(lsty3[(leny3 > 0.88)],bins=24)
    (histy4,binsy4) = np.histogram(lsty4[(leny4 > 0.88)],bins=24)
    (histy5,binsy5) = np.histogram(lsty5[(leny5 > 0.88)],bins=24)
    (histy6,binsy6) = np.histogram(lsty6[(leny6 > 0.88)],bins=24)
    (histfull,binsfull) = np.histogram(lstfull[(lenfull > 0.88)],bins=24)
        
    center = np.arange(0,24) + 0.5

    pergood = .45

    #Read planned file
    plan_tab = read_planned()


    out_tab = Table()
    out_tab['lst'] = center
    out_tab['Plan'] = plan_tab['visits']
    out_tab['yr1'] = np.around(histy1 * pergood,1)
    out_tab['yr2'] = np.around(histy2 * pergood,1)
    out_tab['yr3'] = np.around(histy3 * pergood,1)
    out_tab['yr4'] = np.around(histy4 * pergood,1)
    out_tab['yr5'] = np.around(histy5 * pergood,1)
    out_tab['yr6'] = np.around(histy6 * pergood,1)
    out_tab['Full_Survey'] = np.around(histfull * pergood,1)
    
    out_tab.write('../proj_lst.txt',format='ascii',overwrite=True)

    #LST 1 hour visit by year
    lst_dict = read_lst()
    
    (date,endmjd) = (read_mjd())[1:]
    
    predict_mjd(endmjd, pergood, weng)
    predict_mjd(end6,pergood,weng)
        
    pl.bar(center, lst_dict['Yr1_Visit'],1,color="pink",edgecolor='black')
    pl.bar(center, lst_dict['Yr2_Visit'],1,bottom=lst_dict['Yr1_Visit'],color="#58ACFA",edgecolor='black')
    pl.bar(center, lst_dict['Yr3_Visit'],1,bottom=lst_dict['Yr1_Visit']+lst_dict['Yr2_Visit']
                                                                                      ,color="lightgreen",edgecolor='black')
    pl.xlabel("LST")
    pl.ylabel("Number of visits")
    pl.title("Current ({}) Number of Visits by LST (1 hour bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
    pl.plot(center,out_tab['yr1'],color="r",linewidth=2.0)
    basenum = out_tab['yr1']
    pl.plot(center,out_tab['yr2']+basenum,color="b",linewidth=2.0)
    basenum = basenum + out_tab['yr2']
    pl.plot(center,out_tab['yr3']+basenum,color="g",linewidth=2.0)
    
    pl.plot(plan_tab['mid'],plan_tab['visits'],':',color="r",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*2.0,':',color="b",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*3.0,':',color="g",linewidth=2.0)
    
    pl.legend(('Proj Year1','Proj Year2','Proj Year3','Plan Year1','Plan Year2','Plan Year3','Actual Year1', 'Actual Year2','Actual Year3'),
              ncol=2, fontsize='small')
    pp.savefig()
    pl.clf()
    pp.close()

# -------------
# Main Function
# -------------
def schedule_main():
    #sched_diff()
    predict_visits()

    print("Finished")

if __name__ == '__main__':
    schedule_main()
    
##
#@mainpage
 #@copydetails  schedule
    