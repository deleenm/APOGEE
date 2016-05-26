#!/usr/bin/python
'''
Brief Description

Detailed Description

@package data_plots
@author ndelee
@version \e \$Revision$
@date \e \$Date$

Usage: data_plots.py
'''

# -----------------------------
# Standard library dependencies
# -----------------------------
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
def read_apogee():
    #Create Table
    mytable = Table.read('apogeeled.csv',format='ascii')
    print "APOGEE Table Read Done"
    return(mytable)

def read_plate_visits():
    #Create Table
    mytable = Table.read('plate_visits.csv',format='ascii')
    print "Plate Visits Table Read Done"
    return(mytable)

def read_manga():
    #Create Table
    mytable = Table.read('mangaled.csv',format='ascii')
    print "MaNGA Table Read Done"
    return(mytable)

def read_mjd():
    #Create Table
    colnames = ['mjd','daynum','cum','anc','apok','bulge','clus','disk','halo','goal','sat','sn2']
    mytable = Table.read('mjd.hist',format='ascii',names=colnames)
    print "MJD Table Read Done"

    #Get date from final MJD
    endmjd = (np.sort(mytable['mjd']))[-1]
    print("End MJD: {}".format(endmjd))
    #Deal with shift from MJD to SJD
    jddate = Time(endmjd-1,format='mjd')
    date = (jddate.datetime).strftime('%m/%d/%y')   
    return(mytable,date,endmjd)

def read_proj():
    #Create Table
    mytable = Table.read('proj.txt',format='ascii')
    print "Proj Table Read Done"
    return(mytable)

def read_master():
    #Create Table
    mytable = Table.read('Sch_base.jan16.v2h.txt',format='ascii')
    print "Read MasterTable Read Done"
    return(mytable)

def master_proj(startmjd=0.0,endmjd=None):
    
    mjd_list = list()
    aptime_list = list()
    mantime_list = list()
    
    master_tab = read_master()
    
    #Use specified date range
    if(endmjd != None):
        master_tab = master_tab[(master_tab['MJD']  <= (endmjd + 2400000.0))]
    master_tab = master_tab[(master_tab['MJD'] >= (startmjd + 2400000.0))]
    
    
    for row in range(len(master_tab)):
        mjd_list.append(np.floor(master_tab['MJD'][row]) - 2400000.0)
        #See if apogee time
        if(master_tab['Eng'][row] <= 1):
            bright_time = master_tab['MJD_end_bright'][row] - master_tab['MJD_start_bright'][row]
            good_weather = 0.45
            edge_efficiency = 0.985 #Deals with slack at end of night
            efficiency = good_weather*edge_efficiency
            
            aptime = bright_time*efficiency
        else:
            aptime = 0
        
        aptime_list.append(aptime)
    
    aptime_cum = np.cumsum(aptime_list)
    apvisits_list = np.array(aptime_list) / (87.0 / 60.0 / 24.0)
    apvisits_cum = np.cumsum(apvisits_list)
            
    night_tab = Table()
    
    night_tab['mjd'] = mjd_list
    night_tab['aptime'] = aptime_list
    night_tab['aptime_cum'] = aptime_cum
    night_tab['apvisits'] = apvisits_list
    night_tab['apvisits_cum'] = apvisits_cum
    
    return(night_tab)

def plate_windows(pp):
    plate_tab = read_plate_visits()
    #Remove N188 fields
    plate_tab = plate_tab[(plate_tab['loc_id'] != 5067)]
    
    plate_tab['window'] = (plate_tab['maxha'] - plate_tab['minha'])/15.0
    
    #Look at a single plate
    #row = (plate_tab['plate_id'] == 9006)
    #(plate_tab['plate_id','window'][row]).pprint(show_name=False)
    
    pl.scatter(plate_tab['ha'], plate_tab['dec'], s=25, c=plate_tab['window'])
    
    pl.title('Observing Window Sizes')
    pl.colorbar(label="Obseving Window in Hours")
    pl.xlabel('Drill Hour Angle (Degrees)')
    pl.ylabel("Declination")
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()
    
    pl.hist(plate_tab['window'],20)
    pl.title('Observing Window Sizes')
    pl.xlabel('Obseving Window in Hours')
    pl.ylabel("Number of Plates")
    #pl.savefig('test.png')
    pp.savefig()
    pl.clf()

def combined_plot(pp):
    apogeeled = read_apogee()
    mangaled = read_manga()
    apogeecum = np.cumsum(apogeeled['plates'])
    mangacum  = np.cumsum(mangaled['plates'])
    
    alldates = np.concatenate([apogeeled['mjd'],mangaled['mjd']])
    allplates = np.concatenate([apogeeled['plates'],mangaled['plates']])
    sortindex = np.argsort(alldates)
    
    totalcum = np.cumsum(allplates[sortindex])
    
    pl.plot(apogeeled['mjd'],apogeecum,linewidth=2.0)
    pl.plot(mangaled['mjd'],mangacum,linewidth=2.0)
    pl.plot(alldates[sortindex],totalcum,linewidth=2.0)
    print(apogeecum[-1])
    print(mangacum[-1])
    print(totalcum[-1])
    
    pl.title('Current APOGEE-2 Visits')
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    pl.legend(('APOGEE-led', 'MaNGA-led','Total'),loc=2)
    pl.savefig('test.png')
    pp.savefig()

def apogee_proj(pp):
    (mjdtab,date,endmjd) = read_mjd()
    projtab = master_proj(endmjd=endmjd)
    
    print("Total Vists: {}".format(mjdtab['cum'][-1]))
    print("Total Projected Vists: {}".format(round(projtab['apvisits_cum'][-1])))
    
    pl.plot(mjdtab['mjd'],mjdtab['cum'],linewidth=2.0)
    pl.plot(projtab['mjd'],projtab['apvisits_cum'],linewidth=2.0)
    
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    pl.legend(('Current', 'Projected'),loc=2)
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()
    
def apogee_sn2(pp):
    (mjdtab,date,endmjd) = read_mjd()
    projtab = master_proj(endmjd=endmjd)

    print("Total Vists: {}".format(mjdtab['cum'][-1]))
    print("Total Projected Vists: {}".format(round(projtab['apvisits_cum'][-1])))

    pl.plot(mjdtab['mjd'],mjdtab['cum'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['sn2']/3333.0,linewidth=2.0)
    pl.plot(projtab['mjd'],projtab['apvisits_cum'],linewidth=2.0)
    
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    pl.legend(('Current','SN2 Visits', 'Projected'),loc=2)
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()

    pl.plot(mjdtab['mjd'],mjdtab['cum'] - mjdtab['sn2']/3333.0,linewidth=2.0)
    pl.plot(mjdtab['mjd'],projtab['apvisits_cum'] - mjdtab['sn2']/3333.0,linewidth=2.0)
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    #pl.legend(('Current','SN2 Visits', 'Projected'),loc=2)
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()


def apogee_type(pp):
    (mjdtab,date,endmjd) = read_mjd()
    #projtab = read_proj()
    
    pl.plot(mjdtab['mjd'],mjdtab['cum'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['anc'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['apok'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['bulge'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['clus'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['disk'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['halo'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['goal'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['sat'],linewidth=2.0)
    
    pl.title('Current([]) APOGEE-2 Visits by Type'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    pl.yscale('log')
    #Creates nicer looking log plot tics
    pl.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(
                                                        int(np.maximum(-np.log10(y),0)))).format(y)))
    lgd = pl.legend(('All','Anc','APOKASC','Bulge','Cluster','Disk','Halo','Goal','Satellite')
              ,loc=2,bbox_to_anchor=(1, 1))
    pl.savefig('test.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    pp.savefig(bbox_extra_artists=(lgd,), bbox_inches='tight')
    pl.clf()
    
# -------------
# Main Function
# -------------
def data_plots_main():

    
    pp = PdfPages('progress.pdf')
    #plate_windows(pp)
    #pp = PdfPages('progress.pdf')
    #combined_plot(pp)
    #apogee_proj(pp)
    #apogee_type(pp)
    apogee_sn2(pp)
    pp.close()
    print("Complete!")

if __name__ == '__main__':
    data_plots_main()
    
##
#@mainpage
 #@copydetails  data_plots
    