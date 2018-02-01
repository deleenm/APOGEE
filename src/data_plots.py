#!/usr/bin/env python
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
import argparse
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
    mytable = Table.read('../apogeeled.csv',format='ascii')
    print("APOGEE Table Read Done")
    return(mytable)

def read_plate_visits():
    #Create Table
    mytable = Table.read('../plate_visits.csv',format='ascii')
    print("Plate Visits Table Read Done")
    return(mytable)

def read_manga():
    #Create Table
    mytable = Table.read('../mangaled.csv',format='ascii')
    print("MaNGA Table Read Done")
    return(mytable)

def read_mjd():
    #Create Table
    mytable = Table.read('../mjd.hist',format='ascii')
    print("MJD Table Read Done")

    #Get date from final MJD
    endmjd = (np.sort(mytable['mjd']))[-1]
    print("End MJD: {}".format(endmjd))
    #Deal with shift from MJD to SJD
    jddate = Time(endmjd-1,format='mjd')
    date = (jddate.datetime).strftime('%m/%d/%y')   
    return(mytable,date,endmjd)

def read_proj():
    #Create Table
    mytable = Table.read('../proj.txt',format='ascii')
    print("Proj Table Read Done")
    return(mytable)

def read_master(south=False):
    #Create Table    
    if (south):
        mytable = Table.read('../schedule/Sched_LCO_base18.full.txt',format='ascii')
    else:
        mytable = Table.read('../schedule/Sch_base.Aug16.RM_ELG.txt',format='ascii')
    
    print("Read MasterTable Read Done")
    return(mytable)

def read_weather(south=False):
    #Create Table
    if(south):
        mytable = Table.read('../weather_south.txt',format='ascii')
    else:
        mytable = Table.read('../weather_north.txt',format='ascii')
    print("Read Weather Table Read Done")
    return(mytable)


def master_proj(startmjd=0.0,endmjd=None,withmeng=False,pergood=0.45,south=False):
    
    mjd_list = list()
    aptime_list = list()
    mantime_list = list()
    
    master_tab = read_master(south=south)
    
    #Use specified date range
    if(endmjd != None):
        master_tab = master_tab[(master_tab['MJD']  <= (endmjd + 2400000.0))]
    master_tab = master_tab[(master_tab['MJD'] >= (startmjd + 2400000.0))]
    
    
    for row in range(len(master_tab)):
        mjd_list.append(np.floor(master_tab['MJD'][row]) - 2400000.0)
        #See if apogee time
        if (withmeng == True):
            engint = 2
        else:
            engint = 1
        if(south == True):
            if(master_tab['Lead'][row] != 1):
                aptime_list.append(0)
                continue
        
        if(master_tab['Eng'][row] <= engint):
            bright_time = master_tab['MJD_end_bright'][row] - master_tab['MJD_start_bright'][row]
            good_weather = pergood
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
    pl.clf()

def apogee_proj(pp,south=False):
    (mjdtab,date,endmjd) = read_mjd()
    
    if(south):
        projtab = master_proj(endmjd=endmjd,south=south,pergood=.70)
    else:
        projtab = master_proj(endmjd=endmjd,south=south)
    
    print("Total Visits: {}".format(mjdtab['cum'][-1]))
    print("Total Projected Visits: {}".format(round(projtab['apvisits_cum'][-1])))
    
    pl.plot(mjdtab['mjd'],mjdtab['cum'],linewidth=2.0)
    pl.plot(projtab['mjd'],projtab['apvisits_cum'],linewidth=2.0)
    
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    pl.legend(('Current', 'Projected'),loc=2)
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()
    
def apogee_sn2(pp,south=False):
    (mjdtab,date,endmjd) = read_mjd()
    
    if(south):
        projtab = master_proj(endmjd=endmjd,south=south,pergood=0.70)
    else:
        projtab = master_proj(endmjd=endmjd,south=south)
        mprojtab = master_proj(endmjd=endmjd,withmeng=True,south=south)
        gprojtab = master_proj(endmjd=endmjd,pergood=0.50,south=south)
        gmprojtab = master_proj(endmjd=endmjd,withmeng=True,pergood=0.50,south=south)
    
    print("Total Visits: {}".format(mjdtab['cum'][-1]))
    print("Total SN2 Visits: {}".format(round(mjdtab['sn2'][-1])))
    print("Total SN2corr Visits: {}".format(round(mjdtab['sn2corr'][-1])))
    print("Total Complete: {}".format(round(mjdtab['complete'][-1])))
    print("Total Short Vistis: {}".format(round(mjdtab['2vsn2'][-1])))
    print("Total Projected Visits: {}".format(round(projtab['apvisits_cum'][-1])))
    if(not south):
        print("Total Projected with Eng Visits: {}".format(round(mprojtab['apvisits_cum'][-1])))
        print("Total Projected Visits (50%): {}".format(round(gprojtab['apvisits_cum'][-1])))
        print("Total Projected with Eng Visits (50%): {}".format(round(gmprojtab['apvisits_cum'][-1])))

    pl.plot(mjdtab['mjd'],mjdtab['cum'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['sn2'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['sn2corr'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['complete'],linewidth=2.0)
    pl.plot(projtab['mjd'],projtab['apvisits_cum'],'--',linewidth=2.0)
    if(not south):
        pl.plot(projtab['mjd'],mprojtab['apvisits_cum'],'--',linewidth=2.0)
    
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    if(south):
        pl.legend(('Visits','SN2 Visits','SN2corr Visits','Plate Complete','Projected'),loc=2)
    else:
        pl.legend(('Visits','SN2 Visits','SN2corr Visits','Plate Complete','Projected', 'Projected with Eng'),loc=2)
    #pl.savefig('test.png')
    pp.savefig()
    pl.clf()

    #Compare SN2corr and Plate Complete to Visits
    pl.plot(mjdtab['mjd'],mjdtab['cum'] - mjdtab['sn2corr'],color='r',linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['cum'] - mjdtab['complete'],color='c',linewidth=2.0)
    pl.axvline(57482,color='k',linestyle='--',linewidth=2.0)
    pl.title('Difference from Current({}) APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits - Corrected Visits")
    pl.legend(('SN2corr Visits','Plate Complete',"Twilight Obs Begin"),loc=2)
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()

    pl.subplot(211)
    pl.plot(mjdtab['mjd'],mjdtab['2vsn2'],color='k',linewidth=2.0)
    pl.axvline(57482,color='k',linestyle='--',linewidth=2.0)
    pl.title('Current({}) APOGEE-2 2-Exposure Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    #pl.legend(('SN2corr Visits','Plate Complete'),loc=2)
    pl.subplot(212)
    
    #Prevent divide by zero
    ratio = mjdtab['2vsn2'][(mjdtab['sn2'] != 0)]/mjdtab['sn2'][(mjdtab['sn2'] != 0)]
    pl.plot(mjdtab['mjd'][(mjdtab['sn2'] != 0)],ratio*100,color='k',linewidth=2.0)
    pl.axvline(57482,color='k',linestyle='--',linewidth=2.0)
    pl.xlabel('MJD')
    pl.ylabel("Percent of SN2 Visits")
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()
    
    #Everything combined
    pl.plot(mjdtab['mjd'],mjdtab['cum']-mjdtab['2vsn2'],color='b',linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['sn2corr']-mjdtab['2vsn2'],color='r',linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['complete']-mjdtab['2vsn2'],color='c',linewidth=2.0)
    pl.plot(projtab['mjd'],gprojtab['apvisits_cum'],'--',color='m',linewidth=2.0)
    pl.plot(projtab['mjd'],gmprojtab['apvisits_cum'],'--',color='y',linewidth=2.0)
    pl.axvline(57482,color='k',linestyle='--',linewidth=2.0)
    
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits - 2 Exp. Vists")
    #pl.xlim(57400,57550)
    #pl.ylim(ymin=600)
    pl.legend(('Visits','SN2corr Visits','Plate Complete','Projected (50%)', 
               'Projected with Eng (50%)','Twilight Obs Begin'),loc=2)
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()
    
    pl.plot(mjdtab['mjd'],mjdtab['cum']-mjdtab['2vsn2'],color='b',linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['sn2corr']-mjdtab['2vsn2'],color='r',linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['complete']-mjdtab['2vsn2'],color='c',linewidth=2.0)
    pl.plot(projtab['mjd'],gprojtab['apvisits_cum'],'--',color='m',linewidth=2.0)
    pl.plot(projtab['mjd'],gmprojtab['apvisits_cum'],'--',color='y',linewidth=2.0)
    pl.axvline(57482,color='k',linestyle='--',linewidth=2.0)
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits - 2 Exp. Vists")
    pl.xlim(57400,57550)
    pl.ylim(ymin=600)
    pl.legend(('Visits','SN2corr Visits','Plate Complete','Projected (50%)',
                'Projected with Eng (50%)','Twilight Obs Begin'),loc=2)  
    
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()

def standard_plot(pp,south=False):
    (mjdtab,date,endmjd) = read_mjd()
    if(south):
        projtab = master_proj(endmjd=endmjd,south=south,pergood=0.70)
    else:
        projtab = master_proj(endmjd=endmjd,south=south)
        mprojtab = master_proj(endmjd=endmjd,withmeng=True,south=south)
        gprojtab = master_proj(endmjd=endmjd,pergood=0.50,south=south)
        gmprojtab = master_proj(endmjd=endmjd,withmeng=True,pergood=0.50,south=south)
    
    print("Total Visits: {}".format(mjdtab['cum'][-1]))
    print("Total Complete: {}".format(round(mjdtab['complete'][-1])))
    print("Total Projected Visits: {}".format(round(projtab['apvisits_cum'][-1])))
    if(not south):
        print("Total Projected with Eng Visits: {}".format(round(mprojtab['apvisits_cum'][-1])))
        print("Total Projected Visits (50%): {}".format(round(gprojtab['apvisits_cum'][-1])))
        print("Total Projected with Eng Visits (50%): {}".format(round(gmprojtab['apvisits_cum'][-1])))

    pl.plot(mjdtab['mjd'],mjdtab['cum'],linewidth=2.0)
    pl.plot(mjdtab['mjd'],mjdtab['complete'],linewidth=2.0)
    pl.plot(projtab['mjd'],projtab['apvisits_cum'],'--',linewidth=2.0)
    if(not south):
        pl.plot(projtab['mjd'],mprojtab['apvisits_cum'],'--',linewidth=2.0)
    
    pl.title('Current({}) and Projected APOGEE-2 Visits'.format(date))
    pl.xlabel('MJD')
    pl.ylabel("Number of Visits")
    if(south):
        pl.legend(('Visits','Plate Complete','Projected'),loc=2)
    else:
        pl.legend(('Visits','Plate Complete','Projected', 'Projected with Eng'),loc=2)
    pl.savefig('test.png')
    pp.savefig()
    pl.clf()    

def apogee_type(pp,south=False):
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
    
    pl.title('Current({}) APOGEE-2 Visits by Type'.format(date))
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

def weather_plot(pp):
    wtab = read_weather()
    
    #Let's look at weather
    #https://trac.sdss.org/wiki/APO/Observatory/ObserverForum/SDSSOperationTimeTracking    
    wdate = [datetime.datetime.strptime(x,'%Y-%m-%d') for x in wtab['Date']]
    wmjd = Time(wdate,format='datetime')
    
    med_good = np.median(wtab['per_good'])
    print("Median Percent Good Weather: {}%".format(med_good))
    
    pl.plot(wmjd.mjd,wtab['per_good'],linewidth=2.0)
    
    pl.axhline(45,color='k',linestyle='--',linewidth=2.0)
    pl.axhline(med_good,color='r',linestyle='--',linewidth=2.0)
    pl.title('Percent Good Weather')
    pl.xlabel('MJD')
    pl.ylabel("Percent Good Weather")
    pl.legend(('Data', '45% Good Weather', 'Median Weather'),loc=2)
    pp.savefig()
    pl.clf()


    pl.plot(wdate,wtab['per_good'],linewidth=2.0)
    pl.axhline(45,color='k',linestyle='--',linewidth=2.0)
    pl.axhline(med_good,color='r',linestyle='--',linewidth=2.0)
    pl.gcf().autofmt_xdate()
    pl.title('Percent Good Weather')
    pl.ylabel("Percent Good Weather")
    pl.legend(('Data', '45% Good Weather', 'Median Weather'),loc=2)
    pp.savefig()
    pl.clf()
    
# -------------
# Main Function
# -------------
def data_plots_main(south=False):

    #North or South
    if(south == False):
        print("Running in Northern mode")
    else:
        print("Running in Southern mode")
    
    pp = PdfPages('../progress.pdf')
    #plate_windows(pp)
    #pp = PdfPages('progress.pdf')
    #combined_plot(pp) #Not up to date.
    apogee_proj(pp,south=south)
    apogee_type(pp)
    #apogee_sn2(pp,south=south)
    standard_plot(pp,south=south)
    weather_plot(pp)
    pp.close()
    print("Complete!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates Additional files and plots for tracking APOGEE Progress.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s','--south', action = 'store_true',help='Run this with the Southern options')
    group.add_argument('-n','--north', action = 'store_true',help='Run this with the Northern options')
    args = vars(parser.parse_args())
    data_plots_main(args['south'])
    
##
#@mainpage
 #@copydetails  data_plots
    