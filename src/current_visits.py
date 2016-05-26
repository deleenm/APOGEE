#!/usr/bin/python
'''
Brief Description

Detailed Description

@package current_visits
@author ndelee
@version \e \$Revision$
@date \e \$Date$

Usage: current_visits.py
'''

# -----------------------------
# Standard library dependencies
# -----------------------------
import math as m
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
import os.path
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import astropy as astpy
from scipy.cluster._hierarchy import cluster_dist
from boto.dynamodb.condition import NULL
# -----------------
# Class Definitions
# -----------------

# --------------------
# Function Definitions
# --------------------
def completion(vplan, vdone, sn, cadence_flag):
    ''' Computes completion percentage of APOGEE-II plate '''
    # Something is wrong here...
    if vplan == 0: return 1

    try:
        # 90% of completion percentage is from number of visits
        visit_completion = 0.9 * min([1, vdone / vplan])
        # 10% of completion percentage is from S/N
        sn_completion = 0.1 * calculateSnCompletion(vplan, sn)
        return visit_completion + sn_completion
    except:
        raise RuntimeError("ERROR: unable to calculate completion for vplan: %d, vdone: %d, sn: %d\n%s" %\
            (vplan, vdone, sn, sys.exc_info()))


def calculateSnCompletion(vplan, sn):
    ''' Computes the S/N completion percentage '''
    # Something is wrong here...
    if vplan == 0: return 1

    try:
        sn_completion = min([1, sn / (3136*vplan)])
        return sn_completion
    except:
        raise RuntimeError("ERROR: unable to calculate S/N completion for vplan: %d, sn: %d\n%s" %\
        (vplan, sn, sys.exc_info()))

def read_design():
    #Create Table
    mytable = Table.read('design.csv',delimiter='|')
    print "Design Table Read Done"
    return(mytable)

def read_plates():
    #Create Table
    mytable = Table.read('plates.csv',delimiter='|')
    print "Plates Table Read Done"
    #Note current_survey_mode Null or 1 APOGEE-led 2 or 3 MaNGA-led
    return(mytable)
    
def read_visits():
    #Create Table
    mytable = Table.read('visits.csv',delimiter='|')
    print "Visits Table Read Done"
    return(mytable)

def read_planned():
    #Create Table
    mytable = Table.read('planned_visits.csv',format='ascii')
    print "Planned Visits Table Read Done"
    return(mytable)

def current_planned(plate_tab,visit_tab,design_tab):
    plate_list = list()
    name_list = list()
    design_list = list()
    ra_list = list()
    dec_list = list()
    lst_list = list()
    ha_list = list()
    maxha_list = list()
    minha_list = list()
    loc_list = list()   
    cohort_list = list()
    visits_list = list()
    planned_list = list()
    complete_list = list()
    sn2_list = list()
    type_list = list()
    mjd_list =list()    
    
    for plate in range(len(plate_tab)):
        
        #Get rid of MaNGA-led plates
        if (not (plate_tab['current_survey_mode_pk'][plate] is np.ma.masked 
                 or plate_tab['current_survey_mode_pk'][plate] == 1)):
            continue
        
        plate_list.append(plate_tab['plate_id'][plate])
        name_list.append(plate_tab['name'][plate])
        loc_list.append(plate_tab['location_id'][plate])
        ra_list.append(round(plate_tab['center_ra'][plate]/15.0,4))
        dec_list.append(plate_tab['center_dec'][plate])
        ha_list.append(plate_tab['hour_angle'][plate])
        maxha_list.append(plate_tab['ha_observable_max'][plate])
        minha_list.append(plate_tab['ha_observable_min'][plate])
        
        design_list.append(plate_tab['design_pk'][plate])
        
        #Search for design information
        thisdesign = design_tab[(design_tab['plate_id'] == plate_tab['plate_id'][plate])]

        if (len(thisdesign) == 0):
            print("Plate: {} is missing a design!".format(plate_tab['plate_id'][plate]))
            sys.exit()
        
        #cohort
        tmp = thisdesign['array_to_string'][0].split(',')
        cohort_list.append(100*int(tmp[0]) + 10*int(tmp[2]) + int(tmp[1]))
        planned_list.append(tmp[3])
        type_list.append(tmp[4])
        
        #Search for visit info    
        plate_visits_tab = visit_tab[(visit_tab['plate_id'] == plate_tab['plate_id'][plate])]
        if len(plate_visits_tab) == 0:
            visits_list.append(0)
            sn2_list.append(0.0)
            complete_list.append(0.0)
            mjd_list.append('')
        else:
            #Cycle through visits and see which are good.
            plate_visits = 0
            plate_sn2 = 0
            plate_mjd = list()
            for visit in range(len(plate_visits_tab)):    
                #Determine whether the visit is good
                if (plate_visits_tab['qlcount'][visit] >= 2 or plate_visits_tab['redcount'][visit] >= 2): 
                    plate_visits = plate_visits + 1
                    plate_mjd.append(str(plate_visits_tab['mjd'][visit]))
                    #Add S/N info
                    if(plate_visits_tab['redcount'][visit] >= 2):
                        plate_sn2 = plate_sn2 + plate_visits_tab['redsum'][visit]
                    elif(plate_visits_tab['qlcount'][visit] >=2):
                        plate_sn2 = plate_sn2 + plate_visits_tab['qlsum'][visit]
            visits_list.append(plate_visits)
            sn2_list.append(plate_sn2)
            complete_list.append(0.0)
            mjd_list.append(','.join(plate_mjd))   
        
        
        #Deal with going arcross 360 / 0
        lst = (plate_tab['center_ra'][plate] + plate_tab['hour_angle'][plate])/15.0
        if (lst > 24.0): lst = lst - 24
        if (lst < 0.0): lst = lst + 24
        
        lst_list.append(round(lst,4))
        
    output_tab = Table()
     
    #Make output table
    output_tab['plate_id'] = plate_list
    output_tab['field'] = name_list
    output_tab['loc_id'] = loc_list
    output_tab['cohort'] = cohort_list
    output_tab['design_id'] = design_list
    output_tab['ra'] = ra_list 
    output_tab['dec'] = dec_list
    output_tab['lst'] = lst_list 
    output_tab['ha'] = ha_list 
    output_tab['maxha'] = maxha_list 
    output_tab['minha'] = minha_list 
    output_tab['visits'] = visits_list 
    output_tab['planned'] = planned_list
    output_tab['complete'] = complete_list
    output_tab['sn2'] = sn2_list
    output_tab['type'] = type_list 
    output_tab['mjds'] = mjd_list 
    
    #Now correct visit entries for shared locations+cohorts
    fixed_visit_list = list()
    fixed_sn2_list = list()
    fixed_complete_list = list()
    fixed_mjd_list = list()
    for row in range(len(output_tab)):
        same_loc_tab = output_tab[(output_tab['loc_id'] == output_tab['loc_id'][row]) & 
                                         (output_tab['cohort'] == output_tab['cohort'][row])]
        
        #To check a given loc_id
        #if(output_tab['loc_id'][row] == 5021):
        #    print(same_loc_tab)
        row_visits = 0
        row_sn2 = 0
        row_mjds = list() 
        for plate in range(len(same_loc_tab)):
            row_visits = row_visits + same_loc_tab['visits'][plate]
            row_sn2 = row_sn2 + same_loc_tab['sn2'][plate]
            if(len(same_loc_tab['mjds'][plate]) > 0):
                row_mjds.append(same_loc_tab['mjds'][plate])
        fixed_visit_list.append(row_visits)
        fixed_sn2_list.append(row_sn2)
        fixed_complete_list.append(round(completion(float(same_loc_tab['planned'][plate]),
                                              row_visits,row_sn2,None),4))
        fixed_mjd_list.append(','.join(row_mjds))
        
    output_tab['visits'] = fixed_visit_list
    output_tab['sn2'] = fixed_sn2_list
    output_tab['complete'] = fixed_complete_list
    #I remove the column so that an appropriate column width can be selected
    output_tab.remove_columns('mjds')
    output_tab['mjds'] = fixed_mjd_list
    
    #Write out table
    output_tab.write('plate_visits.csv',format='ascii')   
    
    return(output_tab)

def temporal_stats(visti_tab,comb_tab):
    return(0)       

# -------------
# Main Function
# -------------
def current_visits_main():
    
    plate_list = list()
    mjd_list = list()
    lst_list = list()
    lst_listy1 = list()
    lst_listy2 = list()
    anc_list = list()
    apok_list = list()
    bulge_list = list()
    clus_list = list()
    disk_list = list()
    halo_list = list()
    goal_list = list()
    sat_list = list()
    sn2_list = list()
    
    #Read plates file
    plate_tab = read_plates()
    #Read visits file
    visit_tab = read_visits()
    #Read planned file
    plan_tab = read_planned()
    #Read design file
    design_tab = read_design()
    
    comb_tab = current_planned(plate_tab,visit_tab,design_tab)
    
    for visit in range(len(visit_tab)):
        #Determine whether the visit is good
        if (visit_tab['qlcount'][visit] < 2 and visit_tab['redcount'][visit] < 2): 
            continue 
        
        hold = plate_tab[(plate_tab['plate_id'] == visit_tab['plate_id'][visit])]
        chold = comb_tab[(comb_tab['plate_id'] == visit_tab['plate_id'][visit])]
        
        if len(hold) == 0:
            continue
        #Get rid of MaNGA-led plates
        if (not (hold['current_survey_mode_pk'][0] is np.ma.masked 
                 or hold['current_survey_mode_pk'][0] == 1)):
            continue
            
        plate_list.append(hold['plate_id'][0])
        mjd_list.append(visit_tab['mjd'][visit])
        #Add S/N info
        if(visit_tab['redcount'][visit] >= 2):
            plate_sn2 = visit_tab['redsum'][visit]
        elif(visit_tab['qlcount'][visit] >=2):
            plate_sn2 = visit_tab['qlsum'][visit]
        sn2_list.append(plate_sn2)      
        
        #Determine type
        if (chold['type'] == 'anc'):
            anc_list.append(visit_tab['mjd'][visit])
        if (chold['type'] == 'kep_apokasc'):
            apok_list.append(visit_tab['mjd'][visit])    
        if (chold['type'] == 'bulge' or chold['type'] =='rrlyr'):
            bulge_list.append(visit_tab['mjd'][visit])
        if (chold['type'] == 'cluster_oc' or chold['type'] == 'cluster_gc'):
            clus_list.append(visit_tab['mjd'][visit])
        if (chold['type'] == 'disk' or chold['type'] == 'disk1' or chold['type'] == 'disk2'):
            disk_list.append(visit_tab['mjd'][visit])
        if (chold['type'] == 'halo' or chold['type'] == 'halo_stream'):
            halo_list.append(visit_tab['mjd'][visit])
        if (chold['type'] == 'k2' or chold['type'] == 'kep_koi' or chold['type'] == 'substellar' or chold['type'] == 'yso'):
            goal_list.append(visit_tab['mjd'][visit])
        if (chold['type'] == 'halo_dsph'):
            sat_list.append(visit_tab['mjd'][visit])   
        
        #Deal with going arcross 360 / 0
        lst = (hold['center_ra'][0] + hold['hour_angle'][0])/15.0
        
        if (lst > 24.0): lst = lst - 24
        if (lst < 0.0): lst = lst + 24
        
        lst_list.append(lst)
        if (visit_tab['mjd'][visit] < 57210):
            lst_listy1.append(lst)
        else:
            lst_listy2.append(lst)
    
    output_tab = Table()
     
    output_tab['plate_id'] = plate_list 
    output_tab['mjd'] = mjd_list 
    output_tab['lst'] = lst_list 
    
    #Write out table
    output_tab.write('visit_lst.csv',format='ascii')
    
    #Get the most recent data
    smjd_list = sorted(mjd_list)
    #startmjd = smjd_list[0]
    startmjd = 56840
    endmjd = smjd_list[-1] + 1
    mjddiff = endmjd - startmjd
    
    #Get date from final MJD converted to SJD
    jddate = Time(endmjd-2,format='mjd')
    date = (jddate.datetime).strftime('%m/%d/%y')
       
    #Make Plot
    pp = PdfPages('visit_lst.pdf')
    
    #LST 1 hour visit
    pl.hist(lst_list,bins=24,range=(0,24),color="#58ACFA")
    pl.xlabel("LST")
    pl.ylabel("Number of visits")
    pl.title("Current ({}) Number of Visits by LST (1 hour bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
    pl.plot(plan_tab['mid'],plan_tab['visits'],color="r",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*2.0,color="b",linewidth=2.0)
    pl.legend(('Plan Year1','Plan Year2','Actual'))
    pp.savefig()
    pl.clf()
    
    #LST 1 hour visit by year
    pl.hist([lst_listy1,lst_listy2],bins=24,range=(0,24),stacked=True,color=["pink","#58ACFA"],rwidth=1.0)
    pl.xlabel("LST")
    pl.ylabel("Number of visits")
    pl.title("Current ({}) Number of Visits by LST (1 hour bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
    pl.plot(plan_tab['mid'],plan_tab['visits'],color="r",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*2.0,color="b",linewidth=2.0)
    pl.legend(('Plan Year1','Plan Year2','Actual Year1', 'Actual Year2'))
    pp.savefig()
    pl.clf()
    
    #LST 20 minute bins
    pl.hist(lst_list,bins=72,range=(0,24),color="#58ACFA")
    pl.xlabel("LST")
    pl.ylabel("Number of visits")
    pl.title("Current ({}) Number of Visits by LST (20 min bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
#    pp.savefig()
    pl.clf()
    cbins = m.ceil(mjddiff/15.0)
    pl.hist(mjd_list,bins=cbins,range=(startmjd,endmjd),color="#58ACFA")
    pl.title("Current ({}) Number of Visits by MJD (15 day bins)".format(date))
    pl.xlabel("MJD")
    pl.ylabel("Number of visits")
    pl.xlim(56820,56820+((cbins+5)*15.0))
    pp.savefig()
    pl.clf()
    cbins = m.ceil(mjddiff/5.0)
    pl.hist(mjd_list,bins=cbins,range=(startmjd,endmjd),color="#58ACFA")
    pl.title("Current ({}) Number of Visits by MJD (5 day bins)".format(date))
    pl.xlabel("MJD")
    pl.ylabel("Number of visits")
    pl.xlim(56820,56820+((cbins+15)*5.0))
    pp.savefig()
    pl.clf()
    
    #Plot efficiency
    odate = ['2014-07','2014-09','2014-10','2014-11','2014-12','2015-01','2015-02'
                     ,'2015-03','2015-04','2015-05','2015-06','2015-07','2015-08'
                     ,'2015-09','2015-10','2015-11','2015-12']
                     
    #Create Date Time format
    wdate = [datetime.datetime.strptime(x,'%Y-%m') for x in odate]
    weff = [36.74,25.73,70.81,21.18,43.95,29.75,69.00,49.73,66.95,52.36,24.53,34.72
            ,46.52,53.97,41.73,57.30,27.81]

    fig, ax = pl.subplots()
    
    ax.plot(wdate,weff,'o-')
    mindate = min(wdate) - datetime.timedelta(days=15) 
    maxdate = max(wdate) + datetime.timedelta(days=15) 
    pl.plot((mindate,maxdate), (45,45))
    pl.xlim(mindate,maxdate)
    pl.ylim(0,100)
    pl.title('Percent On-Sky Time by Date')
    pl.ylabel("Percent On-Sky Time")
    pl.legend(('Actual', 'Projected'))
    
    fig.autofmt_xdate()
    
#    pp.savefig()
    
#     pl.clf()
#     pl.hist(lst_list,bins=24,range=(0,24),color="#58ACFA")
#     pl.xlabel("LST")
#     pl.ylabel("Number of visits")
#     pl.title("Current ({}) Number of Visits by LST (1 hour bins)".format(date))
#     #pl.xlim(-1,24)
#     #pl.savefig('lst_visit.png',dpi=400)
#     pl.plot(plan_tab['mid'],plan_tab['visits'],color="r",linewidth=2.0)
#     pl.legend(('Planned','Actual'))
#     pp.savefig()
#     pp.savefig()
    
    pp.close()
    
    #Get visits per day
    (mjd_num,mjd_bin)=np.histogram(mjd_list,range=(startmjd,endmjd),bins=mjddiff)
    (anc_num,anc_bin)=np.histogram(anc_list,range=(startmjd,endmjd),bins=mjddiff)
    (apok_num,apok_bin)=np.histogram(apok_list,range=(startmjd,endmjd),bins=mjddiff)
    (bulge_num,bulge_bin)=np.histogram(bulge_list,range=(startmjd,endmjd),bins=mjddiff)
    (clus_num,clus_bin)=np.histogram(clus_list,range=(startmjd,endmjd),bins=mjddiff)
    (disk_num,disk_bin)=np.histogram(disk_list,range=(startmjd,endmjd),bins=mjddiff)
    (halo_num,halo_bin)=np.histogram(halo_list,range=(startmjd,endmjd),bins=mjddiff)
    (goal_num,goal_bin)=np.histogram(goal_list,range=(startmjd,endmjd),bins=mjddiff)
    (sat_num,sat_bin)=np.histogram(sat_list,range=(startmjd,endmjd),bins=mjddiff)
    
    mjd_cum = np.cumsum(mjd_num)
    
    anc_cum = np.cumsum(anc_num)
    apok_cum = np.cumsum(apok_num)
    bulge_cum = np.cumsum(bulge_num)
    clus_cum = np.cumsum(clus_num)
    disk_cum = np.cumsum(disk_num)
    halo_cum = np.cumsum(halo_num)
    goal_cum = np.cumsum(goal_num)
    sat_cum = np.cumsum(sat_num)
    
    #Get SN2 per day
    sn2_num = list()
    sn2_mjd = list()
    sn2_na = np.array(sn2_list)
    mjd_na = np.array(mjd_list)
    
    for mjd in range(startmjd,endmjd):
        sn2_mjd.append(mjd)
        #Sum up the SN2 on each MJD
        sn2_num.append(np.sum(sn2_na[(mjd_na == mjd)]))

    sn2_cum =np.cumsum(sn2_num)
    
    outputf = open('mjd.hist','w')
    
    for i in range(len(mjd_num)):
        outputf.write("{} {} {} {} {} {} {} {} {} {} {} {}\n".format(mjd_bin[i],mjd_num[i],mjd_cum[i]
                                                      ,anc_cum[i],apok_cum[i],bulge_cum[i],clus_cum[i]
                                                      ,disk_cum[i],halo_cum[i],goal_cum[i],sat_cum[i]
                                                      ,sn2_cum[i]))
        
    outputf.close() 
    print "Total APOGEE-led visits: {}".format(len(mjd_list))
    #Get visits per day per design type
    
    
    
    print "Success!"
    return(0)

if __name__ == '__main__':
    ret = current_visits_main()
    sys.exit(ret)
    
##
#@mainpage
 #@copydetails  current_visits
    