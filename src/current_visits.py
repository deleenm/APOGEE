#!/usr/bin/env python
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
import argparse
import astropy.table
from astropy.table import Table
from astropy.time import Time
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import os.path
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import astropy as astpy
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
    mytable = Table.read('../design.csv',delimiter='|')
    print("Design Table Read Done")
    return(mytable)

def read_plates():
    #Create Table
    mytable = Table.read('../plates.csv',delimiter='|')
    print("Plates Table Read Done")
    #Note current_survey_mode Null or 1 APOGEE-led 2 or 3 MaNGA-led
    return(mytable)
    
def read_visits():
    #Create Table
    mytable = Table.read('../visits.csv',delimiter='|')
    print("Visits Table Read Done")
    return(mytable)

def read_planned():
    #Create Table
    mytable = Table.read('../planned_visits.csv',format='ascii')
    print("Planned Visits Table Read Done")
    return(mytable)

def current_planned(plate_tab,visit_tab,design_tab,south=False):
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
    ploc_list = list()
    temp_list = list()
    exptime_list = list()
    
    for plate in range(len(plate_tab)):
        
        #Get rid of MaNGA-led plates
        if (not (plate_tab['current_survey_mode_pk'][plate] is np.ma.masked 
                 or plate_tab['current_survey_mode_pk'][plate] == 1)):
            continue
        
        if(south):
            #Get rid of commissioning plates
            if (plate_tab['location_id'][plate] == -1):
                print("Plateid {} is a commissioning plate. Ignoring...".format(plate_tab['plate_id'][plate]))
                continue
        
            #Get rid of early plates
            if (plate_tab['plate_id'][plate] < 9280 and plate_tab['plate_id'][plate] >= 9265 ):
                print("Plateid {} has truncated design info. Ignoring...".format(plate_tab['plate_id'][plate]))
                continue
        
        #Search for design information
        thisdesign = design_tab[(design_tab['plate_id'] == plate_tab['plate_id'][plate])]

        if (len(thisdesign) == 0):
            print("Plate: {} is missing a design!".format(plate_tab['plate_id'][plate]))
            sys.exit()
        
        #cohort
        tmp = thisdesign['array_to_string'][0].split(',')
        #Get rid of external plates
        if(south):
            if(int(tmp[0]) < 0):
                continue
            if(len(tmp) == 8):
                if(tmp[7] != 'Main'):
                    continue
            else:
                if(tmp[6] != 'Main'):
                    continue
        
        cohort_list.append(100*int(tmp[0]) + 10*int(tmp[2]) + int(tmp[1]))
        planned_list.append(int(tmp[3]))
        type_list.append(tmp[4])
        if(south):
            if (len(tmp) == 8):
                exptime_list.append(tmp[6])
            else:
                exptime_list.append(500)
        else:
            if (len(tmp) == 7):
                exptime_list.append(tmp[6])
            else:
                exptime_list.append(500)
        
        plate_list.append(plate_tab['plate_id'][plate])
        name_list.append(plate_tab['name'][plate])
        loc_list.append(plate_tab['location_id'][plate])
        ra_list.append(round(plate_tab['center_ra'][plate]/15.0,4))
        dec_list.append(plate_tab['center_dec'][plate])
        ha_list.append(round(plate_tab['hour_angle'][plate]/15.0,4))
        maxha_list.append(plate_tab['ha_observable_max'][plate])
        minha_list.append(plate_tab['ha_observable_min'][plate])
        
        design_list.append(plate_tab['design_pk'][plate])
        ploc_list.append(plate_tab['label'][plate])
        temp_list.append(plate_tab['temperature'][plate])        
        
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
    output_tab['temp'] = temp_list
    output_tab['exptime'] = exptime_list
    output_tab['visits'] = visits_list 
    output_tab['planned'] = planned_list
    output_tab['complete'] = complete_list
    output_tab['sn2'] = sn2_list
    output_tab['type'] = type_list 
    output_tab['location'] = ploc_list
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
    output_tab.write('../plate_visits.csv',format='ascii',overwrite=True)   
    
    return(output_tab)       

def temporal_change(comb_tab,visit_tab):
    
    mjd_dict = dict()
    
    mjd_dict['plate'] = list()
    mjd_dict['mjd'] = list()
    mjd_dict['lst'] = list()
    mjd_dict['lsty1'] = list()
    mjd_dict['lsty2'] = list()
    mjd_dict['lsty3'] = list()
    mjd_dict['anc'] = list()
    mjd_dict['apok'] = list()
    mjd_dict['bulge'] = list()
    mjd_dict['clus'] = list()
    mjd_dict['disk'] = list()
    mjd_dict['halo'] = list()
    mjd_dict['goal'] = list()
    mjd_dict['sat'] = list()
    mjd_dict['sn2'] = list()
    mjd_dict['loc_id'] = list()
    mjd_dict['cohort'] = list()
    mjd_dict['planned'] = list()
    mjd_dict['2vsn2'] = list()
    mjd_dict['2vlst'] = list()
    mjd_dict['2vlsty1'] = list()
    mjd_dict['2vlsty2'] = list()
    mjd_dict['2vlsty3'] = list()
    mjd_dict['lstno2v'] = list()
    mjd_dict['lsty1no2v'] = list()
    mjd_dict['lsty2no2v'] = list()
    mjd_dict['lsty3no2v'] = list()
        
    for visit in range(len(visit_tab)):
        #Determine whether the visit is good
        if (visit_tab['qlcount'][visit] < 2 and visit_tab['redcount'][visit] < 2): 
            continue 
        
        chold = comb_tab[(comb_tab['plate_id'] == visit_tab['plate_id'][visit])]
        
        if len(chold) == 0:
            continue
                    
        mjd_dict['plate'].append(chold['plate_id'][0])
        mjd_dict['loc_id'].append(chold['loc_id'][0])
        mjd_dict['cohort'].append(chold['cohort'][0])
        mjd_dict['planned'].append(chold['planned'][0])
        mjd_dict['mjd'].append(visit_tab['mjd'][visit])
        
        twovisit_sn2 = 0
        
        #Add S/N info
        if(visit_tab['redcount'][visit] >= 2):
            plate_sn2 = visit_tab['redsum'][visit]
            if(visit_tab['redcount'][visit] == 2):
                twovisit_sn2 = visit_tab['redsum'][visit]
        elif(visit_tab['qlcount'][visit] >=2):
            plate_sn2 = visit_tab['qlsum'][visit]
            if(visit_tab['qlcount'][visit] == 2):
                twovisit_sn2 = visit_tab['qlsum'][visit] 
        
        mjd_dict['sn2'].append(plate_sn2)      
        mjd_dict['2vsn2'].append(twovisit_sn2)
        
        #Determine type
        if (chold['type'] == 'anc'):
            mjd_dict['anc'].append(visit_tab['mjd'][visit])
        if (chold['type'] == 'kep_apokasc'):
            mjd_dict['apok'].append(visit_tab['mjd'][visit])    
        if (chold['type'] == 'bulge' or chold['type'] =='rrlyr'):
            mjd_dict['bulge'].append(visit_tab['mjd'][visit])
        if (chold['type'] == 'cluster_oc' or chold['type'] == 'cluster_gc'):
            mjd_dict['clus'].append(visit_tab['mjd'][visit])
        if (chold['type'] == 'disk' or chold['type'] == 'disk1' or chold['type'] == 'disk2'):
            mjd_dict['disk'].append(visit_tab['mjd'][visit])
        if (chold['type'] == 'halo' or chold['type'] == 'halo_stream'):
            mjd_dict['halo'].append(visit_tab['mjd'][visit])
        if (chold['type'] == 'k2' or chold['type'] == 'kep_koi' or chold['type'] == 'substellar' or chold['type'] == 'yso'):
            mjd_dict['goal'].append(visit_tab['mjd'][visit])
        if (chold['type'] == 'halo_dsph'):
            mjd_dict['sat'].append(visit_tab['mjd'][visit])   
        
        #Deal with going arcross 360 / 0
        lst = (chold['ra'][0] + chold['ha'][0])
        
        if (lst > 24.0): lst = lst - 24
        if (lst < 0.0): lst = lst + 24
        
        lst = round(lst,4)
        
        mjd_dict['lst'].append(lst)
        if (visit_tab['mjd'][visit] < 57210):
            mjd_dict['lsty1'].append(lst)
        elif (visit_tab['mjd'][visit] >= 57210 and visit_tab['mjd'][visit] < 57581):
            mjd_dict['lsty2'].append(lst)
        else:
            mjd_dict['lsty3'].append(lst)
            
        #Deal with 2visit LST
        if(visit_tab['redcount'][visit] >= 2):
            if(visit_tab['redcount'][visit] == 2):
                mjd_dict['2vlst'].append(lst)
                if (visit_tab['mjd'][visit] < 57210):
                    mjd_dict['2vlsty1'].append(lst)
                elif (visit_tab['mjd'][visit] >= 57210 and visit_tab['mjd'][visit] < 57581):
                    mjd_dict['2vlsty2'].append(lst)
                else:
                    mjd_dict['2vlsty3'].append(lst)
            else:
                mjd_dict['lstno2v'].append(lst)
                if (visit_tab['mjd'][visit] < 57210):
                    mjd_dict['lsty1no2v'].append(lst)
                elif (visit_tab['mjd'][visit] >= 57210 and visit_tab['mjd'][visit] < 57581):
                    mjd_dict['lsty2no2v'].append(lst)
                else:
                    mjd_dict['lsty3no2v'].append(lst)
        elif(visit_tab['qlcount'][visit] >= 2):
            if(visit_tab['qlcount'][visit] == 2):
                mjd_dict['2vlst'].append(lst)
                if (visit_tab['mjd'][visit] < 57210):
                    mjd_dict['2vlsty1'].append(lst)
                elif (visit_tab['mjd'][visit] >= 57210 and visit_tab['mjd'][visit] < 57581):
                    mjd_dict['2vlsty2'].append(lst)
                else:
                    mjd_dict['2vlsty3'].append(lst)
            else:
                mjd_dict['lstno2v'].append(lst)
                if (visit_tab['mjd'][visit] < 57210):
                    mjd_dict['lsty1no2v'].append(lst)
                elif (visit_tab['mjd'][visit] >= 57210 and visit_tab['mjd'][visit] < 57581):
                    mjd_dict['lsty2no2v'].append(lst)
                else:
                    mjd_dict['lsty3no2v'].append(lst)
        
                
    output_tab = Table()
    
    #Create a single loc_id+cohort
    loc_cohort_list = ["{}_{}".format(a_, b_) 
                       for a_, b_ in zip(mjd_dict['loc_id'],mjd_dict['cohort'])]
     
    output_tab['plate_id'] = mjd_dict['plate'] 
    output_tab['loc_cohort'] = loc_cohort_list
    output_tab['planned'] = mjd_dict['planned']
    output_tab['mjd'] = mjd_dict['mjd'] 
    output_tab['lst'] = mjd_dict['lst']
    output_tab['sn2'] = mjd_dict['sn2'] 
    output_tab['2vsn2'] = mjd_dict['2vsn2']
    
    #Write out table
    output_tab.write('../visit_lst.csv',format='ascii',overwrite=True)
    return (output_tab,mjd_dict)

def lst_plots(mjd_dict,startmjd,endmjd):
    
    #Read planned file
    plan_tab = read_planned()
    
    mjddiff = endmjd - startmjd
    
    #Get date from final MJD converted to SJD
    jddate = Time(endmjd-2,format='mjd')
    date = (jddate.datetime).strftime('%m/%d/%y')
       
    #Make Plot
    pp = PdfPages('../visit_lst.pdf')
    
    #LST 1 hour visit
    pl.hist(mjd_dict['lst'],bins=24,range=(0,24),color="#58ACFA",edgecolor='black')
    pl.xlabel("LST")
    pl.ylabel("Number of visits")
    pl.title("Current ({}) Number of Visits by LST (1 hour bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
    pl.plot(plan_tab['mid'],plan_tab['visits'],color="r",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*2.0,color="b",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*3.0,color="g",linewidth=2.0)
    pl.legend(('Plan Year1','Plan Year2','Plan Year3','Actual'))
    pp.savefig()
    pl.clf()
    
    #LST 1 hour visit by year
    (yrhst1,yrbins1) = np.histogram(mjd_dict['lsty1'],bins=24,range=(0,24))
    (yrhst2,yrbins2) = np.histogram(mjd_dict['lsty2'],bins=24,range=(0,24))
    (yrhst3,yrbins3) = np.histogram(mjd_dict['lsty3'],bins=24,range=(0,24))
    (yrhst1_no,yrbins1) = np.histogram(mjd_dict['lsty1no2v'],bins=24,range=(0,24))
    (yrhst2_no,yrbins2) = np.histogram(mjd_dict['lsty2no2v'],bins=24,range=(0,24))
    (yrhst3_no,yrbins3) = np.histogram(mjd_dict['lsty3no2v'],bins=24,range=(0,24))
    #Write the data out
    lstout = open('../lstdist.txt','w')
    lstout.write("#LST Proj_visit Yr1_Visit Yr2_Visit Yr3_Visit Yr1_no2exp Yr2_no2exp Yr3_visit_no2exp\n")
    for i in range(24):
        lstout.write("{} {} {} {} {} {} {} {}\n".format(plan_tab['mid'][i],plan_tab['visits'][i],yrhst1[i],
                                                  yrhst2[i],yrhst3[i],yrhst1_no[i],yrhst2_no[i],yrhst3_no[i]))
    lstout.close()
    
    pl.hist([mjd_dict['lsty1'],mjd_dict['lsty2'],mjd_dict['lsty3']],bins=24,range=(0,24),
            stacked=True,color=["pink","#58ACFA",'lightgreen'],rwidth=1.0,edgecolor='black')
    pl.xlabel("LST")
    pl.ylabel("Number of visits")
    pl.title("Current ({}) Number of Visits by LST (1 hour bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
    pl.plot(plan_tab['mid'],plan_tab['visits'],color="r",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*2.0,color="b",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*3.0,color="g",linewidth=2.0)
    pl.legend(('Plan Year1','Plan Year2','Plan Year3','Actual Year1', 'Actual Year2','Actual Year3'),
              ncol=2, fontsize='small')
    pp.savefig()
    pl.clf()
    
    #LST 1 hour visit by year - 2exp visits
    pl.hist([mjd_dict['lsty1no2v'],mjd_dict['lsty2no2v'],mjd_dict['lsty3no2v']],bins=24,range=(0,24),
            stacked=True,color=["pink","#58ACFA",'lightgreen'],rwidth=1.0,edgecolor='black')
    pl.xlabel("LST")
    pl.ylabel("Number of visits - 2 expsure visits")
    pl.title("Current ({}) Number of Visits by LST (1 hour bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
    pl.plot(plan_tab['mid'],plan_tab['visits'],color="r",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*2.0,color="b",linewidth=2.0)
    pl.plot(plan_tab['mid'],plan_tab['visits']*3.0,color="g",linewidth=2.0)
    pl.legend(('Plan Year1','Plan Year2','Plan Year3,','Actual Year1', 'Actual Year2', 'Actual Year3'),
              fontsize='small')
    pp.savefig()
    pl.clf()
    
    #LST 20 minute bins
    pl.hist(mjd_dict['lst'],bins=72,range=(0,24),color="#58ACFA",edgecolor='black')
    pl.xlabel("LST")
    pl.ylabel("Number of visits")
    pl.title("Current ({}) Number of Visits by LST (20 min bins)".format(date))
    #pl.xlim(-1,24)
    #pl.savefig('lst_visit.png',dpi=400)
#    pp.savefig()
    pl.clf()
    cbins = m.ceil(mjddiff/15.0)
    pl.hist(mjd_dict['mjd'],bins=int(cbins),range=(startmjd,endmjd),color="#58ACFA",edgecolor='black')
    pl.title("Current ({}) Number of Visits by MJD (15 day bins)".format(date))
    pl.xlabel("MJD")
    pl.ylabel("Number of visits")
    pl.xlim(56820,56820+((cbins+5)*15.0))
    #pp.savefig()
    pl.clf()
    cbins = m.ceil(mjddiff/5.0)
    pl.hist(mjd_dict['mjd'],bins=int(cbins),range=(startmjd,endmjd),color="#58ACFA",edgecolor='black')
    pl.title("Current ({}) Number of Visits by MJD (5 day bins)".format(date))
    pl.xlabel("MJD")
    pl.ylabel("Number of visits")
    pl.xlim(56820,56820+((cbins+15)*5.0))
    #pp.savefig()
    pl.clf()
    
    
    pp.close()

# -------------
# Main Function
# -------------
def current_visits_main(south=False):
    
    #North or South
    if(south == False):
        print("Running in Northern mode")
    else:
        print("Running in Southern mode")
    
    #Read plates file
    plate_tab = read_plates()
    #Read visits file
    visit_tab = read_visits()
    #Read design file
    design_tab = read_design()
    
    comb_tab = current_planned(plate_tab,visit_tab,design_tab,south=south)
    
    (mjd_tab,mjd_dict) = temporal_change(comb_tab,visit_tab)
    
    #Get the most recent data
    smjd_list = sorted(mjd_dict['mjd'])
    if (south):
        #startmjd = 57700
        startmjd = smjd_list[0]
    else:
        startmjd = 56840
    endmjd = smjd_list[-1] + 1
    mjddiff = endmjd - startmjd
    
    lst_plots(mjd_dict, startmjd, endmjd)
    
    #Get visits per day
    (mjd_num,mjd_bin)=np.histogram(mjd_dict['mjd'],range=(startmjd,endmjd),bins=mjddiff)
    (anc_num,anc_bin)=np.histogram(mjd_dict['anc'],range=(startmjd,endmjd),bins=mjddiff)
    (apok_num,apok_bin)=np.histogram(mjd_dict['apok'],range=(startmjd,endmjd),bins=mjddiff)
    (bulge_num,bulge_bin)=np.histogram(mjd_dict['bulge'],range=(startmjd,endmjd),bins=mjddiff)
    (clus_num,clus_bin)=np.histogram(mjd_dict['clus'],range=(startmjd,endmjd),bins=mjddiff)
    (disk_num,disk_bin)=np.histogram(mjd_dict['disk'],range=(startmjd,endmjd),bins=mjddiff)
    (halo_num,halo_bin)=np.histogram(mjd_dict['halo'],range=(startmjd,endmjd),bins=mjddiff)
    (goal_num,goal_bin)=np.histogram(mjd_dict['goal'],range=(startmjd,endmjd),bins=mjddiff)
    (sat_num,sat_bin)=np.histogram(mjd_dict['sat'],range=(startmjd,endmjd),bins=mjddiff)
    
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
    sn2_na = mjd_tab['sn2']
    sn2_num = list()
    twovisitsn2_na = mjd_tab['2vsn2']
    twovisitsn2_num = list()
    
    mjd_na = mjd_tab['mjd']
    
    
    #Results lists for mjd by mjd results
    dbd_mjd = list()
    dbd_visits = list()
    dbd_sn2 = list()
    dbd_sn2corr = list()
    dbd_complete = list()
    
    for mjd in range(startmjd,endmjd):
        
        #Store Date
        dbd_mjd.append(mjd)
        
        #Sum up the SN2 on each MJD
        sn2_num.append(round(np.sum(sn2_na[(mjd_na == mjd)])/3333.0,4))
        twovisitsn2_num.append(round(np.sum(twovisitsn2_na[(mjd_na == mjd)])/3333.0,4))
        
        #Reset values for this mjd
        mjd_visits = 0
        mjd_sn2 = 0
        mjd_sn2corr = 0
        mjd_comp = 0
             
        #Figure out the completeness for each loc_id+cohort on that day
        this_mjd_tab = mjd_tab['loc_cohort','planned'][(mjd_tab['mjd'] <= mjd)]
        
        #Only work on dates that have observations
        if (len(this_mjd_tab) == 0):
            dbd_visits.append(0)
            dbd_sn2.append(0)
            dbd_sn2corr.append(0)
            dbd_complete.append(0)
            continue
        
        #Get rid of duplicate loc_cohort
        this_mjd_tab = astropy.table.unique(this_mjd_tab)
        
        for row in range(len(this_mjd_tab)):
            
            lc_plan = this_mjd_tab['planned'][row]
            loc_cohort = this_mjd_tab['loc_cohort'] [row]                                   
            
            sn2_col = mjd_tab['sn2'][((mjd_tab['loc_cohort'] == loc_cohort) & (mjd_tab['mjd'] <= mjd))]
            lc_visits = len(sn2_col)
            lc_sn2 = np.sum(sn2_col)
            lc_comp = round(completion(float(lc_plan),lc_visits,lc_sn2,None),4)*lc_plan
            lc_sn2corr = round(calculateSnCompletion(lc_plan, lc_sn2),4)*lc_plan
            
            mjd_visits = mjd_visits + lc_visits
            mjd_sn2 = mjd_sn2 + lc_sn2
            mjd_sn2corr = mjd_sn2corr + lc_sn2corr
            mjd_comp = mjd_comp + lc_comp
            
        dbd_visits.append(mjd_visits)        
        dbd_sn2.append(round(mjd_sn2 / 3333.0,4))
        dbd_sn2corr.append(mjd_sn2corr)
        dbd_complete.append(mjd_comp)
        
    sn2_cum = np.cumsum(sn2_num)
    twovisitsn2_cum = np.cumsum(twovisitsn2_num)
    
    cum_tab = Table()
    cum_tab['mjd'] = mjd_bin[0:-1]
    cum_tab['daynum'] = mjd_num
    cum_tab['cum'] = mjd_cum
    cum_tab['anc'] = anc_cum
    cum_tab['apok'] = apok_cum
    cum_tab['bulge'] = bulge_cum
    cum_tab['clus'] = clus_cum
    cum_tab['disk'] = disk_cum
    cum_tab['halo'] = halo_cum
    cum_tab['goal'] = goal_cum
    cum_tab['sat'] = sat_cum
    cum_tab['visits'] = dbd_visits
    cum_tab['sn2'] = dbd_sn2
    cum_tab['sn2corr'] = dbd_sn2corr
    cum_tab['complete'] = dbd_complete
    cum_tab['2vsn2'] = twovisitsn2_cum
   
    cum_tab.write('../mjd.hist',format='ascii',overwrite=True)
    
    print("Total APOGEE-led visits: {}".format(len(mjd_dict['mjd'])))
    
    print("Success!")
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates files and plots for tracking APOGEE Progress.')
    parser.add_argument('-s','--south', action = 'store_true',help='Run this with the Southern options')
    args = vars(parser.parse_args())
    ret = current_visits_main(south=args['south'])
    sys.exit(ret)
    
##
#@mainpage
 #@copydetails  current_visits
    