#!/usr/bin/python
'''
Brief Description

Detailed Description

@package quickred
@author deleenm
@version \e \$Revision$
@date \e \$Date$

Usage: quickred.py
'''

# -----------------------------
# Standard library dependencies
# -----------------------------
import datetime
import sys
# -------------------
# Third-party imports
# -------------------
from astropy.table import Table
from astropy.time import Time
import matplotlib.pyplot as pl
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

def current_sn(plate_tab,visit_tab,design_tab):
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
    qlsn2_list = list()
    redsn2_list = list()
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
        ha_list.append(round(plate_tab['hour_angle'][plate]/15.0,4))
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
        planned_list.append(int(tmp[3]))
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
            qlplate_sn2 = 0
            redplate_sn2 = 0
            plate_mjd = list()
            
            for visit in range(len(plate_visits_tab)):    
                #Determine whether the visit is good
                if (plate_visits_tab['qlcount'][visit] >= 2 or plate_visits_tab['redcount'][visit] >= 2): 
                    plate_visits = plate_visits + 1
                    plate_mjd.append(str(plate_visits_tab['mjd'][visit]))
                    #Add S/N info
                    if(plate_visits_tab['redcount'][visit] >= 2):
                        plate_sn2 = plate_sn2 + plate_visits_tab['redsum'][visit]
                        redplate_sn2 = redplate_sn2 + plate_visits_tab['redsum'][visit]          
                    elif(plate_visits_tab['qlcount'][visit] >=2):
                        plate_sn2 = plate_sn2 + plate_visits_tab['qlsum'][visit]
                    if(plate_visits_tab['qlcount'][visit] >=2):    
                        qlplate_sn2 = qlplate_sn2 + plate_visits_tab['qlsum'][visit]
                             
            visits_list.append(plate_visits)
            sn2_list.append(plate_sn2)
            qlsn2_list.append(qlplate_sn2)
            redsn2_list.append(redplate_sn2)
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
    output_tab['qlsn2'] = qlsn2_list
    output_tab['redsn2'] = redsn2_list
    output_tab['type'] = type_list 
    output_tab['mjds'] = mjd_list 
    
    #Now correct visit entries for shared locations+cohorts
    fixed_visit_list = list()
    fixed_sn2_list = list()
    fixed_qlsn2_list = list()
    fixed_redsn2_list = list()
    fixed_complete_list = list()
    fixed_mjd_list = list()
    for row in range(len(output_tab)):
        same_loc_tab = output_tab[(output_tab['loc_id'] == output_tab['loc_id'][row]) & 
                                         (output_tab['cohort'] == output_tab['cohort'][row])]
        
        row_visits = 0
        row_sn2 = 0
        row_qlsn2 = 0
        row_redsn2 = 0
        row_mjds = list() 
        for plate in range(len(same_loc_tab)):
            row_visits = row_visits + same_loc_tab['visits'][plate]
            row_sn2 = row_sn2 + same_loc_tab['sn2'][plate]
            row_qlsn2 = row_qlsn2 + same_loc_tab['qlsn2'][plate]
            row_redsn2 = row_redsn2 + same_loc_tab['redsn2'][plate]
            if(len(same_loc_tab['mjds'][plate]) > 0):
                row_mjds.append(same_loc_tab['mjds'][plate])
        fixed_visit_list.append(row_visits)
        fixed_sn2_list.append(row_sn2)
        fixed_qlsn2_list.append(row_qlsn2)
        fixed_redsn2_list.append(row_redsn2)
        fixed_complete_list.append(round(completion(float(same_loc_tab['planned'][plate]),
                                              row_visits,row_sn2,None),4))
        fixed_mjd_list.append(','.join(row_mjds))
        
    output_tab['visits'] = fixed_visit_list
    output_tab['sn2'] = fixed_sn2_list
    output_tab['qlsn2'] = fixed_qlsn2_list
    output_tab['redsn2'] = fixed_redsn2_list
    output_tab['complete'] = fixed_complete_list
    #I remove the column so that an appropriate column width can be selected
    output_tab.remove_columns('mjds')
    output_tab['mjds'] = fixed_mjd_list
    
    #Write out table
    output_tab.write('sn_visits.csv',format='ascii',overwrite=True)   
    
    return(output_tab)       

def temporal_sn(comb_tab,visit_tab):
    
    mjd_dict = dict()
    
    mjd_dict['plate'] = list()
    mjd_dict['mjd'] = list()
    mjd_dict['lst'] = list()
    mjd_dict['loc_id'] = list()
    mjd_dict['cohort'] = list()
    mjd_dict['sn2'] = list()
    mjd_dict['qlsn2'] = list()
    mjd_dict['redsn2'] = list() 
        
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
        mjd_dict['mjd'].append(visit_tab['mjd'][visit])
        
        
        #Add S/N info
        if(visit_tab['redcount'][visit] >= 2):
            plate_sn2 = visit_tab['redsum'][visit]
            plate_redsn2 = visit_tab['redsum'][visit]
        elif(visit_tab['qlcount'][visit] >=2):
            plate_sn2 = visit_tab['qlsum'][visit]
        if(visit_tab['qlcount'][visit] >=2):
            plate_qlsn2 = visit_tab['qlsum'][visit]
        
        mjd_dict['sn2'].append(plate_sn2)      

        
# -------------
# Main Function
# -------------
def quickred_main():
    
    #Read plates file
    plate_tab = read_plates()
    #Read visits file
    visit_tab = read_visits()
    #Read design file
    design_tab = read_design()
    
    comb_tab = current_sn(plate_tab,visit_tab,design_tab)
    (mjd_tab,mjd_dict) = temporal_sn(comb_tab,visit_tab)
    
    #Make Plot
    pp = PdfPages('quickred.pdf')
    
    pp.close()

if __name__ == '__main__':
    quickred_main()
    
##
#@mainpage
 #@copydetails  quickred
    