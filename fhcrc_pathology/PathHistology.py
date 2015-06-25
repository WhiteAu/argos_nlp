#
# Copyright (c) 2013-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

'''author@esilgard'''
__version__='PathHistology1.0'

import re
import os
import global_strings
path= os.path.dirname(os.path.realpath(__file__))+'/'


#############################################################################################################################################################

#############################################################################################################################################################

def get(disease_group,dictionary,specimen,module_sections):  
    '''
    extract the histology from the lower cased text of the pathology report       
    return a list of dictionaries of each PathHistology (per specimen) (multiple histologies will be comma seperated and revised to seperate PathFinding)
    '''
    if specimen:
        return_dictionary_list=[]
        ## a list of histologies from the disease relevent histology file
        specimen_histologies=[]
        standardizations={}
        try:
            for line in open(path+'/'+disease_group+'/'+'histologies.txt','r').readlines():
                histologies=line.split(';')
                for each_histology in histologies:
                    each_histology=each_histology.strip().lower()
                    standardizations[each_histology]=histologies[0].strip()
                    specimen_histologies.append(each_histology)
            specimen_histologies=sorted(specimen_histologies,key=lambda x: len(x),reverse=True)
            ## append generic carcinoma histology to end of list as a last resort string match
            specimen_histologies.append('carcinoma');standardizations['carcinoma']='Carcinoma'
        except: return ([{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'ERROR: could not access histology file at '+path+'/'+disease_group+'/'+'histologies.txt -- PathHistology not completed'}],Exception)
       
        specimen_histology_list=[]
        specimen_start_stops_list=[]
        for section in sorted(dictionary):        
            section_specimen=section[3]
            line_onset=section[2]
            header=section[1]                
            if re.search(module_sections[1],header):                    
                for index,results in sorted(dictionary[section].items(),key=lambda x: int(x[0])):               
                   
                    ## meant to weed out references to literature/papers - picking up publication info like this: 2001;30:1-14. ##
                    ## these can contain confusing general statements about the cancer and/or patients in general ##
                    if re.search('[\d]{4}[;,][ ]*[\d]{1,4}:[\d\-]{1,6}',results):pass 
                    else:                          
                        text=results.lower()
                        text=re.sub('[.,:;\\\/\-]',' ',text)                            
                        histology,onset,offset=find_histology(text,specimen_histologies)                           
                        if histology:                               
                            specimen_start_stops_list.append({global_strings.START:line_onset+onset,global_strings.STOP:line_onset+offset})                        
                            already_seen=False
                            for each in specimen_histology_list:
                                if standardizations[histology] in each:
                                    already_seen=True
                            if not already_seen:
                                specimen_histology_list.append(standardizations[histology])
                            
        if specimen_histology_list:
            specimen_histology_set=set(specimen_histology_list)
            confidence=.85
            if len(set(specimen_histology_set))>1:confidence=.7
            return_dictionary_list.append({global_strings.NAME:global_strings.PATHOLOGY_TABLE,global_strings.KEY:specimen,global_strings.TABLE:global_strings.FINDING_TABLE,global_strings.VALUE:';'.join(specimen_histology_set),
                                           global_strings.CONFIDENCE:("%.2f" % confidence),global_strings.VERSION:__version__,global_strings.STARTSTOPS:specimen_start_stops_list})
             
        return (return_dictionary_list,list)        
    else:
        ## for now - do not look for any histologies in the full text - this is too general
        return ([],list)
            
## check for the presence of a non-negated string ##
def find_histology(short_text,specimen_histologies):      
    for histo in specimen_histologies:        
        if re.search(r'([\W]|^)'+histo+r'([\W]|$)',short_text):            
            if not re.search(r'( no |negative |free of |against |(hx|history) of | to rule out|preclud)[\w ]{,50}'+histo+r'([\W]|$)',short_text) and \
               not re.search(r'([\W]|^)'+histo+r'[\w ]{,40}( unlikely| not (likely|identif)| negative)',short_text):                
                return (histo,short_text.find(histo),short_text.find(histo)+len(histo))
    return None,None,None
                      
