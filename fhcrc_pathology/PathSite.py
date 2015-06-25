#
# Copyright (c) 2014-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

'''author@esilgard'''
__version__='PathSite1.0'

import os
import re
import global_strings
import PathFindNumNodes
path= os.path.dirname(os.path.realpath(__file__))
   
def get(disease_group,dictionary,specimen,module_sections):   
    '''
    return a list of dictionaries of each PathSite (location of tumor) per specimen
    from normal cased text of the pathology report
    get specific number of lymph nodes involved if the site is lymph (from imported PathFindNumNodes module)
    '''
   
    ## a list of sites and their standardized forms from the disease relevent sites file and the general sites field##
    def make_lists(disease_group):
        sites=[]
        standardizations={}       
        for line in open(path+'/'+disease_group+'sites.txt','r').readlines():
            site_list=line.split(';')
            for h in site_list:
                h=h.strip().lower()
                standardizations[h]=site_list[0].strip()
                sites.append(h)                
        sites=sorted(sites,key=lambda x: len(x),reverse=True)
       
        return sites,standardizations  
    
    ###############################################################################################################
    
        
    ###############################################################################################################    
    def get_site(site_list,standardizations,specimen):        
        specimen_site_list=[]
        specimen_start_stops_list=[]
        numNodes=None
        numNodesPos=None
        nodes_start_stops=[]
        ######### - don't think we need this not specimen anymore...
        
        for section in dictionary:
            section_specimen=section[3]
            line_onset=section[2]
            header=section[1]            
            if re.search(module_sections[1],header):               
                for index,text in dictionary[section].items():               
                    ## meant to weed out references to literature/papers - picking up publication info like this: 2001;30:1-14. ##
                    ## these can contain confusing general statements about the cancer and/or patients in general ##
                    if re.search('[\d]{4}[;,][ ]*[\d]{1,4}:[\d\-]{1,6}',text):pass               
                    else:
                        text=text.lower()
                        # rudimentry tokenization keeps character count in order to create offsets
                        text=re.sub('[,:;\\\/\-]',' ',text); text=re.sub('[.] ', '  ',text)      ## this should keep decimal places and throw out periods                    
                        for each_site in site_list:                        
                            for each_match in re.finditer('^.*( |^|\")('+each_site+')( |$|\").*',text,re.DOTALL):                                
                                if standardizations[each_site] not in specimen_site_list:                                    
                                    specimen_site_list.append(standardizations[each_site])
                                if 'Lymph' in standardizations[each_site]:
                                    numNodes,numNodesPos=PathFindNumNodes.get(section,text,specimen)
                                specimen_start_stops_list.append({global_strings.START:each_match.start(2)+line_onset,global_strings.STOP:each_match.end(2)+line_onset})                                
        if specimen_site_list:            
            return {global_strings.NAME:"PathSite",global_strings.KEY:specimen,global_strings.TABLE:global_strings.PATHOLOGY_TABLE,global_strings.VALUE:';'.join(set(specimen_site_list)),
                    global_strings.CONFIDENCE:("%.2f" % .85), global_strings.VERSION:__version__,global_strings.STARTSTOPS:specimen_start_stops_list},numNodes,numNodesPos
        else: return None,None,None
                                  
###################################################################################################################
    
    try:
        disease_group_sites,disease_group_standardizations=make_lists(disease_group+'/')
        general_sites,general_standardizations=make_lists('')
    except: return ([{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'ERROR: could not access site file -- PathSite not completed'}],Exception)
    
    
    return_dictionary_list=[]    
    site_list=[]
    start_stops_list=[]              
    specimen_site_dictionary,numNodes,numNodesPos=get_site(disease_group_sites,disease_group_standardizations,specimen)            
    if not specimen_site_dictionary:               
        specimen_site_dictionary,numNodes,numNodesPos=get_site(general_sites,general_standardizations,specimen)    
    if specimen_site_dictionary:
        return_dictionary_list.append(specimen_site_dictionary)
        if numNodes and numNodesPos:                    
            return_dictionary_list.append(numNodes)
            return_dictionary_list.append(numNodesPos)

    
    ## if there were no disease specific sites found for this specimen, look for general sites throughout the specimen report ##
    else:
        overall_site_dictionary,numNodes,numNodesPos=get_site(general_sites,general_standardizations,specimen)
        if overall_site_dictionary:            
            return_dictionary_list.append({global_strings.NAME:"PathSite",global_strings.TABLE:global_strings.PATHOLOGY_TABLE,global_strings.VALUE:overall_site_dictionary[global_strings.VALUE],
                global_strings.CONFIDENCE:0.75,global_strings.VERSION:__version__, global_strings.STARTSTOPS:overall_site_dictionary[global_strings.STARTSTOPS]})
            if numNodes and numNodesPos:                    
                return_dictionary_list.append(numNodes)
                return_dictionary_list.append(numNodesPos)
    return (return_dictionary_list,list)
