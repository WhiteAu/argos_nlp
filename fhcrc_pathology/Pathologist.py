#
# Copyright (c) 2014-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

'''author@esilgard'''
__version__='Pathologist1.0'

import re
import global_strings
def get(disease_group,dictionary,specimen,module_sections):
    '''
    extract the (first) pathologist's name from the end of the report
    ----but --- not picking up PhD's
    '''   
    return_dictionary={global_strings.NAME:"Pathologist",global_strings.VALUE:None,global_strings.VERSION:__version__,
                       global_strings.STARTSTOPS:[]}
    for section in dictionary:
        section_specimen=section[3]
        line_onset=section[2]
        header=section[1]            
        if re.search(module_sections,header):       
            ## this match is non greedy so that the INITIAL pathologist signature is picked out
            name_match=re.match('.*?\n([A-Za-z\'\-,. ]+) MD(, PhD)?[ ]*\n[ ]*Pathologist[ ]*\n.*',dictionary[section],re.DOTALL)            
            if name_match:
                return_dictionary[global_strings.VALUE]=name_match.group(1)
                return_dictionary[global_strings.CONFIDENCE]=1.0        
                return_dictionary[global_strings.STARTSTOPS].append({global_strings.START:name_match.start(1),global_strings.STOP:name_match.end(1)})       
    return ([return_dictionary],list)
    
