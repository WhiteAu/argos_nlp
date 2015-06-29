#
# Copyright (c) 2014-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

'''author@esilgard'''
__version__='PathStageT1.0'
import re,sys
import global_strings

def get(disease_group,dictionary,specimen,module_sections):
    '''
        extract the PathStageT (size/location of tumor)from normal cased text of the pathology report        
    '''
    return_dictionary={global_strings.NAME:"PathStageT",global_strings.VALUE:None,global_strings.CONFIDENCE:0.0,global_strings.VERSION:__version__,
                           global_strings.STARTSTOPS:[]}
    ##currently this will only return one (the last) T value found   
    for section in sorted(dictionary):
        section_specimen=section[3]
        line_onset=section[2]
        header=section[1]
        if re.search(module_sections,header):
            for index,results in sorted(dictionary[section].items(),key=lambda x: int(x[0])):
                t_stage=re.match('.*(pT[012345][abc]?).*',results,re.DOTALL)
                if t_stage:
                    return_dictionary[global_strings.VALUE]=t_stage.group(1)       
                    return_dictionary[global_strings.STARTSTOPS].append({global_strings.START:t_stage.start(),global_strings.STOP:t_stage.end()})
    return ([return_dictionary],list)
    
