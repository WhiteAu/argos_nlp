#
# Copyright (c) 2014-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

'''author@esilgard'''
__version__='PathStageM1.0'
import re
import global_strings

def get(disease_group,dictionary,specimen,module_sections):
   
    '''
    extract the pathological M Stage (evidence of metastasis)from normal cased text of the pathology report        
    '''
    return_dictionary={global_strings.NAME:"PathStageM",global_strings.VALUE:None,global_strings.CONFIDENCE:0.0,global_strings.KEY:specimen,
                       global_strings.VERSION:__version__,global_strings.STARTSTOPS:[],global_strings.TABLE:global_strings.STAGE_GRADE_TABLE}
    ##currently this will only return one (the last) M value found   (but all char offsets)
    for section in dictionary:
        section_specimen=section[3]
        line_onset=section[2]
        header=section[1]
        if re.search(module_sections[1],header):
            for index,results in sorted(dictionary[section].items(),key=lambda x: int(x[0])):              
                m_stage=re.match('.*(pM[012x]).*',results,re.DOTALL)
                if m_stage:
                    return_dictionary[global_strings.VALUE]=m_stage.group(1)   
                    return_dictionary[global_strings.STARTSTOPS].append({global_strings.START:m_stage.start(),global_strings.STOP:m_stage.end()})
    return ([return_dictionary],list)  
  
