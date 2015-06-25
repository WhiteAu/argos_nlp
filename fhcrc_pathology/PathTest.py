#
# Copyright (c) 2014-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

'''author@esilgard'''

__version__='PathTest1.0'

def get(disease_group,dictionary,specimen,module_sections):
    '''
    extract the pathology tests from normal cased text of the pathology report
    
    '''
    return_dictionary={global_strings.NAME:"PathTest",global_strings.VALUE:None,global_strings.CONFIDENCE:0.0,global_strings.VERSION:__version__,
                           global_strings.STARTSTOPS:[],global_strings.KEY:specimen,global_strings.TABLE:global_strings.STAGE_GRADE_TABLE}
    for section in dictionary:
        section_specimen=section[3]
        line_onset=section[2]
        header=section[1]            
        if re.search(module_sections[1],header):
            for index,results in sorted(dictionary[section].items(),key=lambda x: int(x[0])):
                test=re.match('.*...........ppppppppppp.*',results,re.DOTALL)        
                return_dictionary[global_strings.VALUE]=t_stage.group(1)       
                return_dictionary[global_strings.STARTSTOPS].append({global_strings.START:t_stage.start(),global_strings.STOP:t_stage.end()})
                ####do stuff

    return (return_dictionary_list,list) 
