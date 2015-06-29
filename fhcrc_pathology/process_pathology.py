#
# Copyright (c) 2014-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import sys,path_parser,final_logic,re
import os,global_strings
path2= os.path.dirname(os.path.realpath(__file__))+'/'
'''author@esilgard'''
__version__='process_pathology1.0'

#################################################################################################################################################
def return_exec_code(x):
    '''
        helper method to retrieve the returned field value from each module
    '''
    return x

       
### MAIN CLASS ###
def main(arguments,path):
    '''
    current minimum required flags for the pathology parsing in the "arguments" dictionary are:
        -f input pathology file
        -g disease group
        returns a LIST of errors (could be warnings, or fatal exceptions
        and a DICTIONARY of field/value outputs 
    '''
   
    ## parse input documents from OBX table structure into a dictionary
    try:
        pathology_dictionary,return_type=path_parser.parse(arguments.get('-f'))        
        if return_type!=dict: return ({},pathology_dictionary,Exception)
    except:
        return({},[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not parse input pathology file '+arguments.get('-f')+' --- program aborted'+str(sys.exc_info())}],Exception)
    disease_group=arguments.get('-g')


    ## move import statments into these loops
    ## general pathology dictionary of tables, their fields/modules, and their relevant (pipe delimited) section headers ##
    try:
        path_data_dictionary={}
        for lines in open(path2+'/data_dictionary.txt','r').readlines():
            l=lines.strip().split('\t')        
            if l[0][0]!='#':
                path_data_dictionary[l[1]]=path_data_dictionary.get(l[1],{})
                path_data_dictionary[l[1]][l[0]]=l[2]              
                try:
                    if '-a' not in arguments or arguments.get('-a')!='n':                       
                        exec('import '+l[0])                        
                except:        
                    return (return_type,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology.import_modules() -  \
                        no reports completed.  Return error string: '+str(sys.exc_info())}],Exception)       
    except: return ({},[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not access or parse pathology data dictionary at '\
            +path+'/data_dictionary.txt --- program aborted'+str(sys.exc_info())}],Exception)
   
    ## disease group specific dictionary of tables, their fields/modules, and their relevant (pipe delimited) section headers ##
    try:        
        for lines in open(path2+'/'+disease_group+'/data_dictionary.txt','r').readlines():
            l=lines.strip().split('\t')
            if l[0][0]!='#':
                path_data_dictionary[l[1]]=path_data_dictionary.get(l[1],{})
                path_data_dictionary[l[1]][l[0]]=l[2]
                try:
                    if '-a' not in arguments or arguments.get('-a')!='n':
                        exec ('from '+disease_group+' import '+l[0])
                except:        
                    return (return_type,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology.import_modules() -  \
                            no reports completed.  Return error string: '+str(sys.exc_info())}],Exception)
    except: return ({},[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not access or parse disease specific pathology data dictionary at '\
            +path+'/'+disease_group+'/data_dictionary.txt '+str(sys.exc_info())}],Exception)
  
    field_value_output=[]
    error_output=[]    
    ## create a list of output field dictionaries ##
    for mrn in pathology_dictionary:        
        for accession,pathology_report_dictionary in pathology_dictionary[mrn].items():           
            field_value_dictionary={}
            field_value_dictionary[global_strings.REPORT]=accession
            field_value_dictionary[global_strings.MRN]=mrn
            field_value_dictionary[global_strings.SPECS_LOWER]=[]
            return_errors=[]
            ## output cannonical version of text file ##
            try:                
                with open(arguments.get('-f')[:arguments.get('-f').find('.nlp')]+'/'+accession+'.txt','wb') as out:
                          out.write(pathology_dictionary[mrn][accession][(-1,'FullText',0,None)])
            except:
                return (field_value_output,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology attempting to write text to file at'+ \
                        arguments.get('-f')[:arguments.get('-f').find('.nlp')] +'/'+accession+'.txt - unknown number of reports completed. '+str(sys.exc_info())}],list)
            
            if '-a' in arguments and arguments.get('-a')=='n':
                ## the no algorithm flag is on and you've already output text files and metadata, you're done ##
                pass
            else:               
                error_list=[]
                
                ## loop through specimens in dictionary (from the pathology report)  to only send relevant sections to field modules
                specimen_dictionary= pathology_report_dictionary[(0,global_strings.SPECIMEN_SOURCE,0,None)][0]
                try:
                    ## loop through specimens to create specimen level records in the destination tables
                    for specimen,description in specimen_dictionary.items():
                        specimen_pathology_report_dictionary={}
                        table_list=[]                        
                        ## loop through field modules and look for specimen specific values (these are called out explicitly in the "data_dictionary.txt" files
                        for table in path_data_dictionary:                           
                            table_d={global_strings.TABLE:table}
                            field_list=[]                            
                            for field,sections in path_data_dictionary[table].items():
                                if sections=="FullText":
                                    specimen_pathology_report_dictionary[(-1,'FullText',0,None)]=pathology_report_dictionary[(-1,'FullText',0,None)]
                                ## list comprehension is blocked here because of the function scope interfering with the exec commands
                                else:
                                    for section,text in pathology_report_dictionary.items():                                    
                                        ## make a shorter, specimen relevant only, pathology dictionary to pass to the modules
                                        if  section[3] and specimen in section[3]:
                                            specimen_pathology_report_dictionary[section]=text                                        
                                try:                                        
                                    exec("field_value_list,return_type=return_exec_code("+field+".get(disease_group,specimen_pathology_report_dictionary,specimen,sections))")                            
                                    if return_type==list:
                                        for each in field_value_list:
                                            if each[global_strings.VALUE]:
                                                field_list.append(each)                                           
                                    else:                                        
                                        error_list+=field_value_list
                                except:                      
                                    return ({},[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR could not process '+field+\
                                            ' module in general pathology directory --- program aborted. '+str(sys.exc_info())}],Exception)
                            
                            table_d[global_strings.FIELDS]=field_list  
                            table_list.append(table_d)                            
                        field_value_dictionary[global_strings.SPECS_LOWER].append({global_strings.SPEC_LOWER:specimen,global_strings.TABLES:table_list}) 
                            
                except:                
                    return (field_value_output,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology.get_fields() - \
                    unknown number of reports completed. '+str(sys.exc_info())}],list)
                # space HERE for FullText back off model search of any/all elements? (no PathHistology presumably)
            #space HERE for FINAL LOGIC MODULE CALL - BASED ON ALL FIELDS IN ACCESSION/REPORT ...dealing with PathFindings?               
            field_value_output.append(field_value_dictionary)
            return_errors+=error_list
    return (field_value_output,return_errors,list)
