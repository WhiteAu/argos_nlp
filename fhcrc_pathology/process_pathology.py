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

def add_to_report_dictionary(field_value,report_table_d):
    '''
    get returned data elements into output dictionary structure
    specimens: tables : fields
    '''
    for each_field in field_value:
        specimen=each_field.get(global_strings.KEY)
        table=each_field.get(global_strings.TABLE)
        report_table_d[specimen]=report_table_d.get(specimen,{})
        report_table_d[specimen][global_strings.SPECIMEN]=specimen
        report_table_d[specimen][table]=report_table_d.get(table,{})
        report_table_d[specimen][table][global_strings.TABLE]=table
        report_table_d[specimen][table][global_strings.FIELDS]= report_table_d[specimen][table].get(global_strings.FIELDS,[])
        if each_field["value"]: report_table_d[specimen][table][global_strings.FIELDS].append(each_field)
    return report_table_d
            

        
### MAIN CLASS ###
def main(arguments,path):
    '''
    current minimum required flags for the pathology parsing in the "arguments" dictionary are:
        -f input pathology file
        -g disease group
    '''
   
    ## parse input documents from OBX table structure into a dictionary
    try:
        pathology_dictionary,return_type=path_parser.parse(arguments.get('-f'))        
        if return_type!=dict: return ({},pathology_dictionary,Exception)
    except:
        return({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not parse input pathology file '+arguments.get('-f')+' --- program aborted'},Exception)
    disease_group=arguments.get('-g')
    
    ## general pathology dictionary of modules, whether they are a specimen or report level element and their relevant (pipe delimited) section headers ##
    try:    path_data_dictionary=dict((y.split('\t')[0],(y.split('\t')[1],y.split('\t')[2].strip())) for y in open(path2+'/data_dictionary.txt','r').readlines() if y[0]!='#')
    except: return ({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not access or parse pathology data dictionary at '+path+'/data_dictionary.txt --- program aborted'},Exception)
   
    ## disease group specific dictionary of modules, whether they are a specimen or report level element and their relevant (pipe delimited) section headers ##
    try:    disease_group_data_dictionary=dict((y.split('\t')[0],(y.split('\t')[1],y.split('\t')[2].strip())) for y in open(path2+'/'+disease_group+'/data_dictionary.txt','r').readlines() if y[0]!='#')
    except: return ({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not access or parse disease specific pathology data dictionary at '\
                        +path+'/'+disease_group+'/data_dictionary.txt '+str(sys.exc_info())},Exception)

    
    ## import appropriate pathology modules (as long as the "no algorithm" flag isn't on ##
    try:
        if '-a' not in arguments or arguments.get('-a')!='n':            
            '''
            import modules (either general pathology modules, or disease specific depending on parameters)
            disease specific algorithms by the same name will override general ones
            '''    
            for field in path_data_dictionary.keys():                
                ## import general pathology modules
                exec('import '+field)                
            for dz_field in disease_group_data_dictionary.keys():
                ## import disease specific pathology modules  - these will overwrite any general modules of the same name                
                exec ('from '+disease_group+' import '+dz_field)
        else:
            ## the 'no algorithm flag' is set - no importing of modules that aren't needed
            pass
        
    except:        
        return (return_type,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology.import_modules() -  \
                        no reports completed.  Return error string: '+str(sys.exc_info())}],list)

    path_data_dictionary.update(disease_group_data_dictionary)                
    ## seperate out report level data elements/modules from specimen level 
    report_module_dictionary=dict((a,b) for a,b in path_data_dictionary.items() if b[0] =='report_level')
    specimen_module_dictionary=dict((a,b) for a,b in path_data_dictionary.items() if b[0] =='specimen_level')
                
    field_value_output=[]
    error_output=[]    
    ## create a list of output field dictionaries ##
    for mrn in pathology_dictionary:        
        for accession,pathology_report_dictionary in pathology_dictionary[mrn].items():
            print mrn,accession
            field_value_dictionary={}
            field_value_dictionary[global_strings.REPORT]=accession
            field_value_dictionary[global_strings.MRN]=mrn
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
                report_table_d={}
                error_list=[]
                
                ## loop through specimens in dictionary (from the pathology report)  to only send relevant sections to field modules
                specimen_dictionary= pathology_report_dictionary[(0,global_strings.SPECIMEN_SOURCE,0,None)][0]
                try:
                    ## loop through field modules and look for specimen specific values (these are called out explicitly in the "data_dictionary.txt" files                       
                    for field in specimen_module_dictionary:
                        print 'running specimen level ',field
                        value_found=False
                        ## grab values for each of the fields in the disease specific data dictionary
                        for specimen,description in specimen_dictionary.items():
                            print specimen
                            specimen_pathology_report_dictionary={}
                            ## list comprehension is blocked here because of the function scope interfering with the exec commands
                            for section,text in pathology_report_dictionary.items():                            
                                if  section[3] and specimen in section[3]:
                                    specimen_pathology_report_dictionary[section]=text
                                try:            
                                    exec("field_value,return_type=return_exec_code("+field+".get(disease_group,specimen_pathology_report_dictionary,specimen,specimen_module_dictionary[field]))")                            
                                    if return_type==list:                                      
                                    
                                        report_table_d=add_to_report_dictionary(field_value,report_table_d)                                
                                    else:                                        
                                        error_list.append(field_value)
                                except:                      
                                    return ({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR could not process '+field+\
                                            ' module in general pathology directory --- program aborted. '+str(sys.exc_info())},Exception)
                                
                            ## add to the report level dictionary to "back off" to looking for some sort of general value for the report as a whole
                            ## this covers more common formatting edgecases where no individual specimens are listed in the report
                            ## if value_found==False:
                            ##    report_module_dictionary[field]=specimen_module_dictionary[field]                   
                            ##    print 'no spec specific value found for',field
                            for field in report_module_dictionary:
                                print 'running report level ',field
                                exec("field_value,return_type=return_exec_code("+field+".get(disease_group,pathology_report_dictionary,specimen,report_module_dictionary[field]))")
                                if return_type==list:
                                    if field_value:
                                        report_table_d=add_to_report_dictionary(field_value,report_table_d)
                                else:
                                    error_list.append(field_value)
                                    
                        report_table_list=report_table_d.values()              

                except:                
                    return (field_value_output,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology.get_fields() - \
                    unknown number of reports completed. '+str(sys.exc_info())}],list)           
            #space for FINAL LOGIC MODULE CALL - BASED ON ALL FIELDS IN ACCESSION/REPORT                  
            field_value_output.append(report_table_d)
            return_errors+=error_list
    return (field_value_output,return_errors,list)
