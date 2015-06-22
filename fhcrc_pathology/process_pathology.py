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
   
    ## parse input documenbts from OBX table structure into a dictionary
    try:
        pathology_dictionary,return_type=path_parser.parse(arguments.get('-f'))        
        if return_type!=dict: return ({},pathology_dictionary,Exception)
    except:
        return({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not parse input pathology file '+arguments.get('-f')+' --- program aborted'},Exception)
    disease_group=arguments.get('-g')
    
    ## general pathology data dictionary ##
    try:    path_data_dictionary=dict((y.split('\t')[0],y.split('\t')[1].strip()) for y in open(path2+'/data_dictionary.txt','r').readlines())
    except: return ({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not access or parse pathology data dictionary at '+path+'/data_dictionary.txt --- program aborted'},Exception)
   
    ## disease group data dictionary ##
    try:    disease_group_data_dictionary=dict((y.split('\t')[0],y.split('\t')[1].strip()) for y in open(path2+'/'+disease_group+'/data_dictionary.txt','r').readlines())
    except: return ({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR: could not access or parse disease specific pathology data dictionary at '\
                        +path+'/'+disease_group+'/data_dictionary.txt '+str(sys.exc_info())},Exception)

    ## import appropriate pathology modules ##
    try:
        if '-a' not in arguments or arguments.get('-a')!='n':
            #import_modules(disease_group,disease_group_data_dictionary,path_data_dictionary)
            '''
            import modules (either general pathology modules, or disease specific depending on parameters)
            disease specific algorithms by the same name will override general ones
            '''    
            for field in path_data_dictionary.keys():                
                ## import general pathology modules
                exec('import '+field)
                
            for dz_field in disease_group_data_dictionary.keys():
                ## import disease specific pathology modules               
                exec ('from '+disease_group+' import '+dz_field)
    except:        
        return (return_type,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology.import_modules() -  \
                        no reports completed.  Return error string: '+str(sys.exc_info())}],list)
    
    field_value_output=[]
    error_output=[]    
    ## create a list of output field dictionaries ##
    for mrn in pathology_dictionary:        
        for accession in pathology_dictionary[mrn]:           
            field_value_dictionary={}
            field_value_dictionary[global_strings.REPORT]=accession
            field_value_dictionary[global_strings.MRN]=mrn
            return_errors=[]
            ## output cannonical version of text file 
            try:                
                with open(arguments.get('-f')[:arguments.get('-f').find('.nlp')]+'/'+accession+'.txt','wb') as out:
                          out.write(pathology_dictionary[mrn][accession][(-1,'FullText',0,None)])
            except:
                return (field_value_output,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology attempting to write text to file at'+ \
                        arguments.get('-f')[:arguments.get('-f').find('.nlp')] +'/'+accession+'.txt - unknown number of reports completed. '+str(sys.exc_info())}],list)
            #return_fields,return_errors,return_type=get_fields(disease_group,pathology_dictionary[mrn][accession],disease_group_data_dictionary,path_data_dictionary)

            ## if there are no Exceptions and the "no algorithm" isn't there, then run appropriate algorithms
            if return_type!=Exception:
                if arguments.get('-a')=='n':
                    pass
                else:
                    report_table_d={}
                    error_list=[]
                    data_elements=path_data_dictionary.keys()+disease_group_data_dictionary.keys()    
                    for field in data_elements:        
                        ## grab values for each of the fields in the disease specific data dictionary, back off to general module if there is no disease specific version ##
                        try:            
                            exec("field_value,return_type=return_exec_code("+field+".get(disease_group,pathology_dictionary[mrn][accession]))")                            
                            if return_type==list:
                                report_table_d=add_to_report_dictionary(field_value,report_table_d)                                
                            else:
                                error_list+=field_value
                        except:                      
                            return ({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR could not process '+field+\
                                        ' module in general pathology directory --- program aborted. '+str(sys.exc_info())},Exception)
                    
                    report_table_list=report_table_d.values()
                
                field_value_dictionary[global_strings.SPECIMEN]=final_logic.get(report_table_list)                  
                field_value_output.append(field_value_dictionary)
                return_errors+=error_list
            else:                
                return (field_value_output,[{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR in process_pathology.get_fields - unknown number of reports completed.\
                Return error string: '+return_errors[global_strings.ERR_STR]+';'+str(sys.exc_info())}],list)           

    return (field_value_output,return_errors,list)
