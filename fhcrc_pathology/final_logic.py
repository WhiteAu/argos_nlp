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

import sys,path_parser
import os
path2= os.path.dirname(os.path.realpath(__file__))+'/'
'''author@esilgard'''
'''written October 2014'''
__version__='final_pathology_logic1.0'

def get(table_list):
    '''
    use all the extracted pathology elements from a given report to apply final logic
    add/delete values for the final output
    table_list = list of table dictionary where
                table dictionary: tableName and a dictionary
                    
    '''
    return_list=[]
    for table in table_list:
        
        remove_list=[]
        for record,dictionary in table.items():
            record_histology_found=False
            if dictionary['tableName']=='PathologyFinding':
                for each in dictionary['fields']:                    
                    if each['name']=='PathFindHistology':
                        record_histology_found=True
            
                if record_histology_found==False:
                    remove_list.append(record)
        ## take away any PathologyFinding values if there was no cancer in the given specimen ##
        #print 'remove list',remove_list
        table=dict((k,v) for k,v in table.items() if k not in remove_list)
        #print len(table),'appended to return_list'
        return_list.append(table)

    return return_list      