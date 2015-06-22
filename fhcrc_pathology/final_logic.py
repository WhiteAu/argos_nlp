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

import sys,path_parser,global_strings
import os
path2= os.path.dirname(os.path.realpath(__file__))+'/'
'''author@esilgard'''

__version__='final_pathology_logic1.0'

def get(return_list):
    '''
    use all the extracted pathology elements from a given report to apply final logic
    add/delete values for the final output
    return_list = list of specimen dictionaries where
                specimens map to specimen "recordKey"s (A, B, etc)
                and lists of tables where
                    tables map to table names ("Pathology", "PathTest", etc)
                    and dictionaries of field values/data elements               
                    
    '''
    ## TODO - write in "PathFinding" results if there are more than one histology or grade per specimen
    return return_list      
