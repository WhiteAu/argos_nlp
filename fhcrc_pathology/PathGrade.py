#
# Copyright (c) 2014-2015 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

'''author@esilgard'''
__version__='PathGrade1.0'

import re
import global_strings
## mapping of grade numbers to words ##
## ***** NEED TO EXTEND THIS TO A GENERALIZABLE DZ MODEL ***** ##
grades={'high':'high','low':'low','intermediate':'intermediate',
        ' 3 ':'high',' 1 ':'low',' 2 ':'intermediate',' iii ':'high',
        ' i ':'low',' ii ':'intermediate'}
histos=['carcinoma','cancer','sclc']

descriptions=['moderate','poor','well']
## helper method to extract grade from text for each specimen ##
def get_grade(text,line_onset,specimen):
   
    grade_set=set([])
    start_stops_list=[] 

    ## number grades paired with the word "grade" ##
    for g,n in grades.items():
        m=re.match('.*('+g+'.{,15}grade).*',text)
        if m:
            grade_set.add(n)
            start_stops_list.append({global_strings.START:m.start(1)+line_onset,global_strings.STOP:m.end(1)+line_onset})
        m=re.match('.*(grade.{,15'+g+').*',text)
        if m:
            grade_set.add(n)
            start_stops_list.append({global_strings.START:m.start(1)+line_onset,global_strings.STOP:m.end(1)+line_onset})
          
    ## descriptions of cell differentiation ##
    for each_desc in descriptions:
        m= re.match('.*('+each_desc+'.{1,15}differentiated).*',text)
        if m:
            grade_set.add(m.group(1))        
            start_stops_list.append({global_strings.START:m.start(1)+line_onset,global_strings.STOP:m.end(1)+line_onset})
                
    ## specific grading system (FNCLCC)
    m=re.match('.*(([123])[/of ]{1,6}3.{,20}fn[c]?l[c]?c).*',text)    
    if m:
        grade_set.add(grades[' '+m.group(1)+' '])
        start_stops_list.append({global_strings.START:m.start(1)+line_onset,global_strings.STOP:m.end(1)+line_onset})
    else: m=re.match('.*(fn[c]?l[c]?c .{,20}([123])[/of ]{1,6}3).*',text)
    if m:
        grade_set.add(grades[' '+m.group(1)+' '])
        start_stops_list.append({global_strings.START:m.start(1)+line_onset,global_strings.STOP:m.end(1)+line_onset})
    else: m=re.match('.*(fn[c]?l[c]?c .{,20}grade.{,5}([123])).*',text)
    if m:
        grade_set.add(grades[' '+m.group(1)+' '])
        start_stops_list.append({global_strings.START:m.start(1)+line_onset,global_strings.STOP:m.end(1)+line_onset})

    ## discard substrings ##
    grade_list=sorted(grade_set,key=lambda x: len(x))    
    for i in range (len(grade_list)-1):
        if grade_list[i] in grade_list[i+1]:
            grade_set.remove(grade_list[i])
    
    return grade_set,start_stops_list

## main method ##        
def get(disease_group,dictionary,specimen,module_sections):
   
    '''
    extract the grade from the lower cased text of the pathology report       
    '''
    return_dictionary_list=[]    
    whole_start_stops_list=[]
    whole_grade_set=set([])
    grade_set=set([])
    return_dictionary_list=[]
    if specimen:
        for section in dictionary:               
            section_specimen=section[3]
            line_onset=section[2]
            header=section[1]
            if re.search(module_sections, header):                               
                text= dictionary[section].values()[0]                    
                ## meant to weed out references to literature/papers - picking up publication info like this: 2001;30:1-14. ##
                ## these can contain confusing general statements about the cancer and/or patients in general ##
                if re.search('[\d]{4}[;,][ ]*[\d]{1,4}:[\d\-]{1,6}',text):pass               
                else:                        
                    text=text.lower()
                    text=re.sub('[,:;\\\/\-]',' ',text); text=re.sub('[.] ', '  ',text)      ## this should keep decimal places and throw out periods                        
                    specimen_grade,offsets=get_grade(text,line_onset,specimen)                        
                    if specimen_grade:
                        grade_set=grade_set.union(specimen_grade)                            
                        start_stops_list+=offsets  
                                
                if grade_set:
                    confidence=.85
                    if len(grade_set)>1: confidence=.7
                    return_dictionary_list.append({global_strings.NAME:"PathGrade",global_strings.KEY:specimen,global_strings.TABLE:global_strings.PATHOLOGY_TABLE,global_strings.VALUE:';'.join(grade_set),
                        global_strings.CONFIDENCE:("%.2f" % confidence), global_strings.VERSION:__version__,global_strings.STARTSTOPS:start_stops_list} )
                    return_dictionary_list.append({global_strings.NAME:"PathGrade",global_strings.KEY:specimen,global_strings.TABLE:global_strings.STAGE_GRADE_TABLE,global_strings.VALUE:';'.join(grade_set),
                        global_strings.CONFIDENCE:("%.2f" % confidence), global_strings.VERSION:__version__,global_strings.STARTSTOPS:start_stops_list} ) 
    else:
        text=dictionary.lower()
        text=re.sub('[,:;\\\/\-]',' ',text); text=re.sub('[.] ', '  ',text)      ## this should keep decimal places and throw out periods 
        specimen_grade,offsets=get_grade(text,0,'')
        ## clean this up - lot of repeated code
        if specimen_grade:
            confidence=.65
            if len(specimen_grade)>1: confidence=.5
            return_dictionary_list.append({global_strings.NAME:"PathGrade",global_strings.VALUE:';'.join(specimen_grade),
                global_strings.CONFIDENCE:("%.2f" % confidence), global_strings.VERSION:__version__,global_strings.STARTSTOPS:start_stops_list} )
            
    '''
    BRING THIS LOGIC TO THE PROCESS PATHOOLOGY LOOP
    ## if there were no specimens, or no specimen headers in the text - look at the text overall ##
    else:        
        for section in dictionary:               
            section_specimen=section[3]
            line_onset=section[2]
            header=section[1]
            if section_specimen is not None and specimen in section_specimen and ('COMMENT' in header \
                or 'FINAL' in header or 'IMPRESSION' in header or 'SUMMARY' in header):                               
                text= dictionary[section].values()[0]                    
                ## meant to weed out references to literature/papers - picking up publication info like this: 2001;30:1-14. ##
                ## these can contain confusing general statements about the cancer and/or patients in general ##
                if re.search('[\d]{4}[;,][ ]*[\d]{1,4}:[\d\-]{1,6}',text):pass               
                else:                        
                    text=text.lower()
                    text=re.sub('[,:;\\\/\-]',' ',text); text=re.sub('[.] ', '  ',text)      ## this should keep decimal places and throw out periods                        
                    specimen_grade,offsets=get_grade(text,line_onset,specimen)                        
                    if specimen_grade:
                        whole_grade_set=whole_grade_set.union(specimen_grade)                            
                        whole_start_stops_list+=offsets  
                            
        if whole_grade_set:
            return_dictionary_list.append({global_strings.NAME:"PathGrade",global_strings.TABLE:global_strings.PATHOLOGY_TABLE,global_strings.VALUE:';'.join(whole_grade_set),
                global_strings.CONFIDENCE:0.65,global_strings.VERSION:__version__,global_strings.STARTSTOPS:whole_start_stops_list})
            return_dictionary_list.append({global_strings.NAME:"PathGrade",global_strings.TABLE:global_strings.STAGE_GRADE_TABLE,global_strings.VALUE:';'.join(whole_grade_set),
                global_strings.CONFIDENCE:0.65,global_strings.VERSION:__version__,global_strings.STARTSTOPS:whole_start_stops_list}) 
    '''
    return (return_dictionary_list,list)
