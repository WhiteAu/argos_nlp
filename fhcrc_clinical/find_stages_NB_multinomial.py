#Emily Silgard - find prognostic stages using weighted vecotrs (NB) - Multinomial Event Model

import sys, math, global_strings, os
import numpy as np
path= os.path.dirname(os.path.realpath(__file__))+'/'
threshold=.20
classes =['True','False']
prognostic_stages=['I','II','III','IV']
class_labels=['I','II','III','IV','UNKNOWN']


'''normalize log10 values before they're brought back into probability space in order to solve underflow problem'''
def normalize(x):
    a = np.logaddexp.reduce(x)
    return np.exp(x - a)

   
'''iterate through the training_vectors dictionary (label:list of lists [f1,v1],[f2,v2]....)
and output the probabilities of class given the features in the documents'''
def classify_stage_match_vectors(vectors,disease_group):
   
    ## import language models from disease group folder
    try:
        prior_model_file=path+disease_group+'/MASTER_NB_ME_prior_model.txt'
        prob_of_word_model_file=path+disease_group+'/MASTER_NB_ME_prob_of_word_model.txt'
        smoothed_model_file=path+disease_group+'/MASTER_NB_ME_smoothed_model.txt'
    except:
        return ({},{global_strings.ERR_TYPE:'Exception',global_strings.ERR_STR:'FATAL ERROR could not open prognostic staging diseasee specific models --- program aborted. '},Exception)
        
    return_dictionary={}
    pt_dict={}
    pOfClass =  dict((a.split('\t')[0],float(a.split('\t')[1])) for a in open(prior_model_file,'r').readlines())
    smoothed = dict((b.split('\t')[0],float(b.split('\t')[1])) for b in open(smoothed_model_file,'r').readlines())
    wordGivenClass = {}
    for lines in open(prob_of_word_model_file,'r').readlines():
        line=lines.split('\t')
        class_label=line[1]
        wordGivenClass[class_label]=wordGivenClass.get(class_label,{})
        wordGivenClass[class_label][line[0]]=float(line[2])
    
    math_log=math.log10

    for lines in vectors:
        line=lines.split()
        
        instance_data=line[0].split('_')
        stage_match=instance_data[1]        
        #patient=instance_data[0]
        vector=[(line[a],int(line[a+1])) for a in range(2,len(line)-2,2)]
        cProbs=[]
            
        for c in classes:
            num=0.0            
            for each in vector:
                feature=each[0]
                value=each[1]
                if feature in wordGivenClass[c]:
                    try:
                        num+=math_log(wordGivenClass[c][feature])*value
                    except:
                        print 'ERR FINDING LOG OF WORD GIVEN CLASS PROB',wordGivenClass[c][feature]
                else:
                    num+=math_log(smoothed[c])*value                                                    
            cProbs.append(num+math_log(.25))
        ## bring back into probability space ##
        cProbs=normalize(cProbs)        
        cProbDict=dict(zip(classes,cProbs))        
                
        pt_dict[stage_match]= pt_dict.get(stage_match,0)+cProbDict['True']        
       
    return output_pt_stages(sorted(pt_dict.items()))   


def output_pt_stages(candidate_stages):
    
    total=sum([m[1] for m in candidate_stages])
    if total < 0.00000005:  #float value may not match 0.0
        return('UNKNOWN',0.0)
        
    else:
        candidate_probs=[p[1]/total for p in candidate_stages]
        final_probs=dict(zip([q[0] for q in candidate_stages],candidate_probs))                   
        
        final_rank=sorted(final_probs.items(), key=lambda x: x[1], reverse=True)
        initial=final_rank[0][0]
        candidate_list=[(x[0],x[1]) for x in final_rank if abs(final_rank[0][1]-x[1])<threshold and x[1]>0.0]
        candidate_list=sorted(candidate_list)
        return(candidate_list[0][0],candidate_list[0][1])
                    
