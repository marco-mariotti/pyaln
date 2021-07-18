__all__=['sequence_identity', 'weighted_sequence_identity']
import numpy as np


def sequence_identity(a, b, gaps='y'):
    """ 
    gaps:  
       a:  gaps are considered like any other character, may be scored positively
           (in all other cases, match between any character and a gap 
            is never an identity, even gap vs gap)
       y:  count gaps (present in one but not the other sequences)
       n:  don't count gaps (present in either)
       t:  don't count terminal gaps in either sequence

    """
    
    if len(a)!=len(b):
        raise IndexError('sequence_identity ERROR sequences do not have the same length')

    if   gaps=='y':
        pos_to_remove=[i for i in range(len(a)) if a[i]=='-' and b[i]=='-'  ]
    elif gaps=='n':
        pos_to_remove=[i for i in range(len(a)) if a[i]=='-' or  b[i]=='-'  ]        
    elif gaps=='t':     
        pos_to_remove=[i for i in range(len(a)) if a[i]=='-' and b[i]=='-'  ]
        for s in [a,b]:
            for i,c in enumerate(s):
                if c=='-':
                    pos_to_remove.append(i)
                else:
                    break
            for i, c in reversed(list(enumerate(s))):
                if c=='-':
                    pos_to_remove.append(i)
                else:
                    break
    elif gaps=='a':
        count_identical=sum([int(ca == b[i])    for i,ca in enumerate(a)])        
        return count_identical/len(a)        
    else:
        raise Exception('sequence_identity ERROR gaps argument must be one of {a, y, n, t}')

    exclude_pos=set(pos_to_remove)
    count_identical=sum([int(ca == b[i] and ca!='-' )    for i,ca in enumerate(a) if not i in exclude_pos])  
    return count_identical/( len(a) - len(exclude_pos) )

def weighted_sequence_identity(a, b, weights, gaps='y'):
    """ 
    weights:  list or anyhow iterable of weights as long as the alignment
              each position is weighted by this value
              the final score is divided by their sum

    gaps:  
       a:  gaps are considered like any other character, may be scored positively
           (in all other cases, match between any character and a gap 
            is never an identity, even gap vs gap)
       y:  count gaps (present in one but not the other sequences)
       n:  don't count gaps (present in either)
       t:  don't count terminal gaps in either sequence

    """
    
    if len(a)!=len(b):
        raise IndexError('sequence_identity ERROR sequences do not have the same length')

    if len(weights)!=len(a):
        raise IndexError('sequence_identity ERROR weights must be the same length as sequeneces')

    if   gaps=='y':
        pos_to_remove=[i for i in range(len(a)) if a[i]=='-' and b[i]=='-'  ]
    elif gaps=='n':
        pos_to_remove=[i for i in range(len(a)) if a[i]=='-' or  b[i]=='-'  ]        
    elif gaps=='t':     
        pos_to_remove=[i for i in range(len(a)) if a[i]=='-' and b[i]=='-'  ]
        for s in [a,b]:
            for i,c in enumerate(s):
                if c=='-':
                    pos_to_remove.append(i)
                else:
                    break
            for i, c in reversed(list(enumerate(s))):
                if c=='-':
                    pos_to_remove.append(i)
                else:
                    break
    elif gaps=='a':
        total_weight=   sum( weights )    
        count_identical=sum([int(ca == b[i])*weights[i]    for i,ca in enumerate(a)])        
        return count_identical/total_weight
    else:
        raise Exception('sequence_identity ERROR gaps argument must be one of {a, y, n, t}')

    exclude_pos=set(pos_to_remove)
    actual_weights=[w for i,w in enumerate(weights)   if not i in exclude_pos]
    total_weight=   sum( actual_weights )
    
    count_identical=sum([int(ca == b[i] and ca!='-' )*weights[i]    for i,ca in enumerate(a) if not i in exclude_pos])  
    return count_identical/( total_weight )








        
    
                
    
            
            
        
        
        
        
        
        
