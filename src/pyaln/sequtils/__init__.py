__all__=['sequence_identity', 'weighted_sequence_identity', 'sequence_identity_matrix', 'gap_selector']
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
        raise IndexError('sequence_identity ERROR weights must be the same length as sequences')

    total_weight=   sum( weights )    
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
        count_identical=sum([int(ca == b[i])*weights[i]    for i,ca in enumerate(a)])        
        return count_identical/total_weight
    else:
        raise Exception('sequence_identity ERROR gaps argument must be one of {a, y, n, t}')

    exclude_pos=set(pos_to_remove)
    #actual_weights=[w for i,w in enumerate(weights)   if not i in exclude_pos]
    #total_weight=   sum( actual_weights ) #nope
    
    count_identical=sum([int(ca == b[i] and ca!='-' )*weights[i]    for i,ca in enumerate(a) if not i in exclude_pos])  
    return count_identical/( total_weight )


def gap_selector(npt, nps, gaps):
    if gaps=='a':
        selector=np.full(     (npt.shape[0],   *nps.shape) , True, dtype=bool)
    elif gaps=='n':
        ## which positions are taken into account:  those in which neither target and selfseq is gap
        selector=(npt!='-')[:, np.newaxis, :]  & (nps!='-')[np.newaxis, :, :]
    elif gaps=='y':
        selector=(npt!='-')[:, np.newaxis, :]  | (nps!='-')[np.newaxis, :, :]
    elif gaps=='t':
        # code from Alignment.terminal_gap_mask
        terminal_gaps_t=np.apply_along_axis(
            lambda x: ( np.logical_and.accumulate(x)  |  np.logical_and.accumulate( x[::-1] )[::-1] ),
            1,  npt=='-')
        terminal_gaps_s=np.apply_along_axis(
            lambda x: ( np.logical_and.accumulate(x)  |  np.logical_and.accumulate( x[::-1] )[::-1] ),
            1,  nps=='-')        
        selector=(  ~(terminal_gaps_t)[:, np.newaxis, :] &
                    ~(terminal_gaps_s)   [np.newaxis, :, :] &
                    ((npt!='-')[:, np.newaxis, :]  | (nps!='-')[np.newaxis, :, :] ) )
    return selector

def sequence_identity_matrix(npt, nps, selector=None, gaps=None, eq_matrix=None):
    if selector is None:
        if gaps is None:
            raise ValueError("sequence_identity_matrix ERROR you must provide selector or gaps arguments")
        selector=gap_selector(npt, nps, gaps)
    if eq_matrix is None:
        eq_matrix=np.char.equal( nps, npt[:,np.newaxis] )
        
    return ( (eq_matrix & selector).sum(axis=2)  /
             # below: matrix of length of alignments, i.e. the positions actually used for each comparison
             selector.sum(axis=2) )
    




        
    
                
    
            
            
        
        
        
        
        
        
