__all__=['sequence_identity']



def sequence_identity(a, b, gaps='y'):
    """ 
    gaps:  match between any character and a gap is never an identity, even gap vs gap
    
       a:  count all gaps, even those shared 
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
        pos_to_remove=[]
    else:
        raise Exception('sequence_identity ERROR gaps argument must be one of {a, y, n, t}')

    exclude_pos=set(pos_to_remove)
    print( exclude_pos)
    count_identical=sum([int(ca == b[i] and ca!='-')    for i,ca in enumerate(a) if not i in exclude_pos])  
    return count_identical/( len(a) - len(exclude_pos) )
        
        
    
                
    
            
            
        
        
        
        
        
        
