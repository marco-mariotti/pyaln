import os, io
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from functools import lru_cache

__all__=['Alignment']

class AlignmentError(Exception):
    """ """

class Alignment: 
    """This is the mighty alignment class

    Truly, it's very very alignypoweryfully
    - **parameters**, **types**, **return** and **return types**::

          :param arg1: description
          :param arg2: description
          :type arg1: type description
          :type arg1: type description
          :return: return description
          :rtype: the return type description

    - and to provide sections such as **Example** using the double commas syntax::

          :Example:

          followed by a blank line !

      which appears as follow:

      :Example:

      followed by a blank line

    - Finally special sections such as **See Also**, **Warnings**, **Notes**
      use the sphinx syntax (*paragraph directives*)::

          .. seealso:: blabla
          .. warnings also:: blabla
          .. note:: blabla
          .. todo:: blabla

    .. note::
        There are many other Info fields but they may be redundant:
            * param, parameter, arg, argument, key, keyword: Description of a
              parameter.
            * type: Type of a parameter.
            * raises, raise, except, exception: That (and when) a specific
              exception is raised.
            * var, ivar, cvar: Description of a variable.
            * returns, return: Description of the return value.
            * rtype: Return type.

    """
    extension2format={'fa':'fasta', 'fasta':'fasta',
                      'aln':'clustal', 'clustal':'clustal', 'clw':'clustal',
                      'stk':'stockholm', 'stockholm':'stockholm',
                      'phy':'phylip', 'phylip':'phylip'
    }
    max_names_repr=15
    max_sequence_repr=60
    max_cache_size=30  # read only at startup, immutable
        
    def __init__(self, file_or_iter=None, fileformat=None,  tuples=None):
        """And this is how you start.

        Seriously, you start then it gets easier
        and easier

        """
        self._seqs={}
        self._desc={}
        self._ord =[]

        if not file_or_iter is None:            
            if isinstance(file_or_iter, str) or isinstance(file_or_iter, io.IOBase):
                
                # getting file extension for automatic fileformat detection
                if fileformat is None:
                    if isinstance(file_or_iter, io.IOBase):
                        if file_or_iter.name:
                            ext=os.path.splitext(file_or_iter.name)[-1].lstrip('.')
                    else:            
                        ext=os.path.splitext(file_or_iter)[-1].lstrip('.')
                    if ext in self.extension2format:
                        fileformat=self.extension2format[ext]
                        
                if fileformat is None:
                    raise AlignmentError(f"fileformat not specified and not automatically recognized!")
                    
                with (open(file_or_iter) if isinstance(file_or_iter, str)
                      else file_or_iter)      as fh:
                    
                    if fileformat=='fasta':       ## using a faster option for the most common format
                        for title, sequence in SeqIO.FastaIO.SimpleFastaParser(fh):
                            self.add_seq(title, sequence)
                    else:
                        for r in AlignIO.read(fh, fileformat):
                            self.add_seq(r.description, str(r.seq))
                            
            elif isinstance(file_or_iter, Alignment):
                for name, sequence in file_or_iter:
                    self.add_seq(name, sequence, desc=file_or_iter.get_desc(name))
                            
            else:  #an iterable was provided            
                for title, sequence in file_or_iter:
                    self.add_seq(title, sequence)
                    
                                    
    def __iter__(self):
        for name in self._ord:
            yield name, self._seqs[name]

    def __eq__(self, other):
        return self.sequences() == other.sequences() # fastest option

    def __hash__(self):
        return '\n'.join(self.sequences()).__hash__()
        
    def __repr__(self):
        alilen=self.ali_length()
        nseqs=self.n_seqs()
        if not nseqs:
            return '# Empty alignment'
        out=f'# Alignment of {nseqs} sequences and {alilen} positions\n'        
        for name in self._ord[:min(self.max_names_repr-1, nseqs-1)  ]:
            out+=self._oneliner(name)+'\n'
        if nseqs>self.max_names_repr:
            out+='<...>\n'
        if nseqs:
            out+=self._oneliner(self._ord[-1])+'\n'
        return(out)

    def _oneliner(self, name):
        if self.ali_length()<=self.max_sequence_repr:
            seqshown=self.get_seq(name)
        else:                
            seqshown=self.get_seq(name)[:(self.max_sequence_repr-1)] + ' ... ' + self.get_seq(name)[-1]
        return(seqshown+ ' '+name)

    def __getitem__(self, rows_and_cols):
        """ indexed with tuple rows, cols
        rows can be a slice, or a list of names in the alignment"""
        if not type(rows_and_cols) is tuple or len(rows_and_cols)!=2:
            raise AlignmentError('__getitem__ ERROR you must specify an index for rows (names) and one for columns (sequence positions), e.g. ali[5:9,10:20] ali[:,:-5]')
        rows, cols=rows_and_cols
        ali=Alignment()
        names=self._ord[rows] if not type(rows) is list else rows
        for name in names:
            ali.add_seq(name,
                        self.get_seq(name) [cols],
                        desc=self.get_desc(name))
        return(ali)

    def add_seq(self, title, sequence, desc=None):
        """ """
        if desc is None:
            name=title.split(sep=None, maxsplit=1)[0]
            desc=title[len(name)+1:].strip()
        else:
            name, desc=title, desc
        if self.has_name(name):
            raise AlignmentError(f'add_seq ERROR alignment already has {name}. Did you mean to use set_seq() instead?')
        self._seqs [name] = sequence
        self._desc [name] = desc
        self._ord.append(name)

        
    def ali_length(self):
        return len(self.get_seq(self._ord[0])) if self._ord else 0

    def n_seqs(self):
        return len(self._ord)
        
    def set_seq(self, name, sequence):
        if not self.has_name(name):
            raise AlignmentError(f'set_seq ERROR alignment does not have {name}. Did you mean to use add_seq() instead?')
        self._seqs[name]=sequence

    def set_desc(self, name, desc):
        if not self.has_name(name):
            raise AlignmentError(f'set_desc ERROR alignment does not have {name}. Did you mean to use add_seq() instead?')
        self._desc[name]=desc            

    def get_desc(self, name):
        if not name in self._desc:
            raise AlignmentError(f'get_desc ERROR alignment does not have {name}')
        return self._desc[name]

    def get_seq(self, name):
        if not name in self._seqs:
            raise AlignmentError(f'get_seq ERROR alignment does not have {name}')
        return self._seqs[name]

    def names(self):
        return [name for name in self._ord]

    def sequences(self):
        return [self.get_seq(name) for name in self._ord]    
    
    def has_name(self, name):
        return name in self._seqs
            
    def fasta(self, nchar=60):
        out=''
        for name, seq in self:
            desc=self.get_desc(name)
            out+='>' + name + (" "+desc if desc else "") + '\n'
            for i in range(0, len(seq), nchar):
                out+=seq[i:i + nchar]+'\n'
        return out

    def concatenate(self, other):
        if self.names() != other.names():
            raise AlignmentError(f'concatenate ERROR the two alignments must have the same sequence names!')
        a=Alignment()
        for name, seq in self:            
            a.add_seq(name, seq + other.get_seq(name), desc=self.get_desc(name))
        return a

    def copy(self):
        return Alignment(self)

    def remove_seq(self, *names):
        for name in names:
            if not name in self._seqs:
                raise AlignmentError(f'remove_seq ERROR alignment does not have {name}')
            del self._seqs[name]
            del self._desc[name]
        s=set(names)            
        for i in range(len( self._ord )-1, -1, -1):
            if self._ord[i] in s:
                self._ord.pop(i)
                
    def remove_by_index(self, *seqindices):
        for i in sorted(seqindices, reverse=True):
            name=self._ord.pop(i)
            del self._seqs[name]
            del self._desc[name]

    def remove_empty_seqs(self):        
        empty_seq_names=numpy.array(self.names())[  # next: boolean array, True when all gaps
            numpy.char.equal(
                numpy.array( [numpy.array(seq, dtype=numpy.str_) for name, seq in self], dtype=numpy.str_ ),
                '-'*self.ali_length())                ]
        self.remove_seq(*empty_seq_names)        
            
    # def pop(self, index):
    #     name=self._ord.pop(i)
    #     seq=self._seqs[name]
    #     desc=self._desc[name]
    #     return( name, seq, desc )
                                 
    def to_biopython(self):
        return MultipleSeqAlignment(
            [SeqRecord( Seq(seq), id=name, name=name,
                        description=self.get_desc(name) )
             for name, seq in self])

    @lru_cache(maxsize=max_cache_size)
    def to_numpy(self):
        return numpy.array( [numpy.array(list(seq), dtype=numpy.str_) for name, seq in self], dtype=numpy.str_ )

    def to_pandas(self, save=False):
        return pd.DataFrame( self.to_numpy() )

    @lru_cache(maxsize=max_cache_size)
    def conservation_by_column(self):
        return self.to_pandas().apply(pd.value_counts, 0, normalize=True).fillna(0.0)
    
    def write(self, fileformat='fasta'):
        if fileformat=='fasta':
            return self.fasta()
        else:
            return format(self.to_biopython(), fileformat)

    def write_to_file(self, file_or_buff, fileformat='fasta'):
        with (file_or_buff if not isinstance(file_or_buff, str)
              else open(file_or_buff, 'w')) as fh:
            fh.write( self.write(fileformat=fileformat) )

    def convert_sequences(self, seqfn):
        for name, seq in self:
            self.set_seq(name, seqfn(seq))
            
    def trim_gaps(self, pct=1.0, inplace=False):
        ### get positions:
        #  [all([s[pos]=='-'  for n,s in a])   for pos in range(a.ali_length())]   # 65.3 ms
        #  a.conservation_by_column(save=False).loc['-']==1.0                      # 2.6 s
        #  a.conservation_by_column(save=True).loc['-']==1.0                       # 276 µs        
        #  np.apply_along_axis( lambda x: np.char.equal(x, '-').all(), 0,   a.to_numpy() )                 # 74.9 ms
        #  np.apply_along_axis( lambda x: np.char.equal(x, '-').all(), 0,   a.to_numpy(save=True) )        # 30.9 ms
        #  np.vectorize(lambda x:x=='-')(  a.to_numpy(save=False)).all(axis=0)     # 87.7 ms
        #  np.vectorize(lambda x:x=='-')(  a.to_numpy(save=True)).all(axis=0)      # 40.5 ms

        # if not self._cbc is None:   #### If self.conservation_by_column was alread run? couldn't built this check
        #     cbc=self.conservation_by_column()
        #     empty_cols_selector=( cbc.loc['-']==1.0 if '-' in cbc.index
        #                           else pd.Series(False, index=cbc.columns) )
        # else:

        #      .mean(axis=0) below returns 1-D array, size=ali_length, with proportion of gaps                    
        cols_selector=self.gap_mask().mean(axis=0) >= pct
        
        seqs=np.apply_along_axis( ''.join, 1, self.to_numpy()[:,~cols_selector])        
        if inplace:
            for i, (name, _) in enumerate(self):
                self._seqs[name]=seqs[i]  # fast set
                
        else:
            a=Alignment()
            for i, (name, _) in enumerate(self):
                a.add_seq(name, seqs[i], desc=self.get_desc(name))
            return a    

    def consensus(self):
        return ''.join(self.conservation_by_column().apply(pd.Series.idxmax))

    def gap_mask(self):
        return self.to_numpy()=='-'

    def terminal_gap_mask(self):
        return np.apply_along_axis(
            # below: each bool array is as long as sequence;  function is applied to each row
                      #  bool array: is left term gap? # or # bool array: is right term gap?            
            lambda x: ( numpy.logical_and.accumulate(x)  |  numpy.logical_and.accumulate( x[::-1] )[::-1] ),
            1,  self.gap_mask())
        #tested: per row is faster than whole matrix
    
    def _temp_cons(self):

        p=self.to_pandas().rename_axis('seqindex').rename_axis('position', axis=1)
        """position 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
        seqindex
        0         M  -  -  -  -  -  -  -  -  -  -  -  -  -  -  R  C  Y  V  E  -  -  -  -  -  Y  L  P  E  M
        1         M  F  S  -  -  -  -  -  -  -  -  P  G  S  D  R  T  N  L  R  A  E  V  P  P  F  V  P  R  -
        2         M  A  S  E  G  R  R  E  P  A  A  -  -  E  G  -  I  K  L  S  A  D  V  K  P  F  V  P  K  F
        """

        pp=pd.DataFrame( p.unstack().rename('letter') ).reset_index(0)
        """          position letter
        seqindex
        0                0      M
        1                0      M
        2                0      M
        0                1      -
        1                1      F
        ...            ...    ...
        1               28      R
        2               28      K
        0               29      M
        1               29      -
        2               29      F
        [90 rows x 2 columns]"""
        
        q=pd.DataFrame( b.conservation_by_column().unstack() ).set_axis(['conservation'], axis=1).rename_axis( ['position', 'letter'] )
        """                 conservation
        position letter
        0        -                0.0
        A                0.0
        C                0.0
        D                0.0
        E                0.0
        ...                       ...
        29       R                0.0
        S                0.0
        T                0.0
        V                0.0
        Y                0.0        
        [540 rows x 1 columns] """
        
        res=pd.merge( pp, q, left_on=['position', 'letter'], right_index=True)
        """          position letter  conservation
        seqindex
        0                0      M      1.000000
        1                0      M      1.000000
        2                0      M      1.000000
        0                1      -      0.333333
        1                1      F      0.333333
        ...            ...    ...           ...
        1               28      R      0.333333
        2               28      K      0.333333
        0               29      M      0.333333
        1               29      -      0.333333
        2               29      F      0.333333        
        [90 rows x 3 columns] """

        ## average sequence identity per sequence
        res.groupby(res.index).conservation.mean()

