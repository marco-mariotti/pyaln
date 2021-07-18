import os, io
from functools import lru_cache
from typing import Union, TextIO
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO, Seq, SeqRecord, Align

MultipleSeqAlignment=Align.MultipleSeqAlignment
SeqRecord=SeqRecord.SeqRecord
Seq=Seq.Seq

__all__=['Alignment', 'pyaln_folder']

pyaln_folder = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath( __file__)   )))

class AlignmentError(Exception):
    """Exceptions raised by the Alignment class """

    
class Alignment: 
    """Represents a multiple sequence alignment.

    Alignment can contain sequences of any type (nucleotide, protein, or custom).
    Gaps **must** be encoded as `-`.

    Each entry is uniquely identified by a *name*, with an optional *description*.
    When reading alignment files, sequence *titles* are split into the name (the first word)
    and descriptions (the remainder of the title).

    An Alignment can be instanced with a filename (or file buffer), or from any iterable of `[title, sequence]`.
    A variety of file formats are supported, through Bio.AlignIO (see `a full list  <https://biopython.org/wiki/AlignIO>`_).    
    When a filename or buffer is provided but not the file format, Alignment tries to guess it from the extension.

    You can take portions of an Alignment (i.e. take some sequences and/or some columns) by **indexing** it.

    The format is `Alignment[rows_selector, column_selector]`, where: 

        - The `rows_selector` can be an integer (i.e., the vertical position of 
          the sequence in the alignment), or a slice thereof (e.g. `2:5`), or a list of sequence names.
        - The `column_selector` is a integer index (i.e. the horizontal position in the alignment),
          or a slice thereof, or a Numpy boolean array. See examples below. 


    Iterating over an Alignment will yield tuples like `(name, sequence)`. 
    To get the description of a sequence, use `Alignment.get_desc(name)`.

    Parameters
    ----------
    file_or_iter : str | TextIO | iterable
        Filename to sequence file to be loaded, or TextIO buffer already opened on it, 
        or iterable of [title, seq] objects.
    fileformat : str, optional
        When a filename or TextIO is provided, specifies the file format (e.g. fasta, clustal, stockholm ..)

    Examples
    --------

    >>> ali=Alignment(pyaln_folder+'/examples/fep15_protein.fa')

    Default representation (note, it does not contain descriptions): 

    >>> ali
    # Alignment of 6 sequences and 138 positions
    MWLTLVALLALCATGRTAENLSESTTDQDKLVIARGKLVAPSVVGUSIKKMPELYNFLM...L Fep15_danio_rerio
    MWAFLLLTLAFSATGMTEED-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPALYMFLM...L Fep15_S_salar
    MWIFLLLTLAFSATGMTEEN-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPELYTFLM...L Fep15_O_mykiss
    MWAFLVLTFAVAA-GASETV-DNHTAAEEKLLIARGKLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_rubripess
    MWALLVLTFAVTV-GASEEV-KNQTAAEEKLVIARGTLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_nigroviridis
    MWAFVLIAFSVGA---SDS------SNSTAEVIARGKLMAPSVVGUAIKKLPELNRFLM...L Fep15_O_latipes
    <BLANKLINE>

    Many file formats are supported:

    >>> ali=Alignment(pyaln_folder+'/examples/fep15_protein.stockholm', fileformat='stockholm')
    >>> ali
    # Alignment of 6 sequences and 138 positions
    MWLTLVALLALCATGRTAENLSESTTDQDKLVIARGKLVAPSVVGUSIKKMPELYNFLM...L Fep15_danio_rerio
    MWAFLLLTLAFSATGMTEED-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPALYMFLM...L Fep15_S_salar
    MWIFLLLTLAFSATGMTEEN-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPELYTFLM...L Fep15_O_mykiss
    MWAFLVLTFAVAA-GASETV-DNHTAAEEKLLIARGKLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_rubripess
    MWALLVLTFAVTV-GASEEV-KNQTAAEEKLVIARGTLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_nigroviridis
    MWAFVLIAFSVGA---SDS------SNSTAEVIARGKLMAPSVVGUAIKKLPELNRFLM...L Fep15_O_latipes
    <BLANKLINE>

    Initializing from iterable (in this case a list):

    >>> ali=Alignment([ ('seq1', 'ATTCG-'), ('seq2', '--TTGG'), ('seq3', 'ATTCG-')])    
    >>> ali
    # Alignment of 3 sequences and 6 positions
    ATTCG- seq1
    --TTGG seq2
    ATTCG- seq3
    <BLANKLINE>

    To visualize sequence descriptions, use the fasta format:

    >>> ali=Alignment([ ('seq1 this is a seq', 'ATTCG-'), ('seq2 another seq', '--TTGG'), ('seq3', 'ATTCG-')])    
    >>> print(ali.fasta())
    >seq1 this is a seq
    ATTCG-
    >seq2 another seq
    --TTGG
    >seq3
    ATTCG-
    <BLANKLINE>

    **Indexing an alignment**

    Get alignment of first two sequences only:

    >>> ali[:2,:]
    # Alignment of 2 sequences and 6 positions
    ATTCG- seq1
    --TTGG seq2
    <BLANKLINE>

    Trim off the first and last alignment columns:

    >>> ali[:,1:-1]
    # Alignment of 3 sequences and 4 positions
    TTCG seq1
    -TTG seq2
    TTCG seq3
    <BLANKLINE>

    Get subalignment of two sequences, by their name:

    >>> ali[ ['seq1', 'seq3'], : ]
    # Alignment of 2 sequences and 6 positions
    ATTCG- seq1
    ATTCG- seq3
    <BLANKLINE>

    Iterating over an alignment:

    >>> [(name, len(seq))   for name, seq in ali]
    [('seq1', 6), ('seq2', 6), ('seq3', 6)]
    """

    extension2format={'fa':'fasta', 'fasta':'fasta',
                      'aln':'clustal', 'clustal':'clustal', 'clw':'clustal',
                      'stk':'stockholm', 'stockholm':'stockholm',
                      'phy':'phylip', 'phylip':'phylip'
    }
    max_names_repr=15
    max_sequence_repr=60
    max_cache_size=30  # read only at startup, immutable
        
    def __init__(self, file_or_iter=None, fileformat=None):
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
            seqshown=self.get_seq(name)[:(self.max_sequence_repr-1)] + '...' + self.get_seq(name)[-1]
        return(seqshown+ ' '+name)

    def __getitem__(self, rows_and_cols):
        """ indexed with tuple rows, cols
        rows can be a slice, or a list of names in the alignment
        cols con by a slice, or a bool array as long as the alignment, or a list like [ [start, end] ... [start,end]]
        """
        if not type(rows_and_cols) is tuple or len(rows_and_cols)!=2:
            raise AlignmentError('__getitem__ ERROR you must specify an index for rows (names) and one for columns (sequence positions), e.g. ali[5:9,10:20] ali[:,:-5]')
        rows, cols=rows_and_cols

        if type(cols) is np.ndarray:
            ali=self[rows,:].to_numpy()[:,cols]
        else:
            ali=Alignment()
            if type(rows) is slice:                names=self._ord[rows]
            elif type(rows) is list:               names=rows
            elif type(rows) is int:                names=[self._ord[rows]]
            else:
                raise AlignmentError(f'ERROR wrong type for rows! {type(rows)}')

            if type(cols) is slice:                
                for name in names:            
                    ali.add_seq(name,
                                self.get_seq(name) [cols],
                                desc=self.get_desc(name))
            elif type(cols) is list:
                s=self.get_seq(name)
                for name in names:            
                    ali.add_seq(name,
                                ''.join( [s[start:end]  for start, end in cols]),
                                desc=self.get_desc(name))
            else:
                raise AlignmentError(f'ERROR wrong type for cols! {type(cols)}')
            
        return(ali)

    def add_seq(self, title, sequence, desc=None, index=None):
        """Add a sequence to the alignment.

        The sequence name (i.e., its unique id) is derived from title, taking its first word. 
        The rest of title is taken as sequence description. 
        By default, the sequence is added to the bottom of the alignment. 

        Parameters
        ----------
        title : str
            Sequence title, from which name and description are derived
        sequence: str
            Actual sequence, with gaps encoded as "-" characters
        desc : str, optional
            The description can be directly provided here. If so, title is taken as name instead
        index : int, optional
            The position at which the sequence is inserted. If not provided, it goes last

        Returns
        -------
        None
            None


        Examples
        --------

        >>> ali=Alignment()
        >>> ali.add_seq('seq1 custom nt seq', 'ATTCG-')
        >>> ali.add_seq('seq2 another seq',   '--TTGG')
        >>> print(ali.fasta())
        >seq1 custom nt seq
        ATTCG-
        >seq2 another seq
        --TTGG
        <BLANKLINE>

        >>> ali.add_seq('seq3', 'ATT---', desc='some desc')
        >>> ali.add_seq('seq4', 'ATTGG-', index=0)
        >>> print(ali.fasta())
        >seq4
        ATTGG-
        >seq1 custom nt seq
        ATTCG-
        >seq2 another seq
        --TTGG
        >seq3 some desc
        ATT---
        <BLANKLINE>
        """
        
        if desc is None:
            name=title.split(sep=None, maxsplit=1)[0]
            desc=title[len(name)+1:].strip()
        else:
            name, desc=title, desc
        if self.has_name(name):
            raise AlignmentError(f'add_seq ERROR alignment already has {name}. Did you mean to use set_seq() instead?')
        self._seqs [name] = sequence
        self._desc [name] = desc
        if index is None:
            self._ord.append(name)
        else:
            self._ord.insert(index, name)

        
    def ali_length(self):
        """Returns the number of columns in the alignment (i.e., its width) 

        Returns
        -------
        int 
            The number of columns in the alignment


        Examples
        --------
        >>> ali=Alignment([ ('seq1', 'ATTCG-'), ('seq2', '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.ali_length()
        6

        Warning
        -------
        For best performance, the Alignment class does not check that all sequences have the same length.
        This method simply returns the length of the first sequence.
        To check for homogenous sequence length, see :func:`same_length`

        See Also
        --------
        same_length: check that all sequences are truly aligned, i.e. have the same length  

        """                
        return len(self.get_seq(self._ord[0])) if self._ord else 0

    def same_length(self):
        """Check whether sequences are aligned, i.e. they have the same length

        Returns
        -------
        bool 
            Stating if all sequences have the same lengths

        Examples
        --------
        >>> ali=Alignment([ ('seq1', 'ATTCG-'), ('seq2', '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.same_length()
        True

        >>> ali.add_seq('seqX', 'TATTCGGT-')
        >>> ali.same_length()
        False

        See Also
        --------
        ali_length: length of the alignment (i.e. number of columns)

        """                
                
        if not self.n_seqs: return True
        lenfirst=len( self.get_seq(self._ord[0]) )
        for name in self._ord[1:]:
            if lenfirst != len(self.get_seq(name)):
                return False
        return True

    def n_seqs(self):
        """Returns the number of sequences in the alignment (i.e. the number of rows, or alignment height)

        Returns
        -------
        int
            Number of sequences in the alignment

        Examples
        --------
        >>> ali=Alignment([ ('seq1', 'ATTCG-'), ('seq2', '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.n_seqs()
        3

        See Also
        --------
        ali_length: length of the alignment (i.e. number of columns)
        """
        return len(self._ord)
        
    def set_seq(self, name, sequence):
        """Change the sequence of an entry in-place.

        Parameters
        ----------

        name : str
            The name (i.e. identifier) of the sequence to be altered
        sequence : str
            The new sequence to be set

        Returns
        -------
        None
            None

        Examples
        --------
        >>> ali=Alignment([ ('seq1', 'ATTCG-'), ('seq2', '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali
        # Alignment of 3 sequences and 6 positions
        ATTCG- seq1
        --TTGG seq2
        ATTCG- seq3
        <BLANKLINE>

        >>> ali.set_seq('seq1', 'CHANGE')
        >>> ali 
        # Alignment of 3 sequences and 6 positions
        CHANGE seq1
        --TTGG seq2
        ATTCG- seq3
        <BLANKLINE>
        
        See Also
        --------
        add_seq, get_seq
        """       
        if not self.has_name(name):
            raise AlignmentError(f'set_seq ERROR alignment does not have {name}. Did you mean to use add_seq() instead?')
        self._seqs[name]=sequence

    # def convert_sequences(self, seqfn):     
    #     for name, seq in self:
    #         self.set_seq(name, seqfn(seq))
        

    def set_desc(self, name, desc):
        """Change the description of an entry in-place.

        Parameters
        ----------

        name : str
            The name (i.e. identifier) of the entry to be altered
        desc : str
            The new description to be used

        Returns
        -------
        None
            None

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3 not sure',   'ATTCG-')])
        >>> print(ali.fasta())
        >seq1 this is first
        ATTCG-
        >seq2 this is 2nd
        --TTGG
        >seq3 not sure
        ATTCG-
        <BLANKLINE>

        >>> ali.set_desc('seq3', 'obviously third')
        >>> print(ali.fasta())
        >seq1 this is first
        ATTCG-
        >seq2 this is 2nd
        --TTGG
        >seq3 obviously third
        ATTCG-
        <BLANKLINE>
        """       
        if not self.has_name(name):
            raise AlignmentError(f'set_desc ERROR alignment does not have {name}. Did you mean to use add_seq() instead?')
        self._desc[name]=desc            

    def get_desc(self, name):
        """Returns the description for a sequence entry

        If no sequence with that name is present in the alignment, 
        an AlignmentError exception is raised.

        Parameters
        ----------

        name : str
            The name (i.e. identifier) of the sequence

        Returns
        -------
        str
            The description stored for the sequence

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.get_desc('seq1')
        'this is first'
        >>> ali.get_desc('seq3')
        ''
                
        See Also
        --------
        set_desc
        """       
        if not name in self._desc:
            raise AlignmentError(f'get_desc ERROR alignment does not have {name}')
        return self._desc[name]

    def get_seq(self, name):
        """Returns the sequence with the requested name

        If no sequence with that name is present in the alignment, 
        an AlignmentError exception is raised.

        Parameters
        ----------

        name : str
            The name (i.e. identifier) of the sequence requested

        Returns
        -------
        str
            The sequence requested, including any gaps it may contain

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.get_seq('seq1')
        'ATTCG-'
        >>> ali.get_seq('seq2')
        '--TTGG'
                
        See Also
        --------
        set_seq
        """               
        if not name in self._seqs:
            raise AlignmentError(f'get_seq ERROR alignment does not have {name}')
        return self._seqs[name]

    def names(self):
        """Returns a list of all sequence names in the alignment

        The names returned do not include the sequence descriptions.

        Returns
        -------
        list of str
            An ordered list of sequence names (identifiers) in the alignment

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.names()
        ['seq1', 'seq2', 'seq3']
                
        See Also
        --------
        titles: get all sequence titles, including their name and description
        """               
        return [name for name in self._ord]

    def titles(self): 
        """Returns a list of the all sequence titles  in the alignment

        Each title is the concatenation of sequence name and description, separated by a space. 
        If the description is empty for an entry, only the name is returned

        Returns
        -------
        list of str
            An ordered list of sequence titles in the alignment

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.titles()
        ['seq1 this is first', 'seq2 this is 2nd', 'seq3']
                
        See Also
        --------
        names: get all sequence names (their unique identifiers, without description)
        """               
        return [ ' '.join([name, self.get_desc(name)]).rstrip()  for name in self._ord ]
       
        
    def sequences(self):
        """Returns a list of the all sequences in the alignment

        Returns
        -------
        list of str
            An ordered list of sequences in the alignment (without names or descriptions)

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.sequences()
        ['ATTCG-', '--TTGG', 'ATTCG-']
        """               
        return [self.get_seq(name) for name in self._ord]    
    
    def has_name(self, name):
        """Checks whether the alignment contains a sequence with the name provided

        Parameters
        ----------

        name : str
            The name (i.e. identifier) searched in the alignment. 

        Returns
        -------
        bool
            A bool indicating whether the name is present

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.has_name('seq1')
        True
        >>> ali.has_name('seq2 this is 2nd')  # note: the name would be 'seq2' only
        False
        """
        return name in self._seqs
            
    def fasta(self, nchar=60):
        """Returns the alignment in (aligned) fasta format 

        Parameters
        ----------
        nchar : int, default=60
            The number of characters per line for sequences

        Returns
        -------
        str
            A multiline string with the alignment in fasta format, including sequence descriptions

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> print(ali.fasta())
        >seq1 this is first
        ATTCG-
        >seq2 this is 2nd
        --TTGG
        >seq3
        ATTCG-
        <BLANKLINE>

        See also
        --------
        write: generic function supporting many output formats
        """
        out=''
        for name, seq in self:
            desc=self.get_desc(name)
            out+='>' + name + (" "+desc if desc else "") + '\n'
            for i in range(0, len(seq), nchar):
                out+=seq[i:i + nchar]+'\n'
        return out

    def write(self, fileformat='fasta', to_file=None):
        """Returns a string representation of the alignment in a format of choice

        Internally uses Bio.Align to generate output. 
        Supported fileformat arguments include clustal, stockholm, phylip and many others.
        The full list of supported fileformat arguments is  `provided here  <https://biopython.org/wiki/AlignIO>`_).            


        Parameters
        ----------
        fileformat : str, default='fasta'
            text format requested
        to_file : str | TextIO, optional
            filename or buffer to write into. If not specified, the output is returned instead

        Returns
        -------
        str
            String representation of alignment in the requested format

        Examples
        --------
        >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> print(ali.write('phylip'))
         3 6
        seq1       ATTCG-
        seq2       --TTGG
        seq3       ATTCG-
        <BLANKLINE>
        
        """
        if not to_file is None:
            with (to_file if not isinstance(to_file, str)
                  else open(to_file, 'w')) as fh:
                fh.write( self.write(fileformat=fileformat) )    
        else:        
            if fileformat=='fasta':
                return self.fasta()
            else:
                return format(self.to_biopython(), fileformat)

    def __add__(self, other):  return self.concatenate(other)
    def concatenate(self, other):
        """Concatenate two alignments, i.e., add their sequences one next to the other

        The two alignments must have the same names in the same order 
        or an AlignmentError exception is raised. 
        The sequence descriptions in returned alignment are taken from self

        Parameters
        ----------
        other : Alignment
            alignment that will be concatenated to the right of the self 
            in the returned Alignment   

        Returns
        -------
        Alignment
            alignment with same names as inputs, and sequences resulting from their concatenation

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali2=Alignment([ ('seq1 first', 'TTGC-TAG'), ('seq2 this is 2nd'  , '-ATGGGGC'), ('seq3', 'AATCGGCC')])
        >>> ali.concatenate(ali2)
        # Alignment of 3 sequences and 14 positions
        ATTCG-TTGC-TAG seq1
        --TTGG-ATGGGGC seq2
        ATTCG-AATCGGCC seq3
        <BLANKLINE>

        Note that descriptions in the second alignment are ignored:

        >>> ali3= Alignment([ ('seq1 this desc is ignored', 'TTGC-TAG'), ('seq2'  , '-ATGGGGC'), ('seq3 this also', 'AATCGGCC')])
        >>> print( ali.concatenate(ali3).fasta() )
        >seq1 first
        ATTCG-TTGC-TAG
        >seq2 this is 2nd
        --TTGG-ATGGGGC
        >seq3
        ATTCG-AATCGGCC
        <BLANKLINE>
        """
        if self.names() != other.names():
            raise AlignmentError(f'concatenate ERROR the two alignments must have the same sequence names!')
        a=Alignment()
        for name, seq in self:            
            a.add_seq(name, seq + other.get_seq(name), desc=self.get_desc(name))
        return a

    def copy(self):
        """Returns a copy of the alignment

        Returns
        -------
        Alignment
            copy of the self alignment
        """        
        return self.__class__(self)

    def remove_by_name(self, *names):
        """Remove one or more sequences in the alignment by name in-place.
        
        Note that the modification is done in place. To obtain a new object instead, see examples below.

        Parameters
        ----------
        *names : tuple
            name or names of sequences to be removed from the alignment   

        Returns
        -------
        None
            None

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.remove_by_name('seq1')
        >>> ali
        # Alignment of 2 sequences and 6 positions
        --TTGG seq2
        ATTCG- seq3
        <BLANKLINE>
        
        >>> ali.remove_by_name('seq2', 'seq3')
        >>> ali
        # Empty alignment
        
        To return a new alignment without certain names, do not use this function.
        Instead, use indexing by rows:

        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> names_to_omit=set(  ['seq2']  )
        >>> ali[ [n for n in ali.names() if not n in names_to_omit], :]
        # Alignment of 2 sequences and 6 positions
        ATTCG- seq1
        ATTCG- seq3
        <BLANKLINE>

        See Also
        --------
        remove_by_index
        """
        for name in names:
            if not name in self._seqs:
                raise AlignmentError(f'remove_by_seq ERROR alignment does not have {name}')
            del self._seqs[name]
            del self._desc[name]
        s=set(names)            
        for i in range(len( self._ord )-1, -1, -1):
            if self._ord[i] in s:
                self._ord.pop(i)
                
    def remove_by_index(self, *seqindices):
        """Remove one or more sequences in the alignment by their index, in-place.
        
        The input indices refer to the position of the sequence in the alignment, i.e. their row number.
        Note that the modification is done in place. To obtain a new object instead, see examples below.

        Parameters
        ----------
        *seqindices : tuple
            index or indices of sequences to be removed from the alignment   

        Returns
        -------
        None
            None

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> ali.remove_by_index(0)
        >>> ali
        # Alignment of 2 sequences and 6 positions
        --TTGG seq2
        ATTCG- seq3
        <BLANKLINE>
                
        To return a new alignment without certain sequences, do not use this function.
        Instead, use indexing by rows:

        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
        >>> indices_to_omit=set( [0, 2] )
        >>> ali[ [n for i,n in enumerate(ali.names()) if not i in indices_to_omit], :]
        # Alignment of 1 sequences and 6 positions
        --TTGG seq2
        <BLANKLINE>

        See Also
        --------
        remove_by_name
        """        
        for i in sorted(seqindices, reverse=True):
            name=self._ord.pop(i)
            del self._seqs[name]
            del self._desc[name]

    def remove_empty_seqs(self, inplace=True):
        """Remove all sequences which are entirely made of gaps or that are empty.
        
        By default, removal is done in place.

        Parameters
        ----------
        inplace : bool, default:True
            whether the removal should be done in place. If not, a new Alignment is returned instead

        Returns
        -------
        None or Alignment
            If inplace==True, None is returned; otherwise, a new Alignment without empty sequences is returned

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', '------')])
        >>> ali.remove_empty_seqs()
        >>> ali
        # Alignment of 2 sequences and 6 positions
        ATTCG- seq1
        --TTGG seq2
        <BLANKLINE>
                
        See Also
        --------
        trim_gaps, remove_by_name, remove_by_index
        """
        # boolean array, True when all gaps
        selector=np.char.equal(
            np.array( [np.array(seq, dtype=np.str_) for name, seq in self], dtype=np.str_ ),
            '-'*self.ali_length())
        if inplace:
            empty_seq_names=np.array(self.names())[  selector ]
            self.remove_by_name(*empty_seq_names)
        else:
            return self[np.array(self.names())[~selector],:]
                        
    # def pop(self, index):
    #     name=self._ord.pop(i)
    #     seq=self._seqs[name]
    #     desc=self._desc[name]
    #     return( name, seq, desc )
                                 
    def to_biopython(self):
        """Returns a copy of the alignment as a Bio.Align.MultipleSeqAlignment object  
        
        The SeqRecord instances in the returned MultipleSeqAlignment has their id and name 
        attributes set to sequence names, and also possess the description attribute.

        Returns
        -------
        MultipleSeqAlignment
            Alignment in biopython format (Bio.Align.MultipleSeqAlignment)

        See also
        --------
        to_numpy, to_pandas        
        """
        return MultipleSeqAlignment(
            [SeqRecord( Seq(seq), id=name, name=name,
                        description=self.get_desc(name) )
             for name, seq in self])

    @lru_cache(maxsize=max_cache_size)
    def to_numpy(self):
        """Returns a numpy 2-D array representation of the alignment, useful for vectorized sequence methods
        
        Returns
        -------
        np.ndarray
            The returned array has one row per sequence and one column per alignment position.
            Each value is a single character. The dtype is np.str\_
            Note that rows are not indexed by sequence names, just by their order index

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', '--TT--')])
        >>> print(ali.to_numpy())
        [['A' 'T' 'T' 'C' 'G' '-']
         ['-' '-' 'T' 'T' 'G' 'G']
         ['-' '-' 'T' 'T' '-' '-']]

        See Also
        --------
        to_biopython, to_pandas

        Warning
        ------
        This function is cached for best performance. Thus, do not directly modify the returned object.
        The hash key for caching is derived from sequences only: names are not considered.
        """
        return np.array( [np.array(list(seq), dtype=np.str_) for name, seq in self], dtype=np.str_ )

    def to_pandas(self, use_names=False):
        """Returns a pandas DataFrame representation of the alignment

        Parameters
        ----------
        use_names : bool, optional
            Normally, the returned DataFrame has a simply RangeIndex as index.
            Specify this to instead use sequence names as the index.

        Returns
        -------
        pd.DataFrame
            The returned dataframe has one row per sequence and one column per alignment position.
            Each value is a single character. The dtype is object
            Rows are indexed by the sequence names if use_names==True, or by a RangeIndex by default

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', '--TT--')])
        >>> ali.to_pandas()
           0  1  2  3  4  5
        0  A  T  T  C  G  -
        1  -  -  T  T  G  G
        2  -  -  T  T  -  -

        >>> ali.to_pandas(use_names=True)
              0  1  2  3  4  5
        seq1  A  T  T  C  G  -
        seq2  -  -  T  T  G  G
        seq3  -  -  T  T  -  -

        See Also
        --------
        to_biopython, to_numpy
        """
        p=pd.DataFrame( self.to_numpy() )
        if use_names: p.index=self.names()
        return p

    @lru_cache(maxsize=max_cache_size)
    def conservation_by_column(self):
        """Computes the frequency of characters (nucleotides/amino acids) at each column of the alignment

        Gaps are considered as any other character during computation.        
        The returned object reports frequencies at each position, 
        for all characters which are observed at least once in the alignment.
        This may not correspond to the full nucleotide or protein alphabet, if some characters are not present in the alignment.

        Returns
        -------
        pd.DataFrame
            The returned dataframe has one row per observed character (i.e., nucleotide / amino acid)
            and one column per alignment position. Each value is a float ranging from 0 to 1 representing
            the frequency of that character in that alignment column.

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', '--TT--')])
        >>> ali.conservation_by_column()
                  0         1    2         3         4         5
        -  0.666667  0.666667  0.0  0.000000  0.333333  0.666667
        A  0.333333  0.000000  0.0  0.000000  0.000000  0.000000
        C  0.000000  0.000000  0.0  0.333333  0.000000  0.000000
        G  0.000000  0.000000  0.0  0.000000  0.666667  0.333333
        T  0.000000  0.333333  1.0  0.666667  0.000000  0.000000

        Warning
        ------
        This function is cached for best performance. Thus, do not directly modify the returned object.
        The hash key for caching is derived from sequences only: names are not considered.
        """        
        return self.to_pandas(use_names=False).apply(pd.value_counts, 0, normalize=True).fillna(0.0)
                
    def trim_gaps(self, pct=1.0, inplace=False):
        """Removes the alignment columns with more gaps than specified

        By default, a new alignment without the columns identified as 'too gappy' is returned.        
        
        Parameters
        ----------
        pct     : float, default:1.0
            minimal gap frequency (from 0.0 to 1.0) for a column to be removed
        inplace : bool, default:False
            whether the column removal should be done in place. If not, a new Alignment is returned instead

        Returns
        -------
        None or Alignment
            By default, a new Alignment without empty sequences is returned; if inplace==True, None is returned

        Examples
        --------
        >>> ali= Alignment([ ('seq1 first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTG-'), ('seq3', 'ATTCCG')])
        >>> ali
        # Alignment of 3 sequences and 6 positions
        ATTCG- seq1
        --TTG- seq2
        ATTCCG seq3
        <BLANKLINE>        

        >>> ali.trim_gaps(0.5)
        # Alignment of 3 sequences and 5 positions
        ATTCG seq1
        --TTG seq2
        ATTCC seq3
        <BLANKLINE>        
        """
        ### get positions, a bit of benchmark
        #  [all([s[pos]=='-'  for n,s in a])   for pos in range(a.ali_length())]   # 65.3 ms
        #  a.conservation_by_column(save=False).loc['-']==1.0                      # 2.6 s
        #  a.conservation_by_column(save=True).loc['-']==1.0                       # 276 µs        
        #  np.apply_along_axis( lambda x: np.char.equal(x, '-').all(), 0,   a.to_numpy() )                 # 74.9 ms
        #  np.apply_along_axis( lambda x: np.char.equal(x, '-').all(), 0,   a.to_numpy(save=True) )        # 30.9 ms
        #  np.vectorize(lambda x:x=='-')(  a.to_numpy(save=False)).all(axis=0)     # 87.7 ms
        #  np.vectorize(lambda x:x=='-')(  a.to_numpy(save=True)).all(axis=0)      # 40.5 ms

        # if not self._cbc is None:   #### If self.conservation_by_column was alread run? couldn't built this check procedure
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
        
    def consensus(self, ignore_gaps=None ): # to do: what to do with gaps?
        """Compute the consensus sequence, taking the most represented character for each column

        Parameters
        ----------
        ignore_gaps : float, optional
            By default, gaps are treated as any other character, so that they are returned for columns
            in which they are the most common character.
            If you provide ignore_gaps with a float from 0.0 to 1.0, gaps are not present on the output
            except for columns with a frequency equal or greater than the value provided. 
            For example, a value of 1.0 implies gaps are included only if a column is entirely made of gaps

        Returns
        -------
        str
            The consensus sequence 

        Examples
        --------
        >>> ali= Alignment([ ('seq1', 'ATTCG-'), ('seq2'  , '-TTCGT'), ('seq3', 'ACGCG-'),  ('seq4', 'CTTGGT'), ('seq5', '-TGCT-'), ('seq6', '-TGGG-')])
        >>> ali
        # Alignment of 6 sequences and 6 positions
        ATTCG- seq1
        -TTCGT seq2
        ACGCG- seq3
        CTTGGT seq4
        -TGCT- seq5
        -TGGG- seq6
        <BLANKLINE>

        >>> ali.conservation_by_column()
                  0         1    2         3         4         5
        -  0.500000  0.000000  0.0  0.000000  0.000000  0.666667
        A  0.333333  0.000000  0.0  0.000000  0.000000  0.000000
        C  0.166667  0.166667  0.0  0.666667  0.000000  0.000000
        G  0.000000  0.000000  0.5  0.333333  0.833333  0.000000
        T  0.000000  0.833333  0.5  0.000000  0.166667  0.333333

        >>> ali.consensus()
        '-TGCG-'

        >>> ali.consensus(0.6)
        'ATGCG-'
        """
        if ignore_gaps is None:
            return ''.join(self.conservation_by_column().apply(pd.Series.idxmax))
        else:
            cmap=self.conservation_by_column().copy()
            cmap.loc['-'] [ cmap.loc['-']<ignore_gaps ] = np.nan
            return ''.join(cmap.apply(pd.Series.idxmax))

    def gap_mask(self):
        """Returns a boolean numpy array marking gaps
        
        Returns
        -------
        np.ndarray
            The returned array has one row per sequence and one column per alignment position.           
            Each value is a boolean indicating presence of a gap
            Note that rows are indexed by integers, not by sequence names

        Examples
        --------
        >>> ali= Alignment([ ('seq1', 'ATTCG-'), ('seq2'  , '--TTG-'), ('seq3', 'AT--CG')])
        >>> ali.gap_mask()
        array([[False, False, False, False, False,  True],
               [ True,  True, False, False, False,  True],
               [False, False,  True,  True, False, False]])

        See also
        --------
        terminal_gap_mask
        """
        return self.to_numpy()=='-'

    def terminal_gap_mask(self):
        """Returns a boolean numpy array marking terminal gaps
        
        Terminal gaps are those that have only gaps to their left or to their right

        Returns
        -------
        np.ndarray
            The returned array has one row per sequence and one column per alignment position.           
            Each value is a boolean indicating presence of a terminal gap
            Note that rows are indexed by integers, not by sequence names

        Examples
        --------
        >>> ali= Alignment([ ('seq1', 'ATTCG-'), ('seq2'  , '--TTG-'), ('seq3', 'AT--CG')])
        >>> ali.terminal_gap_mask()
        array([[False, False, False, False, False,  True],
               [ True,  True, False, False, False,  True],
               [False, False, False, False, False, False]])

        See also
        --------
        gap_mask
        """        
        return np.apply_along_axis(
            # below: each bool array is as long as sequence;  function is applied to each row
                      #  bool array: is left term gap? # or # bool array: is right term gap?            
            lambda x: ( np.logical_and.accumulate(x)  |  np.logical_and.accumulate( x[::-1] )[::-1] ),
            1,  self.gap_mask())
        #tested: computing per row is faster than whole matrix

    def position_map(self):
        """Compute a numerical matrix instrumental to map alignment positions to sequence positions (and reverse)
        
        Returns
        -------
        pd.DataFrame
            Returns a Pandas DataFrame with one row per sequence (indexed by name), and one column 
            per alignment position. Each value is the index of that particular sequence (without gaps) 
            in that alignment column. All positions are 0-based.
            For gap positions, the position of the closest non-gap to the left is reported.
            For left terminal gaps, the value reported is -1

        Examples
        --------
        >>> ali= Alignment([ ('seq1', 'ATTCG-'), ('seq2'  , '--TTG-'), ('seq3', 'ATTCCG')])
        >>> ali
        # Alignment of 3 sequences and 6 positions
        ATTCG- seq1
        --TTG- seq2
        ATTCCG seq3
        <BLANKLINE>
        
        >>> ali.position_map()
              0  1  2  3  4  5
        seq1  0  1  2  3  4  4
        seq2 -1 -1  0  1  2  2
        seq3  0  1  2  3  4  5

        Note
        ----
        Computing this matrix makes sense only if you will use positions for many or all sequences.
        For the corresponding operations on single sequences, see functions :func:`position_in_seq` 
        and :func:`position_in_ali`
        
        See also
        --------
        position_in_seq: maps an alignment position to sequence position for a single sequence
        position_in_ali: maps a sequence position to alignment position for a single sequence
        """
        m= np.add.accumulate( (self.to_numpy()!='-'), 1) -1
        return pd.DataFrame( m, index=self.names() ) 
    
    def position_in_seq(self, name, pos_in_ali):
        """Maps an alignment column position to the corresponding position in a certain sequence
        
        Parameters
        ----------
        name : str
            the name of the sequence to map to
        pos_in_ali : int
            0-based alignment position, i.e. the column index

        Returns
        -------
        int 
            0-based position in sequence, i.e. the index of the sequence without counting gaps 
            at the requested column.

            **Note** that for gap positions, the position of the closest non-gap to the left is reported.
            For gap positions, the position of the closest non-gap to the left is reported.

        Examples
        --------
        >>> ali= Alignment([ ('seq1', 'ATTCG-'), ('seq2'  , '--TTG-'), ('seq3', 'ATTCCG')])
        >>> ali
        # Alignment of 3 sequences and 6 positions
        ATTCG- seq1
        --TTG- seq2
        ATTCCG seq3
        <BLANKLINE>
        
        >>> ali.position_in_seq('seq1', 4)
        4
        >>> ali.position_in_seq('seq2', 2)
        0

        Checking a left terminal gap returns -1:

        >>> ali.position_in_seq('seq2', 0)
        -1

        Checking a gap position returns the index of the closest non gap char on the left:

        >>> ali.position_in_seq('seq1', 5)
        4
        
        See also
        --------
        position_in_ali: maps a sequence position to alignment position for a single sequence
        position_map: maps all alignment positions for all sequences
        """    
        s=self.get_seq(name)[:pos_in_ali+1] 
        return pos_in_ali-s.count('-') if s else -1
    
    def position_in_ali(self, name, pos_in_seq):
        """Maps an position in a certain sequence (without counting gaps) to its corresponding position in the alignment

        If the requested position is invalid, raise an IndexError
        
        Parameters
        ----------
        name : str
            the name of the sequence
        pos_in_seq : int
            0-based sequence position, i.e. the index mapping to the requested sequence without gaps

        Returns
        -------
        int 
            0-based alignment position , i.e. the column index where the requested sequence position is found

        Examples
        --------
        >>> ali= Alignment([ ('seq1', 'ATTCG-'), ('seq2'  , '--TTG-'), ('seq3', 'AT-CCG')])
        >>> ali
        # Alignment of 3 sequences and 6 positions
        ATTCG- seq1
        --TTG- seq2
        AT-CCG seq3
        <BLANKLINE>
        
        >>> ali.position_in_ali('seq1', 4)
        4
        >>> ali.position_in_ali('seq2', 0)
        2
        >>> ali.position_in_ali('seq3', 2)
        3
        
        See also
        --------
        position_in_seq: maps an alignment position to sequence position for a single sequence
        position_map: maps all alignment positions for all sequences
        """
        if pos_in_seq<0:
            raise IndexError(f'position_in_ali ERROR pos_in_seq requested is negative: {pos_in_seq}')
        p=-1
        for i, c in enumerate( self.get_seq(name) ):
            if c!='-': p+=1
            if  p==pos_in_seq:
                return i
        raise IndexError(f'position_in_ali ERROR pos_in_seq requested ({pos_in_seq}) was greater than sequence length ({self.get_seq(name)})')
    
    def score_similarity(self, targets=None, gaps='y'): # add gaps?

        if not self.same_length():
            raise AlignmentError('score_similarity ERROR sequences are not aligned!')
        
        if targets is None:
            targets=self
        else:
            if not isinstance(targets, Alignment):
                raise AlignmentError(f'score_similarity ERROR "targets" argument must receive an Alignment, instead got {targets.__class__}')
            if not targets.same_length():
                raise AlignmentError('score_similarity ERROR "targets" argument does not have same alignment length as self')
            
        cmap=self.conservation_by_column()
        mc=pd.DataFrame( cmap.unstack() ).set_axis(['conservation'], axis=1).rename_axis( ['position', 'letter'] )

        # removing lines corresponding to gaps: we don't want to score them positively. Later in merging they'll be 0
        mc=mc[  ~mc.index.isin(['-'], level='letter') ]
        """                 conservation
        position letter
        0        A                0.6
                 C                0.0
                 D                0.1
                 E                0.0
        ...                      ...
        """

        mts=( targets.to_pandas(use_names=True).unstack().rename('letter')
          .rename_axis( ['position', 'name']).reset_index(1).set_index('letter', append=True) )
        """              name
        position letter
        0        A       seq1
                 -       seq2
                 A       seq3
        """
        
        
        ## mask by gaps here
        #######
        jt=mts.join(mc, on=None, how='left').rename(columns={'conservation':'av_identity'})    ##join on index, add NAs if missing from mts
        """              name  av_identity
        position letter
        0        -       seq2           NaN
                 A       seq1      0.666667
                 A       seq3      0.666667
        1        -       seq2           NaN
                 T       seq1      0.666667
        """
        # setting score of gaps (and also characters present in targets but not in self) from Na to 0
        jt.fillna(value=0, inplace=True)
        
        # compute weights: max conservation at this position (e.g. if A is most common letter, how common it is)
        we=mc.max(level='position').rename(columns={'conservation':'weight'})
        """       weight
        position
        0         0.666667
        1         0.666667
        """

        ## adding weights per ali position 
        jt=jt.join(we,   on='position',  how='left')
        """              name    av_identity weight
        position letter
        0        -       seq2      0.000000  0.666667
                 A       seq1      0.666667  0.666667
                 A       seq3      0.666667  0.666667
        1        -       seq2      0.000000  0.666667
                 T       seq1      0.666667  0.666667
        """

        jt['w_av_identity']=jt['av_identity']*jt['weight']
        
        
        ### initing scores with:
        # ASI: average sequence identity
        scores=jt.drop('weight', axis=1).groupby('name').mean().rename(columns={'av_identity':'ASI',  'w_av_identity':'AWSI' })
        return scores
        

        
        
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
       ...          ...   ...
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
        ...                      ...
        """
        
        res=pd.merge( pp, q, left_on=['position', 'letter'], right_index=True)
        """          position letter  conservation
        seqindex
        0                0      M      1.000000
        1                0      M      1.000000
        2                0      M      1.000000
        0                1      -      0.333333
        1                1      F      0.333333
       ...          ...  ...          ...
        1               28      R      0.333333
        2               28      K      0.333333
        0               29      M      0.333333
        1               29      -      0.333333
        2               29      F      0.333333        
        [90 rows x 3 columns] """

        ## average sequence identity per sequence
        res.groupby(res.index).conservation.mean()
        

if __name__ == "__main__":
    import doctest, sys
    thepath = os.path.dirname(os.path.dirname( os.path.abspath( __file__)   )) + '/..'
    os.chdir(thepath)
    sys.path.insert(0, thepath)
    
    doctest.testmod()
                    
