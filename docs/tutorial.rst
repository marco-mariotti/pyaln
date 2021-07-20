Tutorial of pyaln
=================

Welcome to the tutorial of pyaln. Here, you will learn how to use the Alignment
class to read, write, process and profile features of  multiple sequence alignment.


.. contents:: Contents of Tutorial
   :depth: 3
   

Installation
------------
We recommend to use conda to install pyaln. With conda installed, run this in a terminal::
  
  conda install -c mmariotti pyaln

Alternatively, you can use pip::
  
  pip install pyaln

To check that installation was successful, you may run::
  
  python -c 'import pyaln'


The Alignment class
-------------------
This class is the core of pyaln, and represents a multiple sequence alignment,
where sequences can be of any type (nucleotide, protein, or custom characters).
Gaps **must** be encoded as `-`.

Each entry is uniquely identified by a *name*, with an optional *description*.

In many file formats (e.g. aligned fasta), an extensive *title* is associated
to each sequence. This will include some form of identifier, plus other information
such as gene/protein name, source etc. When reading alignment files, such *titles*
are split into *name* (the first word, must be unique per alignment) and
*description* (the remainder of the title).

Tutorial set-up
~~~~~~~~~~~~~~~
For the examples below to work correctly, run this in python before anything else::
  
  >>> from pyaln import Alignment, pyaln_folder


Reading alignments
~~~~~~~~~~~~~~~~~~
To use pyaln, you must have already aligned sequences. These may have been generated
by any of the many aligner methods out there (for example:
`ClustalOmega <http://www.clustal.org/omega/>`_,
`Mafft <https://mafft.cbrc.jp/alignment/software/>`_,
`T-coffee <http://tcoffee.crg.cat/>`_).

You load aligned sequences into an Alignment object directly by instantiating it.

In the next few examples, we load alignment files located in pyaln examples folder::
  
  >>> filename=pyaln_folder + '/examples/fep15_protein.fa'
  >>> fep_ali=Alignment(filename, fileformat='fasta')

Many file formats are supported, thanks to Bio.AlignIO (see `a full list here <https://biopython.org/wiki/AlignIO>`_).
For a few common cases, extensions are recognized so it is not compulsory to specify format::
  
  >>> fep_ali2=Alignment(pyaln_folder+'/examples/fep15_protein.stockholm')
  >>> sbp2_ali=Alignment(pyaln_folder+'/examples/SBP2_protein.aln')

Alignments can also be instanced with a IO buffer rather than filename::
  
  >>> fb=open(pyaln_folder+'/examples/SBP2_protein.aln')
  >>> sbp2_ali2=Alignment(fb)  

You may also initialize an alignment manually by providing aligned sequences and
their identifiers. Alignment accepts any iterable of `[title, sequence]`::
  
  >>> ex_ali=Alignment([ ('seq1 description1', 'ATTCG-'), ('seq2 desc2', '--TTGG'), ('seq3', 'ATTCG-')])

Writing alignments
~~~~~~~~~~~~~~~~~~

If you print an Alignment, you will obtain a reduced representation, showing its number of sequences and length::
  
  >>> print(fep_ali)
  # Alignment of 6 sequences and 138 positions
  MWLTLVALLALCATGRTAENLSESTTDQDKLVIARGKLVAPSVVGUSIKKMPELYNFLM...L Fep15_danio_rerio
  MWAFLLLTLAFSATGMTEED-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPALYMFLM...L Fep15_S_salar
  MWIFLLLTLAFSATGMTEEN-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPELYTFLM...L Fep15_O_mykiss
  MWAFLVLTFAVAA-GASETV-DNHTAAEEKLLIARGKLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_rubripess
  MWALLVLTFAVTV-GASEEV-KNQTAAEEKLVIARGTLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_nigroviridis
  MWAFVLIAFSVGA---SDS------SNSTAEVIARGKLMAPSVVGUAIKKLPELNRFLM...L Fep15_O_latipes
  <BLANKLINE>
  
However, note that this representation may not include the full sequence, and omits
descriptions.

The :func:`~pyaln.Alignment.write` method of Alignment offers a variety of output formats (through Bio.AlignIO,
see `a full list here <https://biopython.org/wiki/AlignIO>`_).
The most common, fasta, includes sequence descriptions::
  
  >>> ex_ali=Alignment([ ('seq1 description1', 'ATTCG-'), ('seq2 desc2', '--TTGG'), ('seq3', 'ATTCG-')])
  >>> print( ex_ali.write('fasta') )
  >seq1 description1
  ATTCG-
  >seq2 desc2
  --TTGG
  >seq3
  ATTCG-
  <BLANKLINE>  

The *write* method also accepts a *to_file* argument to write directly to a file::
  
  >>> ex_ali.write('clustal', to_file='ali_file.aln')    # doctest: +SKIP

  
Indexing alignments
~~~~~~~~~~~~~~~~~~~
You can slice portions of an Alignment (i.e. take on some sequences and/or some columns) by **indexing** it.
The format is `Alignment[rows_selector, column_selector]`, where: 

        - The `rows_selector` can be an integer (i.e., the vertical position of 
          the sequence in the alignment), or a slice thereof (e.g. `2:5`), or a list of sequence names.
        - The `column_selector` is a integer index (i.e. the horizontal position in the alignment),
          or a slice thereof, or a list of (start, end) indices, or a Numpy boolean array.


For example, let's create a small alignment::
  
    >>> ali=Alignment([ ('seq1 this is a seq', 'ATTCG-'), ('seq2 another seq', '--TTGG'), ('seq3', 'ATTCG-')])    
    >>> ali
    # Alignment of 3 sequences and 6 positions
    ATTCG- seq1
    --TTGG seq2
    ATTCG- seq3
    <BLANKLINE>

Now, let's get the alignment of first two sequences only::

    >>> ali[:2,:]
    # Alignment of 2 sequences and 6 positions
    ATTCG- seq1
    --TTGG seq2
    <BLANKLINE>

Get the subalignment of two sequences, by their name::

    >>> ali[ ['seq1', 'seq3'], : ]
    # Alignment of 2 sequences and 6 positions
    ATTCG- seq1
    ATTCG- seq3
    <BLANKLINE>

Remove the first and last alignment columns::

    >>> ali[:,1:-1]
    # Alignment of 3 sequences and 4 positions
    TTCG seq1
    -TTG seq2
    TTCG seq3
    <BLANKLINE>

Index columns with a list of (start, end) elements, to get the 1st, 2nd, and 6th position in a single step::

    >>> ali[:, [(0,2), (5, 6)]]
    # Alignment of 3 sequences and 3 positions
    AT- seq1
    --G seq2
    AT- seq3
    <BLANKLINE>
			    
Indexing by row and column at once, to get the 1st character of all sequences except the last::

   >>> ali[:-1, 0:1]
   # Alignment of 2 sequences and 1 positions
   A seq1
   - seq2
   <BLANKLINE> 
 
Complex column selection can be performed by providing a Numpy boolean array.
For example, take all columns except for the 3rd and 4th::

  >>> import numpy as np
  >>> colsel=np.array( [True, True, False, False, True, True] ) 
  >>> ali[:, colsel]
  # Alignment of 3 sequences and 4 positions
  ATG- seq1
  --GG seq2
  ATG- seq3
  <BLANKLINE>
    
Iterating over the alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To go through the sequences in the alignment (i.e. its rows), use a *for* loop.
This will yield tuples like `(name, sequence)`. To get the description of a sequence, use `get_desc(name)`.

For example, here we print the length and description of each sequence, in the same order of the alignment::

  >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
  >>> for name, seq in ali:
  ...   print(  (name, len(seq),     ali.get_desc(name)) )
  ('seq1', 6, 'this is first')
  ('seq2', 6, 'this is 2nd')
  ('seq3', 6, '')

To iterate over alignment positions instead (i.e. its columns) use the :func:`~pyaln.Alignment.positions` method.

For example, here we check at each position whether the two sequences ('seq1' and 'seq2') have the same character::

  >>> for i in ali.positions():
  ...    print(  (i, ali.get_seq('seq1')[i]   ==  ali.get_seq('seq2')[i])  )
  (0, False)
  (1, False)
  (2, True)
  (3, False)
  (4, True)
  (5, False)

  
Sequence composition
~~~~~~~~~~~~~~~~~~~~

An important characteristic of alignments is the composition of each column, meaning the frequencies
of all observed characters at that particular position. Since columns represent homologous positions
and frequencies represent the conservation at those positions, this is referred to as *conservation map*.

Like some other methods of the Alignment class, :func:`~pyaln.Alignment.conservation_map` returns a Pandas DataFrame::

  >>> ali=Alignment([ ('seq1 this is first', 'ATTCG-'), ('seq2 this is 2nd'  , '--TTGG'), ('seq3', 'ATTCG-')])
  >>> ali
  # Alignment of 3 sequences and 6 positions
  ATTCG- seq1
  --TTGG seq2
  ATTCG- seq3
  <BLANKLINE>

  >>> ali.conservation_map()
            0         1    2         3    4         5
  -  0.333333  0.333333  0.0  0.000000  0.0  0.666667
  A  0.666667  0.000000  0.0  0.000000  0.0  0.000000
  C  0.000000  0.000000  0.0  0.666667  0.0  0.000000
  G  0.000000  0.000000  0.0  0.000000  1.0  0.333333
  T  0.000000  0.666667  1.0  0.333333  0.0  0.000000

  
Similarity metrics
~~~~~~~~~~~~~~~~~~

There are various implemented measures to estimate the overall degree of similarity of sequences in the alignment.

In general, they are based on the concept of sequence identity. At first glance, this is very straightforward:
the sequence identity of two sequences is the number of identical positions, divided by their length.
In this example, 4/5  -> 80%

  >>> from pyaln.sequtils import sequence_identity
  >>> sequence_identity('ATGCA',
  ...                   'ATGCC')
  0.8
  
However, when gaps enter in the picture, things get a little more complicated, as you may choose to score them in a few different ways.

...

:func:`~pyaln.sequtils.sequence_identity`

...

Please see function :func:`~pyaln.Alignment.score_similarity` for computation of average sequence identity (ASI)
and average weighted sequence identity (AWSI)




  
