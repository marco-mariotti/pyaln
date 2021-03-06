Tutorial of pyaln
=================

Welcome to the tutorial of pyaln. Here, you will learn how to use the Alignment
class to read, write, process and characterize key features of multiple sequence alignments.


.. contents:: Contents of Tutorial
   :depth: 3
   


Introducing the Alignment class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This class is the core of pyaln, and represents a multiple alignment of homologous sequences.
Sequences can be of any type (nucleotide, protein, or custom characters).
Gaps **must** be encoded as dashes, i.e. ``"-"``.

Each entry is uniquely identified by a *name*, with an optional *description*.

In many file formats (e.g. aligned fasta), an extensive *title* is associated
to each sequence. This will include some form of identifier, plus other information
such as gene/protein name, source etc. When reading alignment files, such *titles*
are split into *name* (the first word, must be unique per alignment) and
*description* (the remainder of the title).

The Alignment class comes into play when you have already aligned sequences.
These may have been generated by any of the numerous aligner methods out there (for example:
`ClustalOmega <http://www.clustal.org/omega/>`_,
`Mafft <https://mafft.cbrc.jp/alignment/software/>`_,
`T-coffee <http://tcoffee.crg.cat/>`_).

The rationale of pyaln Alignment is to provide a convenient and efficient interface for
reading, writing, manipulating, and profiling alignments. Under the hood, pyaln employs
Numpy and Pandas for computationally intensive tasks.


Tutorial set-up
~~~~~~~~~~~~~~~
For the examples below to work correctly, after :doc:`installing pyaln<installation>`,
open python and run this before anything else::
  
  >>> from pyaln import Alignment, pyaln_folder


Reading and writing alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Aligned sequences are **loaded** at the time at the creation of an Alignment object.
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
their identifiers. Alignment accepts any iterable of ``(title, sequence)``::
  
  >>> ex_ali=Alignment([ ('seq1 description1', 'ATTCG-'), ('seq2 desc2', '--TTGG'), ('seq3', 'ATTCG-')])

Various options are available for **writing** alignments. If you print an Alignment,
you will obtain a reduced representation, showing its number of sequences and length::
  
  >>> print(fep_ali)
  # Alignment of 6 sequences and 138 positions
  MWLTLVALLALCATGRTAENLSESTTDQDKLVIARGKLVAPSVVGUSIKKMPELYNFLM...L Fep15_danio_rerio
  MWAFLLLTLAFSATGMTEE-DVTDTAIEERPVIAKGILKAPSVVGUAIKKMPALYMFLM...L Fep15_S_salar
  MWIFLLLTLAFSATGMTEE-NVTDTAIEERPVIAKGILKAPSVVGUAIKKMPELYTFLM...L Fep15_O_mykiss
  MWAFLVLTFAVAA-GASET-VDNHTAAEEKLLIARGKLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_rubripes
  MWALLVLTFAVTV-GASEE-VKNQTAAEEKLVIARGTLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_nigroviridis
  MWAFVLIAFSV---GASDS--SNSTAE----VIARGKLMAPSVVGUAIKKLPELNRFLM...L Fep15_O_latipes
  <BLANKLINE>


However, note that this representation may not include the full sequence, and omits
descriptions.

On the other hand, :func:`~pyaln.Alignment.write` method of Alignment offers a variety of output formats
(again through Bio.AlignIO, `see the full list here <https://biopython.org/wiki/AlignIO>`_ ).
The most common, *fasta*, includes sequence descriptions::
  
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

  
Indexing and transversing alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alignments have two dimensions. By *length* of the alignment, we refer to its width, meaning the number
of alignment columns (aka alignment positions). The other dimension is the *number of sequences* in the alignment
(i.e. its height). These features can be inspected by the methods :func:`~pyaln.Alignment.ali_length`,
and :func:`~pyaln.Alignment.n_seqs`, or at once through the property :func:`~pyaln.Alignment.shape`,
as illustrated below::

  >>> ali=Alignment([ ('seq1 this is a seq', 'ATTCG-'), ('seq2 another seq', '--TTGG'), ('seq3', 'ATTCG-')])
  >>> print( [ali.ali_length(), ali.n_seqs(), ali.shape] )
  [6, 3, (3, 6)]
  
You can slice portions of an Alignment (i.e. take on some sequences and/or some columns) by **indexing** it.
The format is ``Alignment[rows_selector, column_selector]``, where: 

        - The ``rows_selector`` can be an integer (i.e., the vertical position of 
          the sequence in the alignment), or a slice thereof (e.g. ``2:5``), or a list of sequence names.
        - The ``column_selector`` is a integer index (i.e. the horizontal position in the alignment),
          or a slice thereof, or a list of (start, end) indices, or a Numpy boolean array.


.. warning::
   As customary in python, in pyaln all positions are 0-based, and intervals are specified with
   their start included and their end excluded.

For example, we load this small alignment::
  
    >>> ali=Alignment([ ('seq1 this is a seq', 'ATTCG-'), ('seq2 another seq', '--TTGG'), ('seq3', 'ATTCG-')])    
    >>> ali
    # Alignment of 3 sequences and 6 positions
    ATTCG- seq1
    --TTGG seq2
    ATTCG- seq3
    <BLANKLINE>

Let's get the alignment of first two sequences only::

    >>> ali[:2,:]
    # Alignment of 2 sequences and 6 positions
    ATTCG- seq1
    --TTGG seq2
    <BLANKLINE>

We could have done the same by specifying sequences by name::

    >>> ali[ ['seq1', 'seq2'], : ]
    # Alignment of 2 sequences and 6 positions
    ATTCG- seq1
    --TTGG seq2
    <BLANKLINE>

Now let's take the alignment without the first and last columns::

    >>> ali[:,1:-1]
    # Alignment of 3 sequences and 4 positions
    TTCG seq1
    -TTG seq2
    TTCG seq3
    <BLANKLINE>

We can take non-contigous alignment regions by indexing columns with a list of ``(start, end)`` elements.
For example, to get the 1st, 2nd, and 6th position in a single step::

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
    

  
To **iterate** through the sequences in the alignment (i.e. its rows), use a ``for`` loop.
This will yield tuples like ``(name, sequence)``. To get the description of a sequence, use :func:`~pyaln.Alignment.get_desc`.

For example, here we print the name, sequence length, and description of each sequence
(in the same order as they are found in the alignment)::

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

  
Working with alignment columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may want to determine the composition of each column, meaning the frequencies 
of observed characters at each specific  position.
Since alignment columns represent homologous positions in the aligned sequences
and frequencies represent the conservation at those positions, this is referred to as *conservation map*.

Like some other methods of the Alignment class, :func:`~pyaln.Alignment.conservation_map` returns a
`Pandas <https://pandas.pydata.org>`_ DataFrame::

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

Pandas Series (basically the data type of each column of a DataFrame) may be used to index
the columns of pyaln Alignment. This may be convenient, for example, to take all alignment columns
at least one ``"T"`` character::

  >>> ali[:,   ali.conservation_map().loc['T']>0 ]
  # Alignment of 3 sequences and 3 positions
  TTC seq1
  -TT seq2
  TTC seq3
  <BLANKLINE>

A very common operation with alignments involves removing those columns featuring too many gaps.
This is often referred to as **trimming alignments**, and it is achieved through the function
:func:`~pyaln.Alignment.trim_gaps`.

For example, let's remove all columns with at least 50% gaps::
  
  >>> ali.trim_gaps(0.5)
  # Alignment of 3 sequences and 5 positions
  ATTCG seq1
  --TTG seq2
  ATTCG seq3
  <BLANKLINE>
  

Another common operation is **alignment concatenation**: two or more alignments corresponding to different gene families, but
coming from the same set of species, are combined into one. Visually, alignment concatenation corresponds to
stacking one alignment next to the other horizontally. This is achieved in pyaln by adding 
two Alignment instances with the ``+`` operator
(or analogously, calling the :func:`~pyaln.Alignment.concatenate` function).

  >>> ali2=Alignment([ ('seq1', 'AAATAAAA'), ('seq2'  , '-AAGAAAG'), ('seq3', 'ACATAAAC')])
  >>> ali + ali2
  # Alignment of 3 sequences and 14 positions
  ATTCG-AAATAAAA seq1
  --TTGG-AAGAAAG seq2
  ATTCG-ACATAAAC seq3
  <BLANKLINE>

Note that if the two alignments being added do not have exaclty the same *names*, an error occurs.
  
Adding a string to an Alignment is equivalent to adding its content to each sequence of the alignment::
  
  >>> ali + 'NNNN' + ali2
  # Alignment of 3 sequences and 18 positions
  ATTCG-NNNNAAATAAAA seq1
  --TTGGNNNN-AAGAAAG seq2
  ATTCG-NNNNACATAAAC seq3
  <BLANKLINE>
  
Sequence identity
~~~~~~~~~~~~~~~~~  

There are various methods implemented in pyaln to estimate the degree of similarity of sequences in the alignment.
In general, they are based on **sequence identity**. At first glance, this is a very straightforward concept:
the sequence identity of two sequences is the number of identical positions, divided by their length.
In this example, 4/5  -> 80%

  >>> from pyaln.sequtils import sequence_identity
  >>> sequence_identity('ATGCA',
  ...                   'ATGCC')
  0.8
  
However, when gaps come into the picture, things get a little more complicated, as you may choose to score them in a few different ways.
Pyaln offers four options in this regard, each identified by a single letter ``gaps`` code:

#. ``gaps='y'``: gaps are considered and considered mismatches. This is the **default** behaviour.
#. ``gaps='n'``: gaps are ignored
#. ``gaps='t'``: terminal gaps (those at the beginning or the end of sequences) are ignored. Others are considered as in ``'y'``.
#. ``gaps='a'``: gaps are considered as any other character; even gap-to-gap matches are scored as identities

These options can be provided to :func:`~pyaln.sequtils.sequence_identity` and other pyaln methods.
Let's see a few examples of their behavior::

  >>> from pyaln.sequtils import sequence_identity
  >>> seq1='--ATC-GGG-'
  >>> seq2='AAATCGGGGC'
  >>> seq3='--ACC-CCGC'
  >>> ali=Alignment( [('seq1', seq1), ('seq2', seq2), ('seq3', seq3)] )

The first two sequences are identical, but `seq2` has three insertions (i.e. gapped regions) compared to `seq1`.
Comparing them with ``gaps='y'`` will consider all positions (including gaps) as total sequence length,
effectively scoring negatively gaps::
  
  >>> sequence_identity(seq1, seq2, gaps='y')
  0.6

On the other hand, if we ignore gaps with ``gaps='n'``, we obtain 100% sequence identity::
  
  >>> sequence_identity(seq1, seq2, gaps='n')
  1.0

In certain applications, you may want to ignore terminal gaps with ``gaps='t'``.
In this case, this means that the `seq1` subsequence ``ATC-GGG`` is effectively compared to the corresponding
region of `seq2`, resulting in 6/7 --> ~0.86 ::

  >>> sequence_identity(seq1, seq2, gaps='t')
  0.8571428571428571
  
The option ``gaps='a'`` is not recommended for biological alignments. This behaves similarly to ``gaps='y'``,
but with an important difference.
When comparing two sequences coming an alignment that contains many additional ones, it is possible that the two
sequences both have a gap in one or more positions::

  >>> print ( seq1+'\n'+seq3 )
  --ATC-GGG-
  --ACC-CCGC

If we compare them naively, counting all identical characters without differentiating gaps (i.e., the behavior of ``gaps='a'``),
we end up scoring shared gaps positively, with 6/10 matches::

  >>> sequence_identity(seq1, seq3, gaps='a')
  0.6

Shared gaps should be ignored in any pairwise comparison, which is the behavior followed under any other
value of ``gaps`` (``'y', 'n', 't'``)::

  >>> sequence_identity(seq1, seq3, gaps='y')   # 3/7
  0.42857142857142855

  >>> sequence_identity(seq1, seq3, gaps='n')   # 3/6
  0.5

  >>> sequence_identity(seq1, seq3, gaps='t')   # 3/6
  0.5
  
 
  
The function :func:`~pyaln.Alignment.score_similarity` allows to compute
the **Average Sequence Identity (ASI)** of each sequence, when compared to the whole alignment.
This is equivalent to calling the function :func:`~pyaln.sequtils.sequence_identity` introduced above
in all-against-all fashion (but it is implemented differently for better performance).
This measure is instrumental  estimate the overall similarity of sequence in the alignment.

::

  >>> fep_ali=Alignment(pyaln_folder + '/examples/fep15_protein.fa', fileformat='fasta')
  >>> fep_ali.score_similarity()
  metrics                    ASI
  Fep15_danio_rerio     0.777778
  Fep15_S_salar         0.826334
  Fep15_O_mykiss        0.822684
  Fep15_T_rubripes      0.829599
  Fep15_T_nigroviridis  0.815000
  Fep15_O_latipes       0.767438


The :func:`~pyaln.Alignment.score_similarity` method accepts the ``gaps`` parameter to define how to treat gaps.
You may provide a single ``gaps`` argument, or provide multiple ones at once to assess how results would differ::

  >>> fep_ali.score_similarity(gaps=['y', 'n', 't', 'a'])
  gaps                         y         n         t         a
  metrics                    ASI       ASI       ASI       ASI
  Fep15_danio_rerio     0.777778  0.793051  0.777778  0.777778
  Fep15_S_salar         0.826334  0.838283  0.826334  0.827295
  Fep15_O_mykiss        0.822684  0.834522  0.822684  0.823671
  Fep15_T_rubripes      0.829599  0.842566  0.829599  0.830918
  Fep15_T_nigroviridis  0.815000  0.835351  0.815000  0.816425
  Fep15_O_latipes       0.767438  0.805693  0.767438  0.769324

Besides ASI, this method  may also return a variant called
**Average Weighted Sequence Identity (AWSI)**, wherein the most conserved positions in the alignment are given
higher weight. For details, see :func:`~pyaln.Alignment.score_similarity`.
::
   
   >>> fep_ali.score_similarity(metrics=['i', 'w'],  gaps='y')
   metrics                    ASI      AWSI
   Fep15_danio_rerio     0.777778  0.847123
   Fep15_S_salar         0.826334  0.885040
   Fep15_O_mykiss        0.822684  0.882183
   Fep15_T_rubripes      0.829599  0.887255
   Fep15_T_nigroviridis  0.815000  0.874389
   Fep15_O_latipes       0.767438  0.834809
                               
   
These sequence metrics may be employed to assess how some external sequences *fit* in a core alignment.
This may be instrumental to check whether some candidate sequences appear to belong to a certain gene family.
In the following example, we load an alignment containing the same sequences as `fep_ali` above,
with the addition of an extra candidate sequence. We want to test whether this sequence resembles other sequences
in a similar degree as they resemble each other.

::

   >>> cand_ali=Alignment(pyaln_folder + '/examples/fep15_protein.with_candidate.fa', fileformat='fasta')
   >>> cand_ali
   # Alignment of 7 sequences and 163 positions
   MWLTLVALLALCATGRTAENLSESTTDQDKLVIARGKLVAPSVVGUSIKKMPELYNFLM...L Fep15_danio_rerio
   MWAFLLLTLAFSATGMTEE-DVTDTAIEERPVIAKGILKAPSVVGUAIKKMPALYMFLM...L Fep15_S_salar
   MWIFLLLTLAFSATGMTEE-NVTDTAIEERPVIAKGILKAPSVVGUAIKKMPELYTFLM...L Fep15_O_mykiss
   MWAFLVLTFAVAA-GASET-VDNHTAAEEKLLIARGKLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_rubripes
   MWALLVLTFAVTV-GASEE-VKNQTAAEEKLVIARGTLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_nigroviridis
   MWAFVLIAFSV---GASDS--SNSTAE----VIARGKLMAPSVVGUAIKKLPELNRFLM...L Fep15_O_latipes
   ----------------------------------------QSCGGUQLNRLREVKAFVT...L Fep15_candidate
   <BLANKLINE>
   
Let's see the ASI and AWSI metrics for the core alignment (all sequences except the last one):
   
   >>> cand_ali[:-1,:].score_similarity( metrics='iw', gaps='yn' )
   gaps                         y                   n          
   metrics                    ASI      AWSI       ASI      AWSI
   Fep15_danio_rerio     0.777778  0.847123  0.793051  0.856044
   Fep15_S_salar         0.826334  0.885040  0.838283  0.893412
   Fep15_O_mykiss        0.822684  0.882183  0.834522  0.890497
   Fep15_T_rubripes      0.829599  0.887255  0.842566  0.896094
   Fep15_T_nigroviridis  0.815000  0.874389  0.835351  0.891288
   Fep15_O_latipes       0.767438  0.834809  0.805693  0.860639

Now let's see the same metrics but comparing the candidate to the same set of sequences.
This is achieved through the ``targets`` argument of :func:`~pyaln.Alignment.score_similarity`::

  >>> cand_ali[:-1,:].score_similarity( targets=cand_ali[ ['Fep15_candidate'] ,:], metrics='iw', gaps='yn' )
  gaps                    y                   n          
  metrics               ASI      AWSI       ASI      AWSI
  Fep15_candidate  0.213043  0.282854  0.349844  0.362332

We can see that the metrics are well outside the range of the similarity metrics of the core alignments,
indicating that the sequence does not fit in the family just as well. Indeed, this protein is from another family.


Biopython, Numpy, and Pandas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sequences are stored in pyaln Alignment objects in form of built-in string types.
This ensures the most common operations are as fast as possible. For certain procedures, however,
an alternative representation is generated on the fly.

Specifically, pyaln transforms Alignment objects into MultipleSeqAlignments from
`Biopython AlignIO <http://biopython.org/DIST/docs/api/Bio.AlignIO-module.html>`_ to access a variety of
Input / Output capabilities.

On the other hand, fast vectorized operations on alignment columns are performed using alignment representations as
`Numpy <https://numpy.org/>`_ array (one row per sequence, one column per alignment position).
A similar representation, slightly slower but more versatile, is also employed: the
`Pandas <https://pandas.pydata.org/>`_ DataFrame.

Conversions back and forth from these alternative representations of alignments automatically occur
under the hood of pyaln when they are convenient for efficient computation.
If you wish to build on top of pyaln and may find these representation useful, then check the documentation of these methods:

  - :func:`~pyaln.Alignment.to_biopython()`
  - :func:`~pyaln.Alignment.to_numpy()`
  - :func:`~pyaln.Alignment.to_pandas()`
  - :func:`~pyaln.Alignment.from_numpy()`
    





  


