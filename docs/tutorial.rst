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
  
  from pyaln import Alignment, pyaln_folder


Reading alignments
~~~~~~~~~~~~~~~~~~
To use pyaln, you must have already aligned sequences. These may have been generated
by any of the many aligner methods out there (for example:
`ClustalOmega <http://www.clustal.org/omega/>`,
`Mafft <https://mafft.cbrc.jp/alignment/software/>`,
`T-coffee <http://tcoffee.crg.cat/>`).

You load aligned sequences into an Alignment object directly by instantiating it.

In the next few examples, we load alignment files located in pyaln examples folder::
  
  >>> filename=pyaln_folder + '/examples/fep15_protein.fa'
  >>> fep_ali=Alignment(filename, fileformat='fasta')

Many file formats are supported; through Bio.AlignIO (see `a full list here <https://biopython.org/wiki/AlignIO>`_).
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

If you print an Alignment, you will obtain a representation including its size::
  >>> print(ex_ali)
  # Alignment of 6 sequences and 138 positions
  MWLTLVALLALCATGRTAENLSESTTDQDKLVIARGKLVAPSVVGUSIKKMPELYNFLM...L Fep15_danio_rerio
  MWAFLLLTLAFSATGMTEED-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPALYMFLM...L Fep15_S_salar
  MWIFLLLTLAFSATGMTEEN-VTDTAIEERPVIAKGILKAPSVVGUAIKKMPELYTFLM...L Fep15_O_mykiss
  MWAFLVLTFAVAA-GASETV-DNHTAAEEKLLIARGKLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_rubripess
  MWALLVLTFAVTV-GASEEV-KNQTAAEEKLVIARGTLLAPSVVGUGIKKMPELHHFLM...L Fep15_T_nigroviridis
  MWAFVLIAFSVGA---SDS------SNSTAEVIARGKLMAPSVVGUAIKKLPELNRFLM...L Fep15_O_latipes

However, note that this representation may not include the full sequence, and omits
descriptions.

The *write* method of Alignment offers a variety of output formats (through Bio.AlignIO,
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

>>> ali=Alignment('examples/example_ali.fa')

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

>>> ali=Alignment('examples/example_ali.stockholm', fileformat='stockholm')
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
