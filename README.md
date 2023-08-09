# gfflu

[![PyPI - Version](https://img.shields.io/pypi/v/gfflu.svg)](https://pypi.org/project/gfflu)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/gfflu.svg)](https://pypi.org/project/gfflu)
[![CI](https://github.com/CFIA-NCFAD/gfflu/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/CFIA-NCFAD/gfflu/actions/workflows/ci.yml)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gfflu/README.html)
[![](https://img.shields.io/conda/dn/bioconda/gfflu.svg?style=flat)](https://anaconda.org/bioconda/gfflu)

`gfflu` is a Python CLI app to generate annotations of Influenza A virus (IAV) gene segment nucleotide sequences with 
[BLASTX][] and [Miniprot][]  using the same protein sequences as [Influenza Virus Sequence Annotation Tool][] and 
output a [GFF3][] file with the expected genetic features for each of the 8 IAV gene segments.

-----

**Table of Contents**

- [Usage](#usage)
- [Installation](#installation)
  - [Conda](#conda) 
  - [PyPI](#pypi)
  - [From Source](#from-source)
- [Annotation](#annotation)
  - [Segment 1](#segment-1)
  - [Segment 2](#segment-2)
  - [Segment 3](#segment-3)
  - [Segment 4](#segment-4)
  - [Segment 5](#segment-5)
  - [Segment 6](#segment-6)
  - [Segment 7](#segment-7)
  - [Segment 8](#segment-8)
- [License](#license)

## Usage

Below is an example of typical usage with a FASTA nucleotide sequence file ([Segment_4_HA.MH201222.fasta](tests/data/Segment_4_HA.MH201222.fasta)):

```console
gfflu Segment_4_HA.MH201222.fasta
```

Produces an output directory `gfflu-outdir/` by default with the following files:

```console
$ tree gfflu-outdir/
gfflu-outdir/
├── Segment_4_HA.MH201222.blastx.tsv
├── Segment_4_HA.MH201222.faa
├── Segment_4_HA.MH201222.gbk
├── Segment_4_HA.MH201222.gff
└── Segment_4_HA.MH201222.miniprot.gff

1 directory, 4 files
```

> Specify output directory with `-o /path/to/outdir`

Help output:

```console
 Usage: gfflu [OPTIONS] FASTA                                                                                                                                                                                                           
                                                                                                                                                                                                                                        
 Annotate Influenza A virus sequences using Miniprot and BLASTX                                                                                                                                                                         
 The Miniprot GFF for a particular reference sequence gene segment will have multiple annotations for the same gene. This script will select the top scoring annotation for each gene and write out a new GFF file that can be used     
 with SnpEff.                                                                                                                                                                                                                           
                                                                                                                                                                                                                                        
╭─ Arguments ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *    fasta      FILE  Influenza virus nucleotide sequence FASTA file [default: None] [required]                                                                                                                                      │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --outdir              -o      PATH  Output directory [default: gfflu-outdir]                                                                                                                                                         │
│ --force               -f            Overwrite existing files                                                                                                                                                                         │
│ --prefix              -p      TEXT  Output file prefix [default: None]                                                                                                                                                               │
│ --verbose             -v                                                                                                                                                                                                             │
│ --version             -V            Print 'gfflu version 0.0.2' and exit                                                                                                                                                             │
│ --install-completion                Install completion for the current shell.                                                                                                                                                        │
│ --show-completion                   Show completion for the current shell, to copy it or customize the installation.                                                                                                                 │
│ --help                              Show this message and exit.                                                                                                                                                                      │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
                                                                                                                                                                                                                                        
 gfflu version 0.0.2; Python 3.10.5               
```

## Installation

### Conda

This is the recommended installation method.

```console
conda install -c bioconda gfflu
```

### PyPI

```console
pip install gfflu
```

> This install method assumes that you have [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and [Miniprot][]  
> installed and on your `$PATH`.

### From Source

Recommended to use [conda](https://docs.conda.io/en/latest/) to manage the environment from
the provided `environment.yml` file.

```console
git clone https://github.com/CFIA-NCFAD/gfflu.git
cd gfflu
conda env create -f environment.yml
conda activate gfflu
```

## Annotation

`gfflu` outputs a [SnpEff][] compatible GFF with the same features identified as the 
[Influenza Virus Sequence Annotation Tool][].

### Segment 1

[Influenza Virus Sequence Annotation Tool][] output
```
>Feature MH201221
16	2295	gene		
			gene	PB2
16	2295	CDS		
			product	polymerase PB2
			protein_id	MH201221p1
			gene	PB2
    
 INFO: Length: 2316 nucleotides
 INFO: Segment: 1 (PB2)
 INFO: Sequence completeness: protein 1 - complete; nucleotide - complete
 INFO: This sequence (MH201221) contains following signature mutation(s) that might confer high virulence of the virus: (E627K)
 INFO: Virus type: influenza A
```

NCBI Genbank GFF for [MH201221.1](https://www.ncbi.nlm.nih.gov/nuccore/MH201221.1)
```
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
##sequence-region MH201221.1 1 2316
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=11320
MH201221.1	Genbank	region	1	2316	.	+	.	ID=MH201221.1:1..2316;Dbxref=taxon:11320;Name=1;gbkey=Src;isolation-source=embyonated chicken eggs;mol_type=viral cRNA;note=laboratory-derived;segment=1;serotype=H1N1;strain=A/PR/8_RGCDC-4%2C6/34
MH201221.1	Genbank	gene	16	2295	.	+	.	ID=gene-PB2;Name=PB2;gbkey=Gene;gene=PB2;gene_biotype=protein_coding
MH201221.1	Genbank	CDS	16	2295	.	+	0	ID=cds-AVY92608.1;Parent=gene-PB2;Dbxref=NCBI_GP:AVY92608.1;Name=AVY92608.1;gbkey=CDS;gene=PB2;product=polymerase PB2;protein_id=AVY92608.1
```

`gfflu` GFF
```
##gff-version 3
##sequence-region MH201221 1 2295
MH201221	miniprot	gene	16	2295	3747	+	.	ID=gene-PB2;Identity=0.9631;Name=PB2;Positive=0.9842;Rank=1;Target=PB2%7CCDS%7Cpolymerase_PB2%7CSeg1prot1A 1 759;gene=PB2;gene_biotype=protein_coding
MH201221	miniprot	CDS	16	2295	3747	.	0	ID=cds-PB2;Identity=0.9631;Parent=gene-PB2;Rank=1;Target=PB2%7CCDS%7Cpolymerase_PB2%7CSeg1prot1A 1 759;gene=PB2;product=polymerase PB2
```

### Segment 2

[Influenza Virus Sequence Annotation Tool][] output
```
>Feature CY147460
13	2286	gene		
			gene	PB1
13	2286	CDS		
			product	polymerase PB1
			protein_id	CY147460p1
			gene	PB1
107	370	gene		
			gene	PB1-F2
107	370	CDS		
			product	PB1-F2 protein
			protein_id	CY147460p2
			gene	PB1-F2
    
 INFO: Length: 2316 nucleotides
 INFO: Segment: 2 (PB1)
 INFO: Sequence completeness: protein 1 - complete; protein 2 - complete; nucleotide - complete
 INFO: Virus type: influenza A
```

NCBI Genbank GFF for [CY147460.1](https://www.ncbi.nlm.nih.gov/nuccore/CY147460.1)

```
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
##sequence-region CY147460.1 1 2316
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1343803
CY147460.1	Genbank	region	1	2316	.	+	.	ID=CY147460.1:1..2316;Dbxref=taxon:1343803;Name=2;collection-date=1934;country=Puerto Rico;gbkey=Src;lab-host=? + egg2 passage(s);mol_type=viral cRNA;nat-host=human;note=Strain PR8-LVD2 is phenotypically distinct from PR8 molecular clone;segment=2;serotype=H1N1;strain=A/Puerto Rico/8-LVD2/1934
CY147460.1	Genbank	sequence_feature	1	2316	.	+	.	ID=id-CY147460.1:1..2316;Dbxref=IRD:NIGSP_JY2_00027.PB1;gbkey=misc_feature
CY147460.1	Genbank	gene	13	2286	.	+	.	ID=gene-PB1;Name=PB1;gbkey=Gene;gene=PB1;gene_biotype=protein_coding
CY147460.1	Genbank	CDS	13	2286	.	+	0	ID=cds-AGQ47939.1;Parent=gene-PB1;Dbxref=NCBI_GP:AGQ47939.1;Name=AGQ47939.1;gbkey=CDS;gene=PB1;product=polymerase PB1;protein_id=AGQ47939.1
CY147460.1	Genbank	gene	107	370	.	+	.	ID=gene-PB1-F2;Name=PB1-F2;gbkey=Gene;gene=PB1-F2;gene_biotype=protein_coding
CY147460.1	Genbank	CDS	107	370	.	+	0	ID=cds-AGQ47940.1;Parent=gene-PB1-F2;Dbxref=NCBI_GP:AGQ47940.1;Name=AGQ47940.1;gbkey=CDS;gene=PB1-F2;product=PB1-F2 protein;protein_id=AGQ47940.1
```

`gfflu` GFF

```
##gff-version 3
##sequence-region CY147460 1 2286
CY147460        miniprot        gene    13      2286    3892    +       .       ID=gene-PB1;Identity=0.9762;Name=PB1;Positive=0.9974;Rank=1;Target=PB1%7CCDS%7Cpolymerase_PB1%7Cseg2prot1B 1 757;gene=PB1;gene_biotype=protein_coding
CY147460        miniprot        CDS     13      2286    3892    .       0       ID=cds-PB1;Identity=0.9762;Parent=gene-PB1;Rank=1;Target=PB1%7CCDS%7Cpolymerase_PB1%7Cseg2prot1B 1 757;gene=PB1;product=polymerase PB1
CY147460        feature gene    107     370     .       +       .       ID=gene-PB1-F2;Target=PB1-F2%7CCDS%7CPB1-F2_protein%7Cseg2prot2M;gene=PB1-F2;gene_biotype=protein_coding
CY147460        feature CDS     107     370     .       +       0       ID=cds-PB1-F2;Parent=gene-PB1-F2;Target=PB1-F2%7CCDS%7CPB1-F2_protein%7Cseg2prot2M;gene=PB1-F2;product=PB1-F2 protein
```


### Segment 3

[Influenza Virus Sequence Annotation Tool][] output

```
>Feature CY146806
13	2163	gene		
			gene	PA
13	2163	CDS		
			product	polymerase PA
			protein_id	CY146806p1
			gene	PA
13	772	gene		
			gene	PA-X
13	582	CDS		
584	772			
			product	PA-X protein
			protein_id	CY146806p2
			exception	ribosomal slippage
			gene	PA-X
    
 INFO: Length: 2208 nucleotides
 INFO: Segment: 3 (PA)
 INFO: Sequence completeness: protein 1 - complete; protein 2 - complete; nucleotide - complete
 INFO: Virus type: influenza A
```

NCBI Genbank GFF for [CY146806.1](https://www.ncbi.nlm.nih.gov/nuccore/CY146806.1)

```
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
##sequence-region CY146806.1 1 2208
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1346461
CY146806.1	Genbank	region	1	2208	.	+	.	ID=CY146806.1:1..2208;Dbxref=taxon:1346461;Name=3;country=USA: Texas;gbkey=Src;lab-host=? + egg2 passage(s);mol_type=viral cRNA;nat-host=human;segment=3;serotype=H3N2;strain=A/Texas/JY2/unknown
CY146806.1	Genbank	sequence_feature	1	2208	.	+	.	ID=id-CY146806.1:1..2208;Dbxref=IRD:NIGSP_JY2_00014.PA;gbkey=misc_feature
CY146806.1	Genbank	gene	13	2163	.	+	.	ID=gene-PA;Name=PA;gbkey=Gene;gene=PA;gene_biotype=protein_coding
CY146806.1	Genbank	CDS	13	2163	.	+	0	ID=cds-AGO00320.1;Parent=gene-PA;Dbxref=NCBI_GP:AGO00320.1;Name=AGO00320.1;gbkey=CDS;gene=PA;product=polymerase PA;protein_id=AGO00320.1
CY146806.1	Genbank	gene	13	772	.	+	.	ID=gene-PA-X;Name=PA-X;gbkey=Gene;gene=PA-X;gene_biotype=protein_coding
CY146806.1	Genbank	CDS	13	582	.	+	0	ID=cds-AGO00321.1;Parent=gene-PA-X;Dbxref=NCBI_GP:AGO00321.1;Name=AGO00321.1;exception=ribosomal slippage;gbkey=CDS;gene=PA-X;product=PA-X protein;protein_id=AGO00321.1
CY146806.1	Genbank	CDS	584	772	.	+	0	ID=cds-AGO00321.1;Parent=gene-PA-X;Dbxref=NCBI_GP:AGO00321.1;Name=AGO00321.1;exception=ribosomal slippage;gbkey=CDS;gene=PA-X;product=PA-X protein;protein_id=AGO00321.1
```

> TODO: handle/add "exception=ribosomal slippage" to PA-X CDS

`gfflu` GFF
```
##gff-version 3
##sequence-region CY146806 1 2163
CY146806	miniprot	gene	13	2163	3758	+	.	ID=gene-PA;Identity=0.9986;Name=PA;Positive=1.0000;Rank=1;Target=PA%7CCDS%7Cpolymerase_PA%7Cseg3prot 1 716;gene=PA;gene_biotype=protein_coding
CY146806	miniprot	CDS	13	2163	3758	.	0	ID=cds-PA;Identity=0.9986;Parent=gene-PA;Rank=1;Target=PA%7CCDS%7Cpolymerase_PA%7Cseg3prot 1 716;gene=PA;product=polymerase PA
CY146806	miniprot	gene	13	772	1301	+	.	Frameshift=1;ID=gene-PA-X;Identity=0.9987;Name=PA-X;Positive=0.9987;Rank=1;Target=PA-X%7CCDS%7CPA-X_protein%7Cseg3prot2C 1 252;gene=PA-X;gene_biotype=protein_coding
CY146806	miniprot	CDS	13	772	1301	.	0	Frameshift=1;ID=cds-PA-X;Identity=0.9987;Parent=gene-PA-X;Rank=1;Target=PA-X%7CCDS%7CPA-X_protein%7Cseg3prot2C 1 252;gene=PA-X;product=PA-X protein
```

### Segment 4

[Influenza Virus Sequence Annotation Tool][] output

```
>Feature MH201222.1
21	1721	gene		
			gene	HA
21	1721	CDS		
			product	hemagglutinin
			protein_id	MH201222.1p1
			function	receptor binding and fusion protein
			gene	HA
21	71	sig_peptide
72	1052	mat_peptide
			product HA1
1053	1718	mat_peptide
			product HA2
    
 INFO: Length: 1753 nucleotides
 INFO: Segment: 4 (HA)
 INFO: Sequence completeness: protein 1 - complete; nucleotide - complete
 INFO: Serotype: H1
 INFO: Virus type: influenza A
```

NCBI Genbank GFF for [MH201222.1](https://www.ncbi.nlm.nih.gov/nuccore/MH201222.1)

```
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
##sequence-region MH201222.1 1 1753
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=11320
MH201222.1	Genbank	region	1	1753	.	+	.	ID=MH201222.1:1..1753;Dbxref=taxon:11320;Name=4;gbkey=Src;isolation-source=embyonated chicken eggs;mol_type=viral cRNA;note=laboratory-derived;segment=4;serotype=H1N1;strain=A/PR/8_RGCDC-4%2C6/34
MH201222.1	Genbank	gene	21	1721	.	+	.	ID=gene-HA;Name=HA;gbkey=Gene;gene=HA;gene_biotype=protein_coding
MH201222.1	Genbank	CDS	21	1721	.	+	0	ID=cds-AVY92609.1;Parent=gene-HA;Dbxref=NCBI_GP:AVY92609.1;Name=AVY92609.1;gbkey=CDS;gene=HA;product=hemagglutinin;protein_id=AVY92609.1
MH201222.1	Genbank	signal_peptide_region_of_CDS	21	71	.	+	.	ID=id-AVY92609.1:1..17;Parent=cds-AVY92609.1;gbkey=Prot
MH201222.1	Genbank	mature_protein_region_of_CDS	72	1052	.	+	.	ID=id-AVY92609.1:18..344;Parent=cds-AVY92609.1;gbkey=Prot;product=HA1
MH201222.1	Genbank	mature_protein_region_of_CDS	1053	1718	.	+	.	ID=id-AVY92609.1:345..566;Parent=cds-AVY92609.1;gbkey=Prot;product=HA2
```

`gfflu` GFF

```
##gff-version 3
##sequence-region MH201222 1 1721
MH201222	miniprot	gene	21	1721	2545	+	.	ID=gene-HA;Identity=0.8233;Name=HA;Positive=0.8993;Rank=1;Target=HA%7CCDS%7Chemagglutinin%7Cseg4protA 1 566;gene=HA;gene_biotype=protein_coding
MH201222	miniprot	CDS	21	1721	2545	.	0	ID=cds-HA;Identity=0.8233;Parent=gene-HA;Rank=1;Target=HA%7CCDS%7Chemagglutinin%7Cseg4protA 1 566;gene=HA;product=hemagglutinin
MH201222	feature	signal_peptide_region_of_CDS	21	71	.	+	.	ID=signal_peptide-HA;Parent=cds-HA,gene-HA
MH201222	miniprot	mature_protein_region_of_CDS	72	1052	1413	+	0	ID=mature_protein-HA;Identity=0.7737;Parent=cds-HA,gene-HA;Rank=1;Target=HA%7Cmature_protein_region_of_CDS%7CHA1%7Cseg4matureA2 1 327;product=HA1
MH201222	miniprot	mature_protein_region_of_CDS	1053	1718	1109	+	0	ID=mature_protein-HA;Identity=0.9279;Parent=cds-HA,gene-HA;Rank=1;Target=HA%7Cmature_protein_region_of_CDS%7CHA2%7Cseg4matureA3 1 222;product=HA2
```

### Segment 5

[Influenza Virus Sequence Annotation Tool][] output

```
>Feature MH085254
44	1540	gene		
			gene	NP
44	1540	CDS		
			product	nucleocapsid protein
			protein_id	MH085254p1
			gene	NP
    
 INFO: Length: 1561 nucleotides
 INFO: Segment: 5 (NP)
 INFO: Sequence completeness: protein 1 - complete; nucleotide - complete
 INFO: Virus type: influenza A
```

NCBI Genbank GFF for [MH085254.1](https://www.ncbi.nlm.nih.gov/nuccore/MH085254.1)

```
```


`gfflu` GFF
```
##gff-version 3
##sequence-region MH085254 1 1540
MH085254	miniprot	gene	44	1540	2469	+	.	ID=gene-NP;Identity=0.9438;Name=NP;Positive=0.9819;Rank=1;Target=NP%7CCDS%7Cnucleocapsid_protein%7Cseg5prot 1 498;gene=NP;gene_biotype=protein_coding
MH085254	miniprot	CDS	44	1540	2469	.	0	ID=cds-NP;Identity=0.9438;Parent=gene-NP;Rank=1;Target=NP%7CCDS%7Cnucleocapsid_protein%7Cseg5prot 1 498;gene=NP;product=nucleocapsid protein
```

### Segment 6

```
>Feature EF190976
21	1385	gene		
			gene	NA
21	1385	CDS		
			product	neuraminidase
			protein_id	EF190976p1
			gene	NA
    
 INFO: Length: 1413 nucleotides
 INFO: Segment: 6 (NA)
 INFO: Sequence completeness: protein 1 - complete; nucleotide - complete
 INFO: Serotype: N1
 INFO: Virus type: influenza A
```


`gfflu` GFF
```
##gff-version 3
##sequence-region EF190976 1 1385
EF190976	miniprot	gene	21	1385	2231	+	.	ID=gene-NA;Identity=0.8681;Name=NA;Positive=0.9149;Rank=1;Target=NA%7CCDS%7Cneuraminidase%7Cseg6prot1A 1 470;gene=NA;gene_biotype=protein_coding
EF190976	miniprot	CDS	21	1385	2231	.	0	ID=cds-NA;Identity=0.8681;Parent=gene-NA;Rank=1;Target=NA%7CCDS%7Cneuraminidase%7Cseg6prot1A 1 470;gene=NA;product=neuraminidase
```

### Segment 7

[Influenza Virus Sequence Annotation Tool][] output

```
>Feature MH085255
24	782	gene		
			gene	M1
24	782	CDS		
			product	matrix protein 1
			protein_id	MH085255p1
			gene	M1
24	1005	gene		
			gene	M2
24	49	CDS		
738	1005			
			product	matrix protein 2
			protein_id	MH085255p2
			gene	M2
    
 INFO: Length: 1023 nucleotides
 INFO: Segment: 7 (MP)
 INFO: Sequence completeness: protein 1 - complete; protein 2 - complete; nucleotide - complete
 INFO: This sequence (MH085255) contains following signature mutation(s) that might confer amantadine resistance: (V27A) (S31N)
 INFO: Virus type: influenza A
```

`gfflu` GFF

```
##gff-version 3
##sequence-region MH085255 1 1005
MH085255	miniprot	gene	24	782	1238	+	.	ID=gene-M1;Identity=0.9683;Name=M1;Positive=0.9881;Rank=1;Target=M1%7CCDS%7Cmatrix_protein_1%7Cseg7prot1 1 252;gene=M1;gene_biotype=protein_coding
MH085255	miniprot	CDS	24	782	1238	.	0	ID=cds-M1;Identity=0.9683;Parent=gene-M1;Rank=1;Target=M1%7CCDS%7Cmatrix_protein_1%7Cseg7prot1 1 252;gene=M1;product=matrix protein 1
MH085255	miniprot	gene	24	1005	435	+	.	ID=gene-M2;Identity=0.8454;Name=M2;Positive=0.9072;Rank=1;Target=M2%7CCDS%7Cmatrix_protein_2%7Cseg7prot2A 1 97;gene=M2;gene_biotype=protein_coding
MH085255	miniprot	CDS	24	49	41	+	0	ID=cds-M2;Identity=1.0000;Parent=gene-M2;Rank=1;Target=M2%7CCDS%7Cmatrix_protein_2%7Cseg7prot2A 1 8;gene=M2;product=matrix protein 2
MH085255	miniprot	CDS	738	1005	394	.	1	ID=cds-M2;Identity=0.8295;Parent=gene-M2;Rank=1;Target=M2%7CCDS%7Cmatrix_protein_2%7Cseg7prot2A 9 97;gene=M2;product=matrix protein 2
```


### Segment 8

[Influenza Virus Sequence Annotation Tool][] output

```
>Feature MH085256
25	717	gene		
			gene	NS1
25	717	CDS		
			product	nonstructural protein 1
			protein_id	MH085256p1
			gene	NS1
25	862	gene		
			gene	NEP
			gene_syn	NS2
25	54	CDS		
527	862			
			product	nuclear export protein
			note	nonstructural protein 2
			protein_id	MH085256p2
			gene	NEP
    
 INFO: Length: 886 nucleotides
 INFO: Segment: 8 (NS)
 INFO: Sequence completeness: protein 1 - complete; protein 2 - complete; nucleotide - complete
 INFO: Virus type: influenza A
```

`gfflu` GFF

```
##gff-version 3
##sequence-region MH085256 1 862
MH085256	miniprot	gene	25	862	553	+	.	ID=gene-NS2;Identity=0.9174;Name=NS2;Positive=0.9339;Rank=1;Target=NS2%7CCDS%7Cnonstructural_protein_2%7Cseg8prot2A 1 121;gene=NS2;gene_biotype=protein_coding
MH085256	miniprot	CDS	25	54	44	+	0	ID=cds-NS2;Identity=0.9000;Parent=gene-NS2;Rank=1;Target=NS2%7CCDS%7Cnonstructural_protein_2%7Cseg8prot2A 1 10;gene=NS2;product=nonstructural protein 2
MH085256	miniprot	CDS	527	862	509	.	0	ID=cds-NS2;Identity=0.9189;Parent=gene-NS2;Rank=1;Target=NS2%7CCDS%7Cnonstructural_protein_2%7Cseg8prot2A 11 121;gene=NS2;product=nonstructural protein 2
MH085256	miniprot	gene	25	717	1074	+	.	ID=gene-NS1;Identity=0.9130;Name=NS1;Positive=0.9565;Rank=1;Target=NS1%7CCDS%7Cnonstructural_protein_1%7Cseg8prot1J 1 230;gene=NS1;gene_biotype=protein_coding
MH085256	miniprot	CDS	25	717	1074	.	0	ID=cds-NS1;Identity=0.9130;Parent=gene-NS1;Rank=1;Target=NS1%7CCDS%7Cnonstructural_protein_1%7Cseg8prot1J 1 230;gene=NS1;product=nonstructural protein 1
```


## License

`gfflu` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.


## References

- [BLASTX]
- [Miniprot]
- [SnpEff]
- [Influenza Virus Sequence Annotation Tool]
- [GFF3]

[BLASTX]: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
[Miniprot]: https://github.com/lh3/miniprot
[SnpEff]: https://pcingola.github.io/SnpEff/
[Influenza Virus Sequence Annotation Tool]: https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation
[GFF3]: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md