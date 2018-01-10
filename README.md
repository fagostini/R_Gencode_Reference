## Gencode Annotation Processing

If they are not present in the working folder, the pipeline will download the following files:
 - __Comprehensive gene annotation__: It contains the comprehensive gene annotation on the reference chromosomes only.
 - __Long non-coding RNA gene annotation__: It contains the comprehensive gene annotation of lncRNA genes on the reference chromosomes (this is a subset of the main annotation file);
 - __Predicted tRNA genes__: tRNA genes predicted by ENSEMBL on the reference chromosomes using tRNAscan-SE (this dataset does not form part of the main annotation file);
 - __Chromosome sizes__: The file contains the size in nucleorides of each chromosome.

__Note__: the pipeline will search for the genome assembly and Gencode annotation versions specified within the script. Therefore, if you plan to use a different species/assembly/annotation make sure you change the script accordingly. 

The _Comprehensive gene annotation_ file is employed to extract information on genes, transcripts and exons, while the _Long non-coding RNA gene annotation_ is used to include additional non-coding RNAs (assigned as such by Gencode) to those predicted using sequences from [Rfam](http://rfam.xfam.org/) and [miRBase](http://www.mirbase.org/).

By default, the analysis is restricted to standard chromosome only, and to level 1 (validated), 2 (manual annotation) and 3 (automated annotation) genes.

### Output

- __\<genome>\_Gencode\<version>\_annotations.RData__: It contains most of the annotation objects, including genes, transcripts and exons GRanges, genes and transcripts metadata, and longest pre- and mature RNA GRanges.
- __\<genome>\_Gencode\<version>\_annotations.pc.transcript.regions.RData__: It contains protein-coding genes GRanges sub-divided into 5’-UTR, exonic and 3’-UTR regions.
- __\<genome>\_Gencode\<version>\_annotations.all.genes.transcript.regions.RData__: It contains all the annotated features, where each annotated nucleotide is assigned to a transcript biotype using the following hierarchy:
  ncRNA > cds > utr3 > utr5 > intron > other > intergenic
