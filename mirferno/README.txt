miRferno - miRNA-Target prediction module of sPARTA
Updated: version-1.12rc 6/10/2015

======== Description ========
miRferno is the miRNA-target prediction module of small RNA-PARE Target Analyzer (sPARTA).
miRferno utilizes multi-core servers to achieve two-dimensional parallelization
in order to maintain a low memory footprint, imperative to achieve a full genome analysis. 

======== Dependencies ========
miRferno requires bowtie2 in the PATH variable of the user account executing sPARTA
bowtie2 may be downloaded here http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

miRferno requires the following python3 functions to perform properly:
pyfasta - https://pypi.python.org/pypi/pyfasta/

========= Note ===========
miRFerno uses file extensions to identify file types, naming meta-data and selectively 
cleaning up temp files. Therefore, it is recommended to have appropriate file extensions. 
For Ex. a genome/cDNA FASTA file should have '.fa' extension. Please see 'Arguments' section (below) for recommended file extensions.

2.
Make sure that input fasta files do not have integers in name. For ex - test.1.fa or arabidopsis.new.2.4.fa
Files with such names are deleted sometimes while cleanup operation

3.
All the input files 1) miRNAs 2) FASTA file for genome or transcripts and 3) degradome/PARE in tag-count format should be in same directory,including sPARTA script


======== Execution ========
There are command line arguments that are to be used by miRferno for proper
execution. For the first execution, all steps must be performed, but
once this has been completed, provided the genome is same,
the entire analysis will not need to be repeated. Examples of such cases
may be seen below.


======= Arguments ========
-gffFile        GFF3 file for the species being analyzed corresponding
...             to the genome assembly being used. Recommended file
				extension - '.gff' or '.gff3' 

-genomeFile     Genome file in FASTA format that will be used to extract 
...             features (genic or intergenic regions) using GFF3 file.
				Recommended file extension - '.fa'

-featureFile	FASTA file containing sequences of interest (CDS, transcript,
...             intergenic regions etc.) if user already has a set of
...             sequences. This option is mutually exclusive to genome file and
...             gff file. So either genomefile along with gffFile is used or
...             feature set is supplied directly.Recommended file extension - '.fa'

-genomeFeature  0 if prediction is to be done in genic region. 1 if prediction
...             is to be done in intergenic region

-miRNAFile      FASTA format of miRNA sequences. Recommended file extension - '.fa'

-tarPred        Mode of target prediction. H for heuristic. E for exhaustive.
...             H is default if no mode is specified

-misMat		Maximum mismatches allowed for target prediction. 5 is default value
...		with sum of mismatches and wobbles <= 6 

-wob		Maximum wobbles allowed for target prediction. 5 is default value
...		with sum of mismatches and wobbles <= 6

-mirBul		Maximum nucleotides allowed to form symmetrical and asymmetrical bulges
...		in miRNA. 1 is default value

-tarBul		Maximum nucleotides allowed to form symmetrical and asymmetrical bulges
..             in target. 2 is default value

tarScore        Scoring mode for target prediction. S for seedless. N for
...             normal. S is default if no mode is specified 

--repeats       Flag to include PARE reads from repetitive regions

accel           Y to use balanced multiple process scheme or else specify the
...             number of processors to be used. Y is default

======== Genome and Annotation Data ======== 
Both the GFF3 file and corresponding genome FASTA file can be downloaded from
Phytozome [http://www.phytozome.net/]

==============Examples ==================
1. Execution on new genome

This execution should be performed any time a new genome file (along with corresponding GFF file) is being analyzed:

python3 miRferno.py -genomeFile <genomeFile.fa> -gffFile <GFF3file> -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -tarPred -tarScore

OR

a user provided feature set (FASTA file with sequences of interest) is being analyzed:

python3 miRferno.py -featureFile <featureFile.fa> -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -tarPred -tarScore

2. Execution on genome in which genome/feature set has already been processed

This execution should be performed if a genome file has been processed previously but the miRNAs for which targets need to be predicted are new:

python3 miRferno.py -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -tarPred -tarScore


======== Output ==========

Target prediction results can be found in 'predicted' folder under the name 'All.targs.parsed.csv'

======= Contact ===========

Atul Kakrana
kakrana@udel.edu

Reza Hammond
hammond@dbi.udel.edu

===== END of README =======
