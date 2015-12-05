sPARTA: small RNA-PARE Target Analyzer Version 
Updated: version-1.13 11/24/2015

======== Description ========
small RNA-PARE Target Analyzer (sPARTA) is a tool which utilizes
high-throughput sequencing to profile genome-wide cleavage products.
sPARTA begins with a built-in parallelized target prediction module for plant
miRNAs called 'miRferno'. sPARTA as a whole utilizes multi-core servers to
achieve two-dimensional parallelization in order to maintain a low memory
footprint, imperative to achieve a full genome analysis. 

======== Dependencies ========
sPARTA requires bowtie2 in the PATH variable of the user account executing sPARTA
bowtie2 may be downloaded here http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

sPARTA requires the following python3 functions to perform properly:
numpy - http://www.numpy.org/
pyfasta - https://pypi.python.org/pypi/pyfasta/
rpy2 - http://rpy.sourceforge.net/
scipy - http://www.scipy.org/

These may easily be installed using (Python) PIP. Intructions to install PIP - https://pip.pypa.io/en/stable/installing.html

========= Note ===========
1.
sPARTA uses file extensions to identify file types, naming meta-data and selectively 
cleaning up temp files. Therefore, it is recommended to have appropriate file extensions. 
For Ex. a genome/cDNA FASTA file should have '.fa' extension. 
Please see 'Arguments' section (below) for recommended file extensions.

2.
Make sure that input fasta files do not have integers in name. For ex - test.1.fa or arabidopsis.new.2.4.fa
Files with such names are deleted sometimes while cleanup operation

3.
All the input files 1) miRNAs 2) FASTA file for genome or transcripts and 3) degradome/PARE in tag-count format should be in same directory,including sPARTA script

======== Execution ========
There are command line arguments that are to be used by sPARTA for proper
execution. For the first execution, all steps must be performed, but
once this has been completed, provided the miRNAs and genome are the same,
the entire analysis will not need to be repeated. Examples of such cases
may be seen below.



======= Arguments ========
gffFile         GFF3 file for the species being analyzed corresponding
...             to the genome assembly being used. Recommended file
				extension - '.gff' or '.gff3' 
genomeFile      Genome file in FASTA format that will be used to extract 
...             features (genic or intergenic regions) using GFF3 file.
				Recommended file extension - '.fa'

featureFile		FASTA file containing sequences of interest (CDS, transcript,
...             intergenic regions etc.) if user already has a set of
...             sequences. This option is mutually exclusive to genome file and
...             gff file. So either genomefile along with gffFile is used or
...             feature set is supplied directly. Recommended file extension - '.fa'

genomeFeature   0 if prediction is to be done in genic region. 1 if prediction
...             is to be done in intergenic region

miRNAFile       FASTA format of miRNA sequences. Recommended file extension - '.fa'
tarPred         Mode of target prediction. H for heuristic. E for exhaustive.
...             H is default if no mode is specified
tarScore        Scoring mode for target prediction. S for seedless. N for
...             normal. S is default if no mode is specified 

libs            List of PARE library files in tag count format. Data can be
...             easily converted into tag count format using *********
minTagLen       Minimum length of PARE tag. Tags shorter than minTagLen will
...             be discarded. 20 is default
maxTagLen       Maximum length of PARE tag. Tags longer than maxTagLen will be
...             chopped to the specified length. 30 is default
--tag2FASTA     Convert tag count file for PARE libraries to FASTA files for
...             mapping
--map2DD        Map the PARE reads to feature set
--validate      Flag to perform the validation of the potential cleave sites
...             from miRferno
--repeats       Flag to include PARE reads from repetitive regions
--noiseFilter   Flag to include all PARE validations with p-value of <=.5,
...             irrespective of the noise to signal ratio at cleave site and
...             category of PARE read.
accel           Y to use balanced multiple process scheme or else specify the
...             number of processors to be used. Y is default

======== Genome and Annotation Data ======== 
Both the GFF3 file and corresponding genome FASTA file can be downloaded from
Phytozome [http://www.phytozome.net/]

==============Examples ==================
1. Execution on new genome/entirely new dataset
This execution should be performed any time a new genome file (along with corresponding GFF file) is being analyzed:

python3 sPARTA.py -genomeFile <genomeFile.fa> -gffFile <GFF3file> -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -libs <Lib_A.txt Lib_B.txt> -tarPred -tarScore --tag2FASTA --map2DD --validate

OR

a user provided feature set (FASTA file with sequences of interest) is being analyzed:

python3 sPARTA.py -featureFile <featureFile.fa> -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -libs <Lib_A.txt Lib_B.txt> -tarPred -tarScore --tag2FASTA --map2DD --validate

2. Execution on genome in which genome has already been processed
This execution should be performed if a genome file has been processed previously but the miRNAs for which targets need to be predicted are new:

python3 sPARTA.py -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -libs <Lib_A.txt Lib_B.txt> -tarPred -tarScore --tag2FASTA --map2DD --validate

3. Execution on data in which genome and miRNA files have been previously processed
This execution should be performed if targets for a genome file have already been predicted using a miRNA file, but new PARE libraries need to be used for validation of earlier predicted targets:

python3 sPARTA.py -genomeFeature <0/1> -libs <Lib_C.txt Lib_D.txt> --map2DD --validate

4. Execution of 'miRferno', just for target prediction
This execution should be performed in case only predicted targets are required or PARE libraries are not available:
python3 sPARTA.py -genomeFile <genomeFile.fa> -gffFile <GFF3file> -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -tarPred -tarScore

OR

a user provided feature set (FASTA file with sequences of interest) is being analyzed:

python3 sPARTA.py -featureFile <featureFile.fa> -genomeFeature <0/1> -miRNAFile <miRNAFile.fa> -tarPred -tarScore

======== Output ==========
1.  PARE validation results for each library can be found in 'output' folder
    under its corresponding library name. The 'output' folder also contains a combined result file (AllLibValidatedUniq.csv) from all the libraries. 
    Results from all libs were combined by removing redundant miRNA-target interaction with cleavage at same site.

2.  Target prediction results can be found in 'predicted' folder under the name
    'All.targs.parsed.csv'

======= Contact ===========

Atul Kakrana
kakrana@udel.edu

Reza Hammond
hammond@dbi.udel.edu

===== END of README =======
