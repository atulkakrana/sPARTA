## miRferno - miRNA-Target prediction module of sPARTA
## Updated: version-1.12 6/30/2015
## Property of Meyers Lab at University of Delaware
## Author: kakrana@udel.edu

#### PYTHON FUNCTIONS ##############################
import argparse
import sys,os,re,time,glob
import subprocess, multiprocessing
import shutil
import operator
from multiprocessing import Process, Queue, Pool
from operator import itemgetter
import pyfasta
import datetime

#### USER SETTINGS ################################
parser = argparse.ArgumentParser()
parser.add_argument('-gffFile',  default='', help='GFF file for the species '\
    'being analyzed corresponding to the genome assembly being used')
parser.add_argument('-genomeFile', default='', help='Genome file in FASTA '\
    'format')
parser.add_argument('-featureFile', default='', help='Feature file in FASTA '\
    'format')
parser.add_argument('-genomeFeature', required=True, help='0 if prediction is '\
    'to be done in genic region. 1 if prediction is to be done in intergenic '\
    'region')
parser.add_argument('-miRNAFile', default='', help='FASTA format of miRNA '\
    'sequences')
parser.add_argument('-tarPred', nargs='?', const='H', help='Mode of target '\
    'prediction. H for heuristic. E for exhaustive. H is default if no mode '\
    'is specified')
parser.add_argument('-tarScore', nargs='?', const='S', help='Scoring mode '\
    'for target prediction. S for seedless. N for normal. S is default if '\
    'no mode is specified')
parser.add_argument('-tagLen', default=20, help='Minimum length of PARE tag, '\
    'tags longer than tagLen will be chopped to the specified length. 20 is '\
    'default')
parser.add_argument('--repeats', action='store_false', default=True, help=
    'Flag to include PARE reads from repetitive regions')
parser.add_argument('--noiseFilter', action='store_false', default=True,
    help='Flag to include all PARE validations with p-value of <=.5, '\
    'irrespective of the noise to signal ratio at cleave site and category '\
    'of PARE read')
parser.add_argument('-accel', default='Y', help='Y to use '\
    'balanced multiple process scheme or else specify the number of '\
    'processors to be used. Y is default')
parser.add_argument('-misMat', default=5, help='Set the number of '\
    'mismatches allowed for predicting putative targets . 5 is default'\
    'with sum of mismatches and wobbles <= 6')
parser.add_argument('-wob', default=5, help='Set the number of '\
    'wobbles allowed for predicting putative targets . 5 is default'\
    'with sum of mismatches and wobbles <= 6')
parser.add_argument('-tarBul', default=2, help='Set the number bulged '\
    'nuleotides allowed in target for predicting putative targets . 2 is default,'\
    'i.e. max gaps allowed in miRNA')
parser.add_argument('-mirBul', default=1, help='Set the number bulged '\
    'nucleotides allowed in miRNA for predicting putative targets, '\
    'i.e. max gaps allowed in target. 1 is default')

### Developer Options ###
parser.add_argument('--generateFasta', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('--fileFrag', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('--indexStep', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('-splitCutoff', default=20, help=argparse.SUPPRESS)
parser.add_argument('-maxHits', default=30, help=argparse.SUPPRESS)

args = parser.parse_args()

## Various checks for dependencies within command line arguments
# If either gff or genome file is given without the other and featureFile
# is not given, exit.
if(((args.gffFile and not args.genomeFile) or (args.genomeFile and not
        args.gffFile)) and (not args.featureFile)):
    print("gffFile and genomeFile must both be given to create feature set")
    exit()

# If the user input both a genome and feature file, exit as both cannot be
# supplied for proper execution
if(args.genomeFile and args.featureFile):
    print("genomeFile and featureFile cannot both be supplied for execution")
    exit()

# If gffFile and genomeFile are given turn on extraction, frag and index steps
# must be set on
if(args.gffFile and args.genomeFile):
    args.generateFasta = True
    args.fileFrag = True
    args.indexStep = True

# If featureFile is given, frag and index steps must be set on
if(args.featureFile):
    # If featureFile is given and gffFile is given, give a warning letting
    # user know the gffFile will be ignored and the input fasta file may
    # have been intended as a genomeFile
    if(args.gffFile):
        print("Warning: You have input a gffFile but input a FASTA file as "\
        "the featureFile. If you intended for this to be used in conjunction "\
        "with the gff file to create a feature file, please press 'ctrl+c' "\
        "to cancel the execution and rerun with the FASTA file under the "\
        "argument 'genomeFile'. If this is in fact the feature file, allow "\
        "sPARTA to continue its execution.")
        time.sleep(10)
    args.fileFrag = True
    args.indexStep = True

# If indexStep is on and tarPred is off, turn tarPred and tarScore on
if(args.indexStep):
    if(not args.tarPred):
        args.tarPred = 'H'
    if(not args.tarScore):
        args.tarScore = 'S'

# If tarPred is on, then tarScore will default to S
if(args.tarPred and not args.tarScore):
    args.tarScore = 'S'

# If tarPred is on, then miRNAFile must be provided
if(args.tarPred and not args.miRNAFile):
    print("miRNA file must be given to perform target prediction")
    sys.exit()

# genomeFeature must be an integer
args.genomeFeature = int(args.genomeFeature)

#################################################################################################
#### SPARTA FUNCTIONS ###########################

#### Extract coordinates from GFF file
def extractFeatures(genomeFile,gffFile):
    
    fh_in = open(gffFile,'r')
    fh_in.readline() ## GFF version
    gffRead = fh_in.readlines()
    genome_info = [] #
    for i in gffRead:
        ent = i.strip('\n').split('\t')
        #print (ent)
        if i.startswith("#"): ## Fixes Bug: 010115
            pass
        else:
            if ent[2] == 'gene':
                chrID = ent[0]
                strand = ent[6].translate(str.maketrans("+-","WC"))
                geneName = ent[8].split(';')[0].split('=')[1]
                geneType = 'gene'
                #print(chrID,strand,geneName,ent[3],ent[4],geneType)
                genome_info.append((chrID,strand,geneName,int(ent[3]),int(ent[4]),geneType))
        
    #genome_info_inter = genome_info ##
    genome_info_sorted = sorted(genome_info, key=operator.itemgetter(0,1,3)) ## CHR_ID, strand and start 
    genome_info.sort(key=lambda k:(k[0],k[1],k[3])) #
    genome_info_inter = genome_info
    alist = []##
    for i in range(0, int(len(genome_info))+1): #
        #print (i)
        gene1 = (genome_info[i])
        gene2 = (genome_info[i+1])
        gene_type = 'inter' #
        #print(gene1,gene2)
        
        if gene1[3] == gene2[3] and gene1[4] == gene2[4]:
            ##gene is same/overlapping consider next gene
            pass
        else:
            if tuple(gene1[0:2]) not in alist:#
                print ('Caching gene coords for chromosome: %s and strand: %s\n' % (gene1[0], gene1[1]))
                alist.append((gene1[0:2]))
                inter_start1 = 1
                inter_end1 = int(gene1[3])-1#
                
                ## If both the genes belong to same chr and strand i.e. chromosome has atleast two genes
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]: 
                    inter_start2 = int(gene1[4])+1##From end of first gene of chromosome
                    inter_end2 = int(gene2[3])-1###Till start of second gene
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene2[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])

                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))
                    genome_info_inter.append((gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type))

                ## If the first two genes are not from same strand i.e. there is only one gene on this chromosme and strand
                else: ## Intergenic from end of chromosome/scaffold
                    inter_start2 = int(gene1[4])+1##From end of first gene of chromosome
                    inter_end2 = '-' ### Till end of chromosome
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_down' % gene1[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])
                    #if inter_name1 == inter_name2:
                    
                    #print ('\nLoop2 - First gene on this chromosme and strand but also the only one')
                    #print (gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type)
                    #print (gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type)
                    
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))
                    genome_info_inter.append((gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type))

            else:
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]: ### If chr_id and strands are equal than find intergenic.
                    inter_start = int(gene1[4])+1#
                    inter_end = int(gene2[3])-1 #
                    if gene2[1] == 'w': #
                        inter_name = ('%s_up' % (gene2[2]))
                    else:## reverse strand
                        inter_name = ('%s_up' % (gene1[2]))
                    genome_info_inter.append((gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type))
                
                else: #
                    inter_start = int(gene1[4])+1#
                    inter_end = '-' ## Different from MPPP as no table ro query for end of chromosome in public version
                    if gene1[1] == 'w':
                        inter_name = ('%s_down' % (gene1[2]))
                    else: 
                        inter_name = ('%s_up' % (gene1[2]))                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)) ##Chr_id, strand
    
    genome_info_inter_sort = sorted(genome_info_inter, key=operator.itemgetter(0,1))
        
    gene_coords_file = './coords'
    coords_out = open(gene_coords_file, 'w')
    coords = []        
    
    if args.genomeFeature == 2: ## Both gene and inter
        for ent in genome_info_inter_sort:
            print(ent)
            if ent[4] == '-': ## End of chromosome
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                
            elif int(ent[4])-int(ent[3]) > 25:
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
    
    else:
        if args.genomeFeature == 0:
            genomeFilter = 'gene'
        elif args.genomeFeature == 1:
            genomeFilter = 'inter'
        for ent in genome_info_inter_sort:
            if (ent[5] == genomeFilter):
                #print(ent)
                if ent[4] == '-': ## End of chromosome
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                    
                elif int(ent[4])-int(ent[3]) > 25:
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
    
    print ("Number of coords in 'coords' list: %s" % (len(coords)))
    coords_out.close()    
    fh_in.close()
    
    return coords

def getFASTA1(genomeFile,coords):
    fastaOut = './genomic_seq.fa'
    fh_out = open(fastaOut, 'w')

    if os.path.exists("./%s" % genomeFile):
        fh_in = open(genomeFile, 'r')
    else:
        print("\nPlease check the genomeFile - Is it in sPARTA directory? Did you input wrong name?")
        print("If it exists and you input correct name then it could have been deleted while last cleanup operation - Sorry!!")
        print("If deleted then please rename genomeFile (without any integers in name) and re-copy in sPARTA folder\n")
        print("Script will exit for now\n")
        sys.exit()

    genomeFile = fh_in.read()
    genomeList = genomeFile.split('>')[1:] 
    chromoDict = {} 
    for i in genomeList:
        chromoInfo = i.partition('\n') 
        chrid = chromoInfo[0].split()[0] 
        chrSeq = chromoInfo[2].replace("\n", "")
        chromoDict[chrid] = [chrSeq]

    chromo_mem = []
    for i in coords: ## Coords is list from gff parser
        #print (i)
        gene = i[2]
        chr_id = i[0]
        strand = i[1]
        start = i[3]-1#
        end = i[4]
        #print('start:%s End:%s Chr:%s Strand:%s' % (start,end,chr_id,strand))
        
        if tuple(i[0:2]) not in chromo_mem: 
            chromo_mem.append(tuple(i[0:2]))   ##
            print ("Reading chromosome:%s and strand: '%s' into memory to splice genes" % (i[0],i[1]) )
            chrKey = i[0]
            chromo = str(chromoDict[chrKey])
            #print('Chromosome:',chromo)
            
            gene_seq = chromo[start:end]  #
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATCG"))
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
    
        elif end == '-': ##
            print('Fetching gene %s to prepare FASTA file' % (gene))
            gene_seq = chromo[start:]##
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATCG"))
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
            
        else:
            print('Fetching gene %s to prepare FASTA file' % (gene))
            gene_seq = chromo[start:end]##
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATCG"))
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))

    time.sleep(10)
    fh_out.close()
    
    return fastaOut

def fragFASTA(FASTA):
    
    shutil.rmtree('./genome', ignore_errors=True)
    os.mkdir('./genome')
    
    pattern = ".*\.[0-9].*\.fa" ## 
    print ("\n***Purging older files***")
    for file in os.listdir():
        if re.search(pattern,file):
            print (file,'is being deleted')
            os.remove(file)
    
    statInfo = os.stat(FASTA)
    filesize =round(statInfo.st_size/1048576,2)
    print('\n**Size of FASTA file: %sMB**' % (filesize))#
    
    if filesize <= args.splitCutoff: ## 
        fls = []
        fls.append(FASTA)
        print ('No fragmentation performed for file %s' % (fls))
        
    else: #
        #
        fh_in = open(FASTA, 'r')
        seq_count = fh_in.read().count('>')
        print('**Number of headers in file: %s**\n'% (seq_count))
        #if genome == 'N':
        if seq_count >= 30: # 
            
            #
            if filesize <= 3072:
                splitnum = str(args.maxHits)
            elif filesize > 3072 and filesize <= 5120:
                splitnum = str(round(args.maxHits*1.25))
            else:
                splitnum = str(round(args.maxHits*1.5))        
            
            print ("Size based fragmentation in process for '%s' file" % (FASTA))
            retcode = subprocess.call(["pyfasta","split", "-n", splitnum, FASTA])
            fls = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)] ## file list using regex
            #fls = glob.glob(r'%s.[0-9]{1-3}.fa' % (FASTA.split('.')[0])) ## fragmented file list ##
            #print ('The fragments: %s' % (fls))
               
        
        else: ##
            splitnum = str(seq_count) ## 
            if fragFasta == 'Y':
                print ("Header based fragmentation in process for '%s' file" % (FASTA))
                retcode = subprocess.call(["pyfasta","split", "-n", splitnum, FASTA])
            fls = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)]
            #fls = glob.glob(r'%s.[0-9]{1,3}.fa' % (FASTA.split('.')[0])) ## fragmented file list
            #os.chdir("../")
            #print ('The fragments: %s' % (fls))
    
    memFile = open('frag.mem','w')
    memFile.write("fasta=%s\n" % (FASTA))
    memFile.write("genomeFeature=%s\n" % (str(args.genomeFeature)))
    memFile.write("size=%sMB\n" % (filesize))
    memFile.write("frags=%s" % (','.join(fls)))
    memFile.close()

    return fls

def miRinput():
    miRs = [] ## 
    #fh = open(args.miRNAFile)
    
    miRNA_file_clean = CleanHeader(args.miRNAFile)##
    fh_miRNA = open(miRNA_file_clean, 'r')
    fh_out2 = open('miRinput_RevComp.fa', 'w')
    mir_base = fh_miRNA.read()
    mir_blocks= mir_base.split('>')
    for i in mir_blocks[1:]:
        #print (i)
        block = i.strip('\n')##
        ent = block.split('\n')#
        #print (ent)
        #print ('%s,%s,%s,%s' % (ent[0],'None','None',ent[1]))
        miRs.append((ent[0],'None','None',ent[1]))#
        fh_out2.write('>%s\n%s\n' % (ent[0],ent[1].translate(str.maketrans("AUTGC","TAACG"))[::-1]))## make rev comp of miRNA so that it matches the target site in genome rather than mapping miRNA to genome - in target finder file make sure that miRNA is complemented again to main original seq but not direction
    fh_miRNA.close()
    mirTable = 'None'#
    print ('Total number of miRNAs in given file: %s\n' % (len(miRs)))
    
    fh_out2.close()
        
    #for i in miRs:
    #    print (i)
        
    return miRs ##

def tarFind3(frag):

    file_out = './predicted/%s.targ' % (frag.rpartition('.')[0]) ## 

    index = "./index/%s_index" % (frag) ##
    if args.indexStep:
        print('**Creating index of cDNA/genomic sequences:%s\n**' % (index))
        retcode = subprocess.call(["bowtie2-build", frag, index])

    else: 
        if os.path.isfile('%s.1.bt2' % index): ##
            retcode = 0
            print('**Found index of cDNA/genomic sequences:%s\n**' % (index))
        else:
            print('**Could not find index of cDNA/genomic sequences:%s\n**' % (index))
            sys.exit()

    if retcode == 0: ##
        print ('Predicting targets for frag:%s using index:%s' % (frag,index))
        nspread2 = str(nspread)
        if args.tarPred == 'H': 
            intervalFunc = str("L,4,0.1")
            minScoreFunc = str("L,-24,-0.5")
            readGap = str("22,14")
            refGap = str("8,14")
            ### Changed -D 5 to 6, changed -R 1 to 2 | Jan 13 -D 6 -> -D 3
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 3","-R 2","-N 1","-L 8","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--min-score",minScoreFunc,"--norc","--no-unal","--no-hd","-p",nspread2, "-x", index, "-f" ,"miRinput_RevComp.fa","-S", file_out])
        elif args.tarPred == 'E': 
            print ("You chose 'Exhaustive mode' for target identification - Please be patient")
            intervalFunc = str("L,2,0.1")
            minScoreFunc = str("L,-24,-0.5")
            readGap = str("22,14")
            refGap = str("8,14")
            #### Jan 13 -D 7 -> -D 4 |
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 4","-R 2","-N 1","-L 6","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nspread2, "-x", index, "-f", "miRinput_RevComp.fa","-S", file_out])

        else:
            print ('''\nPlease choose correct target prediction mode - Heuristic (H) or Exhaustive (E)\n
                   System will exit now''')
            sys.exit()

    ##
    if retcode2 == 0:#
                    print('\n miRNAs mapped to Fragment: %s' % (frag))
    else:
        print ("There is some problem with miRNA mapping '%s' to cDNA/genomic seq index" % (frag))
        print ("Script exiting.......")
        sys.exit()

## New version added - Apr1/15
def tarFind4(frag):
    
    file_out = './predicted/%s.targ' % (frag.rpartition('.')[0]) ## Result File
    
    ### Make or select index
    index = "./index/%s_index" % (frag)
    if args.indexStep:
        print('**Creating index of cDNA/genomic sequences:%s\n**' % (index))
        retcode = subprocess.call(["bowtie2-build", frag, index])

    else: #
        if os.path.isfile('%s.1.bt2' % index): #
            retcode = 0
            print('**Found index of cDNA/genomic sequences:%s\n**' % (index))
        else:
            print('**Could not find index of cDNA/genomic sequences:%s\n**' % (index))
            sys.exit()
            
    if retcode == 0: ### Index creation sucessful or index already exists
        print ('Predicting targets for frag:%s using index:%s' % (frag,index))
        nspread2 = str(nspread)
        if args.tarPred == 'H': ## Heurustic
            print ("You chose 'Heuristic mode' for target identification")
            intervalFunc = str("L,4,0.1")
            minScoreFunc = str("L,-24,-0.5") ###~34.5 - 35 i.e 34 - Stable v1.08 
            refGap= str("10,8") ## Bulge in miR  + 3MM, 2 seprate or consequite bulges filtered out later - - updated-Mar-23-15 - Can be (10,12) to gain speed - it will effect 1 gap and MM but make 2 gaps impossible
            readGap = str("10,6") ## Bulge in tar + 3MM, 2bulge in tar + 2MM and no bulge 6MM updated-Mar-23-15,
            misPen = str("6,2") ## Mismatch penalty, w/o phred score mx is used i.e. first one - can be changed to (5,2) to improve sensistivity - It will only effect 1 gap scenario and 2 gaps cases will be as is
            matScore = str("0") ## Match score
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 3","-R 2","-N 1","-L 8","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--mp",misPen,"--ma",matScore,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nspread2, "-x", index, "-f" ,"miRinput_RevComp.fa","-S", file_out])
        
        elif args.tarPred == 'E': ##Exhaustive
            print ("You chose 'Exhaustive mode' for target identification - Please be patient")
            intervalFunc = str("L,2,0.1")
            minScoreFunc = str("L,-24,-0.5") ### ~34.5 - 35 i.e 34 - Stable v1.08                      
            refGap= str("10,8") ## Bulge in miR  + 3MM, 2 seprate or consequite bulges filtered out later - updated-Mar-23-15
            readGap = str("10,4") ## Bulge in tar + 4MM, 2bulge in tar + 3MM and no bulge 6MM - updated-Mar-23-15,
            misPen = str("5,2") ## Mismatch penalty, w/o phred score mx is used i.e. first one
            matScore = str("0") ## Match score
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 4","-R 2","-N 1","-L 6","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--mp",misPen,"--ma",matScore,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nspread2, "-x", index,"-f", "miRinput_RevComp.fa","-S", file_out])
        
        else:
            print ('''\nPlease choose correct target prediction mode - Heuristic (H) or Exhaustive (E)\n
                   miRferno will exit now''')
            sys.exit()
    else:
        print("There seems to be a problem with index generation of locating them - miRferno will exit now")
        print("Try reruning the analysis with all steps to generate fresh files")
        sys.exit()
    
    ### Check for proper completion of Target prediction
    if retcode2 == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n miRNAs mapped to Fragment: %s' % (frag))
    else:
        print ("There is some problem with miRNA mapping '%s' to cDNA/genomic seq index" % (frag))
        print ("Script exiting.......")
        sys.exit()

def tarParse3(targComb):
    
    print ('\n**Target prediction results are being generated**')
    ## 
    print ("File for parsing: '%s' in predicted folder\n" % (targComb))
    fh_in = open(targComb,'r')
    TarPred =  './predicted/%s.parsed.csv' % (targComb.rpartition('/')[-1]) ### Similar to parsed target finder format
    fh_out = open(TarPred,'w')
    fh_out.write('miRname,Target,BindSite,miRseq,tarSeq,Score,Mismatch,CIGAR\n')
    
    ###
    acount = 0 ##
    parseCount = 0 #
    for i in fh_in:
        #print(i)
        acount += 1
        ent = i.strip('\n').split('\t')
        #print('\n%s\n' % ent)
        miRrevcomp = ent[9] ### 
        miRrev = miRrevcomp.translate(str.maketrans("TACG","AUGC")) ##      
        tarHash = list(miRrevcomp) ## 
        #tar = miRrev
        #print('Original read mapped i.e miRNA revcomp',miRrevcomp)
        
        #
        gapinfo = ent[5]
        gappos = re.split("[A-Z]",gapinfo) #
        gapNuc = re.findall("[A-Z]",gapinfo)
        posCount = 0
        for x,y in zip(gappos[:-1],gapNuc):##
            #print(x,y)
            if y == 'I':
                #tarHash.insert(posCount,'-') ##
                tarHash[posCount] = '-' #
                posCount += int(x)
            else:
                posCount += int(x)       
        #print('Target seq after gap manipulation: %s' % (''.join(tarHash)))
        #print('This is the mirna in complement',miRrev)
        
        
        misinfo = ent[-2].split(':')[-1] ##
        #print ('This is the mismatch info block:%s' % (misinfo))
        mispos = re.split("[A,T,G,C,N]",misinfo) #
        misposCorrect = [int(x)+1 for x in mispos] #
        misNuc = re.findall("[A,T,G,C,N]",misinfo) #
        posCount = 0
        for x,y in zip(misposCorrect,misNuc):
            #print(x,y)
            posCount += x ## 
            gaps = tarHash[:posCount-1].count('-') #
            #print ('Position of mismatch:%s' % (posCount))
            tarHash[posCount-1+gaps] = y

        tar = ''.join(tarHash).replace("T","U") ###
        bindsite = '%s-%d' % (ent[3],int(ent[3])+(len(miRrev)-1))

        
        gap = [] #
        mis = [] #
        wobble = [] #
        nt_cnt = 1 #
        
    #print('miRNA: %s\n%s' % (miRrevcomp[::-1],miRrevcomp[::-1].replace("T","U") ))

        #for x,y in zip(miRrevcomp[::-1].replace("T","U"),tar[::-1]):## 
        for x,y in zip(miRrevcomp[::-1].replace("T","U"),tar[::-1]):## 
            #print(miRrev[::-1][nt_cnt-1],x,y)##
            if x == '-' or y == '-':
                #print('gap')
                gap.append(nt_cnt)
                if y == '-':
                    nt_cnt+=1
                
            elif x == 'A' and y == 'G': #
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == 'C' and y == 'U': #
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == y:
                #print('match')
                nt_cnt+=1
            else:
                #print('mismatch')
                mis.append(nt_cnt)
                nt_cnt+=1
                
        #print('MimatchList:%s | GapList = %s | WobbleList = %s' % (mis,gap,wobble)) ## Poistion of mismatch gap and wobble

        score = 0   
        #print (mis)
        
        if args.tarScore == 'S': ## Allowed 3 MM, 2 Wob, 1 Gap
            mis2 = list(mis)
            #if set([10,11]).issubset(mis): 
            if 10 in mis and 11 in mis:
                score += 2.5
                #print('Removing 10')
                mis2.remove(10)
                #print ('Removing 11')
                mis2.remove(11) ## 
                
            for i in mis2:
                    score += 1
            for i in gap:
                score += 1.5
            for i in wobble:
                if (i+1 in mis) or (i-1 in mis): ##
                    score += 1.5
                elif (i+1) in mis and (i-1 in mis): ##
                    score += 2
                else:
                    score += 0.5
        else:
            for i in mis:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in gap:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in wobble:
                if i>= 2 and i<=13:
                    score += 1
                    #print ('Wobble pos:%s' % (i))
                else:
                    score += 0.5
                    #print ('Wobble pos:%s' % (i))
        ###################
            
        #print(ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo)## MiRname, Tarname, mirSeq,Taerseq,binding site
        fh_out.write('>%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo))
        parseCount  += 1
    
    
    print("Total number of interactions from 'miRferno':%s AND total interactions scored: %s" % (acount,parseCount))
    fh_in.close()
    fh_out.close()
    
    


    return TarPred

## New Version added - Apr1/15
def tarParse4(targComb):

    ''' Modifying this function is worst nightmare of life - Needs cleaning
    cutoffs w/o bulge or gap - 5MM + 1 wobble 
    bulge in miRNA - 1 bulge + 3MM
    bulge in reference - 1bulge+4mm or 2 consequite bulges +3MM'''
    
    print ('\n**Target prediction results are being generated**')
    ## Input / Output file ######
    print ("File for parsing: '%s' in predicted folder\n" % (targComb))
    fh_in = open(targComb,'r')
    TarPred =  './predicted/%s.parsed.csv' % (targComb.rpartition('/')[-1]) ### Similar to parsed target finder format
    fh_out = open(TarPred,'w')
    fh_out.write('miRname,Target,BindSite,miRseq,tarSeq,Score,Mismatch,CIGAR\n')
    
    #acount = 0
    #for i in fh_in:
    #    print ('Check lines:%s' % (i))
    #    acount +=1
    #print ('Total lines read:',acount)
    #sys.exit()
    
    #### Regenerate Target sequence with all features #####
    acount = 0 ##Total number of interactions from predictions
    parseCount = 0 ## Total number of interactions scores and written to result file
    for i in fh_in:
        # print("\nEntry:",i.strip("\n"))
        acount += 1
        ent = i.strip('\n').split('\t')
        #print('\n%s\n' % ent)
        miRrevcomp = ent[9]                 ## miRNA complemented and reversed to map genome using bowtie. That is target sequence if mismatches and gaps are added
        tarHash = list(miRrevcomp)          ## Strings are immutable convert to list - To rebuilt a traget seq
        # print("\nTarHash",tarHash)
        
        miRrev = miRrevcomp.translate(str.maketrans("TACG","AUGC")) ## Re-translated to get miR but still in reverse orientation - OK
        mirHash = list(miRrev)
        
        #print('Original read mapped i.e miRNA revcomp',miRrevcomp)
        
        ## Gap and Bulges (with reference to miRNA) - Identify gaps and bulges and modify miRNA read used for mapping to regenerate target as well as miRNA
        ## Add gap to target sequence  first to make miR length comparable to target
        gapinfo = ent[5]
        gappos = re.split("[A-Z]",gapinfo)   ## In python format - gap in target seq and bulge in miRNAseq
        gapNuc = re.findall("[A-Z]",gapinfo)
        # print("gappos:",gappos,"| gapNuc:",gapNuc)
        
        ###########################################################################################
        ## SECTION - A - FIND GAPS AND BULGES AND ADD INDICATORS TO MIRNA OR TARGET SEQUENCES
        ###########################################################################################

        posCount = 0
        ## At this point both mirHash and tarHash has same length and perfect complementrity as tarHash is essentially reverse complemented miRNA used for matching
        for x,y in zip(gappos[:-1],gapNuc): ## In gap pos list which is made of alphabet splitting there is always am empty value at end because string has alphabet at last
            # print(x,y)
            if y == 'I':                    ## There was an insertion in miRNA and gap in reference and bulge in miRNA
                for i in range(int(x)):     ## For as many as bulges in miRNA - like 11M 2I 11M              
                    tarHash[posCount] = '-' ## Replace existing nucleotide (from revcomp miRNA) with a gap
                posCount += int(x)          ## This only effects consequitve bulges in miR, if there are multiple insertions like 3I, then "3" needs to be added once and not in every iterneration
            
            elif y == 'D':                  ## There was a deletion in miRNA i.e. gap in miRNA and bulge in reference - In this case length of both miRNA and target will increase
                for i in range(int(x)):             ## For as many as gaps in miRNA                
                    mirHash.insert(posCount,'-')    ## When counted in python insertion will occur after posCount value -  Tested OK
                    tarHash.insert(posCount,'^')    ## Add bulge markers in target sequence as well - Tested OK
                posCount += int(x)                  ## This only effects consequitve gaps in miR, if there are multiple insertions like 3I, then "3" needs to be added once and not in every iterneration
            
            else:
                posCount += int(x)
        
        # print('Target %s seq after manipulation: %s' % (ent[2],''.join(tarHash))) ## Has '-' at gap and '^' at extra nucleotide position (i.e. gap in MiR) -OK
        # print('miRNA %s after after maipulations: %s' % (ent[0],''.join(mirHash))) ## Has '-' at gap pos -OK
        
        #########################################################################################
        ## SECTION -B - GET CORRECT POSITIONS FOR MISMATCHES,GAPS AND BULGES 
        ## AND REGENERATE TARGET BY INSERTING CORRECT NUCLETIDES AT EDITS AND BULGES
        #########################################################################################
        
        ## Mismatches - Identify mismatches and modify miRNA read used for mapping to regenerate target
        misinfoBlock = ent[-2].split(':')[-1] ## Reverse index because XS:i is optional column ## MD:Z:16C3 - these positions are from references - so if there is an insertion/bulge in miRNA i.e. gap that it should be added to these positions before editing miRNA to tar
        # print ('This is the mismatch info block:%s' % (misinfoBlock))
        
        ## Deletion (gap) in miRNAS i.e y= D - which has been added but replace the bulge '^' in target with actual sequence
        ## Should work if there is dletion in miRNA and deletion in traget i.e. two bulges one in miRNA and one in target - NO MM possible
        if '^' in misinfoBlock: 
            misinfo = misinfoBlock.replace('^','')          ## 11^A13 - Here miRNA was 24nt but misinfo shows 25nt as 11+A+13 - Replace the '^' inserted in target with 'A'
            mispos = re.split("[A,T,G,C,N]",misinfo)        ## Found N in one case so included, N confimed in sequence too, will be counted as mismatch
            # print('Mismatch info:%s | Mismatch pos:%s'%(misinfo,mispos))       

            ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
            misposCorrect = []                              ## List hold corrected positions, because edit is a nucleotide next to integers in misPos
            # posIndex = -1                                   ## Index of position for misposCorrect, -1 because after first addition to list it will be incremented to 0
            for x in mispos:
                # print(x)
                ## Assumption: There will be no empty entry in mispos at the begining, because there will be a position to indicate comsequtive edits, others are handled here
                if x:
                    misposCorrect.append(int(x)+1)          ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
                    # posIndex +=1
                else: ## If 'x' is empty like in case of three consequentive gaps in miRNA - Mismatch info:11GTA8 | Mismatch pos:['11', '', '', '8']
                    # y = misposCorrect[posIndex]             ## Get the last corrected position, add one to get position for empty entry
                    # misposCorrect.append(int(y)+1)          ## In consequitve edits this is position just next to last one
                    # posIndex +=1
                    misposCorrect.append(int(0)+1)

            misNuc = re.findall("[A,T,G,C,N]",misinfo)      ## Found N in one case so included, N confirmed in sequence too, will be counted as mismatch
            # print('Misafter:',mispos,'Mispos', misposCorrect,'MisNuc',misNuc)

            ## Replace the nucleotides at bulges(^) and mismatches in tarHash to give actual targets
            ## And also add nucleotide to target (replace ^) if gap in miRNA
            #print('Unedited target:%s-%s' % (''.join(tarHash),len(''.join(tarHash))))
            posCount = 0
            for x,y in zip(misposCorrect,misNuc):
                posCount += x                               ## Convert bowtie positions to python format
                # print(x,y,posCount)
                ## Account for gap before replacing the nucleotide with that in target
                gaps = tarHash[:posCount-1-1].count('-')      ## -1 to convert to python, -1 because - count at positions in target before the current mismatch/bulge position
                tarHash[posCount-1+gaps] = y                  ## Replaced the bulge or mimatch with nucleotide in target - OK
        
        else:   ## Normal i.e y = I - With insertion(bulge) in miRNA and gap in target - OK - What id there is a bulge in miRNA???
                ## In this case miRNA already had inserted nucleltides and target has been added '-' in section-A. Just replace mimatches at correct postions of target
            
            misinfo = misinfoBlock
            mispos = re.split("[A,T,G,C,N]",misinfo)        ## Found N in one case so included, N confimed in sequence too, will be counted as mismatch   
            misposCorrect = [int(x)+1 for x in mispos]      ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
            misNuc = re.findall("[A,T,G,C,N]",misinfo)      ## Found N in one case so included, N confirmed in sequence too, will be counted as mismatch
            # print('Misafter:',mispos,'Mispos', misposCorrect,'MisNuc',misNuc)
            
            ## Count for gaps, since they are added by us and replace MM nucleotides to give actual targets
            posCount = 0
            for x,y in zip(misposCorrect,misNuc):
                posCount += x                                 ## Keep adding the positions, as these are cumulative - MD:Z:2T12C4 - Misafter: ['2', '12', '4'] Mispos [3, 13, 5] MisNuc ['T', 'C'] - Tested OK
                # print(x,y,posCount)
                ## Account for gaps before replacing the nucleotide with that in target
                ## Can give problem if more than one gap? But more than one gap not allowed V07 modification?
                gaps = tarHash[:posCount-1-1].count('-')      ## -1 to convert to python,-1 because - count at positions in target before the current mismatch/bulge position
                tarHash[posCount-1+gaps] = y                  ## TESTED - OK

        tar = ''.join(tarHash).replace("T","U")             ## Target converted to RNA format, will help in catching wobbles ahead
        mir = ''.join(mirHash)
        bindsite = '%s-%d' % (ent[3],int(ent[3])+(len(miRrev)-1))
        
        # print ("Target:%s-%s | miRNA:%s-%s" % (tar,len(tar),mir,len(mir)))
        
        ###################################################################################
        ## SECTION-C - GET POSITIONAL INFORMATION ON GAPS, MISMATCHES, WOBBLES AND BULGES 
        ## AND CHOOSE VALID INTERACTIONS 
        ###################################################################################

        mirGap = [] ## List to hold gaps in miRNA
        tarGap = [] ## List  to hold gap in targets
        mis = []    ## List to hold mismatch position
        wobble = [] ## List to hold Wobble pos
        # print('miRNA: %s\n%s' % (miRrevcomp[::-1],miRrevcomp[::-1].replace("T","U") ))

        ## Read from mapping file -> miRrevComp -> uncomplement -> miRrev -> morhash-> miR 
        valid = 1 ## Validity flag [0] - Invalid and [1] - valid, if more then 1 bulges in miR or 2bulges in tar or bulges in 9,10,11 = invalid
        nt_cnt = 1 ## Keep track of actual position,
        for x,y in zip(mir.translate(str.maketrans("AUGC","UACG",))[::-1],tar[::-1]): ## Orientation changed to read from 5' miRNA - OK
            
            #print(miRrev[::-1][nt_cnt-1],x,y)## Print miRNA, rev complemmnetry miRNA used for matching, target
            # if x == '-' or y == '-' or x == '^' or y == '^':
            #     #print('gap')
            #     gap.append(nt_cnt)
            #     if y == '-':
            #         nt_cnt+=1

            if x == '-' or x == '^': ## Don't think '^' would be here, since it has been replaced with nucleotide
                #print('miRNA gap')
                mirGap.append(nt_cnt)
                if nt_cnt == 9 or nt_cnt == 10 or nt_cnt == 11:
                    # print("@miRNA has gap/bulge in 9th and 10th position")
                    valid = 0 ## This is invalid as it has bulge in miRNA at 9th and 10th pos, this interaction will be skipped

            elif y == '-' or y == '^':
                #print('target gap')
                tarGap.append(nt_cnt)
                if y == '-': ## Don't think '^' would be here, since it has been replaced with nucleotide
                    nt_cnt+=1
                    if nt_cnt == 9 or nt_cnt == 10 or nt_cnt == 11:
                        # print("@target has gap/bulge in 9th and 10th position")
                        valid = 0 ## This is invalid as it has bulge in miRNA at 9th and 10th pos, this interaction will be skipped

            elif x == 'A' and y == 'G': ### If in reference its 'G' than miRNA should have 'U' i.e. T but this is revcomplememnt of miRNA so complement of 'U' is A - Tested OK - v08 modifcation
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == 'C' and y == 'U': ### If in reference its 'U' than miRNA should have 'G' but this is rev complememnt of miRNA so complement of 'G' is C - Tested OK - v08 modification
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == y:
                #print('match')
                nt_cnt+=1
            else:
                #print('mismatch')
                mis.append(nt_cnt)
                nt_cnt+=1
        # print('MismatchList:%s | miRGapList = %s | tarGapList = %s | WobbleList = %s' % (mis, mirGap, tarGap, wobble)) ## Poistion of mismatch gap and wobble

        ## Check if there are more then permitted gaps in miRNA and target
        # if len(mirGap) > 2:
        if len(mirGap) > int(args.tarBul):
            # print("@miRNA %s has more then two gaps" % (ent[0]))
            valid = 0
        # elif len(tarGap) > 1:
        elif len(tarGap) > int(args.mirBul):
            # print("@target has more then one gap")
            valid = 0
        # elif len(mis) > 5:
        elif len(mis) > int(args.misMat):
            # print("More then 5 mismatches not allowed")
            valid = 0
        elif len(wobble) > int(args.wob):
            # print("More then 3 wobbles not allowed")
            valid = 0
        elif (len(mis) + len(wobble)) > 6: 
            # print(" More then six edits are not allowed - It's too much")
            valid = 0
        else:
            pass

        ## Decide to report this interaction if its valid
        ## Validity flag, if more then 1 bulges in miR or 2bulges in tar or bulges in 9,10,11 = invalid
        if valid == 0: 
            # print("Skipping this entry - It's biologically invalid\n")
            continue ## Go to main for loop
        else:
            ## Go for scoring
            pass

        gap = mirGap+tarGap ## Combine gaps or scoring
        
        #####################################################################################
        ## SECTION-D - SCORE THE INTERACTIONS 
        #####################################################################################

        score = 0   ## Initialize
        #print (mis)
        if args.tarScore == 'S': ## Allowed 3 MM, 2 Wob, 1 Gap
            mis2 = list(mis)
            #if set([10,11]).issubset(mis): ## Works well but took 1 sec more than below in Rice timed test
            if 10 in mis and 11 in mis: ## Check for sunsequent mismatch at 10 and 11 if yes than strict penalty ## if set(['a','b']).issubset( ['b','a','foo','bar'] )
                score += 2.5
                #print('Removing 10')
                mis2.remove(10)
                #print ('Removing 11')
                mis2.remove(11) ## So that they are not counted again
                
            for i in mis2:
                    score += 1
            for i in gap:
                score += 1.5
            for i in wobble:
                if (i+1 in mis) or (i-1 in mis): ## Mismatches around wobble - Strong penalty
                    score += 1.5
                elif (i+1) in mis and (i-1 in mis): ## Mismatches on both sides - Stronger penalty
                    score += 2
                else:
                    score += 0.5
        else:
            ##Heuristic and Exhaustive
            for i in mis:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in gap:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in wobble:
                if i>= 2 and i<=13:
                    score += 1
                    #print ('Wobble pos:%s' % (i))
                else:
                    score += 0.5
                    #print ('Wobble pos:%s' % (i))
        
        ## Correctly output mismatches - If there is no mismatch then misinfo is just the length of match like - 15G2A1T1 (if there is mismatch) and 21(if no mismatch)
        ## Since the second ouput which originally is part of mapping file confussing in mismatch column - take mismatch info from 'mis' list
        
        if mis or gap: ## If list of mismathes has positions of mismatches, then output the block from mapping file
            mismatches = misinfo
        else:
            mismatches = 0

        #print(ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo)## MiRname, Tarname, mirSeq,Taerseq,binding site
        fh_out.write('>%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[2],bindsite,mir.replace("U","T"),tar.replace("U","T"),score,mismatches,gapinfo))
        parseCount  += 1
    
    
    print("Total number of interactions from 'miRferno':%s AND total interactions scored: %s" % (acount,parseCount))
    fh_in.close()
    fh_out.close()

    return TarPred

def tag2FASTA2(lib):
    print("'%s' tag count file being converted to FASTA format" % (lib))
    fh_in = open(lib,'r') #
    fh_out = open('./PARE/%s_PARE_tags.fa' % (lib), 'w')#
    tag_num = 1 #
    for tag in fh_in:#
        #print(tag.strip('\n').split('\t'))
        ent = tag.strip('\n').split('\t')
        tag = ent[0]
        if len(tag) >= args.tagLen: #
            fh_out.write('>%s\n%s\n' % (tag_num, tag[:args.tagLen]))
            tag_num += 1
        else:
            #print ('Length is not 20nt')
            pass
    fh_out.close()

def mapdd2trans(anIndex):#
    mismatch = str(0) ##
    nspread2 = str(nspread)
    index = anIndex.rsplit('.', 2)[0]
    indexLoc = './index/%s' % index
    #for lib in libs:
    dd_file = ('./PARE/%s_PARE_tags.fa' % (templib))
    map_out = ('./dd_map/%s_%s_map' % (templib,index))
    print ('\n**The library %s is being mapped to transcriptome index file: %s**\n' % (dd_file,indexLoc))
    
    retcode2 = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nspread2, "-x", indexLoc, "-f", dd_file,"-S",map_out]) ###
    ## Optimized | no-unal was the culprit | 
    #retcode2 = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "-t","-p",nspread2, "-f", indexLoc, dd_file,"-S",map_out]) ##
    
    if retcode2 == 0:#
        print('\nDegradome from PARE lib: %s mapped to cDNA/Trascript file' % (templib)) 
    else:
        print ("There is some problem with mapping of PARE lib: %s to cDNA/genomic seq index" % (templib))
        print ("Script exiting.......")
        sys.exit()

def FileCombine():

    print('\n****************************************')
    targ_fls = [file for file in os.listdir('./predicted') if file.endswith ('.targ')]
    print ('Target files:',targ_fls)
    print ('\nCombining all the target prediction files for parsing and scoring\n')
    
    targComb = './predicted/All.targs'
    targ_out = open(targComb ,'w')
    
    for x in targ_fls:
        print (x)
        targfile = open('./predicted/%s' % (x), 'r')
        #targfile.readline()
        data = targfile.read()
        targfile.close()
        targ_out.write(data)
    
    targ_out.close()
        
    return targComb

def fileDelete():
    rm_fls = [file for file in os.listdir('./predicted') if file.endswith ('.targ')]
    rm_alltargs = [file for file in os.listdir('./predicted') if file.endswith ('.targs')]
    rm_fls.append(rm_alltargs[0])
    print ('Files for cleanup',rm_fls)
    for file in rm_fls:
        print (file)
        os.remove('./predicted/%s' % (file))

def CleanHeader(filename):
    #read file
    fh_in=open(filename, 'r')
    #write file
    out_file = ('%s_new_head.fa' % (filename))
    fh_out =open(out_file, 'w')
    
    print ('\nProcessing "%s" file to clean FASTA headers\n' % (filename))
    
    acount = 0 ## count the number of entries
    for i in fh_in:
        if re.match('>', i):
            header = i.split()#
            new_head = header[0].split('|')[0]#
            fh_out.write('%s\n' % new_head)
            acount+=1
    #        print(i)
    #        print(new_head)
        else:
            fh_out.write('%s' % i)
        
    fh_in.close()
    fh_out.close()
    return out_file

    print('The fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file, acount))

def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    nprocPP = round((args.accel/int(nspread))+1) 
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)
    
def PPmultiple(module,alist1,alist2):
    start = time.time()
    npool = Pool(int(args.accel))
    npool.map(lambda args: module(*args), alist2)

def PPResults(module,alist):
    npool = Pool(int(args.accel))    
    res = npool.map_async(module, alist)
    results = (res.get())
    npool.close()
    return results
    
def feed(queue, parlist):
    print ('Feeder function started')
    for par in parlist:
        print ('Echo from Feeder: %s' % (par))
        queue.put(par)
    print ('**Feeder finished queing**')

def calc(queueIn, queueOut):
    print ('Worker function started')
    while True:
        try:
            par = queueIn.get(block = False)
            print ('Echo from Worker \n Dealing with:', par)
            res = function(par)
            queueOut.put((par,res))
        except:
            break
    print ('**Worker finished **')

def write(queue, fname):
    print ('Writer function started')
    fhandle = open(fname, "w")
    while True:
        
        try:
            par, res = queue.get(block = False)
            print >>fhandle, par, res
        except:
            break
    fhandle.close()

#### MAIN FUNCTION ################################
def main():
    if args.generateFasta:
        coords = extractFeatures(args.genomeFile,args.gffFile) ## 
        fastaOut = getFASTA1(args.genomeFile,coords) ##
        print('This is the extracted file: %s' % (fastaOut))
    
    elif args.featureFile:
        print("\nThe input FASTA file is considered 'as is' for analysis\n")
        fastaOut = args.featureFile ### Make it better
        
    
    runLog = 'runtime_%s' % datetime.datetime.now().strftime("%m_%d_%H_%M")
    fh_run = open(runLog, 'w')
    print('tarPred: %s | tarScore: %s | Uniq filter: %s' % (args.tarPred,args.tarScore,args.repeats))
    fh_run.write('tarPred:%s | tarScore: %s | Uniq filter:%s\nGenomeFile:%s | GenomeFeature:%s' % (args.tarPred,args.tarScore,args.repeats,args.genomeFile,args.genomeFeature))
    #fh_run.write ('\nLibs: %s' % (','.join(args.libs)))
    FragStart = time.time()
    
    if args.fileFrag:
        start = time.time()##
        fragList = fragFASTA(fastaOut)##
        end = time.time()
        print ('fileFrag time: %s' % (round(end-start,2)))
    else:
        fragMem = open('frag.mem','r')
        fragRead = fragMem.readlines()
        print("Frags from earlier processed file: '%s' will be used in this run" % (fragRead[0].strip('\n').split('=')[1]))
        for i in fragRead:
            akey,aval = i.strip('\n').split('=')
            if akey == "genomeFeature":
                if int(aval) != args.genomeFeature:
                    print('\n------ERROR-------')
                    print("The earlier processed genome or feature set belongs to genomeFeature: %s" % (aval))
                    print("Your current genomeFeature input: %s" % (str(args.genomeFeature)))
                    print("Please input either correct genomeFeature value or re-run the analysis with '-featurefile' or '-genomeFile and -gffFile'")
                    print("\nSystem will exit now")
                    sys.exit(0)
                else:
                    pass
                
            if akey =="frags":
                fragList =  i.split('=')[1].split(',')
                #print('fragList from fragMem:', fragList)

        fragMem.close()       
        #fragList = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)] #deprecated
        print ('The fragments: %s' % (fragList))
        
    FragEnd = time.time()
    #print ('\n\nThe fragmentaton time is %s\n\n' % (round(FragEnd-FragStart,2)))
    fh_run.write('Fragmentation time is : %s\n' % (round(FragEnd-FragStart,2)))

    #####################################
    
    ## TARGET PREDICTION ###################
    
    
    TFStart = time.time()
    
    if args.indexStep:
        shutil.rmtree('./index', ignore_errors=True)
        os.mkdir('./index')
    if args.tarPred and args.tarScore:
        shutil.rmtree('./predicted', ignore_errors=True)
        os.mkdir('./predicted')
        print('\nFragments to be indexed and used for TP: %s' % (fragList))
        
        miRs = miRinput()
        start = time.time()###time start
        #for i in fragList:
        #    tarFind4(i)
        ## Parallel mode
        PP(tarFind4,fragList)
        end = time.time()
        print ('Target Prediction time: %s' % (round(end-start,2)))
        
        
        targComb = FileCombine()
        #targComb = 'All.targs' ## Test - open line above when real
        
        start = time.time()###time start
        predTargets = tarParse4(targComb)
        end = time.time()
    
        #print ('Target Prediction time: %s' % (round(end-start,2)))
        
    elif not args.tarPred and args.tarScore:
        targComb = FileCombine()
        #targComb = 'All.targs' ## Test - open line above when real
        
        start = time.time()###time start
        predTargets = tarParse4(targComb)
        end = time.time()
        
        #print ('Target Scoring time: %s' % (round(end-start,2)))
    
    else: ## Target prediction is OFF
        print("!!Target prediction is OFF - Files in 'predicted' folder might be old!!")
        predTargets = './predicted/All.targs.parsed.csv'
        

    fileDelete()
    
    ###Timer
    TFEnd = time.time()
    print ('\n\nTarget prediction time is %s seconds\n\n' % (round(TFEnd-TFStart,2)))
    fh_run.write('Target prediction time is : %s seconds\n' % (round(TFEnd-TFStart,2)))
    fh_run.close()

#### RUN ##########################################

if __name__ == '__main__':
    nspread = 6
    if args.accel == 'Y':
        args.accel = int(multiprocessing.cpu_count()*0.85)
    else:
        args.accel = int(args.accel)

    
    start = time.time()
    main()
    end = time.time()
    print ("The complete 'miRferno' run time is %s seconds" % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting..\n')
    sys.exit(0)

## v1.02 -> v1.03
## Reverse fix, miRNA input moved to TF function, implemented frag memory files

## v1.03 -> 1.04
## TarFind tweak

## v1.04 -> 1.06
## Version jump to match sPARTA
## Fixed compatibility to lBowtie 2.4.4 - Tested OK

##v1,06 -> v1.07
## Fix: Fixes Bug: 010115 - GFF file lines start with ##
## extractFeature bug updated from MPPP - scaffold fix
## File check and message display to rename genome file without integers

## v1.07 -> 1.10 [Major update]
## Two modules updated for improved target prediction and postprocessing to capture more secondary structures
## Sensitivity improved for both Heuristic and Exhaustive modes [But not implemented in this script - Need to be tested for sconfidence score]
## Added a filter whch throws our biologically invalid interactions

## v1.10 -> v1.11rc [stable]
## Fixed a bug introduced in v1.10 - 1) Fixed bugs with args in tarfind4 - Fixed 2) nspread got changed to nthred by mistake - Fixed

## v1.11 -> v1.12
## Fixed a bug introduced due to regresssin, effected fetching seqeunces for reverse strand