import pysam, re, sys

def getContig(x, feature):
    ret= re.sub ("transcript:","", re.sub(".t1", "", x.reference_name))
    if feature == "genes":
        ret = re.sub("_\w\d\Z","",ret)
    return(ret)

mylog = open(snakemake.log.log,'wt')
sys.stdout=mylog
sys.stderr=mylog


### main ###
# Iterate the BAM file to set the gene name as the transcriptome's entry
infile = pysam.AlignmentFile(snakemake.input.bam, "rb")
outfile = open(snakemake.output.counts, 'wt')
dbReadPass = {}
count= 0
name = None
contig = None
score = 0
nreads = 0
ignored = 0 # TODO: Capturing reasons failed reads could be made more fine-grained
for read in infile.fetch(until_eof=True):
    if name == None:
        name = read.query_name
        if read.is_unmapped or \
           read.is_reverse or \
           read.template_length < 0 or \
           read.mapping_quality < snakemake.params.minMapq or \
           (snakemake.params.multimappers == "ignore" and read.get_tag("NH") > 1):
            ignored += 1
            name = "ignore"
            continue
        else:
            contig = getContig(read, snakemake.params.feature) #read.reference_name
            score = 1.0 #/read.get_tag("NH")
            if score > 1.0:
                print("ERROR: Read pair matesScore should not be higher than 1.0")
                sys.exit(-1)
    else:
        if name != "ignore":
            if read.is_unmapped or \
               not read.is_reverse or \
               read.mapping_quality < snakemake.params.minMapq or \
               (snakemake.params.multimappers == "ignore" and read.get_tag("NH") > 1):
                ignored +=1
            else:
                if name != read.query_name:
                    print("ERROR: Read pair mates not in order")
                    sys.exit(-1)
                if  contig != getContig(read, snakemake.params.feature): 
                    print("ERROR: Read pair contig not matching")
                    sys.exit(-1)
                if read.mapping_quality > snakemake.params.minMapq:
                    if snakemake.params.multimappers != "ignore" and score != 1.0: 
                        print("ERROR: Read pair mates multimaps differently")
                        sys.exit(-1)
                    else:
                        if contig not in dbReadPass.keys():
                            dbReadPass[contig] = 0
                        dbReadPass[contig] += score
                else:
                    ignored += 1
        nreads += 1
        name = None
        contig = None
print("{} read pairs read, {} of these were ignored".format(nreads, ignored))
for contig in dbReadPass.keys():
    outfile.write("{}\t{}\n".format(contig, int(dbReadPass[contig])))
                
infile.close()
outfile.close()
mylog.close()
