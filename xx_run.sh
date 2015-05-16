#!/bin/bash

#
# Overview
# ==========
# This script is prepared to reproduce hierarchical
# clustering of CAGE and RNA-seq data described in
# Kawaji H, 2014. Genome Res 24: 708–717
# 
# This is a supplemental material for a chaptor of
# upcomming book (in Japanese).
#
# You can run this script as below:
#     $ mkdir NEW_DIR
#     $ cd NEW_DIR
#     $ git clone ...
#     $ sh xx_run.sh
# on the computer where the software below is installed 
#
#
#
# requirements 
# ============
# * general utilities
# - lftp (http://lftp.yar.ru/)
# - wget (http://www.gnu.org/software/wget/)
# - Ruby/RubyOnRails/activesupport (https://github.com/rails/rails/tree/master/activesupport)
# - jq (http://stedolan.github.io/jq/)
#
# * utilities for biology
# - bedtools (https://github.com/arq5x/bedtools2)
# - jksrc (http://hgdownload.soe.ucsc.edu/admin/exe/)
# - R/bioconductor (and edgeR)
#
#


obtain_TSS_counts_CAGE() {
    lftp -c "open http://fantom.gsc.riken.jp/5/suppl/Kawaji_et_al_2013/data/HeliScopeCAGE/;mget *rep0[1-3]*" 
    lftp -c "open http://fantom.gsc.riken.jp/5/suppl/Kawaji_et_al_2013/data/IlluminaCAGE/; mget *rep0[1-3]*"
}

obtain_RefSeq() {
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
}

prepare_promoter_and_exon() {
  obtain_RefSeq

    # prep promoter regions
    gunzip -c refFlat.txt.gz \
    | awk 'BEGIN{OFS="\t"}{$2=$1","$2","$3","$4","$5","$6;print}' \
    | cut -f 2-11 \
    | genePredToBed stdin stdout \
    | awk 'BEGIN{OFS="\t"}{if ($6 == "+"){$2=$2-500}else{$2=$3-500};$3=$2+1000;print $1,$2,$3,$4,$5,$6}' \
    | sort -k1,1 -k2,2n \
    | uniq \
    | gzip -c > refFlat.promoter.bed.gz

    # prep exon regions
    gunzip -c refFlat.txt.gz \
    | awk 'BEGIN{OFS="\t"}{$2=$1","$2","$3","$4","$5","$6;print}' \
    | cut -f 2-11 \
    | genePredToBed stdin stdout \
    | bed12ToBed6 -i stdin \
    | sort -k1,1 -k2,2n \
    | uniq \
    | gzip -c > refFlat.exon.bed.gz
}


count_promoter_read ()
{
    infileFwd=$1
    outfile=$2

    # set reverse strand file name
    infileRev=$(echo $1 | sed -e s'/fwd.bw/rev.bw/')

    # get total number of reads
    sum=$( echo "$(bigWigToBedGraph $infileFwd /dev/stdout \
           | awk '{ sum = sum + $4 * ($3 - $2) }END{print sum}' ) \
           +  \
           $(bigWigToBedGraph $infileRev /dev/stdout \
           | awk '{ sum = sum + $4 * ($3 - $2) }END{print sum}' ) " | bc )
    echo "01STAT:MAPPED\t$sum" > .tmp$$

    # count on foward strand
    gunzip -c refFlat.promoter.bed.gz \
    | awk '{if($6 == "+"){print}}'\
    | bigWigAverageOverBed $infileFwd /dev/stdin /dev/stdout \
    | cut -f 1,4 \
    >> .tmp$$

    # count on reverse  strand
    gunzip -c refFlat.promoter.bed.gz \
    | awk '{if($6 == "-"){print}}'\
    | bigWigAverageOverBed $infileRev /dev/stdin /dev/stdout \
    | cut -f 1,4 \
    >> .tmp$$

    # post processing
    sort .tmp$$ | gzip -c > $outfile
    rm -f .tmp$$
}


obtain_RNAseq_alignment ()
{
    # get metadata of DRA records
    wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA001/DRA001100/DRA001100.analysis.xml

    # download aligned files.
    # DRA XML file is convered into JSON, and 
    # processed by 'jq' to get filenames to be downloaded
    for X in $( cat DRA001100.analysis.xml  \
      | ruby -ractive_support -ractive_support/core_ext -e "puts Hash.from_xml(  STDIN.read ).to_json" \
      | jq -r  ' .ANALYSIS_SET.ANALYSIS[] | [.TITLE, .accession, .DATA_BLOCK.FILES.FILE.filename ] | @csv ' \
      | grep "R[A-F]" \
      | grep "rep0[1-3]" \
      | grep -e "RNig" \
      | sed -e 's/"//g' )
    do
      local_file=$(echo $X | cut -f 1 -d ',').bam
      repo_acc=$(echo $X | cut -f 2 -d ',')
      repo_file=$(echo $X | cut -f 3 -d ',')
      url=ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA001/DRA001100/${repo_acc}/provisional/${repo_file}
      wget -O $local_file $url
    done
}


count_exon_read ()
{
    infile=$1
    outfile=$2

    # set library size as total counts of mapped reads
    sum=$(samtools view -q 20 ${infile}  | wc -l)
    echo "01STAT:MAPPED\t$sum" > .tmp$$

    # counts in each exon (htseq can be used nowadays)
    samtools view -uq 20 ${infile}  \
    | bamToBed -i stdin \
    | annotateBed -s -counts -i refFlat.exon.bed.gz -files stdin \
    | cut -f 4,7 \
    | sort \
    | groupBy -grp 1  -c 2 -o sum \
    >> .tmp$$

    gzip -c .tmp$$ > $outfile
    rm -f .tmp$$
}


prep_dgeLset ()
{
  outfile=$1

  cat << 'EOF' | R --slave

    library(edgeR)
    # --------------
    #   subroutine
    # --------------

    # get read count of individual files, and set as a single table
    read_count_into_table <- function(infiles)
    {
      sapply(
        infiles,
        function(infile)
        {
          tmp1 = read.table(infile,sep="\t",as.is=T,nrow=-1)
          tmp2 = tmp1[,2]
          names(tmp2) = tmp1[,1]
          tmp2
        }
      )
    }

    # --------
    #   main
    # --------

    # specify data files for individual platforms
    infile_pattern = c()
    infile_pattern["IlluminaCAGE"] = 'R.*J.*-DA.*count.txt.gz'
    infile_pattern["HeliScopeCAGE"] = 'R.*CNhs.*count.txt.gz'
    infile_pattern["RNAseq"] = 'R.*RNig.*count.txt.gz'

    # get the read counts
    dgeLset = lapply(
      names(infile_pattern),
      function(n)
      {
        # read counts
        infiles = dir( path=".",  pattern=infile_pattern[n] , full.names=T )
        counts = read_count_into_table(infiles)

        # define groups, based on the file names.
        group =  colnames(counts)
        group = sub("/","", group)
        group = sapply( group , function(str) strsplit(str, "\\.")[[1]][2] )
        group[ group == "RA" ] = "100:0"
        group[ group == "RB" ] = "99:1"
        group[ group == "RC" ] = "95:5"
        group[ group == "RD" ] = "90:10"
        group[ group == "RE" ] = "50:50"
        group[ group == "RF" ] = "0:100"

        # set the 1st row as library sizes, and the remaining rows as read counts
        lib.size = counts[1,]
        counts = counts[-1,]
        dgeL = DGEList(counts = counts, lib.size = lib.size, group = group)
        dgeL = calcNormFactors(dgeL, method="RLE")
        dgeL
      }
    )
    names(dgeLset) = names(infile_pattern)
    save(dgeLset, file=".tmp.RData")
EOF
  mv -f .tmp.RData $outfile
}

plotClust ()
{
  infile=$1
  outfile=$2

    mv -f $infile .tmp.RData
    cat << 'EOF' | R --slave

    library(edgeR)
    # -------------
    #   subroutine
    # -------------
    dgeLdendrogram <- function(dgeL)
    {
      log2cpm = cpm(dgeL, log=T)
      grp = dgeL$samples$group
      grp.u = unique(grp)

      idx = which( apply( log2cpm , 1, max) > log2(10) )
      log2cpm = log2cpm[idx,]

      log2cpm.grp = sapply(
        grp.u,
        function(g)
        {
          idx = which( grp == g )
          rowMeans( log2cpm[,idx] )
        }
      )
      colnames( log2cpm.grp  ) = grp.u

      dist = as.dist(1 - cor(log2cpm.grp, method="spearman"))
       hclust(dist,method="average")
    }

    # --------
    #   main
    # --------
    infile = ".tmp.RData"
    load(infile)
    pdf(".tmp.pdf", width=15,height=5)
    layout(t(c(1,2,3)))
    for ( n in names( dgeLset ) )
    {
      hcl = dgeLdendrogram(dgeLset[[n]])
      plot(hcl, main=n)
    }
    dev.off()
EOF
    mv -f .tmp.RData $infile
    mv -f .tmp.pdf $outfile
}


###  <main> ###
if false ; then

  # prepareation of reference data
  prepare_promoter_and_exon

  # download data
  obtain_RefSeq
  obtain_TSS_counts_CAGE
  obtain_RNAseq_alignment

  # counting - CAGE
  for X in *.ctss.fwd.bw
  do
    count_promoter_read $X $(basename $X .ctss.fwd.bw).count.txt.gz
  done

  # counting - RNAseq
  for X in *.align.bam
  do
    count_exon_read $X $(basename $X .bam).count.txt.gz

  done

  # store the data as DGEList object of edgeR (R/bioconductor)
  prep_dgeLset xxA_dgeLset.RData

  # plot
  plotClust xxA_dgeLset.RData xx_plot.pdf
  rm -f xxA_dgeLset.RData

fi
