#!/bin/bash

write_usage() {
  echo ""
  echo "Usage: `basename $0` <path to the target tumor bam file> <path to the target normal bam file> <path to the output directory> <path to the text file of the list for normal reference samples> [<path to the config.sh>]"
  echo ""
}

INPUTBAM_TUM=$1
INPUTBAM_NOR=$2
OUTPUTPATH=$3
REFERENCELIST=$4
config=$5

readonly DIR=`dirname ${0}`

if [ $# -le 3 -o $# -ge 6 ]; then
  echo "wrong number of arguments"
  write_usage
  exit 1
fi

if [ $# -eq 4 ]; then
  config=${DIR}/config.sh
fi

if [ $# -eq 5 ]; then
  if [ ! -f ${config} ]; then
    echo "${config} dose not exists"
    write_usage
    exit 1
  fi
fi

source ${config}
source ${DIR}/utility.sh


check_file_exists ${INPUTBAM_TUM}
check_file_exists ${INPUTBAM_NOR}
check_file_exists ${REFERENCELIST}

targetfile=${REFERENCELIST}
cat $targetfile | while read line; do
  if [ ! -f $line ]; then
    echo "$line does not exists."
    echo "Please check $targetfile"
    exit 1
  fi
done
if [ $? -ne 0 ]; then
  exit 1
fi

check_file_exists ${PATH_TO_REF}
check_mkdir ${OUTPUTPATH}/tmp

# exec > >(tee ${OUTPUTPATH}/tmp/command.log)


# pileup the tumor bam file
##########
echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_TUM} > ${OUTPUTPATH}/tmp/temp.tumor.pileup"
${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_TUM} > ${OUTPUTPATH}/tmp/temp.tumor.pileup
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_NOR} > ${OUTPUTPATH}/tmp/temp.normal.pileup"
${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_NOR} > ${OUTPUTPATH}/tmp/temp.normal.pileup
check_error $?
##########

# make count files for mismatches, insertions and deletions
# mismatch count is performed considering bases whose quality is more than ${TH_BASE_QUAL}.
##########
echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "${DIR}/subscript/Pileup2Base.o ${TH_BASE_QUAL} ${MIN_TUMOR_DEPTH} ${OUTPUTPATH}/tmp/temp.tumor.pileup ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.tumor.depth"
${DIR}/subscript/Pileup2Base.o ${TH_BASE_QUAL} ${MIN_TUMOR_DEPTH} ${OUTPUTPATH}/tmp/temp.tumor.pileup ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.tumor.depth
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "${DIR}/subscript/Pileup2Base.o ${TH_BASE_QUAL} ${MIN_NORMAL_DEPTH} ${OUTPUTPATH}/tmp/temp.normal.pileup ${OUTPUTPATH}/tmp/temp.normal.base ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.normal.depth"
${DIR}/subscript/Pileup2Base.o ${TH_BASE_QUAL} ${MIN_NORMAL_DEPTH} ${OUTPUTPATH}/tmp/temp.normal.pileup ${OUTPUTPATH}/tmp/temp.normal.base ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.normal.depth
check_error $?

# extract putative germline mutations
##########
echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/filterBase.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.base ${MIN_NORMAL_DEPTH} > ${OUTPUTPATH}/tmp/temp.normal.base.filt"
perl ${DIR}/subscript/filterBase.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.base ${MIN_NORMAL_DEPTH} > ${OUTPUTPATH}/tmp/temp.normal.base.filt
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/filterBase_del.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.del ${MIN_NORMAL_DEPTH} > ${OUTPUTPATH}/tmp/temp.normal.del.filt"
perl ${DIR}/subscript/filterBase_del.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.del ${MIN_NORMAL_DEPTH} > ${OUTPUTPATH}/tmp/temp.normal.del.filt
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/filterBase_ins.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.ins ${MIN_NORMAL_DEPTH} > ${OUTPUTPATH}/tmp/temp.normal.ins.filt"
perl ${DIR}/subscript/filterBase_ins.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.ins ${MIN_NORMAL_DEPTH} > ${OUTPUTPATH}/tmp/temp.normal.ins.filt
check_error $?


# filter candidate of variation between tumor and normal
##########
echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "${DIR}/subscript/CompBase.o ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.normal.base ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt"
${DIR}/subscript/CompBase.o ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.normal.base ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt"
perl ${DIR}/subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt"
perl ${DIR}/subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt
check_error $?


# add information of normal reference
##########
echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/getRefNor_base.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt ${REFERENCELIST} ${TH_BASE_QUAL_REF} ${TH_MAPPING_QUAL_REF} ${OUTPUTPATH}/tmp ${PATH_TO_SAMTOOLS} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref"
perl ${DIR}/subscript/getRefNor_base.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt ${REFERENCELIST} ${TH_BASE_QUAL_REF} ${TH_MAPPING_QUAL_REF} ${OUTPUTPATH}/tmp ${PATH_TO_SAMTOOLS} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt ${REFERENCELIST} 1 ${TH_MAPPING_QUAL_REF} ${OUTPUTPATH}/tmp ${PATH_TO_SAMTOOLS} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref"
perl ${DIR}/subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt ${REFERENCELIST} 1 ${TH_MAPPING_QUAL_REF} ${OUTPUTPATH}/tmp ${PATH_TO_SAMTOOLS} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt ${REFERENCELIST} 2 ${TH_MAPPING_QUAL_REF} ${OUTPUTPATH}/tmp ${PATH_TO_SAMTOOLS} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref"
perl ${DIR}/subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt ${REFERENCELIST} 2 ${TH_MAPPING_QUAL_REF} ${OUTPUTPATH}/tmp ${PATH_TO_SAMTOOLS} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref
check_error $?


# caluculate p-value for each candidate by the Beta-binominal sequencing error model 
##########
echo "`date '+%Y-%m-%d %H:%M:%S'`"
if [ ! -s ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref ]; then
    echo "make empty file : ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb"
    rm -f ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb
    touch ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb
else
    echo "${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb < ${DIR}/subscript/proc_EBcall.R"
    ${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb < ${DIR}/subscript/proc_EBcall.R
    check_error $?
fi

echo "`date '+%Y-%m-%d %H:%M:%S'`"
if [ ! -s ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref ]; then
    echo "make empty file : ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb"
    rm -f ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb
    touch ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb
else
    echo "${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb < ${DIR}/subscript/proc_EBcall.R"
    ${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb < ${DIR}/subscript/proc_EBcall.R
    check_error $?
fi

echo "`date '+%Y-%m-%d %H:%M:%S'`"
if [ ! -s ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref ]; then
    echo "make empty file : ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb"
    rm -f ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb
    touch ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb
else
    echo "${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb < ${DIR}/subscript/proc_EBcall.R"
    ${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb < ${DIR}/subscript/proc_EBcall.R
    check_error $?
fi

# merge and convert the format of three variation files (mutation, insertion and deletion) for Annovar
##########
echo "`date '+%Y-%m-%d %H:%M:%S'`"
echo "perl ${DIR}/subscript/procForAnnovar.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb > ${OUTPUTPATH}/output.txt"
perl ${DIR}/subscript/procForAnnovar.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb > ${OUTPUTPATH}/output.txt
check_error $?

echo "`date '+%Y-%m-%d %H:%M:%S'`"

: <<'#__COMMENT_OUT__'
#__COMMENT_OUT__
