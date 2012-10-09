#!/bin/bash
#$ -S /bin/bash
#$ -cwd

INPUTBAM_TUM=$1
INPUTBAM_NOR=$2
OUTPUTPATH=$3
REFERENCELIST=$4


source ./config.sh
source ./utility.sh


# :<<_COMMENT_OUT_

# pileup the tumor bam file
##########
echo "${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_TUM} > ${OUTPUTPATH}/tmp/temp.tumor.pileup"
${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_TUM} > ${OUTPUTPATH}/tmp/temp.tumor.pileup
check_error $?

echo "${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_NOR} > ${OUTPUTPATH}/tmp/temp.normal.pileup"
${PATH_TO_SAMTOOLS}/samtools mpileup -BQ0 -d10000000 -f ${PATH_TO_REF} -q ${TH_MAPPING_QUAL} ${INPUTBAM_NOR} > ${OUTPUTPATH}/tmp/temp.normal.pileup
check_error $?
##########




# make count files for mismatches, insertions and deletions
# mismatch count is performed considering bases whose quality is more than ${TH_BASE_QUAL}.
##########
echo "perl subscript/pileup2base.pl ${TH_BASE_QUAL} ${OUTPUTPATH}/tmp/temp.tumor.pileup ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.tumor.depth"
perl subscript/pileup2base.pl ${TH_BASE_QUAL} ${OUTPUTPATH}/tmp/temp.tumor.pileup ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.tumor.depth
check_error $?

echo "perl subscript/pileup2base.pl ${TH_BASE_QUAL} ${OUTPUTPATH}/tmp/temp.normal.pileup ${OUTPUTPATH}/tmp/temp.normal.base ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.normal.depth"
perl subscript/pileup2base.pl ${TH_BASE_QUAL} ${OUTPUTPATH}/tmp/temp.normal.pileup ${OUTPUTPATH}/tmp/temp.normal.base ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.normal.depth
check_error $?
##########


# extract putative germline mutations
##########
# echo "perl filterBase.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.base > ${OUTPUTPATH}/tmp/temp.normal.base.filt"
# perl filterBase.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.base > ${OUTPUTPATH}/tmp/temp.normal.base.filt 
# check_error $?

echo "perl subscript/filterBase_del.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.del > ${OUTPUTPATH}/tmp/temp.normal.del.filt"
perl subscript/filterBase_del.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.del > ${OUTPUTPATH}/tmp/temp.normal.del.filt
check_error $?

echo "perl subscript/filterBase_ins.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.ins > ${OUTPUTPATH}/tmp/temp.normal.ins.filt"
perl subscript/filterBase_ins.barcode.pl ${OUTPUTPATH}/tmp/temp.normal.ins > ${OUTPUTPATH}/tmp/temp.normal.ins.filt
check_error $?
##########

# _COMMENT_OUT_

# filter candidate of variation between tumor and normal
##########
echo "perl subscript/compBase.pl ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.normal.base ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt"
perl subscript/compBase.pl ${OUTPUTPATH}/tmp/temp.tumor.base ${OUTPUTPATH}/tmp/temp.normal.base ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt
check_error $?

echo "perl subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt"
perl subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.ins ${OUTPUTPATH}/tmp/temp.normal.ins ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt
check_error $?

echo "perl subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt"
perl subscript/compInsDel.pl ${OUTPUTPATH}/tmp/temp.tumor.del ${OUTPUTPATH}/tmp/temp.normal.del ${OUTPUTPATH}/tmp/temp.tumor.depth ${OUTPUTPATH}/tmp/temp.normal.depth ${MIN_TUMOR_DEPTH} ${MIN_NORMAL_DEPTH} ${MIN_TUMOR_VARIANT_READ} ${MIN_TUMOR_ALLELE_FREQ} ${MAX_NORMAL_ALLELE_FREQ} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt
check_error $?
##########

# _COMMENT_OUT_

# add information of normal reference
##########
echo "perl subscript/getRefNor_base.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt ${REFERENCELIST} ${TH_BASE_QUAL_REF} ${TH_MAPPING_QUAL_REF} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref"
perl subscript/getRefNor_base.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt ${REFERENCELIST} ${TH_BASE_QUAL_REF} ${TH_MAPPING_QUAL_REF} > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref 

echo "perl subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt ${REFERENCELIST} 1 ${TH_MAPPING_QUAL_REF} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref"
perl subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt ${REFERENCELIST} 1 ${TH_MAPPING_QUAL_REF} > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref

echo "perl subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt ${REFERENCELIST} 2 ${TH_MAPPING_QUAL_REF} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref"
perl subscript/getRefNor_insdel.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt ${REFERENCELIST} 2 ${TH_MAPPING_QUAL_REF} > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref
##########


# caluculate p-value for each candidate by the Beta-binominal sequencing error model 
##########
if [ ! -s ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref ]; then
    echo "make empty file : ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref"
    echo -n > ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb
else
    echo "${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref ${OUTPUTPATH}/tmp/temp$.tumor_normal.base.eb < subscript/proc_EBcall.R"
    ${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.base.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb < subscript/proc_EBcall.R
    check_error $?
fi

if [ ! -s ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref ]; then
    echo "make empty file : ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb"
    echo -n > ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb
else
    echo "${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb < subscript/proc_EBcall.R"
    ${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb < subscript/proc_EBcall.R
    check_error $?
fi

if [ ! -s ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref ]; then
    echo "make empty file : ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb"
    echo -n > ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb
else
    echo "${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb < subscript/proc_EBcall.R"
    ${PATH_TO_R}/R --vanilla --slave --args ${OUTPUTPATH}/tmp/temp.tumor_normal.del.filt.ref ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb < subscript/proc_EBcall.R
    check_error $?
fi
##########

# _COMMENT_OUT_

# merge and convert the format of three variation files (mutation, insertion and deletion) for Annovar
##########
echo "perl subscript/procForAnnovar.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb > ${OUTPUTPATH}/output.txt"
perl subscript/procForAnnovar.pl ${OUTPUTPATH}/tmp/temp.tumor_normal.base.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.ins.eb ${OUTPUTPATH}/tmp/temp.tumor_normal.del.eb > ${OUTPUTPATH}/output.txt
check_error $?
##########


