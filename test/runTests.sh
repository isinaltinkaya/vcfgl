#!/usr/bin/env bash
####################################################################################################
# runTests.sh
#
set -uo pipefail


SCRIPTPATH=$(realpath "$0")
SCRIPTDIR=$(dirname "$SCRIPTPATH")
TESTWD=$(realpath "$SCRIPTDIR/testwd")
EXEC=$(realpath "$SCRIPTDIR/../vcfgl")

rm -rfv ${TESTWD}/
mkdir -pv ${TESTWD}/
echo ${TESTWD}



initMainLog(){
	local id=${1}
	printf "\n\n"
	printf "###############################################\n"
	printf "# RUNNING TEST: ${id}\n"
	printf "###############################################\n"
	printf "\n\n"
}


printMainLog(){
	local msg=${@}
	printf "\n# ${msg}\n"
}

testSuccess(){
	local id=${1}
	printf "\n\n"
	printf "###############################################\n"
	printf "# TEST ${id}: SUCCESS\n"
	printf "###############################################\n"
	printf "\n\n"
}

testFail(){
	local id=${1}
	printf "\n\n"
	printf "###############################################\n"
	printf "# TEST ${id}: FAILED\n"
	printMainLog "Output file:\n${outFile}"
	printMainLog "Reference file:\n${refFile}"
	printMainLog "Log file:\n${logFile}"
	printMainLog "Diff file:\n${diffFile}"
	printf "###############################################\n"
	printf "\n\n"
	exit 1
}


runTestDiff(){
	diff -s -I '^##' ${outFile} ${refFile} > ${diffFile} 2>&1

	if [ $? -eq 0 ]; then
		testSuccess ${id}
	else
		testFail ${id}
	fi
}

runTest(){

	local id=${1}
	local inFileName=${2}
	local args=${3}
	local inFile=${SCRIPTDIR}/${inFileName}
	local outPref=${TESTWD}/${id}
	local outFile=${outPref}.vcf
	local refFile=${SCRIPTDIR}/reference/${id}.vcf
	local diffFile=${outPref}.diff
	local logFile=${outPref}.log
	local cmd="${EXEC} -i ${inFile} -o ${outPref} ${args}"
	initMainLog ${id}
	printMainLog "Command:\n${cmd}"

	${cmd} > ${logFile} 2>&1

	runTestDiff ${id} ${outFile} ${refFile} ${diffFile}
}



# ###############################################
# defaults:
# ARGS=" -O v \
# --depth 1 \
# --error-rate 0.01 \
# --error-qs 0 \
# --beta-mean -1.0 \
# --beta-variance -1.0 \
# --precise-gl 1 \
# --pos0 0 \
# --trim-alt-alleles 0 \
# --use-unknown-allele 0 \
# --platform 0 \
# --rm-invar-sites 0 \
# -addGL 1 \
# -addGP 0 \
# -addPL 0 \
# -addI16 0 \
# -addQS 0 \
# -addFormatDP 1 \
# -addInfoDP 0 \
# -addFormatAD 0 \
# -addFormatADF 0 \
# -addFormatADR 0 \
# -addInfoAD 0 \
# -addInfoADF 0 \
# -addInfoADR 0 \
# -explode 0 \
# --threads 1 \
# --verbose 0"
# # --seed <UNSET> \
# ###############################################





###############################################
# Test 1: test basic

ID="test1"
INFILENAME="t1.vcf"
ARGS=" -O v \
--depth 1 \
--error-rate 0.01 \
--error-qs 0 \
--precise-gl 1 \
-addGL 1 \
-addFormatDP 1 \
--seed 42 \
--pos0 0 \
-explode 0"

runTest ${ID} ${INFILENAME} "${ARGS}"

###############################################
# Test 2: test explode
# 0 pos0
# 1 explode **

ID="test2"
INFILENAME="t1.vcf"
ARGS=" -O v \
--depth 1 \
--error-rate 0.01 \
--error-qs 0 \
--precise-gl 1 \
-addGL 1 \
-addFormatDP 1 \
--seed 42 \
-explode 1 \
--pos0 0"

runTest ${ID} ${INFILENAME} "${ARGS}"

###############################################
# Test 3: test pos0, 2 sample vcf(t2.vcf)
# 1 pos0 **
# 0 explode

ID="test3"
INFILENAME="t2.vcf"

ARGS=" -O v \
--depth 1 \
--error-rate 0.01 \
--error-qs 0 \
--precise-gl 1 \
-addGL 1 \
-addFormatDP 1 \
--seed 42 \
--pos0 1 \
-explode 0"

runTest ${ID} ${INFILENAME} "${ARGS}"


###############################################
# Test 5: test addGP, addPL (interacts with pos0,explode)
# 1 pos0 *
# 1 explode *
# 1 addGL (default)
# 1 addGP **
# 1 addPL **

ID="test5"
INFILENAME="t2.vcf"

ARGS=" -O v \
--depth 1 \
--error-rate 0.01 \
--error-qs 0 \
--precise-gl 1 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addFormatDP 1 \
--seed 42 \
--pos0 1 \
-explode 1"

runTest ${ID} ${INFILENAME} "${ARGS}"

###############################################
# Test 6:
# 1 addQS**
# 5 depth* for comparison with test7,8
# 			which are equivalent with only diff
#			being --error-qs

ID="test6"
INFILENAME="t2.vcf"


ARGS=" -O v \
--depth 5 \
--error-rate 0.01 \
--error-qs 0 \
--precise-gl 1 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addQS 1 \
-addFormatAD 1 \
-addFormatDP 1 \
--seed 42 \
--pos0 1 \
-explode 1"

runTest ${ID} ${INFILENAME} "${ARGS}"


# ###############################################
# # Test 7:
# 5		depth* for comparison with test6,8
# 0.01	beta-mean*
# 1e-5 	beta-variance*
# 1		error-qs**

ID="test7"
INFILENAME="t2.vcf"

ARGS=" -O v \
--depth 5 \
--error-rate 0.01 \
--beta-variance 1e-4 \
--error-qs 1 \
--precise-gl 1 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addQS 1 \
-addFormatDP 1 \
-addFormatAD 1 \
--seed 42 \
--pos0 1 \
-explode 1"

runTest ${ID} ${INFILENAME} "${ARGS}"


###############################################
# Test 8:
# 5		depth* for comparison with test6,7
# 0.01	beta-mean*
# 1e-5 	beta-variance*
# 2		error-qs**

ID="test8"
INFILENAME="t2.vcf"


ARGS=" -O v \
--depth 5 \
--error-rate 0.01 \
--beta-variance 1e-4 \
--error-qs 2 \
--precise-gl 1 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addQS 1 \
-addFormatDP 1 \
--seed 42 \
--pos0 1 \
-explode 1"

runTest ${ID} ${INFILENAME} "${ARGS}"

###############################################
# Test 9
# 5		depth* for comparison with test6,7
# 0.01	beta-mean*
# 1e-5 	beta-variance*
# 2		error-qs**
# 0		precise-gl*

ID="test9"
INFILENAME="t2.vcf"


ARGS=" -O v \
--depth 5 \
--error-rate 0.01 \
--beta-variance 1e-4 \
--error-qs 2 \
--precise-gl 0 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addQS 1 \
-addFormatDP 1 \
--seed 42 \
--pos0 1 \
-explode 1"

runTest ${ID} ${INFILENAME} "${ARGS}"


# #TODO does not work since trim alts does not work right now
#named: old_test8
# 	# ./vcfgl -i test/t2.vcf -o test/t2_pos01_explode1_gp_gl_qs_trim1_test -O v --seed 42 --trim-alt-alleles 0 --rm-invar-sites 0 --use-unknown-allele 0 -d 1 -e 0.01 --pos0 1 -explode 1 -addQS 1 -addGP 1 -addPL 1 --trim-alt-alleles 0 --rm-invar-sites 0 --use-unknown-allele 0;
# 	# bash -c "diff -I '^##'  test/t2_pos01_explode1_gp_gl_qs_trim1_test.vcf test/reference/t2_pos01_explode1_gp_gl_qs_trim1.vcf";

# # ID="test8"
# # INFILENAME="t2.vcf"
# # ARGS="-O v --seed 42 --trim-alt-alleles 0 --rm-invar-sites 0 --use-unknown-allele 0 -d 1 -e 0.01 --pos0 1 -explode 1 -addQS 1 -addGP 1 -addPL 1 --trim-alt-alleles 0 --rm-invar-sites 0 --use-unknown-allele 0"

# # runTest ${ID} ${INFILENAME} "${ARGS}"
# # REFFILE="reference/t2_pos01_explode1_gp_gl_qs_trim1.vcf"
# # bash ${SCRIPTDIR}/get_gls.sh ${OUTFILE} | LC_ALL=C LC_COLLATE=C sort -t"," -k1,1 -k2,2 > ${TESTWD}/${ID}_cmp_gls_1
# # bash ${SCRIPTDIR}/get_gls.sh ${REFFILE} | LC_ALL=C LC_COLLATE=C sort -t"," -k1,1 -k2,2 > ${TESTWD}/${ID}_cmp_gls_2
# # diff -s ${TESTWD}/${ID}_cmp_gls_1 ${TESTWD}/${ID}_cmp_gls_2

# # bash ${SCRIPTDIR}/get_gts.sh ${OUTFILE} | LC_ALL=C LC_COLLATE=C sort -t"," -k1,1 -k2,2 > ${TESTWD}/${ID}_cmp_gts_1
# # bash ${SCRIPTDIR}/get_gts.sh ${REFFILE} | LC_ALL=C LC_COLLATE=C sort -t"," -k1,1 -k2,2 > ${TESTWD}/${ID}_cmp_gts_2
# # diff -s ${TESTWD}/${ID}_cmp_gts_1 ${TESTWD}/${ID}_cmp_gts_2

# # cp ${OUTFILE} ${REFFILE}
# # diff -s -I '^##' ${OUTFILE} ${REFFILE}

# ###############################################
# # Test 10:
# vcfgl --input test/t3.vcf --output-mode v --error-rate 0.010000 --depth 0.010000 --pos0 0 --seed 42 --error-qs 0 --beta-variance 0.000000e+00 --rm-invar-sites 0 --trim-alt-alleles 0 --platform 0 -explode 1 -addGL 1 -addGP 1 -addPL 1 -addI16 1 -addQS 1  -addFormatDP 1 -addFormatAD 1 -addFormatADF 0 -addFormatADR 0 -addInfoAD 1 -addInfoADF 0 -addInfoADR 0
# vcfgl -i /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_dev_231008/test/t3.vcf -o /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_dev_231008/test/testwd/test10 --output-mode v --error-rate 0.010000 --depth 1 --pos0 0 --seed 42 --error-qs 0 --beta-variance 0.000000e+00 --rm-invar-sites 0 --trim-alt-alleles 0 --platform 0 -explode 1 -addGL 1 -addGP 1 -addPL 1 -addI16 1 -addQS 1  -addFormatDP 1 -addFormatAD 1 -addFormatADF 0 -addFormatADR 0 -addInfoAD 1 -addInfoADF 0 -addInfoADR 0 --use-unknown-allele 0
#
# ID="test10"
# INFILENAME="t3.vcf"
# ARGS="--output-mode v --error-rate 0.010000 --depth 1 --pos0 0 --seed 42 --error-qs 0 --beta-variance 0.000000e+00 --rm-invar-sites 0 --trim-alt-alleles 0 --platform 0 -explode 1 -addGL 1 -addGP 1 -addPL 1 -addI16 1 -addQS 1  -addFormatDP 1 -addFormatAD 1 -addFormatADF 0 -addFormatADR 0 -addInfoAD 1 -addInfoADF 0 -addInfoADR 0 --use-unknown-allele 0"
#
# runTest ${ID} ${INFILENAME} "${ARGS}"



${EXEC}


echo "__"
# check if logfile looks ok 
wc -l ${TESTWD}/${ID}.log
echo "__"

wait

printf "\n\n"
echo -e "                                                                           
 -------------------
< ALL TESTS PASSED! >
 -------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/
                ||----w |
                ||     ||
"

printf "\n\n"




