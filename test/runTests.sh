#!/usr/bin/env bash
################################################################################
# runTests.sh
#
set -uo pipefail

SCRIPTPATH=$(realpath "$0")

TESTTYPE=${1:-"regular"}

SCRIPTDIR=$(dirname "$SCRIPTPATH")
DATADIR=$(realpath "$SCRIPTDIR/data")
TESTWD=$(realpath "$SCRIPTDIR/testwd")
EXEC=$(realpath "$SCRIPTDIR/../vcfgl")

GREEN='\033[0;32m'
RED='\033[0;31m'
NOCOLOR='\033[0m'

rm -rfv ${TESTWD}/
mkdir -pv ${TESTWD}/
echo ${TESTWD}


testExec(){
	if ! command -v ${EXEC} &> /dev/null; then
		printf "${RED}\n\n"
		printf "###############################################################################\n"
		printf "# ERROR\n"
		printf "# Executable could not be found at:\n"
		printf "${EXEC}\n"
		printf "###############################################################################\n"
		printf "\n\n"
		printf ${NOCOLOR}
		exit 1;
	fi
}

testValgrind(){
	if ! command -v valgrind &> /dev/null; then
		echo "valgrind could not be found";
		exit 1;
	fi
}


if [ ${TESTTYPE} == "regular" ]; then
	testExec
elif [ ${TESTTYPE} == "vg" ]; then
	testExec
	testValgrind
elif [ ${TESTTYPE} == "all" ]; then
	testExec
	testValgrind
else
	echo "Unknown test type: ${TESTTYPE}"
	exit 1;
fi

printf "\n\n"
printf "###############################################################################\n"
printf "# Starting ${TESTTYPE} tests\n"
printf "\n# Script path:\n"
printf "${SCRIPTPATH}\n"
printf "\n# Script directory:\n"
printf "${SCRIPTDIR}\n"
printf "\n# Data directory:\n"
printf "${DATADIR}\n"
printf "\n# Test working directory:\n"
printf "${TESTWD}\n"
printf "\n# Executable:\n"
printf "${EXEC}\n"
printf "\n# Test data directory:\n"
printf "${DATADIR}\n"
printf "###############################################################################\n"
printf "\n\n"


runTestDiff(){
	if [ ${TESTTYPE} == "regular" ]  || [ ${TESTTYPE} == "all" ]; then
		local id=${1}
		local outFile=${2}
		local refFile=${3}
		local outPref=${TESTWD}/${id}
		local logFile=${outPref}.log
		local diffFile=${outPref}.diff

		local diffcmd="diff -s ${outFile} ${refFile} > ${diffFile} 2>&1"

		printf "# ${id} -> Running diff\n"
		printf "# Command:\n${diffcmd}\n"

		eval ${diffcmd}

		if [ $? -eq 0 ]; then
			printf "${GREEN}"
			printf "# ${id} -> Diff: OK\n"
			printf "${NOCOLOR}"
			printf "\n\n"
		else
			printf "\n\n"
			printf "${RED}"
			printf "###############################################################################\n"
			printf "# ${id} FAILED\n"

			printf "\n# Command:\n${diffcmd}\n"
			printf "\n# Output file:\n${outFile}\n"
			printf "\n# Reference file:\n${refFile}\n"
			printf "\n# Log file:\n${logFile}\n"
			printf "\n# Diff file:\n${diffFile}\n"
			printf "###############################################################################\n"
			printf "${NOCOLOR}"
			printf "\n\n"
			exit 1
		fi
	fi
}

runTestDiffPileup(){
        if [ ${TESTTYPE} == "regular" ]  || [ ${TESTTYPE} == "all" ]; then
                local id=${1}
                local outFile=${2}
                local refFile=${3}
                local outPref=${TESTWD}/${id}
                local logFile=${outPref}.log
                local diffFile=${outPref}.diff

                local diffcmd="diff -s <(zcat ${outFile}) <(zcat ${refFile}) > ${diffFile} 2>&1"

                printf "# ${id} -> Running diff\n"
                printf "# Command:\n${diffcmd}\n"

                eval ${diffcmd}

                if [ $? -eq 0 ]; then
                        printf "${GREEN}"
                        printf "# ${id} -> Diff: OK\n"
                        printf "${NOCOLOR}"
                        printf "\n\n"
                else
                        printf "\n\n"
                        printf "${RED}"
                        printf "###############################################################################\n"
                        printf "# ${id} FAILED\n"

                        printf "\n# Command:\n${diffcmd}\n"
                        printf "\n# Output file:\n${outFile}\n"
                        printf "\n# Reference file:\n${refFile}\n"
                        printf "\n# Log file:\n${logFile}\n"
                        printf "\n# Diff file:\n${diffFile}\n"
                        printf "###############################################################################\n"
                        printf "${NOCOLOR}"
                        printf "\n\n"
                        exit 1
                fi
        fi
}


runTestDiffVcf(){
	if [ ${TESTTYPE} == "regular" ]  || [ ${TESTTYPE} == "all" ]; then
		local id=${1}
		local outFile=${2}
		local refFile=${3}
		local outPref=${TESTWD}/${id}
		local logFile=${outPref}.log
		local diffFile=${outPref}.diff

		local diffcmd="diff -s -I '^##' ${outFile} ${refFile} > ${diffFile} 2>&1"

		printf "# ${id} -> Running diff\n"
		printf "# Command:\n${diffcmd}\n"

		eval ${diffcmd}

		if [ $? -eq 0 ]; then
			printf "${GREEN}"
			printf "# ${id} -> Diff: OK\n"
			printf "${NOCOLOR}"
			printf "\n\n"
		else
			printf "\n\n"
			printf "${RED}"
			printf "###############################################################################\n"
			printf "# ${id} FAILED\n"

			printf "\n# Command:\n${diffcmd}\n"
			printf "\n# Output file:\n${outFile}\n"
			printf "\n# Reference file:\n${refFile}\n"
			printf "\n# Log file:\n${logFile}\n"
			printf "\n# Diff file:\n${diffFile}\n"
			printf "###############################################################################\n"
			printf "${NOCOLOR}"
			printf "\n\n"
			exit 1
		fi
	fi
}

runTest(){

	local id=${1}
	local inFileName=${2}
	local inOpt=${3}
	local args=${4}
	local outPref=${TESTWD}/${id}
	local outFile=${outPref}.vcf
	local refFile=${SCRIPTDIR}/reference/${id}/${id}.vcf
	local diffFile=${outPref}.diff
	local logFile=${outPref}.log

	local cmd;

	if [ ${TESTTYPE} == "regular" ]; then
		cmd="${EXEC} ${inOpt} ${inFileName} -o ${outPref} ${args} 2> ${logFile}"
	elif [ ${TESTTYPE} == "vg" ]; then
		cmd="valgrind --leak-check=full -q --error-exitcode=1 --log-fd=9 9>> ${outPref}.vg.log ${EXEC} ${inOpt} ${inFileName} -o ${outPref} ${args} 2> ${logFile}"
	elif [ ${TESTTYPE} == "all" ]; then
		cmd="valgrind --leak-check=full -q --error-exitcode=1 --log-fd=9 9>> ${outPref}.vg.log ${EXEC} ${inOpt} ${inFileName} -o ${outPref} ${args} 2> ${logFile}"
	fi

	printf "\n\n"
	printf "###############################################################################\n"
	printf "# ${id} -> Running ${TESTTYPE} test \n"
	printf "\n\n"

	printf "# Command:\n${cmd}\n"

	
	if [ ${TESTTYPE} == "vg" ] || [ ${TESTTYPE} == "all" ]; then

		local vgFile=${outPref}.vg.log
		local vgrun;
		eval ${cmd}
		local exitCode=$?

		if [ ${exitCode} -eq 0 ]; then
			printf "${GREEN}"
			printf "# ${id} -> Valgrind test: OK\n"
			printf "${NOCOLOR}"
			printf "\n\n"
		else
			printf "${RED}\n\n"
			printf "# ${id} -> Valgrind test: FAILED\n"
			printf "# Valgrind output:\n"
			printf "${vgFile}\n"

			printf "\n\n"
			printf "${NOCOLOR}"
			exit 1

		fi

	elif [ ${TESTTYPE} == "regular" ]; then
		eval ${cmd}
		local exitCode=$?

		if [ ${exitCode} -eq 0 ]; then
			printf "${GREEN}"
			printf "# ${id} -> Run: OK\n"
			printf "${NOCOLOR}"
			printf "\n\n"
		else
			printf "${RED}\n\n"
			printf "# ${id} -> Run: FAILED\n"
			printf "# Log file:\n"
			printf "${logFile}\n"

			printf "\n\n"
			printf "${NOCOLOR}"
			exit 1

		fi

	fi

}






################################################################################
# TEST1
# test all tags together, gl 1
ID="test1"

INFILENAME=${DATADIR}/"data1.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 1 \
--error-rate 0.2 \
--gl-model 1 \
--precise-gl 0 \
--adjust-qs 3 \
-explode 1 \
-doUnobserved 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# TEST2
# # ** 2 sample input file

ID="test2"

INFILENAME=${DATADIR}/"data2.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 2 \
--error-rate 0.01 \
--gl-model 2 \
--precise-gl 0 \
-explode 0 \
-doUnobserved 2 \
-printTruth 1 \
--adjust-qs 3 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf  ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf


# ################################################################################
# TEST3
# data3
# ** err 0
# ** --rm-empty-sites 1 // removes site 6 and 8
# ** --rm-invar-sites 3 // 1+2, 1:removes sites 1,2,3,7,10 2:removes site 9
# ** --doUnobserved 3 // explode A,C,G,T
# ** -explode 1 // explode to unobserved invar sites

ID="test3"

INFILENAME=${DATADIR}/"data3.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 1 \
--error-rate 0 \
--gl-model 1 \
--precise-gl 0 \
-explode 1 \
--rm-invar-sites 3 \
--rm-empty-sites 1 \
--adjust-qs 3 \
-doUnobserved 3 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

# ################################################################################
# TEST4
# test --depth inf
ID="test4"

INFILENAME=${DATADIR}/"data3.vcf"


ARGS="--seed 42 \
--output-mode v \
--depth inf \
--error-rate 0 \
--gl-model 1 \
--precise-gl 0 \
-explode 1 \
--rm-empty-sites 1 \
--adjust-qs 1 \
-doUnobserved 1 \
-printTruth 1 \
-addGP 1 \
-addPL 1 \
-addI16 0 \
-addQS 0 \
-addFormatDP 0
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf  ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf


################################################################################
# TEST5
ID="test5"
# 0.01	error-rate
# 0		error-qs **
# ** -doUnobserved 4 // <*> + explode A,C,G,T
INFILENAME=${DATADIR}/"data3.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 2 \
--error-rate 0.01 \
--gl-model 2 \
--precise-gl 0 \
--adjust-qs 3 \
-explode 1 \
-doUnobserved 4 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# TEST6
ID="test6"

INFILENAME=${DATADIR}/"data3.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 2 \
--error-rate 0.01 \
--error-qs 2 \
--beta-variance 1e-5 \
--gl-model 2 \
--precise-gl 0 \
--adjust-qs 3 \
-explode 1 \
-doUnobserved 4 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# TEST7
ID="test7"
# ** gvcf with --gvcf-dps 1
# ** -doUnobserved 2 // <NON_REF>
INFILENAME=${DATADIR}/"data1.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 1 \
--error-rate 0.2 \
--gl-model 1 \
--gl1-theta 0.83 \
--precise-gl 0 \
-explode 1 \
--rm-empty-sites 1 \
--adjust-qs 3 \
-doUnobserved 2 \
-addGP 1 \
-addPL 1 \
-addI16 0 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 \
-doGVCF 1 \
--gvcf-dps 1
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# TEST8
ID="test8"

INFILENAME=${DATADIR}/"data2.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 4 \
--error-rate 0.001 \
--gl-model 2 \
--precise-gl 0 \
-explode 1 \
--rm-empty-sites 1 \
-doUnobserved 2 \
--adjust-qs 3 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addI16 0 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 \
-doGVCF 1 \
--gvcf-dps 1
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# TEST9
ID="test9"

INFILENAME=${DATADIR}/"data3.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 2 \
--error-rate 0.024 \
--error-qs 2 \
--beta-variance 1e-5 \
--gl-model 2 \
--precise-gl 0 \
--adjust-qs 3 \
-explode 1 \
-doUnobserved 4 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# TEST10
ID="test10"

# same as test9 but --precise-gl 1
INFILENAME=${DATADIR}/"data3.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 2 \
--error-rate 0.024 \
--error-qs 2 \
--beta-variance 1e-5 \
--gl-model 2 \
--precise-gl 1 \
--adjust-qs 2 \
-explode 1 \
-doUnobserved 4 \
-printPileup 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffPileup ${ID} ${TESTWD}/${ID}.pileup.gz ${SCRIPTDIR}/reference/${ID}/${ID}.pileup.gz



################################################################################
# TEST11
ID="test11"

# same as test9 but use depths file
INFILENAME=${DATADIR}/"data3.vcf"
depthsFile=${DATADIR}/"data3.depthsfile.txt"

ARGS="--seed 42 \
--output-mode v \
--depths-file ${depthsFile} \
--error-rate 0.024 \
--error-qs 2 \
--beta-variance 1e-5 \
--gl-model 2 \
--precise-gl 0 \
--adjust-qs 3 \
-explode 1 \
-doUnobserved 4 \
-addGL 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-addFormatADF 1 \
-addInfoADF 1 \
-addFormatADR 1 \
-addInfoADR 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

###############################################################################
# TEST12
ID="test12"

INFILENAME=${DATADIR}/"data4_acgt_biallelic.vcf"

ARGS="--seed 42 \
--depth 2 \
--output-mode v \
--error-rate 0.024 \
--gl-model 2 \
--source 1 \
--adjust-qs 1 \
-explode 1 \
-doUnobserved 1 \
-addPL 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf



###############################################################################
# TEST13
ID="test13"
#
# diff -s <(./misc/fetchGl -i test/testwd/test12.vcf -gt AA 2> /dev/null) <(./misc/fetchGl -i test/testwd/test13.vcf -gt GG 2> /dev/null)
# diff -s <(./misc/fetchGl -i test/testwd/test12.vcf -gt CC 2> /dev/null) <(./misc/fetchGl -i test/testwd/test13.vcf -gt TT 2> /dev/null)

INFILENAME=${DATADIR}/"data4_acgt_biallelic_a2g_c2t.vcf"

ARGS="--seed 42 \
--depth 2 \
--output-mode v \
--error-rate 0.024 \
--gl-model 2 \
--source 1 \
--adjust-qs 1 \
-explode 1 \
-doUnobserved 2 \
-addPL 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


###############################################################################
# TEST14
ID="test14"
#
# diff -s <(./misc/fetchGl -i test/testwd/test12.vcf -gt AA 2>&1 ) <(./misc/fetchGl -i test/testwd/test14.vcf -gt AA 2>&1 )
# diff -s <(./misc/fetchGl -i test/testwd/test12.vcf -gt TT 2>&1 ) <(./misc/fetchGl -i test/testwd/test14.vcf -gt TT 2>&1 )
# diff -s <(./misc/fetchGl -i test/testwd/test12.vcf -gt GA 2>&1 ) <(./misc/fetchGl -i test/testwd/test14.vcf -gt GA 2>&1 )
# diff -s <(./misc/fetchGl -i test/testwd/test12.vcf -gt AC 2>&1 ) <(./misc/fetchGl -i test/testwd/test14.vcf -gt AC 2>&1 )


INFILENAME=${DATADIR}/"data4_binary_biallelic.vcf"

ARGS="--seed 42 \
--depth 2 \
--output-mode v \
--error-rate 0.024 \
--gl-model 2 \
--source 0 \
--adjust-qs 3 \
-explode 1 \
-doUnobserved 0 \
-addGL 1 \
-addQS 1
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


###############################################################################
# TEST15
ID="test15"
# same as test2, test multithreading

INFILENAME=${DATADIR}/"data2.vcf"

ARGS="--seed 42 \
--threads 4 \
--output-mode b \
--depth 2 \
--error-rate 0.01 \
--gl-model 2 \
--precise-gl 0 \
-explode 0 \
-doUnobserved 2 \
--adjust-qs 3 \
-printTruth 1 \
-printPileup 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"


###############################################################################
# TEST16
ID="test16"

# same as test2, test qs binning

INFILENAME=${DATADIR}/"data2.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 2 \
--error-rate 0.01 \
--gl-model 2 \
--precise-gl 0 \
-explode 0 \
-doUnobserved 2 \
-printTruth 1 \
--qs-bins test/data/rta3_qs_bins.csv \
--adjust-qs 3 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf ${SCRIPTDIR}/reference/test2/test2.truth.vcf

###############################################################################
# TEST17
ID="test17"

INFILENAME="${DATADIR}/data6.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 2 \
--error-rate 0 \
--gl-model 2 \
--precise-gl 0 \
--source 1 \
-explode 0 \
-doUnobserved 5 \
-printTruth 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf

###############################################################################
# TEST18
# test multiple contigs
ID="test18"

INFILENAME="${DATADIR}/data7.vcf"

ARGS="--seed 42 \
--output-mode v \
--depth 100 \
--error-rate 0 \
--gl-model 2 \
--precise-gl 0 \
--source 1 \
-explode 1 \
-doUnobserved 5 \
-printTruth 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf


###############################################################################
# TEST19
# test multiple contigs with gvcf
ID="test19"

INFILENAME="${DATADIR}/data7.vcf"


ARGS="--seed 42 \
--output-mode v \
--depth 100 \
--error-rate 0 \
--gl-model 2 \
--precise-gl 0 \
--source 1 \
-explode 1 \
-doUnobserved 5 \
-printTruth 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
-doGVCF 1 \
--gvcf-dps 1
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf



###############################################################################
# TEST20
# test --retain-refalt
ID="test20"

INFILENAME="${DATADIR}/data8.vcf"


ARGS="--seed 42 \
--output-mode v \
--depth 1 \
--error-rate 0 \
--gl-model 2 \
--precise-gl 0 \
--source 1 \
-explode 1 \
-doUnobserved 5 \
-printTruth 1 \
-addGP 1 \
-addPL 1 \
-addI16 1 \
-addQS 1 \
-addInfoDP 1 \
-addFormatAD 1 \
-addInfoAD 1 \
--retain-refalt 1
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

###############################################################################


${EXEC}


echo "__"
# check if logfile looks ok 
wc -l ${TESTWD}/test1.log
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




