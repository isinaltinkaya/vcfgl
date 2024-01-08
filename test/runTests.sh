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
		exitCode=$?

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
		exitCode=$?

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

threads=1
seed=42
depth=1
depthsFile=""
err=0.2
qs=0
betavar=""
GL=1
gl1theta="--gl1-theta 0.83"
platform=0
usePreciseGlError=0
explode=1
rmInvarSites=0
rmEmptySites=0
doUnobserved=1 # **
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=0
gvcfDps=""


ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# TEST2
# # ** 2 sample input file

ID="test2"

INFILENAME=${DATADIR}/"data2.vcf"

threads=1
seed=42
depth=2
depthsFile=""
err=0.01
qs=0
betavar=""
GL=2 # **
gl1theta=""
platform=0
usePreciseGlError=0
explode=0 # **
rmInvarSites=0
rmEmptySites=0
doUnobserved=2 # **
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=0
addInfoADF=0
addFormatADR=0
addInfoADR=0
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf  ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf


# ################################################################################
# # TEST3
# data3
# # ** err 0
# # ** --rm-empty-sites 1 // removes site 6 and 8
# # ** --rm-invar-sites 3 // 1+2, 1:removes sites 1,2,3,7,10 2:removes site 9

ID="test3"

INFILENAME=${DATADIR}/"data3.vcf"

threads=1
seed=42
depth=1
depthsFile=""
err=0
qs=0
betavar=""
GL=1
gl1theta="--gl1-theta 0.83"
platform=0
usePreciseGlError=0
explode=1 # **
rmInvarSites=3 # **
rmEmptySites=1 # **
doUnobserved=3 # **
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

# ################################################################################
# # TEST4
# test --depth inf
ID="test4"

INFILENAME=${DATADIR}/"data3.vcf"

threads=1
seed=42
depth="inf"
depthsFile=""
err=0
qs=0
betavar=""
GL=1
gl1theta="--gl1-theta 0.83"
platform=0
usePreciseGlError=0
explode=1 
rmInvarSites=0 
rmEmptySites=1 
doUnobserved=1
doGVCF=0
printPileup=0
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=0
addQS=0
addFormatDP=0
addInfoDP=0
addFormatAD=0
addInfoAD=0
addFormatADF=0
addInfoADF=0
addFormatADR=0
addInfoADR=0
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf  ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf


################################################################################
# # TEST5
ID="test5"
# 0.01	error-rate
# 0		error-qs **
INFILENAME=${DATADIR}/"data3.vcf"

threads=1
seed=42
depth=2
depthsFile=""
err=0.01
qs=0
betavar=""
GL=2
gl1theta=""
platform=0
usePreciseGlError=0
explode=1 
rmInvarSites=0
rmEmptySites=0
doUnobserved=4 # **
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# # TEST6
ID="test6"

INFILENAME=${DATADIR}/"data3.vcf"

threads=1
seed=42
depth=2
depthsFile=""
err=0.01
qs=2
betavar="--beta-variance 1e-5"
GL=2
gl1theta=""
platform=0
usePreciseGlError=0
explode=1 
rmInvarSites=0
rmEmptySites=0
doUnobserved=4
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# # TEST7
ID="test7"

# gvcf with --gvcf-dps 1
INFILENAME=${DATADIR}/"data1.vcf"

threads=1
seed=42
depth=1
depthsFile=""
err=0.2
qs=0
betavar=""
GL=1
gl1theta="--gl1-theta 0.83"
platform=0
usePreciseGlError=0
explode=1 
rmInvarSites=0
rmEmptySites=1
doUnobserved=2
doGVCF=1
printPileup=0
printTruth=0
addGL=1
addGP=1
addPL=1
addI16=0
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=1
gvcfDps="--gvcf-dps 1"

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf



################################################################################
# # TEST8
ID="test8"

# gvcf with --gvcf-dps 
INFILENAME=${DATADIR}/"data2.vcf"

threads=1
seed=42
depth=4
depthsFile=""
err=0.001
qs=0
betavar=""
GL=2
gl1theta="--gl1-theta 0.83"
platform=0
usePreciseGlError=0
explode=1 
rmInvarSites=0
rmEmptySites=1
doUnobserved=2
doGVCF=1
printPileup=0
printTruth=0
addGL=1
addGP=1
addPL=1
addI16=0
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=1
gvcfDps="--gvcf-dps 1"

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# # TEST9
ID="test9"

INFILENAME=${DATADIR}/"data3.vcf"

threads=1
seed=42
depth=2
depthsFile=""
err=0.024
qs=2
betavar="--beta-variance 1e-5"
GL=2
gl1theta=""
platform=0
usePreciseGlError=0
explode=1 
rmInvarSites=0
rmEmptySites=0
doUnobserved=4
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# # TEST10
ID="test10"

# same as test9 but --precise-gl 1
INFILENAME=${DATADIR}/"data3.vcf"

threads=1
seed=42
depth=2
depthsFile=""
err=0.024
qs=2
betavar="--beta-variance 1e-5"
GL=2
gl1theta=""
platform=0
usePreciseGlError=1
explode=1 
rmInvarSites=0
rmEmptySites=0
doUnobserved=4
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depth ${depth} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
"

runTest ${ID} ${INFILENAME} "-i" "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# # TEST11
ID="test11"

# same as test9 but use depths file
INFILENAME=${DATADIR}/"data3.vcf"
depthsFile=${DATADIR}/"data3.depthsfile.txt"

threads=1
seed=42
err=0.024
qs=2
betavar="--beta-variance 1e-5"
GL=2
gl1theta=""
platform=0
usePreciseGlError=0
explode=1 
rmInvarSites=0
rmEmptySites=0
doUnobserved=4
doGVCF=0
printPileup=1
printTruth=1
addGL=1
addGP=1
addPL=1
addI16=1
addQS=1
addFormatDP=1
addInfoDP=1
addFormatAD=1
addInfoAD=1
addFormatADF=1
addInfoADF=1
addFormatADR=1
addInfoADR=1
doGVCF=0
gvcfDps=""

ARGS="--verbose 0 \
--threads ${threads} \
--seed ${seed} \
--output-mode v \
--depths-file ${depthsFile} \
--error-rate ${err} \
--error-qs ${qs} \
${betavar} \
--gl-model ${GL} \
${gl1theta} \
--platform ${platform} \
--precise-gl ${usePreciseGlError} \
-explode ${explode} \
--rm-invar-sites ${rmInvarSites} \
--rm-empty-sites ${rmEmptySites} \
-doUnobserved ${doUnobserved} \
-doGVCF ${doGVCF} \
-printPileup ${printPileup} \
-printTruth ${printTruth} \
-addGL ${addGL} \
-addGP ${addGP} \
-addPL ${addPL} \
-addI16 ${addI16} \
-addQS ${addQS} \
-addFormatDP ${addFormatDP} \
-addInfoDP ${addInfoDP} \
-addFormatAD ${addFormatAD} \
-addInfoAD ${addInfoAD} \
-addFormatADF ${addFormatADF} \
-addInfoADF ${addInfoADF} \
-addFormatADR ${addFormatADR} \
-addInfoADR ${addInfoADR} \
-doGVCF ${doGVCF} \
${gvcfDps}
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




