#!/usr/bin/env bash
################################################################################
# runTests.sh
#
set -uo pipefail


SCRIPTPATH=$(realpath "$0")
SCRIPTDIR=$(dirname "$SCRIPTPATH")
TESTWD=$(realpath "$SCRIPTDIR/testwd")
EXEC=$(realpath "$SCRIPTDIR/../vcfgl")

GREEN='\033[0;32m'
RED='\033[0;31m'
NOCOLOR='\033[0m'

rm -rfv ${TESTWD}/
mkdir -pv ${TESTWD}/
echo ${TESTWD}



initMainLog(){
	local id=${1}
	printf "\n\n"
	printf "###############################################################################\n"
	printf "# RUNNING TEST: ${id}\n"
	printf "\n\n"
}


printMainLog(){
	local msg=${@}
	printf "\n# ${msg}\n"
}

testSuccess(){
	local id=${1}
	printf "${GREEN}\n\n"
	printf "# FINISHED ${id} -> OK\n"
	printf "${NOCOLOR}"
	printf "###############################################################################\n"
	printf "\n\n"
}

testFail(){
	local id=${1}
	local outFile=${2}
	local refFile=${3}
	local logFile=${4}
	local diffFile=${5}

	printf "\n\n"
	printf "${RED}"
	printf "###############################################################################\n"
	printf "# TEST ${id}: FAILED\n"
	printMainLog "Output file:\n${outFile}"
	printMainLog "Reference file:\n${refFile}"
	printMainLog "Log file:\n${logFile}"
	printMainLog "Diff file:\n${diffFile}"
	printf "###############################################################################\n"
	printf "${NOCOLOR}"
	printf "\n\n"
	exit 1
}

runTestDiff(){
	local id=${1}
	local outFile=${2}
	local refFile=${3}
	local outPref=${TESTWD}/${id}
	local logFile=${outPref}.log
	local diffFile=${outPref}.diff

	diff -s ${outFile} ${refFile} > ${diffFile} 2>&1

	if [ $? -eq 0 ]; then
		testSuccess ${id}
	else
		testFail ${id} ${outFile} ${refFile} ${logFile} ${diffFile}
	fi
}


runTestDiffVcf(){
	local id=${1}
	local outFile=${2}
	local refFile=${3}
	local outPref=${TESTWD}/${id}
	local logFile=${outPref}.log
	local diffFile=${outPref}.diff

	diff -s -I '^##' ${outFile} ${refFile} > ${diffFile} 2>&1

	if [ $? -eq 0 ]; then
		testSuccess ${id}
	else
		testFail ${id} ${outFile} ${refFile} ${logFile} ${diffFile}
	fi
}


runTest(){

	local id=${1}
	local inFileName=${2}
	local args=${3}
	local inFile=${SCRIPTDIR}/data/${inFileName}
	local outPref=${TESTWD}/${id}
	local outFile=${outPref}.vcf
	local refFile=${SCRIPTDIR}/reference/${id}/${id}.vcf
	local diffFile=${outPref}.diff
	local logFile=${outPref}.log
	local cmd="${EXEC} -i ${inFile} -o ${outPref} ${args}"
	initMainLog ${id}
	printMainLog "Command:\n${cmd}"

	${cmd} > ${logFile} 2>&1

}



################################################################################
# TEST1
# test all tags together, gl 1
ID="test1"


INFILENAME="data1.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# TEST2
# # ** 2 sample input file

ID="test2"

INFILENAME="data2.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf  ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf


# ################################################################################
# # TEST3
# data3
# # ** err 0
# # ** --rm-empty-sites 1 // removes site 6 and 8
# # ** --rm-invar-sites 3 // 1+2, 1:removes sites 1,2,3,7,10 2:removes site 9

ID="test3"

INFILENAME="data3.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

# ################################################################################
# # TEST4
# test --depth inf
ID="test4"

INFILENAME="data3.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf
runTestDiffVcf ${ID} ${TESTWD}/${ID}.truth.vcf  ${SCRIPTDIR}/reference/${ID}/${ID}.truth.vcf

# //TODO

################################################################################
# # TEST5
ID="test5"
# 0.01	error-rate
# 0		error-qs **
INFILENAME="data3.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# # TEST6
ID="test6"

INFILENAME="data3.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf


################################################################################
# # TEST7
ID="test7"

# gvcf with --gvcf-dps 1
INFILENAME="data1.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf



################################################################################
# # TEST8
ID="test8"

# gvcf with --gvcf-dps 
INFILENAME="data2.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# # TEST9
ID="test9"

INFILENAME="data3.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
runTestDiffVcf ${ID} ${TESTWD}/${ID}.vcf ${SCRIPTDIR}/reference/${ID}/${ID}.vcf

################################################################################
# # TEST10
ID="test10"

# same as test9 but --precise-gl 1
INFILENAME="data3.vcf"

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

runTest ${ID} ${INFILENAME} "${ARGS}"
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




