# isinaltinkaya 230920
#
# script for testing gp genotype order match


vcf=${1}


readarray -t data < <(bcftools query -f "%POS %REF,%ALT[ %GP]\n" ${vcf})

declare -A lut=(["CA"]="AC" ["GA"]="AG" ["TA"]="AT" ["CG"]="GC" ["TG"]="GT" ["TC"]="CT" ["AC"]="AC" ["AG"]="AG" ["AT"]="AT" ["GC"]="GC" ["GT"]="GT" ["CT"]="CT" ["AA"]="AA" ["CC"]="CC" ["GG"]="GG" ["TT"]="TT")



len="${#data[@]}"

for i in $(seq 0 $((${len} - 1)));do
	line=(${data[${i}]})
	alleles=($(echo ${line[1]}|tr ',' ' '))
	gps=($(echo ${line[2]}|tr ',' ' '))
	site=$(echo ${line[0]})

	nAlleles=${#alleles[@]}

	gts=()

	for x in $(seq 0 $((${nAlleles}-1)) );do
		for y in $(seq 0 ${x} );do
			# printf "${site}"
			# gts+=("${alleles[${y}]}${alleles[${x}]}")
			gtraw=${alleles[${y}]}${alleles[${x}]}
			gt=${lut[${gtraw}]}

			gts+=(${gt})
			# echo ${y} ${x}
		done
	done

	nGts=${#gts[@]}

	for i in $(seq 0 $((${nGts}-1)));do
		printf "${site},${gts[i]},${gps[i]}\n"
	done
done
