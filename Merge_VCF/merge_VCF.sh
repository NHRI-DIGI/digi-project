#!/bin/bash
filename='CKD_target.list'
current_dir='/NovaSeq_128/Kenny/CKD/test/Input/Targetdata/VCF/'
tool_dir='/NovaSeq_128/Kenny/tool/bin/'

cd ${current_dir}

## Check for dir, if not found create it using the mkdir ##
create_dir="/NovaSeq_128/Kenny/CKD/test/Input/Targetdata/recode_vcf/"
[ ! -d "$create_dir" ] && mkdir -p "$create_dir"
output_path=${create_dir}


while IFS='' read -r line || [[ -n "$line" ]]; do
	# Check file exist or not
        echo "==> $line"
        echo -e "`ls ${current_dir} | grep ${line}`"
	# Modify : --vcf
	${tool_dir}vcftools --vcf "${current_dir}${line}".vcf --recode --recode-INFO-all \
		--out ${output_path}${line}_DP10_MAF21.vcf \
		--chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr X --chr Y --chr M \
		--min-meanDP 10 \
		--non-ref-af-any 0.21
done < $filename


#Merge all VCF
merge_filename='merge_CKD.list'

## Check for dir, if not found create it using the mkdir ##
merge_dir='/NovaSeq_128/Kenny/CKD/test/Input/Targetdata/VCF_merge/'
[ ! -d "$merge_dir" ] && mkdir -p "$merge_dir"
merge_output_path=${merge_dir}

disease='CKD'
Amount_of_Samples='12'

# move.log file to log file
cd ${create_dir}
mkdir ../log
mv *.log ../log
mv ${merge_filename} ${merge_dir}

# Build 'merge_HBOC_48.list'
ls ${create_dir} > ${merge_dir}${merge_filename}
mv ${merge_dir}${merge_filename} ${create_dir}

while IFS='' read -r line || [[ -n "$line" ]]; do
	# Check file exist or not
        echo "==> $line"
        echo -e "`cat ${create_dir}${merge_filename} | grep ${line}`"
	# Make VCF.gz and index
	${tool_dir}bgzip -c "${create_dir}${line}" > "${create_dir}${line}".gz
  ${tool_dir}tabix -p vcf "${create_dir}${line}".gz
done < ${merge_filename}

# Build list for merging, with route and .gz in the end
sed -e 's/^AA01/\/NovaSeq_128\/Kenny\/CKD\/test\/Input\/Targetdata\/recode_vcf\/AA01/' -i ${create_dir}${merge_filename}
sed -e 's/$/.gz/' -i ${create_dir}${merge_filename}

#Merge VCF.gz
${tool_dir}bcftools merge -0 -l ${create_dir}${merge_filename} --threads 50 --force-samples -Ov -o ${merge_output_path}${disease}_${Amount_of_Samples}.merge.vcf
