Input_dir='/staging2/reserve/flagship/u3121714/kenny/Input/BaseData/hg38/HBOC_200/VCF/'
tool='/staging2/reserve/flagship/u3121714/kenny/Aging_script/Team1/Script/tool/bin/vcftools_v0.1.16/bin/'
filename='HBOC_base.list'

## Check for dir, if not found create it using the mkdir ##
convert_dir='/staging2/reserve/flagship/u3121714/kenny/Input/BaseData/hg38/HBOC_200/convert_vcf/'
[ ! -d "$convert_dir" ] && mkdir -p "$convert_dir"
output_path=${convert_dir}


while IFS='' read -r line || [[ -n "$line" ]]; do
	# Check file exist or not
        echo "==> $line"
        echo -e "`ls ${Input_dir} | grep ${line}`"
	# Remove indels
	${tool}vcftools --vcf "${Input_dir}${line}".vcf --remove-indels --recode --recode-INFO-all --out ${Input_dir}${line}
  # Convert VCF4.1 to VCF4.2
  cat "${Input_dir}${line}".recode.vcf | ${tool}vcf-convert -v 4.2 > ${convert_dir}${line}.vcf
done < ${Input_dir}${filename}
