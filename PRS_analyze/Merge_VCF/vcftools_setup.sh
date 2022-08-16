###run code once per row

git clone https://github.com/vcftools/vcftools.git
# Start from here ------------------------------------
cd /staging2/reserve/flagship/u3121714/kenny/Aging_script/Team1/Script/tool/bin/vcftools/
#
export PERL5LIB=/staging2/reserve/flagship/u3121714/kenny/Aging_script/Team1/Script/tool/bin/vcftools/src/perl/
#
./autogen.sh
#
./configure --prefix=/staging2/reserve/flagship/u3121714/kenny/Aging_script/Team1/Script/tool/bin/vcftools_v0.1.16/
#
make
#
make install
#
cd ../
#
vcftools_v0.1.16/bin/vcf-convert -h
#------------------------------------------------------
