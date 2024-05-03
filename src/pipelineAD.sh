for case in "$@"
do
  mkdir $case
  cd $case
  
  echo
  echo "========================="
  echo "    Bowtie alignement    "
  echo "========================="
  echo
  
  bowtie2 -U /home/BCG2024_genomics_exam/case${case}_father.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o case${case}_father.bam
  bowtie2 -U /home/BCG2024_genomics_exam/case${case}_child.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o case${case}_child.bam
  bowtie2 -U /home/BCG2024_genomics_exam/case${case}_mother.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o case${case}_mother.bam

  samtools index ./${case}_mother.bam 
  samtools index ./${case}_father.bam 
  samtools index ./${case}_child.bam 
  echo
  echo "========================="
  echo "  Bedgraph Generation    "
  echo "========================="
  echo

  bedtools genomecov -ibam case${case}_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > fatherCov.bg
  bedtools genomecov -ibam case${case}_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > motherCov.bg
  bedtools genomecov -ibam case${case}_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > childCov.bg

  echo
  echo "========================="
  echo "Freebayes variant calling"
  echo "========================="
  echo

  freebayes -f /home/BCG2024_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 -q 10 --min-coverage 10 case${case}_mother.bam case${case}_child.bam case${case}_father.bam  > case${case}.vcf

  bcftools query -l case${case}.vcf | sort > samples.txt
  bcftools view -S samples.txt case${case}.vcf > case${case}.sorted.vcf 

  grep "#" case${case}.sorted.vcf > candilist${case}.vcf
  grep "0/1.*0/0.*0/0" case${case}.sorted.vcf  >> candilist${case}.vcf

  bedtools intersect -a candilist${case}.vcf -b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u > ${case}candilistTG.vcf

  echo
  echo "========================="
  echo "     Qualimap Reports    "
  echo "========================="
  echo

  qualimap bamqc -bam case${case}_father.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case${case}_father
  qualimap bamqc -bam case${case}_mother.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case${case}_mother
  qualimap bamqc -bam case${case}_child.bam -gff  /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case${case}_child

  echo
  echo "========================="
  echo "     FastQC Reports      "
  echo "========================="
  echo

  mkdir ./fastqc_${case}/
  fastqc /home/BCG2024_genomics_exam/case${case}_father.fq.gz -o ./fastqc_${case}/
  fastqc /home/BCG2024_genomics_exam/case${case}_child.fq.gz -o ./fastqc_${case}/
  fastqc /home/BCG2024_genomics_exam/case${case}_mother.fq.gz -o ./fastqc_${case}/

  echo
  echo "========================="
  echo "     MultiQC Reports     "
  echo "========================="
  echo

  multiqc . -o ./multiqc_${case}


  cd ..
done

echo "FINISHED" 

#command to consider "sort -n -k 6 vcffile | tail -n20 | less"