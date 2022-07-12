#conda activate cutadaptenv

for line in ./*_1.fastq.gz
do 
FORWARD=$line
REVERSE="${line/_1.fastq.gz/_2.fastq.gz}"
MERGE="${line/_1.fastq.gz/}"
Merge_output="${line/_1.fastq.gz/.assembled.fastq}"
Filter_output="${line/_1.fastq.gz/_merge_filter.fastq}"
CAF="${line/_1.fastq.gz/.unassembled.forward.fastq}"
CAR="${line/_1.fastq.gz/.unassembled.reverse.fastq}"
CAD="${line/_1.fastq.gz/.discarded.fastq}"

#echo $FORWARD
#echo $REVERSE
echo $MERGE
#echo $Merge_output
#echo $Filter_output

PEAR -f $FORWARD -r $REVERSE -o $MERGE >> PEAR_output.txt
cutadapt -g NNNNCCTACGGGAGGCAGCAG -e 9 $Merge_output > $Filter_output --quiet
cutadapt -a ATTAGATACCCSBGTAGTCCCC -e 11 $Filter_output > $Merge_output --quiet
mv $Merge_output ./Merge_Filter
mv $Filter_output ./Cut_Adapt_output
mv $CAF ./Cut_Adapt_output
mv $CAR ./Cut_Adapt_output 
mv $CAD ./Cut_Adapt_output 

done
