#conda activate cutadaptenv

for line in ../*_1.fastq.gz
do
Forward=$line
Forward_output="${line/.fastq.gz/_noPrimer.fastq}"

echo $Forward
#echo $Forward_output

cutadapt -g NNNNCCTACGGGAGGCAGCAG -e 9 $Forward > $Forward_output --quiet
done



for line in ../*_2.fastq.gz
do
Reverse=$line
Reverse_output="${line/.fastq.gz/_noPrimer.fastq}"

echo $Reverse

cutadapt -g GGGGACTACVSGGGTATCTAAT -e 11 $Reverse > $Reverse_output --quiet
done
