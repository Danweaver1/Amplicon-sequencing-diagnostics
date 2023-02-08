workingdir=$( pwd )
#NB - folder 'database' containing bowtie2 indexes is required in working directory
#NB BC_file containing forward index list and R_BC_file containing Reverse index list are required in working dir
analysisdir="TEF_merge_analysis"



for file in *R1_001.fastq.gz
do
R2=$( echo $file | sed 's/R1/R2/' )
PREFIX=$( echo $file | sed 's/_001.fastq.gz//' )

bbmerge.sh in1=$file in2=$R2 ecco mix merge=t out="$PREFIX"_bbmerge.fastq outu=unmerged 1>>bbmerge.out 2>>bbmerge.out
done


fastx_barcode_splitter.pl --bcfile BC_file_v2 --prefix FDM_ --bol --exact <"$PREFIX"_bbmerge.fastq 1>>barcode_splitter_F.out 2>>barcode_splitter_F.out &&

for file in FDM_*
do
fastx_barcode_splitter.pl --bcfile R_BC_file_v2 --prefix "$file"_RDM_ --eol --exact <$file 1>>barcode_splitter_R.out 2>>barcode_splitter_R.out
done

#file_name_key should be a csv (comma seperated) with column 1 containing barcode splitter filenames and col2 containing new name to change to (eg. TEF1_MA1.1)
while read p
do
  name=$( echo $p | cut -d"," -f1)
  newname=$( echo $p | cut -d"," -f2)
  mv "$name" "$newname".fastq &&
  echo "$name" renamed to "$newname".fastq >> file_rename_info
done <file_name_key*.csv


mkdir $analysisdir
mkdir $analysisdir/unmatched_fastq
mv *unmatched $analysisdir/unmatched_fastq
mkdir $analysisdir/demulti_fastq
mv TEF*.fastq $analysisdir/demulti_fastq &&

rm FDM_*
cd $analysisdir/demulti_fastq

mkdir Primer_untrimmed
mkdir Primer_QC

for file in TEF*.fastq
  do
  Prefix=$( echo $file | sed 's/.fastq//' )
  cutadapt -q 30 --minimum-length 150 --untrimmed-output $Prefix.Primer_untrimmed.fastq --overlap 10 -g YGGYGARTTCGARGCYGGTAT -o $Prefix.temp.fastq $file 1>>cutadaptF.out
  cutadapt -a GGHTKCMAYGGHGAYAAYATGHTNT --untrimmed-output $Prefix.Primer_untrimmed2.fastq --overlap 10 -o $Prefix.QC.fastq $Prefix.temp.fastq 1>>cutadaptR.out
  done
mv *QC.fastq Primer_QC
mv *untrimmed*.fastq Primer_untrimmed
rm *temp.fastq

cd Primer_QC


source activate samtools1_9
#above line necessary on Dan desktop


val=0.06
db="TEF_v12.5"
dir="$db"_"$val"

mkdir $dir
mkdir $dir/unaligned
mkdir $dir/BAM
mkdir $dir/Counts
mkdir $dir/Counts/Seqs
mkdir $dir/Counts/Species
mkdir $dir/Counts/Genera
mkdir $dir/Counts/Reads


for file in TEF*QC.fastq
do
        if [ -e $file ]
        then
		            PREFIX=`echo $file | sed -e 's/.QC.fastq//'`
		            SEQ="$PREFIX.seq_count"
		            UNALIGNED="$PREFIX.unaligned.fastq"
                SPECIES="$PREFIX.unique_species"
                GENUS="$PREFIX.unique_genera"
                READS="$PREFIX"_reads.csv
                FvHREADS="$PREFIX"_FvHreads.csv
                echo $PREFIX >> $dir/Merge_"$db"_"$val"_bowtie.out
		bowtie2 --end-to-end -D 25 -R 5 -N 0 -L 25 -i S,1,2.5 --rdg 3,3 --score-min L,-0.6,-$val --un $dir/unaligned/$UNALIGNED -U $file -x $workingdir/database/$db 2>>$dir/Merge_"$db"_"$val"_bowtie.out | samtools sort -o $dir/BAM/"$PREFIX".BAM  &&

                samtools view $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep -v "*" | sort | uniq -c | sort -nr >> $dir/Counts/Seqs/$SEQ
                samtools view $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep -v "*" | sort | cut -d"_" -f1-2 | sort | uniq -c | sort -nr > $dir/Counts/Species/$SPECIES
                #create species counts and csv compatible with R plotting script
                cat $dir/Counts/Species/$SPECIES | sed -e 's/^ *//;s/ /,/' | awk -v where=`basename $PREFIX` '{print $1","where}' > $dir/Counts/Species/$PREFIX.Rspecies_counts.csv
                samtools view $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep -v "*" | sort | cut -d"_" -f1 | sort | uniq -c | sort -nr > $dir/Counts/Genera/$GENUS
                cat $dir/Counts/Genera/$GENUS | sed -e 's/^ *//;s/ /,/' | awk -v where=`basename $PREFIX` '{print $1","where}' > $dir/Counts/Genera/$PREFIX.Rgenera_counts.csv
                #generate overview of unaligned and aligned
                samtools view -f 4 $dir/BAM/"$PREFIX".BAM | grep -vE '(SO|PN|LN)' | wc -l | awk -v where=`basename $PREFIX` '{print $1",unaligned,"where}' >> $dir/Counts/Reads/$READS
                samtools view -F 4 $dir/BAM/"$PREFIX".BAM | grep -vE '(SO|PN|LN)' | wc -l | awk -v where=`basename $PREFIX` '{print $1",aligned,"where}' >> $dir/Counts/Reads/$READS

                ##generate overview of unaligned, human and fungal read numbers
                samtools view -f 4 $dir/BAM/"$PREFIX".BAM | grep -vE '(SO|PN|LN)' | wc -l | awk -v where=`basename $PREFIX` '{print $1",unaligned,"where}' >> $dir/Counts/Reads/$FvHREADS
                samtools view -F 4 $dir/BAM/"$PREFIX".BAM | cut -f3 | grep -vE '(SO|PN|LN)' | grep -v "*" | grep -v "Homo_sapiens" | wc -l | awk -v where=`basename $PREFIX` '{print $1",fungal,"where}' >> $dir/Counts/Reads/$FvHREADS
                samtools view -F 4 $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep "Homo_sapiens" | wc -l | awk -v where=`basename $PREFIX` '{print $1",human,"where}' >> $dir/Counts/Reads/$FvHREADS
        else
                echo "Can't find fastq files!"
        fi
done
