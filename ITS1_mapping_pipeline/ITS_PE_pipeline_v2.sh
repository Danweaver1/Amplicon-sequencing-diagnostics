workingdir=$( pwd )
#NB - folder 'database' containing bowtie2 indexes (create using bowtie2-build <db file> <db_name>) is required in working directory
#set val for bowtie stringency
#set db for bowtie indexes
analysisdir="ITS_PE_analysis"

##bbduk QC
for file in *R1_001.fastq.gz
do
READ2=$( echo $file | sed 's/R1/R2/' )
Prefix=$( echo $file | cut -d"_" -f 1,4 )
R2pre=$( echo $Prefix | sed 's/R1/R2/' )
bbduk.sh in1=$file in2=$READ2 out1=BB_"$Prefix".fastq out2=BB_"$R2pre".fastq qtrim=r trimq=30 minlen=75 2>"$Prefix"_bbduk.out
done

mkdir $analysisdir
mv BB* $analysisdir
mkdir $analysisdir/bbduk_info
mv *bbduk.out $analysisdir/bbduk_info
cd $analysisdir

mkdir Primer_untrimmed
mkdir Primer_QC

#remove primers
for file in BB*R1.fastq
do
  READ2=$( echo $file | sed 's/R1/R2/' )
  Prefix=$( echo $file | sed 's/BB_//' | sed 's/.fastq//' )
  R2pre=$( echo $Prefix | sed 's/R1/R2/' )
  cutadapt -q 10 --minimum-length 75 --untrimmed-output $Prefix.Primer_untrimmed.fastq --untrimmed-paired-output $R2pre.Primer_untrimmed.fastq --overlap 6 -g ^TCCGTAGGTGAACCTGCGG -G ^GCYRCGTTCTTCATCGAYGC -o "$Prefix"_QC.fastq -p "$R2pre"_QC.fastq $file $READ2 1>"$Prefix"_cutadapt.out
done

mkdir cutadapt_info
mv *_cutadapt.out cutadapt_info

#file_name_key should be a csv (comma seperated) with column 1 containing raw fastq prefix (ie. filename text before first _) and col2 containing new name to change to
while read p
do
    name=$( echo $p | cut -d"," -f1)
    newname=$( echo $p | cut -d"," -f2)
    if [ -e "$name"_R1_QC.fastq ]
    then
    mv "$name"_R1_QC.fastq "$newname"_R1_QC.fastq &&
    mv "$name"_R2_QC.fastq "$newname"_R2_QC.fastq &&
    echo "$name" renamed to "$newname" >> file_rename_info
    else
      echo "Sample "$name" not found"
    fi
done <../file_name_key.csv


mv *QC.fastq Primer_QC
mv *untrimmed.fastq Primer_untrimmed

cd Primer_QC
########bowtie2 alignment

source activate samtools1_9
#above line necessary on desktop

val=0.06
db="v2_ISHAM_ITS"
dir="$db"_"$val"

mkdir $dir
mkdir $dir/unaligned
mkdir $dir/BAM
mkdir $dir/Counts
mkdir $dir/Counts/Seqs
mkdir $dir/Counts/Species
mkdir $dir/Counts/Genera
mkdir $dir/Counts/Reads

for file in *R1_QC.fastq
do
        if [ -e $file ]
        then
		            PREFIX=$( echo $file | sed 's/_R1_QC.fastq//' )
                READ2=$( echo $file | sed 's/R1/R2/' )
		            SEQ="$PREFIX.seq_count"
		            UNALIGNED=""$PREFIX"_unaligned.fastq"
                SPECIES=""$PREFIX"_unique_species"
                GENUS=""$PREFIX"_unique_genera"
                READS="$PREFIX"_reads.csv
                FvHREADS="$PREFIX"_FvHreads.csv

                echo $PREFIX >> $dir/PE_"$db"_"$val"_bowtie.out

		            bowtie2 --end-to-end -D 25 -R 5 -N 0 -L 25 -i S,1,2.5 --rdg 3,3 --score-min L,-0.6,-$val --un-conc $dir/unaligned/$UNALIGNED -1 $file -2 $READ2 -x $workingdir/database/$db 2>>$dir/PE_"$db"_"$val"_bowtie.out | samtools sort -o $dir/BAM/"$PREFIX".BAM &&
                samtools view -f 2 $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep -v "*" | sort | uniq -c | sort -nr >> $dir/Counts/Seqs/$SEQ
                samtools view -f 2 $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)'  | grep -v "*" | sort | cut -d"_" -f2-3 | sort | uniq -c | sort -nr > $dir/Counts/Species/$SPECIES
                cat $dir/Counts/Species/$SPECIES | sed -e 's/^ *//;s/ /,/' | awk -v where=`basename $PREFIX` '{print $1","where}' > $dir/Counts/Species/$PREFIX.Rspecies_counts.csv
                #create genera counts and csv compatible with R plotting script
                samtools view -f 2 $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep -v "*" | sort | cut -d"_" -f2 | sort | uniq -c | sort -nr > $dir/Counts/Genera/$GENUS &&
                cat $dir/Counts/Genera/$GENUS | sed -e 's/^ *//;s/ /,/' | awk -v where=`basename $PREFIX` '{print $1","where}' > $dir/Counts/Genera/$PREFIX.Rgenera_counts.csv
                #generate overview of unaligned and aligned
                samtools view -F 2 $dir/BAM/"$PREFIX".BAM | grep -vE '(SO|PN|LN)' | wc -l | awk -v where=`basename $PREFIX` '{print $1",unaligned,"where}' >> $dir/Counts/Reads/$READS
                samtools view -f 2 $dir/BAM/"$PREFIX".BAM | cut -f3 | grep -vE '(SO|PN|LN)' | wc -l | awk -v where=`basename $PREFIX` '{print $1",aligned,"where}' >> $dir/Counts/Reads/$READS

                ##generate overview of unaligned, human and fungal read numbers
                samtools view -F 2 $dir/BAM/"$PREFIX".BAM | grep -vE '(SO|PN|LN)' | wc -l | awk -v where=`basename $PREFIX` '{print $1",unaligned,"where}' >> $dir/Counts/Reads/$FvHREADS
                samtools view -f 2 $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep -v "*" | grep -v "Homo_sapiens" | wc -l | awk -v where=`basename $PREFIX` '{print $1",fungal,"where}' >> $dir/Counts/Reads/$FvHREADS
                samtools view -f 2 $dir/BAM/"$PREFIX".BAM | cut -f3 | sort | grep -vE '(SO|PN|LN)' | grep -v "*" | grep "Homo_sapiens" | wc -l | awk -v where=`basename $PREFIX` '{print $1",human,"where}' >> $dir/Counts/Reads/$FvHREADS
        else
                echo "Can't find fastq files!"
        fi
done

#########
