# ASTool
ASTool: Accurate Identification of Alternative Splicing Events from Plant RNA-Seq Data<br><br>
Step1: Download and installation:<br>
Dependence:<br>
```
Perl,threads(Perl module),SRA Toolkit, STAR
```
Download:
```
wget http://zzdlab.com/ASTool/ASTool_v1/ASTool.zip

Genome FASTA and gene annotation file can be downloaded from Ensemble Plants: https://plants.ensembl.org/index.html

For example:

wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget http://ftp.ensemblgenomes.org/pub/plants/release-52/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz
```
Installation:<br>
```
unzip ASTools_v1.zip
cd ASTools_v1
```
Step2: Preparing SAM file for ASTool<br>
1.Raw RNA-Seq data is downloaded from NCBI SRA and .sra file is converted to .fastq file using fastq-dump of SRA Toolkit<br>
```
e.g.
fastq-dump SRR4048211.sra
```
2. Adaptor removal and quanlity filtering are performed by Trimmomatic<br>
```
java -jar trimmomatic-0.33.jar SE -phred33 -threads 10 SRR4048211.fastq SRR4048211_trim.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36&>SRR4048211_trim_log.txt
```
3. Sequencing reads are mapped to reference genome by STAR<br>
```
1. Build index

mkdir genome_index_length100
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir genome_index_length100 --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --sjdbGTFfile Arabidopsis_thaliana.TAIR10.52.gtf --sjdbOverhang 99

2. Alignment

STAR --runThreadN 20 --genomeDir genome_index_length100/ --outFileNamePrefix SRR4048211. --readFilesIn SRR4048211_trim.fasta -outSJfilterReads Unique --outFilterMismatchNmax 4 --quantMode GeneCounts
```
Step3: Running ASTool<br>
Calculating PSI of four events: IR, ES,A5SS and A3SS.<br>
```
1 Calculate the junction count (Reference genome notes are available)
perl junction_count.pl --gtf [gene annotation file] --sam [SAM file] --thread [thread] --readlength [readlength] --m [m value] --outdir [outdir]
e.g.
perl junction_count.pl --gtf Arabidopsis_thaliana.TAIR10.51.gtf --sam SRR4048211.Aligned.out.sam --thread 10 --readlength 100 --m 8 --outdir SRR4048211_junction_ref

2. Calculate the PSI
2.1 IR events:
perl IR_PSI.pl --junction [the file of junction count] --gtf [gene annotation file] --m [m value] --outdir [outdir]
(m value: minimum length between read and exon or intron,the same below)
e.g.
perl IR_PSI.pl --junction ./SRR4048211_junction_ref/junction_count.txt --gtf Arabidopsis_thaliana.TAIR10.52.gtf --m 8 --outdir SRR4048211_IR_PSI.txt

2.2 ES events:
perl ES_PSI.pl --junction [the file of junction count] --gtf [gene annotation file] --m [m value] --outdir [outdir]
e.g.
perl ES_PSI.pl SRR4048211_junction_ref/junction_count.txt Arabidopsis_thaliana.TAIR10.52.gtf 8 SRR4048211_ES_PSI.txt

2.3 A5SS/A3SS events:
perl A5SS_A3SS_PSI.pl --junction [the file of junction count] --gtf [gene annotation] --m [m value] --A5SS_outdir [A5SS events outdir] --A3SS_outdir [A3SS events outdir]
e.g.
perl A5SS_A3SS_PSI.pl --junction SRR4048211_junction_ref/junction_count.txt --gtf Arabidopsis_thaliana.TAIR10.52.gtf --m 8 --A5SS_outdir SRR4048211_A5SS_PSI.txt --A3SS_outdir SRR4048211_A3SS_PSI.txt
Main outfile
1. junction_count.txt
Column1: ch

Column2: position of junction

Column3: effect length of junction

Column3: reads mapping to junction
2. IR_PSI.txt
Column1: Gene

Column2: ch

Column3: position of intron (I)

Column4: position of flanking exon 1(E1)

Column5: position of flanking exon 2(E2)

Column6: intron type ("Known" or "Unknown")

Column7: intron type ("Clean" of "Not Clean")

Column8: reads mapping to junction E1_I

Column9: reads mapping to junction I_E2

Column10: reads mapping to junction E1_E2

Column11: warning 1 ("Low Count" or "-")

Column12: warning 2 ("Imblance" or "-")

Column13: PSI
2. ES_PSI.txt
Column1: Gene

Column2: ch

Column3: position of exon (E)

Column4: position of flanking intron 1 (intron between exon E1 and E)

Column5: position of flanking intron 2 (intron between exon E1 and E)

Column6: intron type ("Known" or "Unknown")

Column7: exon type ("Clean" of "Not Clean")

Column8: reads mapping to junction E1_E

Column9: reads mapping to junction E_E2

Column10: reads mapping to junction E1_E2

Column11: warning 1 ("Low Count" or "-")

Column12: warning 2 ("Imblance" or "-")

Column13: PSI
3. A5SS_PSI.txt and A3SS_PSI.txt
Column1: Gene

Column2: ch

Column3: strand

Column4: shoter intron

Column5: longer intron

Column6: -

Column7: -

Column8: reads mapping to junction span shorter intron

Column9: eads mapping to junction span longer intron

Column10: -

Column11: warning 1 ("Low Count" or "-")

Column12: -

Column13: PSI
```
Visualization<br>
Visualize introns of interest<br>
```
4 Intron visualization
4.1.1 Calculate the junction count (Reference genome notes are available) See section 1.1 for details

4.1.2 Count junction count without reference comments
perl junction_count_no_ref.pl --sam [SAM file] --thread [thread] --m [m value] --outdir [outdir]
e.g.
perl junction_count_no_ref.pl --sam SRR4048211.Aligned.out.sam --thread 10 --m 8 --outdir SRR4048211_junction_no_ref

4.2 Visualization
perl ASTools_IR_view2.pl --gtf [gene annotation file] --bam [BAM file] --intron [Introns of interest(.txt)] --outdir [outdir] --junction [Junction count with reference comments] --psi [the file of IR PSI] --no_junction [Junction count with no reference comments]
e.g.
The format of interest intron(intron.txt):AT1G01520_1_160303_160417
perl ASTools_IR_view2.pl --gtf Arabidopsis_thaliana.TAIR10.52.gtf --bam SRR5197909.Aligned.out.bam --intron intron.txt --outdir ./SRR4048211_plot --junction ./SRR4048211_junction_ref/junction_count.txt --psi SRR4048211_IR_PSI.txt --no_junction SRR4048211_junction_no_ref/junction_count.txt
```
For more details, see [ASTool](http://zzdlab.com/ASTool/manual.php)
