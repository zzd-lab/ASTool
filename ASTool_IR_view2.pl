#!/usr/bin/perl -w
use strict;

use Getopt::Long;
my $usage="$0 [options]
--gtf           <str>           gene annotation file
--bam		<str>		BAM file
--intron	<str>		The format of interest intron(intron.txt):AT1G01520_1_160303_160417
--outdir	<str>		the path of outdir
--junction      <str>           the file of junction count
--psi		<str>		the file of IR PSI
--no_junction	<str>		Junction count with no reference comments
        ";
my($gtf,$bam,$intron,$outdir,$junction,$psi,$no_junction);
GetOptions("gtf:s" =>\$gtf,
        "bam:s" =>\$bam,
        "intron:s" =>\$intron,
	"outdir:s" =>\$outdir,
	"junction:s" =>\$junction,
	"psi:s"=>\$psi,
	"no_junction:s" =>\$no_junction


);

die $usage if !defined $gtf;



my $gtf_file=$gtf;
my $sam_file=$bam;
my $intron_file=$intron;
my $out_dir=$outdir;
my $junction_count_file=$junction;
my $psi_file=$psi;
my $junction_no_ref=$no_junction;
if (-e $out_dir) {
	`rm -rf $out_dir`;
    `mkdir $out_dir`;
}else{
	`mkdir $out_dir`;
}
my $dir = `pwd`;
chomp($dir);
$out_dir=$dir."\/".$out_dir;

`dos2unix $intron_file`;
open(IN,"$intron_file") or die "can't open $intron_file";
while (<IN>) {
    chomp(my $line=$_);
    my ($gene,$ch,$intron_start,$intron_end)=split /_/,$line;
    my $intron=$ch."_".$intron_start."_".$intron_end;
    my $edge_file=$out_dir."\/"."edge_".$line.".txt";
    `perl junction_edge2.pl $junction_count_file $psi_file $intron $edge_file $junction_no_ref`;
    my $gene_bam=$out_dir."\/".$gene.".bam";
    my $gene_sam=$out_dir."\/".$gene.".sam";
	my $gene_sort_bam=$out_dir."\/".$gene.".sorted.bam";
	$intron=$gene."_".$ch."_".$intron_start."_".$intron_end;
	`perl specific_region_sam.pl $sam_file $intron $gene_sam`;
	`samtools view -S $gene_sam -b > $gene_bam`;
	`samtools sort $gene_bam -o $gene_sort_bam`;
	`samtools index $gene_sort_bam`;
	my $out_tmp_depth=$out_dir."\/".$gene."_depth.txt";
	`perl sam2depth.pl $gene_sam $out_tmp_depth`;
	my $out_pdf=$out_dir."\/".$line.".pdf";
	`Rscript plot_view.r $edge_file $out_tmp_depth $out_pdf`
}
close IN;

