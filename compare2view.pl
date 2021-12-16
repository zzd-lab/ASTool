#!/usr/bin/perl -w
use strict;
open(IN, $ARGV[0]) or die "can't open $ARGV[0]";
my @para;
my $para_i=0;
while (<IN>) {
    chomp (my $line=$_);
    if ($line=~m/^#/) {
        next
    }
    $para[$para_i]=(split /\s+/,$line)[1];
    $para_i+=1
}
my ($gtf_file,$sam_files,$psi_file,$out_dir,$sample_name,$sample_group)=@para;
my @sample_names=split /,/,$sample_name;
my @psi_file=split /,/,$psi_file;

close IN;
my $intron=$ARGV[1];
my $dir = `pwd`;
chomp($dir);
$out_dir=$dir."\/".$out_dir;
if (-e $out_dir) {
    `rm -rf $out_dir`;
	`mkdir $out_dir`;
}else{
    `mkdir $out_dir`;
}

my %gene_pos;
open(IN,"$gtf_file") or die "can,t open $gtf_file";
while (<IN>) {
    chomp(my $line=$_);
    unless ($line=~m/^#/){
        my ($tag,$start,$end,$descript)=(split /\t/,$line)[2,3,4,8];
        if ($tag eq "gene") {
            $descript=~s/"//gi;
            my $gene=(split /; /,$descript)[0];
            $gene=(split / /,$gene)[1];
            my $pos=$start."_".$end;
            $gene_pos{$gene}=$pos;
        } 
    }

}
close IN;


my $intron_pos;
my @sam_files=split /,/,$sam_files;
my @depth_file;
my @psi;
my $i=0;
    my ($gene,$ch,$intron_start,$intron_end)=split /_/,$intron;
    my $pos=$gene_pos{$gene};
    my ($start,$end)=split /_/,$pos;
    foreach my $sam_file(@sam_files){
        my $gene_sort_bam=$out_dir."\/".$gene.".sorted.bam".$i;
		my $gene_bam=$out_dir."\/".$gene.".bam".$i;
		my $gene_sam=$out_dir."\/".$gene.".sam".$i;
		`perl specific_region_sam.pl $sam_file $intron $gene_sam`;
		`samtools view -S $gene_sam -b > $gene_bam`;
		`samtools sort $gene_bam -o $gene_sort_bam`;
        `samtools index $gene_sort_bam`;
         my $out_tmp_depth=$out_dir."\/".$gene."_".$sample_names[$i]."_depth.txt";
        `perl sam2depth.pl $gene_sam $out_tmp_depth`;
        `rm $gene_sam`;
        `rm $gene_sort_bam`;
        my $out_tmp_bai=$gene_sort_bam.".bai";
        `rm $out_tmp_bai`;
        $depth_file[$i]=$out_tmp_depth;
        open(IN,"$psi_file[$i]") or die "can't open psi_file[$i]";
        #1_704495_704587
        while (<IN>) {
            chomp (my $line=$_);
            my ($gene,$ch,$intron_start_end,$exon1,$exon2,$psi)=(split /\t/, $line)[0,1,2,3,4,12];
            unless ($psi eq "-"){
                $psi=sprintf "%.2f",$psi;
            }else {
                $psi="NA"
            }
            my $intron2=$gene."_".$ch."_".$intron_start_end;
            if ($intron2 eq $intron) {
                $intron_pos=$exon1."_".$exon2;
                $psi[$i]=$psi;
                next;
            }
        }
        close IN;
        $i+=1;
    }
close IN;
my $gene_depth_file= join ",",@depth_file;
my $out_pdf=$out_dir."\/".$intron.".pdf";
my $psi= join ",",@psi;
$intron=~s/_/,/;
$intron=~s/_/:,/;
$intron=~s/_/-/;
`Rscript compare2view.r $gene_depth_file $psi $sample_group $out_pdf $intron_pos $intron $sample_name`;

