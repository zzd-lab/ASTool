#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage="$0 [options]
--junction      <str>           the file of junction count
--gtf           <str>           gene annotation file
--m             <int>           minimum length between read and exon or intron(default:8)
--outdir        <str>           the path of outdir
        ";
my($junction,$gtf,$m,$outdir);
GetOptions("junction:s" =>\$junction,
	"gtf:s" =>\$gtf,
        "m:i" =>\$m,
        "outdir:s" =>\$outdir


);

die $usage if !defined $junction;

$m ||=8;







my $junction_file=$junction;
my $Annotation_file=$gtf;
my $min_overlap=$m;
my $psi_file=$outdir;
open(IN,"$junction_file") or die "can't open $junction_file";
my %junction_length_count;
my %junction_length_count_other;
while (<IN>) {
    chomp(my $line=$_);
    my ($ch,$junction,$length,$count,$count_other)=split /\t/,$line;
    unless ($length>0){
        next;
    }
    $junction_length_count{$ch}{$junction}="$length\t$count";
    $junction_length_count_other{$ch}{$junction}="$length\t$count_other";
}
close IN;

my $file_type= (split /\./,$Annotation_file)[-1];
if (($file_type eq "gtf")||($file_type eq "gff3") ) {
    my %versions;
    my %transcript_gene;
    my %exon_pos;
    my %gene_transcripts;
    my %ch_transcript;
    my %exon_gene;
    my %exons;
    if ($file_type eq "gtf") {
        open(IN,"$Annotation_file") or die "can't open $Annotation_file";
        while (<IN>) {
            chomp(my $line=$_);
            unless ($line=~m/^#/){
                $line=~s/;$//;
                my ($ch,$version,$tag,$start,$end,$strand)=(split /\t/,$line)[0,1,2,3,4,6];
                $versions{$version}++;
                if ($tag=~m/exon/gi) {
                    #$number=0;
                    my $exon=$ch."\t".$strand."\t".$start."\t".$end;
                    my $descrip=(split /\t/,$line)[8];
                    $descrip=~s/"//gi;
                    my @descript=split /; /,$descrip;
                     my ($transcript,$gene);
                    foreach my $des(@descript){
                        if ($des=~m/^transcript_id/gi) {
                            $transcript=(split / /,$des)[1];
                        }elsif($des=~m/gene_id/gi){
                            $gene=(split / /,$des)[1];
                        }
                    }
                    $gene_transcripts{$gene}{$transcript}++;
                    $transcript_gene{$transcript}=$gene;
                    $exons{$transcript}{$start}=$exon;
                    $exon=$start."_".$end;
                    my $pos1=int($start/1000);
                    my $pos2=int($end/1000);
                    for(my $pos=$pos1;$pos<=$pos2;$pos++){
                             $exon_pos{$ch}{$pos}{$exon}=$gene;
                    }
                    $exon=$ch."_".$start."_".$end;
                    $exon_gene{$exon}{$gene}++;
                }
            }
        }
        close IN;
    }else{
        open(IN,"$Annotation_file") or die "can't open $Annotation_file";
        while (<IN>) {
            chomp(my $line=$_);
            unless($line=~/^#/) {
                $line=~s/;$//;
                my ($ch,$version,$type,$start,$end,$strand,$descript)=(split /\t/,$line)[0,1,2,3,4,6,8];
                $versions{$version}++;
                if (($descript=~m/^ID/)&&($descript=~m/Parent/)) {
                    my ($transcript,$gene);
                    my @descript=split /;/,$descript;
                    foreach my $des(@descript){
                        if ($des=~m/^ID/gi) {
                            $transcript=(split /=|:/,$des)[-1];
                        }elsif($des=~m/Parent/gi){
                            $gene=(split /=|:/,$des)[-1];
                        }
                    }
                    if ($type=~m/gene|RNA|transcript/gi) {
                        $gene_transcripts{$gene}{$transcript}++;
                        $transcript_gene{$transcript}=$gene;
                    }
                }
                if ($type eq "exon") {
                     my $exon=$ch."\t".$strand."\t".$start."\t".$end;
                    my ($transcript,$gene);
                    my @descript=split /;/,$descript;
                    foreach my $des(@descript){
                        if ($des=~m/^Parent/gi) {
                            $transcript=(split /=|:/,$des)[-1];
                            
                        }
                    }
                    $gene=$transcript_gene{$transcript};
                    $gene_transcripts{$gene}{$transcript}++;
                    $exons{$transcript}{$start}=$exon;
                    $exon=$start."_".$end;
                    my $pos1=int($start/1000);
                    my $pos2=int($end/1000);
                    for(my $pos=$pos1;$pos<=$pos2;$pos++){
                             $exon_pos{$ch}{$pos}{$exon}=$gene;
                    }
                    $exon=$ch."_".$start."_".$end;
                    $exon_gene{$exon}{$gene}++;
                }   
            }
        }
        close IN;
    }
    my %exons2;
    my $end0;
    my %intron_numner;
    my %introns;
    my %intron_pos;
    foreach my $transcript(sort keys %exons){
        my $gene=$transcript_gene{$transcript};
        my $temp=0;
        foreach my $number(sort {$a<=>$b} keys %{$exons{$transcript}}){
            $temp+=1;
            my $exon=$exons{$transcript}{$number};
            $exons2{$transcript}{$temp}=$exon;
            my ($ch,$strand,$start,$end )=split /\t/,$exon;
            if ( $temp>1) {
                    my $intron_start=$end0+1;
                    my $intron_end=$start-1;
                    my $intron=$ch."\t".$strand."\t".$intron_start."\t".$intron_end;
                    $intron_numner{$transcript}++;
                    my $temp2=$temp-1;
                    $introns{$transcript}{$temp2}=$intron;
                    my $intron_length=$intron_end-$intron_start+1;
                    my $pos1 =int($intron_start/1000);
                    my $pos2 =int($intron_end/1000);
                    $intron=$intron_start."_".$intron_end;
                    for(my $pos=$pos1;$pos<=$pos2;$pos++){
                        $intron_pos{$ch}{$pos}{$intron}=$gene;
                    }
                }
            $end0=$end;
        }
    }
    open OUT, ">$psi_file";
    my %repeat_out;
    foreach my $transcript(sort keys %introns){
        my $intron_number;
        my $gene=$transcript_gene{$transcript};
        my @numbers=sort {$a<=>$b} keys %{$introns{$transcript}};
        my $total_number=$#numbers+1;
        foreach my $number(@numbers){
            my $intron=$introns{$transcript}{$number};
            my ($ch,$strand,$start,$end )=split /\t/,$intron;
            my $exon1=$exons2{$transcript}{$number};
            my ($exon1_start,$exon1_end )=(split /\t/,$exon1)[2,3];
            my $exon2=$exons2{$transcript}{$number+1};
            my ($exon2_start,$exon2_end )=(split /\t/,$exon2)[2,3];
            $intron=$start."_".$end;
            $exon1=$exon1_start."_".$exon1_end;
            $exon2=$exon2_start."_".$exon2_end;
            my $intron_length=$end-$start+1;
            my $exon1_length=$exon1_end-$exon1_start+1;
            my $exon2_length=$exon2_end-$exon2_start+1;
            unless(($intron_length>=$min_overlap)&&($exon1_length>=$min_overlap)&&($exon2_length>=$min_overlap)){
                next;  
            }
            if (exists $repeat_out{$gene}{$intron}) {
                next;
            }
            $repeat_out{$gene}{$intron}++;
            my $e1_i_junction=$exon1_start."_".$exon1_end."_".$start."_".$end;
            my $e1_e2_junction=$exon1_start."_".$exon1_end."_".$exon2_start."_".$exon2_end;
            my $i_e2_junction=$start."_".$end."_".$exon2_start."_".$exon2_end;
            my ($e1_i_count,$e1_e2_count,$i_e2_count,$e1_e2_count_other)=(0,0,0,0);
            my ($normalized_e1_i_count,$normalized_e1_e2_count,$normalized_i_e2_count,$normalized_e1_e2_count_other)=(0,0,0,0);
            if ($strand eq "-") {
                $intron_number=$total_number-$number+1;
            }else{
                $intron_number=$number;
            }
            my $type="Novel";
            foreach my $transcript2(keys %{$gene_transcripts{$gene}}){ 
                unless ($transcript eq $transcript2){
                    foreach my $number(sort {$a<=>$b} keys %{$exons2{$transcript2}}){
                        my $exon=$exons2{$transcript2}{$number};
                        my ($ch,$strand,$exon_start,$exon_end)=split /\t/,$exon;
                        if (($start>=$exon_start)&&($end<=$exon_end)) {
                            $type="Known";
                            last;
                        }
                    }
                }
            }
            my $clean1="Clean";
            my $pos1=int($start/1000);
            my $pos2=int($end/1000);
            my %exon_temp;
            for(my $pos=$pos1;$pos<=$pos2;$pos++){
                for my $exon(keys %{$exon_pos{$ch}{$pos}}){
                    if (exists $exon_temp{$exon}) {
                        next
                    }
                    my ($exon_start,$exon_end)=split /_/,$exon;
                    if (($start>=$exon_start)&&($end<=$exon_end)) {
                        if ($type eq "Novel") {
                            $clean1="Not_Clean";
                            last;
                        }  
                    } 
                    if(($exon_start >= $start) && ($exon_start< $end)){
                        $clean1="Not_Clean";
                        last;
                    }	
                    if(($exon_end>$start)&&($exon_end <= $end)){
                        $clean1="Not_Clean";
                        last;
                    }
                    # 1 bp overlap
                    if(($exon_start == $end) and ($exon_end>= $end)){
                        $clean1="Not_Clean";
                        last;
                    }
                     # 1 bp overlap
                    if(($exon_end == $start) and ($exon_start<= $start)){
                        $clean1="Not_Clean";
                        #print "$exon_start\t$exon_end\t$start\t$end\n";
                        last;
                    }
                    $exon_temp{$exon}++
                }
                if ( $clean1 eq "Not_Clean") {
                    last;
                }
            }
            if ($clean1 eq "Clean") {
                my %intron_temp;
                for(my $pos=$pos1;$pos<=$pos2;$pos++){
                    for my $intron(keys %{$intron_pos{$ch}{$pos}}){
                        if (exists $intron_temp{$intron}) {
                            next
                        }
                        my ($intron_start,$intron_end)=split /_/,$intron;
                        if(($intron_start== $start) && ($intron_end== $end)){
                            next
                        }
                        if(($intron_start< $start) && ($intron_end > $start)){
                            $clean1="Not_Clean";
                            next
                        }
                        if( ($intron_start> $start) &&  ($intron_start < $end)){
                            $clean1="Not_Clean";
                            next
                        } 
                        if(($intron_start== $start) || ($intron_end== $end)|| ($intron_end== $start)||($intron_start== $end)){
                            $clean1="Not_Clean";
                            next
                        }
                                  
                        $intron_temp{$intron}++;
                    }
                    if ( $clean1 eq "Not_Clean") {
                        last;
                    }
                }
            }
            my $length;
            if (exists $junction_length_count{$ch}{$e1_i_junction}) {
                my $length_count=$junction_length_count{$ch}{$e1_i_junction};
                ($length,$e1_i_count)=split /\t/,$length_count;
                $normalized_e1_i_count=$e1_i_count/$length;
            }
            if (exists $junction_length_count{$ch}{$e1_e2_junction}) {
                my $length_count=$junction_length_count{$ch}{$e1_e2_junction};;
                ($length,$e1_e2_count)=split /\t/,$length_count;
                $normalized_e1_e2_count=$e1_e2_count/$length;
            }
            if (exists $junction_length_count{$ch}{$i_e2_junction}) {
                my $length_count=$junction_length_count{$ch}{$i_e2_junction};
                ($length,$i_e2_count)=split /\t/,$length_count;
                $normalized_i_e2_count=$i_e2_count/$length;
            }
            if (exists $junction_length_count_other{$ch}{$e1_e2_junction}) {
                my $length_count=$junction_length_count_other{$ch}{$e1_e2_junction};;
                ($length,$e1_e2_count_other)=split /\t/,$length_count;
            }
            my $lowCount="Low_Count";
            my $imbalance="-";
            if ((($e1_i_count>0)&&($i_e2_count==0))||(($e1_i_count==0)&&($i_e2_count>0))) {
               $imbalance="Imbalance";
            }
            if (min($e1_i_count,$i_e2_count)>0) {
                my $fold=max($e1_i_count,$i_e2_count)/min($e1_i_count,$i_e2_count);
                if ($type eq "Novel") {
                    if (($e1_e2_count<=1)||($fold>3)) {
                        $imbalance="Imbalance";
                    }   
                }
                    if (($e1_e2_count_other>2)&&(($e1_e2_count/$e1_e2_count_other)<10)) {
                        $imbalance="NotClean";
                    }
            }
            if ((mean($e1_i_count,$i_e2_count)+$e1_e2_count)>0) {
                if ((mean($e1_i_count,$i_e2_count)+$e1_e2_count)>=10) {
                    $lowCount="-";
                }
                my $psi=mean($normalized_e1_i_count,$normalized_i_e2_count)/(mean($normalized_e1_i_count,$normalized_i_e2_count)+$normalized_e1_e2_count);
                print OUT "$gene\t$ch\t$intron\t$exon1\t$exon2\t$type\t$clean1\t$e1_i_count\t$i_e2_count\t$e1_e2_count\t$lowCount\t$imbalance\t$psi\n";
            }else{
                print OUT "$gene\t$ch\t$intron\t$exon1\t$exon2\t$type\t$clean1\t$e1_i_count\t$i_e2_count\t$e1_e2_count\t$lowCount\t$imbalance\t-\n";
            } 
        }
    }
    close OUT;
    
}else{
    print "Notice: The format of the input file should be gtf or gff3!\n";
}



sub min{
    my @numbers=@_;
    my $min=$numbers[0];
    foreach my $number (@numbers){
        if ($number < $min) {
            $min=$number
        }
    }
    return $min;
}

sub max{
    my @numbers=@_;
    my $max=$numbers[0];
    foreach my $number (@numbers){
        if ($number > $max) {
            $max=$number
        }
    }
    return $max;
}
sub mean{
    my @numbers=@_;
    my $num=$#numbers+1;
    my $sum=0;
    foreach my $value(@numbers){
        $sum+=$value
    }
    my $mean_value=$sum/$num;
    return $mean_value;
}

