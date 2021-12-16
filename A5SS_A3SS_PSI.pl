  #!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $usage="$0 [options]
--junction      <str>           the file of junction count
--gtf           <str>           gene annotation file
--m             <int>           minimum length between read and exon or intron(default:8)
--A5SS_outdir   <str>           the path of A5SS outdir
--A3SS_outdir	<str>		the path of A3SS outdir
        ";
my($junction,$gtf,$m,$A5SS_outdir,$A3SS_outdir);
GetOptions("junction:s" =>\$junction,
        "gtf:s" =>\$gtf,
        "m:i" =>\$m,
        "A5SS_outdir:s" =>\$A5SS_outdir,
	"A3SS_outdir:s" =>\$A3SS_outdir


);

die $usage if !defined $junction;

$m ||=8;




my $junction_file=$junction;
my $Annotation_file=$gtf;
my $min_overlap=$m;
my $A5SS_psi_file=$A5SS_outdir;
my $A3SS_psi_file=$A3SS_outdir;
#junction_count
open(IN,"$junction_file") or die "can't open $junction_file";
my %junction_length_count;
while (<IN>) {
    chomp(my $line=$_);
    my ($ch,$junction,$length,$count)=split /\t/,$line;
    unless ($length>0){
        next;
    }
    $junction_length_count{$ch}{$junction}="$length\t$count";
}
close IN;

my $file_type= (split /\./,$Annotation_file)[-1];
if (($file_type eq "gtf")||($file_type eq "gff3") ) {
    my %versions;
    my %exon_site;
    my %exon_strand;
    my %exon_ch;
    my %transcript_gene;
    my %gene_transcripts;
    my %exon_pos;
    my %intron_pos;
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
                    my $exon=$ch."\t".$strand."\t".$start."\t".$end;
                    $gene_transcripts{$gene}{$transcript}++;
                    $exon_site{$transcript}{$exon}++;
                    $exon_strand{$transcript}=$strand;
                    $exon_ch{$transcript}=$ch;
                    my $pos1=int($start/1000);
                    my $pos2=int($end/1000);
                    for(my $pos=$pos1;$pos<=$pos2;$pos++){
                        $exon_pos{$ch}{$pos}{$exon}=$gene;
                    }
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
                    my ($transcript,$gene);
                    my @descript=split /;/,$descript;
                    foreach my $des(@descript){
                        if ($des=~m/^Parent/gi) {
                            $transcript=(split /=|:/,$des)[-1];
                            
                        }
                    }
                    $gene=$transcript_gene{$transcript};
                    my $exon=$ch."\t".$strand."\t".$start."\t".$end;
                    $gene_transcripts{$gene}{$transcript}++;
                    $exon_site{$transcript}{$exon}++;
                    $exon_strand{$transcript}=$strand;
                    $exon_ch{$transcript}=$ch;
                    my $pos1=int($start/1000);
                    my $pos2=int($end/1000);
                    for(my $pos=$pos1;$pos<=$pos2;$pos++){
                             $exon_pos{$ch}{$pos}{$exon}=$gene;
                    }
                }   
            }
        }
        close IN;
    }
    
    
    
    my %A5SS_pair;
    foreach my $gene(keys %gene_transcripts){
        my @transcripts=keys %{$gene_transcripts{$gene}};
        my $tr_number=$#transcripts+1;
        if ($tr_number>1) {
            for (my $i=0;$i<$#transcripts;$i++){
                for (my $j=$i+1;$j<=$#transcripts;$j++){
                    my $tr1=$transcripts[$i];
                    my $tr2=$transcripts[$j];
                    my @exon_sites1;
                    my $exon_site_index1=-1;
                    foreach my $exon1(keys %{$exon_site{$tr1}}){
                        my ($site1,$site2)=(split /\t/,$exon1)[2,3];
                        $exon_site_index1+=1;
                        $exon_sites1[$exon_site_index1]=$site1;
                        $exon_site_index1+=1;
                        $exon_sites1[$exon_site_index1]=$site2;   
                    }
                    my @sites1=sort { $a <=> $b } @exon_sites1;
                    my $site_numbers1=$#sites1+1;
                    my @exon_sites2;
                    my $exon_site_index2=-1;
                    foreach my $exon2(keys %{$exon_site{$tr2}}){
                        my ($site1,$site2)=(split /\t/,$exon2)[2,3];
                        $exon_site_index2+=1;
                        $exon_sites2[$exon_site_index2]=$site1;
                        $exon_site_index2+=1;
                        $exon_sites2[$exon_site_index2]=$site2;   
                    }
                    my @sites2=sort { $a <=> $b } @exon_sites2;
                    my $site_numbers2=$#sites2+1;
                    if (($site_numbers1>=4)&&($site_numbers2>=4)) {
                        for(my $i=0;$i<=$site_numbers1-4;$i+=2){
                            for(my $j=0;$j<=$site_numbers2-4;$j+=2){
                                my $strand=$exon_strand{$tr1};
                                my $ch=$exon_ch{$tr1};
                                my $ch_gene=$gene."\t".$ch."\t".$strand;
                                if ($strand eq "+") {
                                    if (($sites1[$i+2]==$sites2[$j+2])&&($sites1[$i+1]!=$sites2[$j+1])) {
                                        my $a5ss;
                                        if ($sites1[$i+1]<$sites2[$j+1]) {
                                            if ($sites2[$j]<=$sites1[$i+1]) {
                                                $a5ss=$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3]."_".$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3];
                                                $A5SS_pair{$ch_gene}{$a5ss}++;   
                                            }
                                        }else{
                                            if ($sites1[$i]<=$sites2[$j+1]) {
                                                $a5ss=$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3]."_".$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3];
                                                $A5SS_pair{$ch_gene}{$a5ss}++;   
                                            }
                                        }
                                    } 
                                }else{
                                    if (($sites1[$i+2]!=$sites2[$j+2])&&($sites1[$i+1]==$sites2[$j+1])) {
                                        my $a5ss;
                                        if ($sites1[$i+2]<$sites2[$j+2]) {
                                           if($sites1[$i+3]>=$sites2[$j+2]){
                                                $a5ss=$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3]."_".$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3];
                                                $A5SS_pair{$ch_gene}{$a5ss}++;
                                           }
                                        }else{
                                            if ($sites2[$j+3]>=$sites1[$i+2]) {
                                                $a5ss=$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3]."_".$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3];
                                                $A5SS_pair{$ch_gene}{$a5ss}++;
                                            }
                                        }
                                    }
                                }  
                            }
                        }
                    } 
                }
            }
        }
    }
    
    
    
    
    my %a3ss_pair;
    foreach my $gene(keys %gene_transcripts){
        my @transcripts=keys %{$gene_transcripts{$gene}};
        my $tr_number=$#transcripts+1;
        if ($tr_number>1) {
            for (my $i=0;$i<$#transcripts;$i++){
                for (my $j=$i+1;$j<=$#transcripts;$j++){
                    my $tr1=$transcripts[$i];
                    my $tr2=$transcripts[$j];
                    my @exon_sites1;
                    my $exon_site_index1=-1;
                    foreach my $exon1(keys %{$exon_site{$tr1}}){
                        my ($site1,$site2)=(split /\t/,$exon1)[2,3];
                        $exon_site_index1+=1;
                        $exon_sites1[$exon_site_index1]=$site1;
                        $exon_site_index1+=1;
                        $exon_sites1[$exon_site_index1]=$site2;   
                    }
                    my @sites1=sort { $a <=> $b } @exon_sites1;
                    my $site_numbers1=$#sites1+1;
                    my @exon_sites2;
                    my $exon_site_index2=-1;
                    foreach my $exon2(keys %{$exon_site{$tr2}}){
                        my ($site1,$site2)=(split /\t/,$exon2)[2,3];
                        $exon_site_index2+=1;
                        $exon_sites2[$exon_site_index2]=$site1;
                        $exon_site_index2+=1;
                        $exon_sites2[$exon_site_index2]=$site2;   
                    }
                    my @sites2=sort { $a <=> $b } @exon_sites2;
                    my $site_numbers2=$#sites2+1;
                    if (($site_numbers1>=4)&&($site_numbers2>=4)) {
                        for(my $i=0;$i<=$site_numbers1-4;$i+=2){
                            for(my $j=0;$j<=$site_numbers2-4;$j+=2){
                                my $strand=$exon_strand{$tr1};
                                my $ch=$exon_ch{$tr1};
                                my $ch_gene=$gene."\t".$ch."\t".$strand;
                                if ($strand eq "-") {
                                    if (($sites1[$i+2]==$sites2[$j+2])&&($sites1[$i+1]!=$sites2[$j+1])) { 
                                        my $a3ss;
                                        if ($sites1[$i+1]<$sites2[$j+1]) {
                                            if ($sites2[$j]<=$sites1[$i+1]) {
                                                $a3ss=$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3]."_".$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3];
                                                $a3ss_pair{$ch_gene}{$a3ss}++;   
                                            }
                                        }else{
                                            if ($sites1[$i]<=$sites2[$j+1]) {
                                                $a3ss=$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3]."_".$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3];
                                                $a3ss_pair{$ch_gene}{$a3ss}++;   
                                            }
                                        }
                                    } 
                                }else{
                                    if (($sites1[$i+2]!=$sites2[$j+2])&&($sites1[$i+1]==$sites2[$j+1])) {
                                        my $a3ss;
                                        if ($sites1[$i+2]<$sites2[$j+2]) {
                                            if($sites1[$i+3]>=$sites2[$j+2]){
                                                $a3ss=$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3]."_".$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3];
                                                $a3ss_pair{$ch_gene}{$a3ss}++;
                                            }
                                        }else{
                                            if ($sites2[$j+3]>=$sites1[$i+2]) {
                                                $a3ss=$sites1[$i]."_".$sites1[$i+1]."_".$sites1[$i+2]."_".$sites1[$i+3]."_".$sites2[$j]."_".$sites2[$j+1]."_".$sites2[$j+2]."_".$sites2[$j+3];
                                                $a3ss_pair{$ch_gene}{$a3ss}++;
                                            }
                                        }
                                    }
                                }  
                            }
                        }
                    }   
                }
            }
        }
    }

    my %repeat_out;
    open OUT,">$A5SS_psi_file";
    foreach my $ch_gene(keys %A5SS_pair){
        my ($gene,$ch,$strand)=split /\t/,$ch_gene;
        foreach my $a5ss(keys %{$A5SS_pair{$ch_gene}}){
            my ($ss1,$ss2,$ss3,$ss4,$ss5,$ss6,$ss7,$ss8)=split /_/,$a5ss;
            my $s_junction=$ss1."_".$ss2."_".$ss3."_".$ss4;
            my $i_junction=$ss5."_".$ss6."_".$ss7."_".$ss8;
            my ($s_count,$i_count)=(0,0);
            my ($normalized_s_count,$normalized_i_count)=(0,0);
            my $exon1_length=$ss2-$ss1+1;
            my $exon2_length=$ss4-$ss3+1;
            my $exon3_length=$ss6-$ss5+1;
            my $exon4_length=$ss8-$ss7+1;
            unless(($exon1_length>=$min_overlap)&&($exon2_length>=$min_overlap)&&($exon3_length>=$min_overlap)&&($exon4_length>=$min_overlap)){
                next;  
            }
            my $temp=$ss2."_".$ss3."_".$ss6."_".$ss7;
            if (exists $repeat_out{$ch_gene}{$temp}) {
                next;  
            }
            $repeat_out{$ch_gene}{$temp}++;
            my $lowCount="Low_Count";
            
            my $length;
            if (exists $junction_length_count{$ch}{$i_junction}) {
                my $length_count=$junction_length_count{$ch}{$i_junction};
                ($length,$i_count)=split /\t/,$length_count;
                $normalized_i_count=$i_count/$length;
            }
            if (exists $junction_length_count{$ch}{$s_junction}) {
                my $length_count=$junction_length_count{$ch}{$s_junction};
                ($length,$s_count)=split /\t/,$length_count;
                $normalized_s_count=$s_count/$length;
            }
            if (($i_count+$s_count)>0) {
                if (($i_count+$s_count)>=10) {
                    $lowCount="-";
                }
                my $psi=$i_count/($i_count+$s_count);
                print OUT "$gene\t$ch\t$strand\t$ss6\_$ss7\t$ss2\_$ss3\t-\t-\t$i_count\t$s_count\t-\t$lowCount\t-\t$psi\n";
            }else{
                print OUT "$gene\t$ch\t$strand\t$ss6\_$ss7\t$ss2\_$ss3\t-\t-\t$i_count\t$s_count\t-\t$lowCount\t-\t-\n";
            } 
            
            
        }
        
    }
    close OUT;
   
    %repeat_out=();
    open OUT,">$A3SS_psi_file";
    foreach my $ch_gene(keys %a3ss_pair){
        my ($gene,$ch,$strand)=split /\t/,$ch_gene;
        foreach my $a3ss(keys %{$a3ss_pair{$ch_gene}}){
            my ($ss1,$ss2,$ss3,$ss4,$ss5,$ss6,$ss7,$ss8)=split /_/,$a3ss;
            my $s_junction=$ss1."_".$ss2."_".$ss3."_".$ss4;
            my $i_junction=$ss5."_".$ss6."_".$ss7."_".$ss8;
            my ($s_count,$i_count)=(0,0);
            my ($normalized_s_count,$normalized_i_count)=(0,0);
            my $exon1_length=$ss2-$ss1+1;
            my $exon2_length=$ss4-$ss3+1;
            my $exon3_length=$ss6-$ss5+1;
            my $exon4_length=$ss8-$ss7+1;
            unless(($exon1_length>=$min_overlap)&&($exon2_length>=$min_overlap)&&($exon3_length>=$min_overlap)&&($exon4_length>=$min_overlap)){
                next;  
            }
            my $temp=$ss2."_".$ss3."_".$ss6."_".$ss7;
            if (exists $repeat_out{$ch_gene}{$temp}) {
                next;  
            }
            $repeat_out{$ch_gene}{$temp}++;
            my $lowCount="Low_Count";
            
            my $length;
            if (exists $junction_length_count{$ch}{$i_junction}) {
                my $length_count=$junction_length_count{$ch}{$i_junction};
                ($length,$i_count)=split /\t/,$length_count;
                $normalized_i_count=$i_count/$length;
            }
            if (exists $junction_length_count{$ch}{$s_junction}) {
                my $length_count=$junction_length_count{$ch}{$s_junction};
                ($length,$s_count)=split /\t/,$length_count;
                $normalized_s_count=$s_count/$length;
            }
            if (($i_count+$s_count)>0) {
                if (($i_count+$s_count)>=10) {
                    $lowCount="-";
                }
                my $psi=$i_count/($i_count+$s_count);
                print OUT "$gene\t$ch\t$strand\t$ss6\_$ss7\t$ss2\_$ss3\t-\t-\t$i_count\t$s_count\t-\t$lowCount\t-\t$psi\n";
            }else{
                print OUT "$gene\t$ch\t$strand\t$ss6\_$ss7\t$ss2\_$ss3\t-\t-\t$i_count\t$s_count\t-\t$lowCount\t-\t-\n";
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
