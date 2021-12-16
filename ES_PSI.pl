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
    my $version=(sort {$versions{$b}<=>$versions{$a}} keys %versions)[0];
    my %es_pair;
    my %uniq_out;
    #判断是否是已知的ES
    foreach my $gene(keys %gene_transcripts){
        my @transcripts=sort keys %{$gene_transcripts{$gene}};
        my $tr_number=$#transcripts+1;
        if ($tr_number>1) {
            for (my $i=0;$i<=$#transcripts;$i++){
                for (my $j=0;$j<=$#transcripts;$j++){
                    unless ($i==$j){
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
                        if (($site_numbers1>=6)&&($site_numbers2>=4)) {
                            for(my $i=0;$i<=$site_numbers1-6;$i+=2){
                                for(my $j=0;$j<=$site_numbers2-4;$j+=2){
                                    if (($sites1[$i+1]==$sites2[$j+1])&&($sites1[$i+4]==$sites2[$j+2])) {
                                        my $strand=$exon_strand{$tr1};
                                        my $ch=$exon_ch{$tr1};
                                        $sites1[$i+1]+=1;
                                        $sites1[$i+2]-=1;
                                        $sites1[$i+3]+=1;
                                        $sites1[$i+4]-=1;
                                        #position of intron
                                        my $es=$ch."\t".$gene."\t".$strand."\t".$sites1[$i+1]."\t".$sites1[$i+2]."\t".$sites1[$i+3]."\t".$sites1[$i+4];
                                        $es_pair{$es}++;
                                        $sites1[$i+1]-=1;
                                        $sites1[$i+2]+=1;
                                        $sites1[$i+3]-=1;
                                        $sites1[$i+4]+=1;
                                    }
                                }
                            }
                        } 
                        
                    }
                }
            }
        }
    }
    my $file_tmp=(split /\//,$ARGV[0])[-1];
    my $out_file="Intron_pair"."_".$file_tmp;
    my %flank_intron;
    foreach my $gene(keys %gene_transcripts){
        foreach my $tr (keys %{$gene_transcripts{$gene}}){
            my @exon_sites;
            my $exon_site_index=-1;
            my $ch=$exon_ch{$tr};
            my $strand=$exon_strand{$tr};
            foreach my $exon(keys %{$exon_site{$tr}}){
               
                my ($site1,$site2)=(split /\t/,$exon)[2,3];
                $exon_site_index+=1;
                $exon_sites[$exon_site_index]=$site1;
                $exon_site_index+=1;
                $exon_sites[$exon_site_index]=$site2;
                 
            }
            my @sites=sort { $a <=> $b } @exon_sites;
            my $site_numbers=$#sites+1;
            if ($site_numbers>=4) {
                for(my $i=0;$i<=$site_numbers-2;$i+=2){
                    my $start=$sites[$i];
                    my $end=$sites[$i+1];
                    my $exon=$ch."\t".$strand."\t".$start."\t".$end;
                    my ($left1,$left2,$right1,$right2);
                    if ($i<($site_numbers-2)) {
                        $right1=$sites[$i+1]+1;
                        $right2=$sites[$i+2]-1;
                    }else{
                        $right1=0;
                        $right2=0
                    }
                    if ($i>0) {
                        $left1=$sites[$i-1]+1;
                        $left2=$sites[$i]-1;
                    }else{
                        $left1=0;
                        $left2=0;
                    }
                    my $flank=$left1."_".$left2."_".$right1."_".$right2;
                    my $pos1=int($left1/1000);
                    my $pos2=int($right2/1000);
                    if ($left1==0) {
                        $pos1=int($start/1000);
                    }
                    if ($right2==0) {
                        $pos2=int($end/1000);
                    }
                    for(my $pos=$pos1;$pos<=$pos2;$pos++){
                        $flank_intron{$ch}{$pos}{$flank}++;
                    }
                }
            }
            
            
            
        }    
    }
    
    open OUT,">$psi_file";
    my %overlap_s;
    my %overlap_i;
    my %overlap_e;
    foreach my $gene(keys %gene_transcripts){
        foreach my $tr (keys %{$gene_transcripts{$gene}}){
            my @exon_sites;
            my $exon_site_index=-1;
            my $ch=$exon_ch{$tr};
            my $strand=$exon_strand{$tr};
            foreach my $exon_q(keys %{$exon_site{$tr}}){
                my ($site1,$site2)=(split /\t/,$exon_q)[2,3];
                $exon_site_index+=1;
                $exon_sites[$exon_site_index]=$site1;
                $exon_site_index+=1;
                $exon_sites[$exon_site_index]=$site2;   
            }
            my @sites=sort { $a <=> $b } @exon_sites;
            my $site_numbers=$#sites+1;
            if ($site_numbers>=6) {
                for(my $i=0;$i<=$site_numbers-6;$i+=2){
                    my $start=$sites[$i+2];
                    my $end=$sites[$i+3];
                    
                    $sites[$i+1]+=1;
                    $sites[$i+2]-=1;
                    $sites[$i+3]+=1;
                    $sites[$i+4]-=1;
                    my $exon_q=$ch."\t".$strand."\t".$start."\t".$end;
                    my $pos1=int($sites[$i+1]/1000);
                    my $pos2=int($sites[$i+4]/1000);
                    my %flank_temp;
                    my $es=$ch."\t".$gene."\t".$strand."\t".$sites[$i]."\t".$sites[$i+1]."\t".$sites[$i+2]."\t".$sites[$i+3]."\t".$sites[$i+4]."\t".$sites[$i+5];
                    for(my $pos=$pos1;$pos<=$pos2;$pos++){
                        for my $flank(keys %{$flank_intron{$ch}{$pos}}){
                            if (exists $flank_temp{$flank}) {
                                next
                            }
                            my ($intron1_start,$intron1_end,$intron2_start,$intron2_end)=split /_/,$flank;
                            my $exon_start=$intron1_end+1;
                            my $exon_end=$intron2_start-1;
                            my $exon_db=$ch."\t".$strand."\t".$exon_start."\t".$exon_end;
                            unless (($exon_start == $start) && ($exon_end== $end)) {
                                if(($exon_start<$sites[$i+1])||($exon_end>$sites[$i+4])){
                                    if(($exon_start >= $start) && ($exon_start< $end)){
                                        $overlap_e{$exon_q}++;
                                        $overlap_e{$exon_db}++;
                                    }	
                                    if(($exon_end>$start)&&($exon_end <= $end)){
                                        $overlap_e{$exon_q}++;
                                        $overlap_e{$exon_db}++;
                                    }
                                    if(($exon_start == $end) and ($exon_end>= $end)){
                                        $overlap_e{$exon_q}++;
                                        $overlap_e{$exon_db}++;
                                    }
                                    if(($exon_end == $start) and ($exon_start<= $start)){
                                        $overlap_e{$exon_q}++;
                                        $overlap_e{$exon_db}++;
                                    }
                                }
                            }
                            unless (($sites[$i+1] == $intron1_start)&&($sites[$i+2]== $intron1_end)&&($sites[$i+3]== $intron2_start)&&($sites[$i+4]== $intron2_end)){  
                                if ((($sites[$i+1] == $intron1_start)&&($sites[$i+2]== $intron1_end))||(($sites[$i+3]== $intron2_start)&&($sites[$i+4]== $intron2_end))){
                                    $overlap_i{$es}++;
                                }
                                if (($sites[$i+1] == $intron1_start)&&($sites[$i+4]==$intron2_end)) {
                                    $overlap_s{$es}++;
                                }
                                
                            }
                            $flank_temp{$flank}++;
                        }
                    } 
                    $uniq_out{$es}++;
                    $sites[$i+1]-=1;
                    $sites[$i+2]+=1;
                    $sites[$i+3]-=1;
                    $sites[$i+4]+=1;
                }
            } 
        }
    }
    my %repeat_out;
    foreach my $es(sort keys %uniq_out){
        my ($ch,$gene,$strand,$exon1_start,$intron1_start,$intron1_end,$intron2_start,$intron2_end,$exon2_end)=split /\t/,$es;
        my $exon1_end=$intron1_start-1;
        my $exon2_start=$intron2_end+1;
        my $exon_start=$intron1_end+1;
        my $exon_end=$intron2_start-1;
        my $e1_e_junction=$exon1_start."_".$exon1_end."_".$exon_start."_".$exon_end;
        my $e_e2_junction=$exon_start."_".$exon_end."_".$exon2_start."_".$exon2_end;
        my $e1_e2_junction=$exon1_start."_".$exon1_end."_".$exon2_start."_".$exon2_end;
        my ($e1_e_count,$e1_e2_count,$e_e2_count)=(0,0,0);
        my ($normalized_e1_e_count,$normalized_e1_e2_count,$normalized_e_e2_count)=(0,0,0);
        my $exon_length=$exon_end-$exon_start+1;
        my $exon1_length=$exon1_end-$exon1_start+1;
        my $exon2_length=$exon2_end-$exon2_start+1;
        unless(($exon_length>=$min_overlap)&&($exon1_length>=$min_overlap)&&($exon2_length>=$min_overlap)){
            next;  
        }
        my $exon=$ch."\t".$strand."\t".$exon_start."\t".$exon_end;
        my $o_e="Yes";
        if (exists $overlap_e{$exon}){
            $o_e="No"
        }
        my $o_S="Yes";
        my $o_I="Yes";
        
        if (exists $overlap_i{$es}){
            $o_I="No";
        }
        if (exists $overlap_s{$es}){
            $o_S="No";
        }
        my $annotated_ES="Novel";
        $es=$ch."\t".$gene."\t".$strand."\t".$intron1_start."\t".$intron1_end."\t".$intron2_start."\t".$intron2_end;
        if (exists $repeat_out{$es}) {
            next;  
        }
        $repeat_out{$es}++;
        if (exists $es_pair{$es}) {
            $annotated_ES="Known";
        }
        my $intron1=$intron1_start."_".$intron1_end;
        my $intron2=$intron2_start."_".$intron2_end;
        my $clean="Not_Clean";
        if (($o_e eq "Yes")&&($o_S eq "Yes")&&($o_I eq "Yes")) {
           $clean="Clean";
        }
        #if ($gene eq "AT4G05020") {
        #    print "$exon_start\_$exon_end\t$o_e\t$o_I\t$o_S\n";
        #}
        
        my $length;
        if (exists $junction_length_count{$ch}{$e1_e_junction}) {
            my $length_count=$junction_length_count{$ch}{$e1_e_junction};
            ($length,$e1_e_count)=split /\t/,$length_count;
            $normalized_e1_e_count=$e1_e_count/$length;
        }
        if (exists $junction_length_count{$ch}{$e1_e2_junction}) {
            my $length_count=$junction_length_count{$ch}{$e1_e2_junction};;
            ($length,$e1_e2_count)=split /\t/,$length_count;
            $normalized_e1_e2_count=$e1_e2_count/$length;
        }
        if (exists $junction_length_count{$ch}{$e_e2_junction}) {
            my $length_count=$junction_length_count{$ch}{$e_e2_junction};
            ($length,$e_e2_count)=split /\t/,$length_count;
            $normalized_e_e2_count=$e_e2_count/$length;
        }
        my $lowCount="Low_Count";
        my $imbalance="-";
        if ((($e1_e_count>0)&&($e_e2_count==0))||(($e1_e_count==0)&&($e_e2_count>0))) {
            $imbalance="Imbalance";
        }
        if ((mean($e1_e_count,$e_e2_count)+$e1_e2_count)>0) {
            if ((mean($e1_e_count,$e_e2_count)+$e1_e2_count)>=10) {
                $lowCount="-";
            }
            my $psi=mean($normalized_e1_e_count,$normalized_e_e2_count)/(mean($normalized_e1_e_count,$normalized_e_e2_count)+$normalized_e1_e2_count);
            print OUT "$gene\t$ch\t$exon_start\_$exon_end\t$intron1_start\_$intron1_end\t$intron2_start\_$intron2_end\t$annotated_ES\t$clean\t$e1_e_count\t$e_e2_count\t$e1_e2_count\t$lowCount\t$imbalance\t$psi\n";
        }else{
            print OUT "$gene\t$ch\t$exon_start\_$exon_end\t$intron1_start\_$intron1_end\t$intron2_start\_$intron2_end\t$annotated_ES\t$clean\t$e1_e_count\t$e_e2_count\t$e1_e2_count\t$lowCount\t$imbalance\t-\n"
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
