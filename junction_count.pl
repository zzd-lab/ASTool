#!/usr/bin/perl -w
use strict;
use Time::HiRes qw( time );
use Getopt::Long;

my $usage="$0 [options]
--gtf		<str>		the gene annotation file
--sam		<str>		SAM file
--thread	<int>		the number of thread(default:1)
--readlength	<int>		the length of read
--m		<int>		minimum length between read and exon or intron(default:8)
--outdir	<str>		the path of outdir
	";
my($gtf,$sam,$thread,$readlength,$m,$outdir);
GetOptions("gtf:s" =>\$gtf,
	"sam:s" =>\$sam,
	"thread:i" =>\$thread,
	"readlength:i" =>\$readlength,
	"m:i" =>\$m,
	"outdir:s" =>\$outdir


);

die $usage if !defined $gtf;

$thread ||=1;
$m ||=8;




my $start = time();
###input parameter
#my $input_gtf=$ARGV[0];
#my $read_lengh=$ARGV[1];#read_length
my $Annotation_file=$gtf;
my $sam_file=$sam;
my $th_p=$thread;
my $read_length=$readlength;
my $min_overlap=$m;
my $out_dir=$outdir;
#my $min_overlap=$ARGV[5];
###creat output file
$out_dir=~s/\/$//gi;
unless (-e $out_dir) {
    `mkdir $out_dir`;
}
###junction_position

my %junction_pos;
my %junction_length;
my $file_type= (split /\./,$Annotation_file)[-1];
if (($file_type eq "gtf")||($file_type eq "gff3") ) {
    my %versions;
    my %transcript_gene;
    my %exon_pos;
    my %gene_transcripts;
    my %ch_transcript;
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
                    $exon_pos{$transcript}{$start}=$end;
                    $gene_transcripts{$gene}{$transcript}++;
                    $transcript_gene{$transcript}=$gene;
                    $ch_transcript{$transcript}=$ch;
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
                    $exon_pos{$transcript}{$start}=$end;
                    $gene_transcripts{$gene}{$transcript}++;
                    $ch_transcript{$transcript}=$ch;  
                }   
            }
        }
        close IN;
    }
    my %repeat_junction;
    foreach my $transcript(sort keys %exon_pos){
        my $ch=$ch_transcript{$transcript};
        my @exon_starts=sort {$a<=>$b} keys  %{$exon_pos{$transcript}};
        if ($#exon_starts>0) {
            for (my $i=0;$i<=$#exon_starts-1;$i+=1){
                my $exon1_start=$exon_starts[$i];
                my $exon1_end=$exon_pos{$transcript}{$exon1_start};
                my $exon2_start=$exon_starts[$i+1];
                my $exon2_end=$exon_pos{$transcript}{$exon2_start};
                my $intron_start= $exon1_end+1;
                my $intron_end= $exon2_start-1;
                my $exon1_length=$exon1_end-$exon1_start+1;
                my $exon2_length=$exon2_end-$exon2_start+1;
                my $intron_length=$intron_end-$intron_start+1;
                my $e1_i_junction=$exon1_start."_".$exon1_end."_".$intron_start."_".$intron_end;
                my $i_e2_junction=$intron_start."_".$intron_end."_".$exon2_start."_".$exon2_end;
                my $e1_e2_junction=$exon1_start."_".$exon1_end."_".$exon2_start."_".$exon2_end;
                my $e1_i_length=$read_length-2*$min_overlap+1;
                my $i_e2_length=$read_length-2*$min_overlap+1;
                my $e1_e2_length=$read_length-2*$min_overlap+1;
                my $i_junction=$intron_start."_".$intron_end;
                my $i_length=$intron_end-$intron_start+1;
                if ($i==0) {
                    $e1_i_length=min($exon1_length,$read_length-$min_overlap)-$min_overlap+1;
                    if ($i==($#exon_starts-1)) {
                        $i_e2_length=min($exon2_length,$read_length-$min_overlap)-$min_overlap+1;
                        $e1_e2_length=min($exon1_length,$read_length-$min_overlap)+min($exon2_length,$read_length-$min_overlap)-$read_length+1;
                    }else{
                        $e1_e2_length=min($exon1_length,$read_length-$min_overlap)-$min_overlap+1;
                    }
                    
                }else{
                    if ($i==($#exon_starts-1)) {
                        $i_e2_length=min($exon2_length,$read_length-$min_overlap)-$min_overlap+1;
                        $e1_e2_length=min($exon2_length,$read_length-$min_overlap)-$min_overlap+1;
                    }
                }
                $junction_length{$ch}{$e1_i_junction}=$e1_i_length;
                $junction_length{$ch}{$i_e2_junction}=$i_e2_length;
                $junction_length{$ch}{$e1_e2_junction}=$e1_e2_length;
                $junction_length{$ch}{$i_junction}=$i_length;
                unless (exists $repeat_junction{$ch}{$e1_i_junction}){
                    #print OUT "$ch\t$e1_i_junction\t$e1_i_length\tFlank\n";
                    my $type="Flank";
                    my ($junc1,$junc2,$junc3,$junc4)=split /_/,$e1_i_junction;
                    my $important1=int($junc1/1000); 
                    my $important2=int($junc2/1000);
                    my $important3=int($junc3/1000); 
                    my $important4=int($junc4/1000); 
                    for (my $important=$important1;$important<=$important2;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$e1_i_junction}++;
                    }
                    for (my $important=$important3;$important<=$important4;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$e1_i_junction}++;
                    }
                    $repeat_junction{$ch}{$e1_i_junction}++;
                }
                unless (exists $repeat_junction{$ch}{$i_e2_junction}){
                    my $type="Flank";
                    my ($junc1,$junc2,$junc3,$junc4)=split /_/,$i_e2_junction;
                    my $important1=int($junc1/1000); 
                    my $important2=int($junc2/1000);
                    my $important3=int($junc3/1000); 
                    my $important4=int($junc4/1000); 
                    for (my $important=$important1;$important<=$important2;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$i_e2_junction}++;
                    }
                    for (my $important=$important3;$important<=$important4;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$i_e2_junction}++;
                    }
                    $repeat_junction{$ch}{$i_e2_junction}++;
                }
                unless (exists $repeat_junction{$ch}{$e1_e2_junction}){
                    my $type="Span";
                    my ($junc1,$junc2,$junc3,$junc4)=split /_/,$e1_e2_junction;
                    my $important1=int($junc1/1000); 
                    my $important2=int($junc2/1000);
                    my $important3=int($junc3/1000); 
                    my $important4=int($junc4/1000); 
                    for (my $important=$important1;$important<=$important2;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$e1_e2_junction}++;
                    }
                    for (my $important=$important3;$important<=$important4;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$e1_e2_junction}++;
                    }
                    $repeat_junction{$ch}{$e1_e2_junction}++;
                }  
            }
        }
        if ($#exon_starts>1) {
            for (my $i=0;$i<=$#exon_starts-2;$i+=1){
                my $exon1_start=$exon_starts[$i];
                my $exon1_end=$exon_pos{$transcript}{$exon1_start};
                my $exon3_start=$exon_starts[$i+2];
                my $exon3_end=$exon_pos{$transcript}{$exon3_start};
                my $e1_e3_junction=$exon1_start."_".$exon1_end."_".$exon3_start."_".$exon3_end;
                my $e1_e3_length=$read_length-2*$min_overlap+1;
                my $exon1_length=$exon1_end-$exon1_start+1;
                my $exon3_length=$exon3_end-$exon3_start+1;
                if ($i==0) {
                    if ($i==($#exon_starts-2)) {
                        $e1_e3_length=min($exon1_length,$read_length-$min_overlap)+min($exon3_length,$read_length-$min_overlap)-$read_length+1;
                    }else{
                        $e1_e3_length=min($exon1_length,$read_length-$min_overlap)-$min_overlap+1;
                    }
                }else{
                    if ($i==($#exon_starts-2)) {
                        $e1_e3_length=min($exon3_length,$read_length-$min_overlap)-$min_overlap+1;
                    }
                }
                $junction_length{$ch}{$e1_e3_junction}=$e1_e3_length;
                unless (exists $repeat_junction{$ch}{$e1_e3_junction}){
                    my $type="Span";
                    my ($junc1,$junc2,$junc3,$junc4)=split /_/,$e1_e3_junction;
                    my $important1=int($junc1/1000); 
                    my $important2=int($junc2/1000);
                    my $important3=int($junc3/1000); 
                    my $important4=int($junc4/1000); 
                    for (my $important=$important1;$important<=$important2;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$e1_e3_junction}++;
                    }
                    for (my $important=$important3;$important<=$important4;$important++){
                        my $tag=$ch.",".$important;
                        $junction_pos{$tag}{$type}{$e1_e3_junction}++;
                    }
                    $repeat_junction{$ch}{$e1_e3_junction}++;
                }  
            } 
        }
    }
    close OUT;
   
    
}else{
    print "Notice: The format of the input file should be gtf or gff3!\n";
}

###split same file
 my $sam_file_name;
if ($th_p>1) {
    use threads; 
    use Config;
    my $num=`ls -l $sam_file| awk '{ print \$5}'`;
    chomp $num;
    my $number=(split /\s+/,$num)[0];
    my $size=int($number/$thread/1024/1024)+1;
    $size=$size."m";
    $sam_file_name=(split /\//,$sam_file)[-1];
    $sam_file_name=$out_dir."/".$sam_file_name;
    `split -C $size $sam_file -d -a 3 $sam_file_name`;
}


sub read_count
{
    #open TEMP, ">/home/qihuan/PlantIR/try01/temp.txt";
    my  $all_count_temp=0;
    my %junction_count;
    my %junction_count_other;
    my %exon_count;
    my %intron_count;
    my @canshu=@_;
    open(IN,"$canshu[0]") or die "can't open $canshu[0]";
    while (<IN>){
        chomp(my $line=$_);
        
        if($line=~m/^\@/){
            next
        }
        unless(($line=~m/NH:i:1\t/)||($line=~m/NH:i:1$/)){
            next
        }
        #print TEMP "$line\n";
        my ($read_name,$ch,$read_start,$code)=(split /\t/, $line)[0,2,3,5];
        if($code=~m/D|H|S|I/){
            next
        }
        $all_count_temp++;
       
        if ($code=~m/N/){
            my @n_cishu=($code=~m/([0-9]+N[0-9]+)/g);
            my $n_cishu=$#n_cishu+1;
            if ($n_cishu==1) {
                my @n_cishu=($code=~m/([0-9]+[A-Z])/gi);
                my $read_e1_length=$n_cishu[0];
                my $read_intron_length=$n_cishu[1];
                my $read_e2_length=$n_cishu[2];
                $read_e1_length=~s/[A-Z]//gi;
                $read_intron_length=~s/[A-Z]//gi;
                $read_e2_length=~s/[A-Z]//gi;
                my $read_exon1_start=$read_start;
                #my $read_intron_start=$read_start+$read_e1_length;
                #my $read_exon1_end=$read_intron_start-1;
                my $read_exon1_end=$read_start+$read_e1_length-1;
                #my $read_intron_end=$read_start+$read_e1_length+$read_intron_length-1;
                my $read_exon2_start=$read_start+$read_e1_length+$read_intron_length;
                my $read_exon2_end=$read_start+$read_e1_length+$read_intron_length-1+$read_e2_length;
                my $min_value=int($read_exon1_start/1000);
                my $max_value=int($read_exon2_end/1000)+1;
                my %junction_temp;
                for (my $important=$min_value;$important<=$max_value;$important++){
                    my $tag=$ch.",".$important;
                    if (exists $junction_pos{$tag}{"Span"}){
                        foreach my $junction(keys %{$junction_pos{$tag}{"Span"}}){
                            if (exists $junction_temp{$ch}{$junction}) {
                                next;
                            }
                            $junction_temp{$ch}{$junction}++;
                            my ($exon1_end,$exon2_start)=(split /_/, $junction)[1,2];
                            if (($read_exon1_start<=($exon1_end-$min_overlap+1))&&($exon1_end==$read_exon1_end)&&($exon2_start==$read_exon2_start)&&($read_exon2_end>=($read_exon2_start+$min_overlap-1))) {
                                $junction_count{$ch}{$junction}++;
                            }elsif((($read_exon1_start<=($exon1_end-$min_overlap+1))&&($exon1_end==$read_exon1_end))||(($exon2_start==$read_exon2_start)&&($read_exon2_end>=($read_exon2_start+$min_overlap-1)))){
                                unless(($exon1_end==$read_exon1_end)&&($exon2_start==$read_exon2_start)){
                                    $junction_count_other{$ch}{$junction}++;
                                }
                                
                            }
                        }
                    }
                    if (exists $junction_pos{$tag}{"Flank"}){
                        foreach my $junction(keys %{$junction_pos{$tag}{"Flank"}}){
                            if (exists $junction_temp{$ch}{$junction}) {
                                next;
                            }
                            $junction_temp{$ch}{$junction}++;
                            my ($exon1_end,$exon2_start)=(split /_/, $junction)[2,3];
                            my ($junc2,$junc3)=(split /_/, $junction)[1,2];
                            if ((($read_exon1_start<=($junc2-$min_overlap+1))&&($read_exon1_end>=($junc3+$min_overlap-1)))||($read_exon2_start<=($junc2-$min_overlap+1))&&($read_exon2_end>=($junc3+$min_overlap-1)))  {
                                $junction_count{$ch}{$junction}++;
                            }
                        }
                    }
                    
                    ################
                                  
                }
            }else{
                my @zuhe;
                my $zuhe_type;
                my $zuhe_type_exon;
                my $zuhe_index=-1;
                my $code_temp=$code;
                for(my $i=0;$i<=10;$i++){
                    if ($code_temp=~m/N/) {
                        $code_temp=~m/([0-9]+M[0-9]+N[0-9]+M)/;
                        my $pattern=$1;
                        $zuhe_index++;
                        $zuhe[$zuhe_index]=$pattern;
                        $pattern=~m/([0-9]+M[0-9]+N)/;
                        $code_temp=~s/^$1//;
                    }else{
                        last
                    }
                }
                my $read_start1=$read_start;
                my $zuhe_num=$#zuhe+1;
                my $zuhe_order=0;
                foreach my $code1(@zuhe){
                    $zuhe_order+=1;
                    if ($zuhe_order==1) {
                        $zuhe_type="First"
                    }elsif($zuhe_order==$zuhe_num){
                        $zuhe_type="Last"
                    }else{
                        $zuhe_type="Other"
                    }
                    my @n_cishu=($code1=~m/([0-9]+[A-Z])/gi);
                    my $read_e1_length=$n_cishu[0];
                    my $read_intron_length=$n_cishu[1];
                    my $read_e2_length=$n_cishu[2];
                    $read_e1_length=~s/[A-Z]//gi;
                    $read_intron_length=~s/[A-Z]//gi;
                    $read_e2_length=~s/[A-Z]//gi;
                    #my $read_intron_start=$read_start1+$read_e1_length;
                    my $read_exon1_start=$read_start1;
                    #my $read_exon1_end=$read_intron_start-1;
                    my $read_exon1_end=$read_start1+$read_e1_length-1;
                    #my $read_intron_end=$read_start1+$read_e1_length+$read_intron_length-1;
                    my $read_exon2_start=$read_start1+$read_e1_length+$read_intron_length;
                    my $read_exon2_end=$read_start1+$read_e1_length+$read_intron_length-1+$read_e2_length;
                    my $read_length=$read_e1_length+$read_e2_length;
                    $read_start1=$read_start1+$read_e1_length+$read_intron_length;
                    #my $read_exon2_start=$read_intron_end+1;
                    if ($zuhe_type=~m/First|Other/){
                        my $min_value=int($read_exon1_start/1000);
                        my $max_value=int($read_exon1_start/1000)+1;
                        my %junction_temp;
                        for (my $important=$min_value;$important<=$max_value;$important++){
                            my $tag=$ch.",".$important;
                            if (exists $junction_pos{$tag}{"Span"}){
                                foreach my $junction(keys %{$junction_pos{$tag}{"Span"}}){
                                    if (exists $junction_temp{$ch}{$junction}) {
                                        next;
                                    }
                                    $junction_temp{$ch}{$junction}++;
                                    my ($exon1_end,$exon2_start)=(split /_/, $junction)[1,2];
                                    if (($read_exon1_start<=($exon1_end-$min_overlap+1))&&($exon1_end==$read_exon1_end)&&($exon2_start==$read_exon2_start)&&($read_exon2_end>=($read_exon2_start+$min_overlap-1))) {
                                        $junction_count{$ch}{$junction}++; 
                                    }elsif((($read_exon1_start<=($exon1_end-$min_overlap+1))&&($exon1_end==$read_exon1_end))||(($exon2_start==$read_exon2_start)&&($read_exon2_end>=($read_exon2_start+$min_overlap-1)))){
                                        unless(($exon1_end==$read_exon1_end)&&($exon2_start==$read_exon2_start)){
                                            $junction_count_other{$ch}{$junction}++;
                                        }
                                    }
                                }
                            }
                            if (exists $junction_pos{$tag}{"Flank"}){
                                foreach my $junction(keys %{$junction_pos{$tag}{"Flank"}}){
                                    if (exists $junction_temp{$ch}{$junction}) {
                                        next;
                                    }
                                    $junction_temp{$ch}{$junction}++;
                                    my ($exon1_end,$exon2_start)=(split /_/, $junction)[2,3];
                                    my ($junc2,$junc3)=(split /_/, $junction)[1,2];
                                    if (($read_exon1_start<=($junc2-$min_overlap+1))&&($read_exon1_end>=($junc3+$min_overlap-1)))  {
                                        $junction_count{$ch}{$junction}++;
                                    }
                                }
                            }
                            
                        }
   
                    }else{
                        my $min_value=int($read_exon1_start/1000);
                        my $max_value=int($read_exon2_end/1000)+1;
                        my %junction_temp;
                        for (my $important=$min_value;$important<=$max_value;$important++){
                            my $tag=$ch.",".$important;
                            if (exists $junction_pos{$tag}{"Span"}){
                                foreach my $junction(keys %{$junction_pos{$tag}{"Span"}}){
                                    if (exists $junction_temp{$ch}{$junction}) {
                                        next;
                                    }
                                    $junction_temp{$ch}{$junction}++;
                                    my ($exon1_end,$exon2_start)=(split /_/, $junction)[1,2];
                                    if (($read_exon1_start<=($exon1_end-$min_overlap+1))&&($exon1_end==$read_exon1_end)&&($exon2_start==$read_exon2_start)&&($read_exon2_end>=($read_exon2_start+$min_overlap-1))) {
                                        $junction_count{$ch}{$junction}++; 
                                    }elsif((($read_exon1_start<=($exon1_end-$min_overlap+1))&&($exon1_end==$read_exon1_end))||(($exon2_start==$read_exon2_start)&&($read_exon2_end>=($read_exon2_start+$min_overlap-1)))){
                                        unless(($exon1_end==$read_exon1_end)&&($exon2_start==$read_exon2_start)){
                                            $junction_count_other{$ch}{$junction}++;
                                        }
                                    }
                                }
                            }
                            if (exists $junction_pos{$tag}{"Flank"}){
                                foreach my $junction(keys %{$junction_pos{$tag}{"Flank"}}){
                                if (exists $junction_temp{$ch}{$junction}) {
                                    next;
                                }
                                $junction_temp{$ch}{$junction}++;
                                    my ($exon1_end,$exon2_start)=(split /_/, $junction)[2,3];
                                    my ($junc2,$junc3)=(split /_/, $junction)[1,2];
                                    if ((($read_exon1_start<=($junc2-$min_overlap+1))&&($read_exon1_end>=($junc3+$min_overlap-1)))||($read_exon2_start<=($junc2-$min_overlap+1))&&($read_exon2_end>=($junc3+$min_overlap-1)))  {
                                        $junction_count{$ch}{$junction}++;
                                    }
                                }
                            }
                            
                            
                        }
                    }
                }
            }
        }else{
            unless($code=~m/M/){
                next;
            }
            $code=~m/([0-9]+)M/;
            my $read_length=$1;
            my $read_end=$read_start+$read_length-1;
            my $min_value=int($read_start/1000);
            my $max_value=int($read_start/1000)+1;
            my %junction_temp;
            for (my $important=$min_value;$important<=$max_value;$important++){
                my $tag=$ch.",".$important;
                if (exists $junction_pos{$tag}{"Flank"}){
                    foreach my $junction(keys %{$junction_pos{$tag}{"Flank"}}){
                        if (exists $junction_temp{$ch}{$junction}) {
                            next;
                        }
                        $junction_temp{$ch}{$junction}++;
                        my ($exon1_end,$exon2_start)=(split /_/, $junction)[2,3];
                        my ($junc2,$junc3)=(split /_/, $junction)[1,2];
                        if (($read_start<=($junc2-$min_overlap+1))&&($read_end>=($junc3+$min_overlap-1)))  {
                            $junction_count{$ch}{$junction}++;
                        }
                     }
                }  
            }
        }
    }
    open OUT,">$canshu[0]_junction_count_tmp.txt";
   
    foreach my $ch(keys %junction_count){
        foreach my $junction(keys %{$junction_count{$ch}}){
            my $count=$junction_count{$ch}{$junction};
            my $length=$junction_length{$ch}{$junction};
            if (exists $junction_count_other{$ch}{$junction}) {
                my $count_other=$junction_count_other{$ch}{$junction};
                print OUT "$ch\t$junction\t$length\t$count\t$count_other\n";
            }else{
                print OUT "$ch\t$junction\t$length\t$count\t0\n";
            }
            
            
        }
    }
    foreach my $ch(keys %junction_count_other){
        foreach my $junction(keys %{$junction_count_other{$ch}}){
            my $count_other=$junction_count_other{$ch}{$junction};
            my $length=$junction_length{$ch}{$junction};
            unless (exists $junction_count{$ch}{$junction}){
                my $count_other=$junction_count_other{$ch}{$junction};
                print OUT "$ch\t$junction\t$length\t0\t$count_other\n"; 
            }
            
            
        }
    }
    close OUT;
    #close TEMP;
}
# Simultaneous computation
my %junction_count_all;
my %junction_count_all_other;
if ($th_p>1) {
    my @thread_all=qw /1...$th_p/;
    for(my $i=0;$i<$th_p;$i++ ){
        my $num;
        if ($i<10) {
            $num="00".$i;
        }else{
            $num="0".$i
        }
        my $sam="$sam_file_name$num";
        $thread_all[$i] = threads->create( \&read_count,$sam);
    }    
    for(my $i=0;$i<$th_p;$i++ ){
        $thread_all[$i] ->join();
    }
    opendir (DIR, $out_dir) or die "can't open the directory!";
    my @dir = readdir DIR;
    foreach my $file (@dir){
        my $query=(split /\//, $sam_file_name)[-1];
        if ($file=~m/$query[0-9]+_junction_count_tmp\.txt/gi) {
            open(IN,"$out_dir/$file") or die "";
            while (<IN>) {
                chomp (my $line=$_);
                my ($ch,$junction,$count,$count_other)=(split /\t/,$line)[0,1,3,4];
                $junction_count_all{$ch}{$junction}+=$count;
                $junction_count_all_other{$ch}{$junction}+=$count_other;
            }
            close IN;
        }
    }
    my $file=$out_dir."/"."junction_count\.txt";
    open OUT, ">$file";
    foreach my $ch(keys %junction_count_all){
        foreach my $junction(keys %{$junction_count_all{$ch}}){
            my $count=$junction_count_all{$ch}{$junction};
            my $count_other=$junction_count_all_other{$ch}{$junction};
            my $length=$junction_length{$ch}{$junction};
            print OUT "$ch\t$junction\t$length\t$count\t$count_other\n";
        }
        
    }
    close OUT;
    for (my $i=0;$i<$th_p;$i++){
        my $j;
        if ($i<10) {
            $j="00".$i;
        }elsif($i<100){
            $j="0".$i;
        }
        my $file="$sam_file_name$j";
        my $file2="$sam_file_name$j"."_junction_count_tmp.txt";
        `rm $file`;
        `rm $file2`;
    }
}else{
    
    read_count($sam_file);
    my $file_temp=$sam_file."_junction_count_tmp.txt";
    my $file=$out_dir."/"."junction_count.txt";
    `mv $file_temp $file`;
    
}
my $end = time();
printf("Execution Time: %0.02f s\n", $end - $start);

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
    return $max
}
sub median{
    #仅限于3个值
    my @numbers=@_;
    my @number_new=sort {$a<=>$b} @numbers;
    my $median=$number_new[1];
    return $median;
}
sub median_true{
    my @numbers=@_;
    my @number_new=sort {$a<=>$b} @numbers;
    my $num=$#number_new+1;
    my $median;
    if (($num%2)>0) {
        $median=$number_new[int($num/2)]
    }else{
        $median=($number_new[int($num/2)]+$number_new[int($num/2)-1])/2
    }
    return $median;
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
