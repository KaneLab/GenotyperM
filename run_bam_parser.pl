
$sample_list = shift;
$bam_dir = shift;
$vcf = shift;
$target_dir = shift;

die "Usage: $0 sample_list bam_dir vcf_file target_dir\n" unless $target_dir;


$i=0;
open(FI, $sample_list);
while(<FI>)
{
	chomp;
	$i ++;
	&run ( "bam_parser $bam_dir/$_.sorted.bam $vcf $target_dir/$_ $i");
}
close FI;



sub run
{
    $_ = shift;
    print $_, "\n";
    system $_;
}



