# script to translate mach output format to vcf
# now only works for one chromosome, but can easily be updated

use Getopt::Long;

GetOptions(
	   "sites=s" => \$sites_file,
	   # "refalt=s" => \$refalt_file,
	   "individual=s" => \$individual_file,
	   "mach=s"  => \$mach,
	  );

# die "Usage: $0 -m MACH_outpt -s sites_file -i indivdual -r refalt > out\n" unless $mach;
die "Usage: $0 -m MACH_outpt -s sites_file -i indivdual > out\n" unless $mach;

open(SI, $sites_file);
$m = 0;
while(<SI>)
  {
    @li = split /\s+/;
    push @markers, $li[0];
    $marker_index{$li[0]} = scalar @markers;
    $mar[$m] = $li[0];
    $ref[$m] = $li[1];
    $alt[$m] = $li[2];
    $m ++;
    # $ref{$li[0]} = $li[1];
    # $alt{$li[0]} = $li[2];
  }
close SI;

open(SA, $individual_file);
while(<SA>)
  {
    @li = split /\s+/;
    push @samples, $li[0];
    $sample{$li[0]} = scalar @samples;
  }
close SA;

# now the header
print "##fileformat=VCFv4.0\n";
print join "\t", qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
print "\t";
print join "\t", @samples;
print "\n";


$CHR = 20;
$ID = '.';
$QUAL = 1;
$INFO = '.';
$FILTER = ".";
$FORMAT = 'GT';
@sa = @samples;

open(MA, $mach);
while(<MA>)
  {
    ($sam, $hap, $seq) = split /\s+/;
    if ($sam =~ /(\d+)->(\d+)/)
    { $sample = $1; }
    foreach $m (0..(length($seq))) {
      $b = substr $seq, $m, 1;
      $h = (lc($b) eq lc($ref[$m]))?0:1;
      $hap{$CHR}{$markers[$m]}{$sample}{$hap} = $h;
    }
  }
close MA;

# now final output
foreach $marker (@markers)
  {
    printf "$CHR\t%d\t$ID\t%s\t%s\t", $marker, $ref[$marker_index{$marker}-1], $alt[$marker_index{$marker}-1];
    print "$QUAL\t$FILTER\t$INFO\t$FORMAT";
    foreach $sid (1..(scalar @samples))
      {
	print "\t", $hap{$CHR}{$marker}{$sid}{'HAPLO1'};
	print "|",  $hap{$CHR}{$marker}{$sid}{'HAPLO2'};
      }
    print "\n";
  }

