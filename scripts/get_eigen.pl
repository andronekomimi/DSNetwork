#!/usr/bin/perl 

# replace the file name in the following tabix command
my $tabix_cmd = 'tabix Eigen_hg19_0916_chr22.tab.bgz';


#list contains a list of files, each containing a list of variants
open(LST, "./list.txt");
my @files = <LST>;
chomp(@files);
close LST;

foreach $f (@files){
  print "processing $f ...\n";
  my $file = './'.$f;
  open(DATA, "$file") or die $!;
  my @cor_list = <DATA>;
  chomp(@cor_list);
  close DATA;

  my $out = $file.'.out';
  foreach $l (@cor_list){
    if($l =~ /chr/){
    }
    elsif($l !~ /^$/){
	my($chr,$start) = split(/\s+/,$l);
	my $corr = $chr.':'.$start.'-'.$start;
	print "$tabix_cmd $corr  >> $out\n";
  	system("$tabix_cmd $corr  >> $out");
    }
  }
}
