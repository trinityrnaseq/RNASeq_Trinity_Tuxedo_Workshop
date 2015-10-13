#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $FULL_CLEAN = $ARGV[0] || 0;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";




my @files_to_keep = qw (
Sp_ds.left.fq.gz
Sp_ds.right.fq.gz
Sp_hs.left.fq.gz
Sp_hs.right.fq.gz
Sp_log.left.fq.gz
Sp_log.right.fq.gz
Sp_plat.left.fq.gz
Sp_plat.right.fq.gz
cleanme.pl
genes.bed
genes.gff3
genome.fa
runTrinityDemo.pl
runTuxedoDemo.pl
docs
cummeRbund.demo.R
);                      


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		
        if (-f $file) {
            unlink($file);
        }
        elsif (-d $file) {
            if ($file =~ /trinity_out_dir/) {
                if ($FULL_CLEAN) {
                `rm -rf $file`;
                }
                else {
                    unlink "trinity_out_dir/Trinity.fasta";
                }
            }
            else {
                `rm -rf $file`;
            }
        }
    }
}



exit(0);
