#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use Cwd;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;


#################################################################################
#
#  --autopilot         automatically run the pipeline end-to-end
#
#################################################################################



__EOUSAGE__

    ;


my $help_flag = 0;
my $AUTO_MODE = 0;


&GetOptions( 'help|h' => \$help_flag,
             'autopilot' => \$AUTO_MODE,
             
    );


if ($help_flag) {
    die $usage;
}

## first check for tools needed.

my @tools = qw (
    bowtie
    tophat
    cufflinks
    cuffcompare
    cuffmerge
    cuffdiff
    samtools
    R
    
);

{
    my $missing_tool_flag = 0;
    foreach my $tool (@tools) {
        my $path = `which $tool`;
        unless ($path =~ /\w/) {
            print STDERR "Error, cannot find path to: $tool\n";
            $missing_tool_flag = 1;
        }
    }
    
    if ($missing_tool_flag) {
        die "\n\nTools must be in PATH setting before proceeding.\n\n";
    }

}

## check for env vars
unless (defined $ENV{IGV}) {
    die "Error, must set env var for IGV path";
}

{ 
    ## unzip the gzipped files.
    foreach my $file (<*.gz>) {
        my $unzipped_file = $file;
        $unzipped_file =~ s/\.gz//;
        unless (-s $unzipped_file) {
            my $ret = system("gunzip -c $file > $unzipped_file");
            if ($ret) {
                die "Error, guznip of $file died ";
            }
        }
    }
}


main: {

    unless (-e "genome.1.ebwt") {
        &process_cmd("bowtie-build genome.fa genome");
    }
    
    my @samples = ("Sp_ds", "Sp_hs", "Sp_log", "Sp_plat");
    
    my @tophat_bam_files;
    my @cufflinks_gtf_files;

    foreach my $sample (@samples) {

        ## Tophat alignments of reads to the genome:
        
        my $tophat_out_dir = "tophat.$sample.dir";
        my $tophat_bam = "$tophat_out_dir/$sample.bam";
        
        unless (-s $tophat_bam) {
            
            &process_cmd("tophat -I 1000 -i 20 --bowtie1 --library-type fr-firststrand -o $tophat_out_dir genome $sample.left.fq $sample.right.fq");
            &process_cmd("mv $tophat_out_dir/accepted_hits.bam $tophat_bam");
            
        }
        
        &process_cmd("samtools index $tophat_bam") unless (-s "$tophat_bam.bai");
            
        push (@tophat_bam_files, $tophat_bam);

        ## Run Cufflinks to generate transcript structures.
        
        my $cuff_out_dir = "cufflinks.$sample.dir";
        my $cuff_gtf_file = "$cuff_out_dir/$sample.transcripts.gtf";
        
        unless (-s $cuff_gtf_file) {
            &process_cmd("cufflinks --overlap-radius 1 --library-type fr-firststrand -o $cuff_out_dir $tophat_bam");
        
            &process_cmd("mv $cuff_out_dir/transcripts.gtf $cuff_gtf_file");
        }
        
        push (@cufflinks_gtf_files, $cuff_gtf_file);
        
        if ($sample eq "Sp.ds") {
        
            eval { 

                my $cmd = "java -Xmx2G -jar $ENV{IGV}/igv.jar -g `pwd`/genome.fa `pwd`/$cuff_gtf_file,`pwd`/genes.bed,`pwd`/$tophat_bam";
                if ($AUTO_MODE) {
                    $cmd .= " & ";
                }
                &process_cmd($cmd);
            };
        }
        
    }
    

    
    {  ## Cuffmerge the individual cufflinks gtf files:
        
        unless (-s "assemblies.txt") {
            foreach my $cuff_file (@cufflinks_gtf_files) {
                &process_cmd("echo $cuff_file >> assemblies.txt");
            }
        }
        
        &process_cmd("cat assemblies.txt");
        
        &process_cmd("cuffmerge -s genome.fa assemblies.txt") unless (-s "merged_asm/merged.gtf");
        
        
        eval { 
            
            ## show merged cuff transcripts in IGV along with the reads for all conditions.
            
            my $cmd = "java -Xmx2G -jar $ENV{IGV}/igv.jar -g `pwd`/genome.fa `pwd`/merged_asm/merged.gtf,`pwd`/genes.bed";
            foreach my $bam (@tophat_bam_files) {
                $cmd .= ",`pwd`/$bam";
            }

            if ($AUTO_MODE) {
                $cmd .= " & ";
            }
            
            &process_cmd($cmd);
            

        };
        

        ################
        ## Run Cuffdiff
        ################
        
        &process_cmd("cuffdiff --library-type fr-firststrand  -o diff_out -b genome.fa -L " . join(",", @samples) . " -u merged_asm/merged.gtf " . join(" ", @tophat_bam_files)) unless (-e "diff_out/gene_exp.diff");
                     
        &process_cmd("head diff_out/gene_exp.diff");
        
    }
    

    exit(0);


}



####
sub process_cmd {
    my ($cmd) = @_;
    
    unless ($AUTO_MODE) {
        
        my $response = "";
        while ($response !~ /^[YN]/i) {
            print STDERR "\n\n"
                . "###############################################\n"
                . "CMD: $cmd\n"
                . "###############################################\n\n"
                . "Execute (Y/N)? ";

            $response = <STDIN>;
        }

        if ($response =~ /^N/i) {
            print STDERR "\t *** Exiting on demand. ****\n\n"
                . "Goodbye. \n\n\n";
            exit(0);
        }
    }
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    return;
}

