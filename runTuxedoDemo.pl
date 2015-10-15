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
    bowtie2
    tophat2
    cufflinks
    cuffcompare
    cuffmerge
    cuffdiff
    samtools
    R
    igv.sh
    
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


if (0) { 
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


my $BASEDIR = $FindBin::Bin;
my $RNASEQ_DATA_DIR = "$BASEDIR/RNASEQ_data";


main: {

    my $checkpoints_dir = cwd() . "/__Tuxedo_checkpoints";
    unless (-d $checkpoints_dir) {
        mkdir $checkpoints_dir or die "Error, cannot mkdir $checkpoints_dir";
    }
        
    unless (-e "genome.1.bt2") {
        &process_cmd("bowtie2-build GENOME_data/genome.fa genome", "$checkpoints_dir/bowtie2.genome.idx.build.ok");
    }
    
    

    my %RNASEQ_DATASETS = ( 'Sp_ds' => [ "$RNASEQ_DATA_DIR/Sp_ds.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_ds.right.fq.gz" ],
                            'Sp_hs' =>  [ "$RNASEQ_DATA_DIR/Sp_hs.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_hs.right.fq.gz" ],
                            'Sp_log' => [ "$RNASEQ_DATA_DIR/Sp_log.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_log.right.fq.gz" ],
                            'Sp_plat' => [ "$RNASEQ_DATA_DIR/Sp_plat.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_plat.right.fq.gz" ],
        );
    
    my @samples = sort keys %RNASEQ_DATASETS;
    
    my @tophat_bam_files;
    my @cufflinks_gtf_files;

    foreach my $sample (@samples) {

        ## Tophat alignments of reads to the genome:
        
        my $tophat_out_dir = "tophat.$sample.dir";
        my $tophat_bam = "$tophat_out_dir/$sample.bam";

        my ($left_fq, $right_fq) = @{$RNASEQ_DATASETS{$sample}};
        
        &process_cmd("tophat2 -I 1000 -i 20 --library-type fr-firststrand -o $tophat_out_dir genome $left_fq $right_fq", "$checkpoints_dir/tophat2.$sample.genome.align.ok");
        &process_cmd("mv $tophat_out_dir/accepted_hits.bam $tophat_bam", "$checkpoints_dir/rename.$sample.tophat2.bam.ok");
            
        
        
        &process_cmd("samtools index $tophat_bam", "$checkpoints_dir/idx.tophat.$sample.bam.ok");
            
        push (@tophat_bam_files, $tophat_bam);

        ## Run Cufflinks to generate transcript structures.
        
        my $cuff_out_dir = "cufflinks.$sample.dir";
        my $cuff_gtf_file = "$cuff_out_dir/$sample.transcripts.gtf";
        
        &process_cmd("cufflinks --no-update-check --overlap-radius 1 --library-type fr-firststrand -o $cuff_out_dir $tophat_bam", "$checkpoints_dir/cufflinks.$sample.ok");
        
        &process_cmd("mv $cuff_out_dir/transcripts.gtf $cuff_gtf_file", "$checkpoints_dir/cuff.rename.$sample.ok");
        
        
        push (@cufflinks_gtf_files, $cuff_gtf_file);
        
        if ($sample eq "Sp.ds") {
        
            eval { 

                my $cmd = "igv.sh -g `pwd`/genome.fa `pwd`/$cuff_gtf_file,`pwd`/genes.bed,`pwd`/$tophat_bam";
                if ($AUTO_MODE) {
                    $cmd .= " & ";
                }
                &process_cmd($cmd);
            };
        }
        
    }
    

    
    {  ## Cuffmerge the individual cufflinks gtf files:
        
        my $counter = 0;
        foreach my $cuff_file (@cufflinks_gtf_files) {
            if ($counter == 0) {
                &process_cmd("echo $cuff_file > assemblies.txt", "$checkpoints_dir/echo." . basename($cuff_file) . ".ok");
            }
            else {
                &process_cmd("echo $cuff_file >> assemblies.txt", "$checkpoints_dir/echo." . basename($cuff_file) . ".ok");
            }
            $counter++;
        }
                
        &process_cmd("cat assemblies.txt", "$checkpoints_dir/cat.assemblies.ok");
        
        &process_cmd("cuffmerge -s GENOME_data/genome.fa assemblies.txt", "$checkpoints_dir/cuffmerge.ok");
        
        
        ## show merged cuff transcripts in IGV along with the reads for all conditions.
        
        my $cmd = "igv.sh -g `pwd`/GENOME_data/genome.fa `pwd`/merged_asm/merged.gtf,`pwd`/GENOME_data/genes.bed";
        foreach my $bam (@tophat_bam_files) {
            $cmd .= ",`pwd`/$bam";
        }
        
        if ($AUTO_MODE) {
            $cmd .= " & ";
        }
        
        &process_cmd($cmd, "$checkpoints_dir/igv.view.cuff.ok");
                

        ################
        ## Run Cuffdiff
        ################

        my $cuffdiff_cmd = "cuffdiff  --no-update-check --library-type fr-firststrand  -o diff_out -b GENOME_data/genome.fa -L " . join(",", @samples) . " -u merged_asm/merged.gtf " . join(" ", @tophat_bam_files);
        
        &process_cmd($cuffdiff_cmd, "$checkpoints_dir/cuffdiff.ok");
                     
        &process_cmd("head diff_out/gene_exp.diff", "$checkpoints_dir/head.gene_exp.diff.ok");
        
    }


    print STDERR "\n\n\n\tDone with Tuxedo data processing.  Now explore data analysis using CummeRbund.\n\n\n";
    

    exit(0);


}



####
sub process_cmd {
    my ($cmd, $checkpoint) = @_;

    unless ($checkpoint) {
        die "Error, need checkpoint file defined";
    }
    
    if (-e $checkpoint) { return; }
        
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

    system("touch $checkpoint");
    
    return;
}

