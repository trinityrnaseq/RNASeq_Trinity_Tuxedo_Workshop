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

my $BASEDIR = $FindBin::Bin;
my $RNASEQ_DATA_DIR = "$BASEDIR/RNASEQ_data";

my %RNASEQ_DATASETS = ( 'Sp_ds' => [ "$RNASEQ_DATA_DIR/Sp_ds.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_ds.right.fq.gz" ],
                        'Sp_hs' =>  [ "$RNASEQ_DATA_DIR/Sp_hs.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_hs.right.fq.gz" ],
                        'Sp_log' => [ "$RNASEQ_DATA_DIR/Sp_log.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_log.right.fq.gz" ],
                        'Sp_plat' => [ "$RNASEQ_DATA_DIR/Sp_plat.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_plat.right.fq.gz" ],
    );

my $trinity_dir = $ENV{TRINITY_HOME} or die "Error, need env var TRINITY_HOME set to Trinity installation directory";
$ENV{PATH} .= ":$trinity_dir";  ## adding it to our PATH setting.



my $OS_type = `uname`;

## first check for tools needed.

my @tools = qw (Trinity
    bowtie
    bowtie2
    tophat2
    samtools
    gmap_build
    gmap
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
                die "Error, could not gunzip file $file";
            }
        }
    }
}


# Run Trinity.
my @left_fqs;
my @right_fqs;
foreach my $sample_type (keys %RNASEQ_DATASETS) {
    my ($left_fq, $right_fq) = @{$RNASEQ_DATASETS{$sample_type}};
    push (@left_fqs, $left_fq);
    push (@right_fqs, $right_fq);
}

my $checkpoints_dir = $FindBin::Bin . "/__TrinDemo_checkpoints_dir";
unless (-d $checkpoints_dir) {
    mkdir $checkpoints_dir or die "Error, cannot mkdir $checkpoints_dir";
}


my $run_Trinity_cmd = "Trinity --seqType fq --SS_lib_type RF "
    . " --left " . join(",", @left_fqs)
    . " --right " . join(",", @right_fqs)
    . " --CPU 2 --max_memory 1G";
&process_cmd($run_Trinity_cmd, "$checkpoints_dir/trinity.ok");

# Examine top of Trinity.fasta file
&process_cmd("head trinity_out_dir/Trinity.fasta", "$checkpoints_dir/head_trinity.ok");

# Get Trinity stats:
&process_cmd("$trinity_dir/util/TrinityStats.pl trinity_out_dir/Trinity.fasta", "$checkpoints_dir/trin_stats.ok");

## run gmap
&process_cmd("gmap_build -d genome -D . -k 13 GENOME_data/genome.fa", "$checkpoints_dir/gmap_build.ok");
&process_cmd("gmap -n 0 -D . -d genome trinity_out_dir/Trinity.fasta -f samse > trinity_gmap.sam", "$checkpoints_dir/gmap_align.ok");

## convert to bam file format
&process_cmd("samtools view -Sb trinity_gmap.sam > trinity_gmap.bam", "$checkpoints_dir/gmap_sam_to_bam.ok");
&process_cmd("samtools sort trinity_gmap.bam trinity_gmap", "$checkpoints_dir/gmap_sort_bam.ok");
&process_cmd("samtools index trinity_gmap.bam", "$checkpoints_dir/gmap_index_bam.ok");



## TODO:  examine the trinity gmap alignments.


## align the rna-seq reads against the genome, too, for comparison 
&process_cmd("bowtie2-build GENOME_data/genome.fa genome", "$checkpoints_dir/bowtie2_build_genome.ok");

my $tophat_align_cmd = "tophat2 -I 300 -i 20 genome " 
    . join(",", @left_fqs) . " " 
    . join(",", @right_fqs);

&process_cmd($tophat_align_cmd, "$checkpoints_dir/tophat_align_genome.ok");
&process_cmd("samtools index tophat_out/accepted_hits.bam", "$checkpoints_dir/index_tophat_bam.ok");


# use IGV

my $cmd = "igv.sh -g `pwd`/GENOME_data/genome.fa `pwd`/GENOME_data/genes.bed,`pwd`/tophat_out/accepted_hits.bam,`pwd`/trinity_gmap.bam";
if ($AUTO_MODE) {
    $cmd .= " & ";
}
&process_cmd($cmd, "$checkpoints_dir/igv.view_all.ok");



###################################
## Abundance estimation using RSEM
###################################

my @rsem_result_files;

foreach my $sample (sort keys %RNASEQ_DATASETS) {
    
    my ($left_fq, $right_fq) = @{$RNASEQ_DATASETS{$sample}};
    
    my $output_dir = "$sample.RSEM";
    
    my $rsem_result_file = "$output_dir/$sample.isoforms.results";
    push (@rsem_result_files, $rsem_result_file);
    
            
    my $align_estimate_command = "$trinity_dir/util/align_and_estimate_abundance.pl --seqType fq "
        . " --left $left_fq --right $right_fq "
        . " --transcripts trinity_out_dir/Trinity.fasta "
        . " --output_prefix $sample --est_method RSEM "
        . " --aln_method bowtie --trinity_mode --prep_reference"
        . " --output_dir $output_dir";
    
    &process_cmd($align_estimate_command, "$checkpoints_dir/$sample.align_estimate.ok");
        
    
    # look at the output
    &process_cmd("head $rsem_result_file", "$checkpoints_dir/head.$sample.rsem.ok");
    
    
    
}


## generate matrix of counts and perform TMM normalization
&process_cmd("$trinity_dir/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix Trinity_trans @rsem_result_files", "$checkpoints_dir/counts_matrix.ok");

## Look at the matrix
&process_cmd("head -n20 Trinity_trans.counts.matrix", "$checkpoints_dir/head.counts.matrix.ok");

## run edgeR
&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Trinity_trans.counts.matrix --method edgeR --dispersion 0.1 --output edgeR", "$checkpoints_dir/run.edgeR.ok");

# take a look at what edgeR generated:
&process_cmd("ls -ltr edgeR/", "$checkpoints_dir/ls.edgeR.dir.ok");


&process_cmd("head edgeR/Trinity_trans.counts.matrix.Sp_log_vs_Sp_plat.edgeR.DE_results", "$checkpoints_dir/head.edgeR.DE_results.ok");

&show("edgeR/Trinity_trans.counts.matrix.Sp_log_vs_Sp_plat.edgeR.DE_results.MA_n_Volcano.pdf");

&process_cmd("sed '1,1d' edgeR/Trinity_trans.counts.matrix.Sp_log_vs_Sp_plat.edgeR.DE_results | awk '{ if (\$5 <= 0.05) print;}' | wc -l", "$checkpoints_dir/count_signif_DE_trans.ok");



&process_cmd("cd edgeR", "$checkpoints_dir/cd.edgeR.ok");
# now do it in the script. :)
chdir("edgeR") or die "Error, could not cd to edgeR/"; 


&process_cmd("$trinity_dir/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../Trinity_trans.TMM.EXPR.matrix -P 1e-3 -C 2",
             "$checkpoints_dir/analyze_diff_expr.ok");


&process_cmd("wc -l diffExpr.P1e-3_C2.matrix", "$checkpoints_dir/wc_diff_expr_matrix.ok"); # number of DE transcripts + 1

&show("diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData",
    "$checkpoints_dir/cut_clusters_tree.ok");
&show("diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60/my_cluster_plots.pdf");

print STDERR "\n\n\tDemo complete.  Congratulations!  :)\n\n\n\n";


exit(0);

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
    
    print STDERR "\tNow running:\n\t\t$cmd\n\n\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    system("touch $checkpoint");
    
    return;
}


sub show {
    my ($image) = @_;

    my $cmd;

    if ($OS_type =~ /linux/i) {
        ## use evince
        $cmd = "evince $image";
    }
    else {
        ## rely on ImageMagick:
        $cmd = "open $image";
    }
    
    if ($AUTO_MODE) {
        $cmd .= " & ";
    }
    
    &process_cmd($cmd, "$checkpoints_dir/view." . basename($image) . ".ok");

    return;
}
