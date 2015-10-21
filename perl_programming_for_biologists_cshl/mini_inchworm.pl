#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 kmer_counts.txt min_contig_length=100\n\n";

my $kmer_counts_file = $ARGV[0] or die $usage;
my $MIN_CONTIG_LENGTH = $ARGV[1] || 100;

my $VERBOSE = 1; # print lots of messages to track run along the way.

my $KMER_SIZE = 25;


# parse the count(tab)kmer file into  kmer => count  hashtable.
my %KMER_COUNTS = &parse_KMER_COUNTS($kmer_counts_file);
print STDERR "-done parsing kmer counts\n";

# sort the kmers in descending order of abundance
my @sorted_kmers_desc = reverse sort {$KMER_COUNTS{$a}<=>$KMER_COUNTS{$b}} keys %KMER_COUNTS;
print STDERR "-done sorting kmers descendingly by count.\n";

my $contig_counter = 0;

while (@sorted_kmers_desc) {

    # take the seed as the most abundant kmer
    my $seed_kmer = shift @sorted_kmers_desc; # remove kmer at index [0], then shift indices to left.


    # use it to seed an extension (as long as we haven't 
    # already used it from a previous iteration)
    
    if ($KMER_COUNTS{$seed_kmer} > 0) {
        
        my $contig = &build_contig($seed_kmer);

        if (length($contig) >=  $MIN_CONTIG_LENGTH) {
            $contig_counter++;
            print "Contig [$contig_counter]: $contig\n";
        }
    
        # Now discard those kmers that exist in our contig.
        # We can only use each kmer once.
        
        &remove_contig_kmers($contig);
        
    }
        
}

exit(0);
    

####
sub build_contig {
    my ($seed_kmer) = @_;

    ##############################################################
    ## Greedily extend contig using the highest scoring extension.
    ##############################################################
    

    # start the contig as the seed kmer itself
    my $contig = $seed_kmer;

    my %seen_kmer;
    
    my $have_extension_flag = 1;
    while ($have_extension_flag) {

        $have_extension_flag = 0; # by default, let's assume no extension possible.
        
        # get the prefix for the next kmer (it's the k-1 last characters of the seed:
        my $next_kmer_prefix = substr($contig, -1 * ($KMER_SIZE-1));
                
        
        my $best_kmer_count = 0;
        my $best_char = "";

        foreach my $nuc_char ('G', 'A', 'T', 'C') {

            ### <----  Write code here to find best extension kmer

            ## Todo:
            #     -contruct the next possible kmer based on the above nucleotide character extension
            #     -check - is it the best extension so far?
            #            - set $best_kmer_count and $best_char accordingly.

            





            
            ###   ----->
            
        }

        
        ### <----  Write code to extend contig by $best_char (if we have a $best_char)
        
        ## Todo:
        #       -if we have a best char, use it to extend the contig.
        #       -if we do not have an extension, must stop (Hint: see while loop condition)
        #       -beware of loops...  ie. if seed kmer is 'AAAAAAAAAAAAAA', you might keep
        #           extending it by A for an eternity.
        #           Hint: if you've used an extension kmer in this contig, don't use it again.
        
        

        


        

        ###  -----> 


    }
        
    return($contig);

}


####
sub remove_contig_kmers {
    my ($contig) = @_;

    print STDERR "-removing kmers from contig: $contig\n" if $VERBOSE;
    
    for (my $i = 0; $i <= length($contig) - $KMER_SIZE; $i++) {

        my $kmer = substr($contig, $i, $KMER_SIZE);
        
        $KMER_COUNTS{$kmer} = 0;
    }

    return;
}
    

####
sub parse_KMER_COUNTS {
    my ($kmers_file) = @_;

    print STDERR "-parsing kmer counts from file: $kmers_file\n" if $VERBOSE;
    
    my $count;
    my %kmer_counts;
    
    open (my $fh, $kmers_file) or die "Error, cannot open file $kmers_file";
    while (<$fh>) {
        chomp;
        my ($count, $kmer) = split(/\t/);
        $kmer_counts{$kmer} = $count;
    }
    close $fh;


    return(%kmer_counts);
    
}
