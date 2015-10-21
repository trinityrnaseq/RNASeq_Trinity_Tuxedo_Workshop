# Perl Programming Assignment - Complete the mini-inchworm program

Mini-inchworm is a simplified version of the Inchworm utility of Trinity.

The Mini-inchworm algorithm assembles contigs like so:

1.  reads in a list of kmers and their abundance values and create a dictionary of kmer counts.
2.  sorts the kmers by abundance value (descendingly), prioritizing them by count.
3.  taking a seed kmer as the most abundant kmer.
4.  greedily extend the seed kmer to the right, generating a contig.
5.  remove contig-containing kmers from our dictionary
6.  goto step 3, continue looping until the kmer catalog is depleted.


You'll find in the file 'kmer_counts.txt' and has format like so:

     1946    AAAAAAAAAAAAAAAAAAAAAAAAA
     2       GTGGCGGTGCTTGATGTAGATTTTA
     8       GCTGGCCCTTGGCCCCCACCCCCCC
     2       CTTCCTAAACTTGAAATCCGGTACT
     2       TATACCTTTTCTAAAAGAAGACAAC
     4       TCATGGGACGTTCATTTCGGGATGG
     2       ACCAATGAATAGGCAGGCGCCCGAT
     2       CTCTGTATAAGATTGACCCCAAAAA
     4       TGGTACCAGGCCTATAGTGTGTTCC
     234     GGGAAACCATTCTGAGATTTTAAAC
     2       CAGTGTTTCCTTCCTCTCACTTGCC
     42      GGGTTGCTTTTTCATCATGCAGTGG
     4       GCTCGGACTGCATTGTCTGGCTTTA



The Perl script 'mini_inchworm.pl' implements several of the steps above, but is missing the critical part of the code that involves assembling the contigs via the greedy extension mechanism.  That is left for you!


The command-line usage for mini_inchworm is shown below:

     %   ./mini_inchworm.pl


	usage: ./mini_inchworm.pl kmer_counts.txt min_contig_length=100


So, it takes two parameters: 1. the 'kmer_counts.txt' file described above, and then a setting for the minimum contig length (default set to 100).  

View the script and locate the subroutine

```Perl
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

```

It's your job to complete the code in the designated areas, using the variables that are defined in that subroutine.

## Things to consider

For testing purposes, try running the min_inchworm.pl with a min_contig_length set to 1.  You'll see a VERBOSE variable in the script that is set to 1 (true).  Simply add print statements within the script to track its execution, like so:

     print STDERR "printing something...."  if $VERBOSE;


Good luck!!!

