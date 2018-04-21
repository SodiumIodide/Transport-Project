#!/usr/bin/env perl
use warnings;
use strict;

# Execute the compiled Fortran files and produce a plot

my $EXEC_NAME = "mc_slab_non_mark.exe";

sub main() {
    print "Running $EXEC_NAME...\n";
    system "./bin/$EXEC_NAME";
    system "python ./src_py/mc_plot_non_mark.py";
    print "Done\n";
}

main() unless caller;

__END__
