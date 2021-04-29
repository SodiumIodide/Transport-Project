#!/usr/bin/env perl
use warnings;
use strict;

# Execute the compiled Fortran files and produce a plot

my $EXEC_NAME = "steady_state_slab_non_mark_closure.exe";

sub main() {
    print "Running $EXEC_NAME...\n";
    system "./bin/$EXEC_NAME";
    system "python ./src_py/non_mark_closure_plot.py";
    system "python ./src_py/non_mark_closure_plot_overlap.py";
    print "Done\n";
}

main() unless caller;

__END__
