#!/usr/bin/env perl
use warnings;
use strict;

# Execute the compiled Fortran files and produce a plot

my $EXEC_NAME = "geometry_test.exe";

sub main() {
    print "Running $EXEC_NAME...\n";
    system "./bin/$EXEC_NAME";
    system "python ./src_py/geometry_test_plot.py";
    print "Done\n";
}

main() unless caller;

__END__
