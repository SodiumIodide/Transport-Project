#!/usr/bin/env perl
use warnings;
use strict;

# Compile Fortran code into executables

my $SRC = "./src_f";
my $BIN = "./bin";
my $COMPILER = "gfortran";
my @FILES = qw/
    self_library.f90
    mcnp_random.f90
    dynamic_arrays.f90
    geometry_gen.f90
    mesh_map.f90
    steady_state_slab_non_mark.f90
/;
my $EXEC_NAME = "steady_state_slab.exe";

sub main {
    print "Compiling using $COMPILER...\n";
    # Build overall object list for executable compilation

    my $obj_string = "";
    # Compile individual objects
    foreach my $src_name (@FILES) {
        my $obj_name = $src_name =~ s/\.f90/\.o/r;
        comp("$COMPILER -c $SRC/$src_name -o $BIN/$obj_name");
        $obj_string = $obj_string . " $BIN/" . $obj_name;
    }

    # Compile executable
    comp("$COMPILER -o $BIN/$EXEC_NAME" . $obj_string);
    print "Compilation completed\n";
    clean();
}

sub comp {
    my $command = shift;
    open my $output, "-|", $command;
    while (<$output>) {
        die "Compilation failed\n" if length($_) > 0;
    }
    close $output;
}

sub clean {
    opendir(my $dir, $SRC) or die "Couldn't opn $SRC: $!";
    while (my $filename = readdir($dir)) {
        unlink "$SRC/$filename" if $filename =~ /.*[\.o|\.exe]$/;
    }
    print "$SRC cleaned of binary\n";
}

main() unless caller;

__END__
