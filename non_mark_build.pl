#!/usr/bin/env perl
use warnings;
use strict;

# Compile Fortran code into executables

my $SRC = "./src_f";
my $BIN = "./bin";
my $COMPILER = "gfortran";
my @FILES = qw/
    self_library.f90
    steady_state_slab_non_mark.f90
/;
my $EXEC_NAME = "steady_state_slab_non_mark.exe";
my $OPT_LEVEL = 3;

sub main {
    print "Compiling using $COMPILER...\n";
    # Build overall object list for executable compilation

    unlink "$SRC/$EXEC_NAME" if -e "$SRC/$EXEC_NAME";
    my $obj_string = "";
    # Compile individual objects
    foreach my $src_name (@FILES) {
        my $obj_name = $src_name =~ s/\.f90/\.o/r;
        comp("$COMPILER -Wall -c $SRC/$src_name -o $BIN/$obj_name -O$OPT_LEVEL");
        $obj_string = $obj_string . " $BIN/" . $obj_name;
    }

    # Compile executable
    comp("$COMPILER -Wall -O$OPT_LEVEL -o $BIN/$EXEC_NAME" . $obj_string);
    die "Compilation failed\n" unless -e "$BIN/$EXEC_NAME";
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
    opendir(my $dir, $SRC) or die "Couldn't open $SRC: $!";
    while (my $filename = readdir($dir)) {
        unlink "$SRC/$filename" if $filename =~ /.*[\.o|\.exe]$/;
    }
    closedir $dir;
    print "$SRC cleaned of binary\n";
    opendir(my $curdir, '.') or die "Couldn't open current directory: $!";
    while (my $filename = readdir($curdir)) {
        unlink "./$filename" if $filename =~ /.*\.mod$/;
    }
    closedir $curdir;
    print "Working directory cleaned of compiler artifacts\n";
}

main() unless caller;

__END__
