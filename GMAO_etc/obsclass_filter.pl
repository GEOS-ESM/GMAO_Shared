#!/usr/bin/env perl
use strict;
use warnings;
#=======================================================================
# name - obsclass_filter.pl
# purpose - filter an obsclass list and return the subset of data sets in the
#           obsys rcfile that are available for a specified date/time range.
#
# Notes:
# 1. See usage() information for calling parameters and options.
# 2. The data must be available at the beginning date/time and through the
#    number of hours specified in order for a "true" value to be returned.
# 3. The script does not check whether the observation data files actually
#    exist but only whether the obsys.rc file says that they are available.
#
# output (printed to STDOUT)
# => $obsStr (see usage info)
#
# revision history
# 20171223  jstassi  modified version of obsys_check.pl script 
#=======================================================================
use File::Basename;
use FindBin qw($Bin);
use lib ("$Bin");
use Getopt::Long;
use Manipulate_time ("num_days_in_month", "tick");

# global variables
#-----------------
my ($rcfile, $sortFLG, %obsclass, $qm);
my ($ymd1, $ymd2, $hr1, $hr2, $numhrs);

# main program
#-------------
{
    my ($found, $label, @obslist, $obsStr);
    init();

    # check availability of each member of obsclass
    #----------------------------------------------
    @obslist = ();

    open RC, "< $rcfile" or die "Error opening rcfile: $rcfile\n";
    while (<RC>) {
        $label = "";
        if (/^\s*BEGIN\s*(\S+)/) {
            $label = $1;
            next unless $obsclass{$label} or $obsclass{"all"};
        } else { next }

        $found = 0;
        if ($label) {
            while (<RC>) {
                chomp;
                next if /^\s*\#/;
                last if /^\s*END\b/ or /^\s*BEGIN\b/;
                last if ($found = search_date($_));
            }
            push @obslist, $label if $found;
            $obsclass{$label} = 0;  # include each label only once
        }
    }

    # print results
    #--------------
    @obslist = sort(@obslist) if $sortFLG;
    if (@obslist) { $obsStr = $qm .join(",", @obslist) .$qm }
    else          { $obsStr = 0 }

    print "$obsStr\n";
}

#=======================================================================
# name - init
# purpose - get runtime parameters and flags
#=======================================================================
sub init {
    my ($obsclass_str, $help, $year, $month, $day, $lastday, $hhmmss2);

    # get runtime options
    #--------------------
    GetOptions( "rc=s"   => \$rcfile,
                "sort"   => \$sortFLG,
                "h|help" => \$help );
    usage() if $help;

    $rcfile = "obsys.rc" unless $rcfile;
    die "Error. Cannot find file, $rcfile;" unless -e $rcfile;
    
    # get runtime parameters
    #-----------------------
    usage() unless scalar(@ARGV) == 4;
    ($obsclass_str, $ymd1, $hr1, $numhrs) = @ARGV;

    # remove leading and trailing quotes
    #-----------------------------------
    $qm = ""; 
    if ($obsclass_str =~ m/^('|")/ and $obsclass_str =~ m/($1)$/) {
        $qm = $1;
        $obsclass_str = substr($obsclass_str, 1, -1);
    }

    # extract classes from class string
    #----------------------------------
    %obsclass = ();
    foreach (split /[,]/, $obsclass_str) { $obsclass{$_} = 1 if $_ }

    # check inputs
    #-------------
    die "Error. Undecipherable date, $ymd1;" unless $ymd1 =~ m/\b\d{8}\b/;
    die "Error. Undecipherable hour, $hr1;"  unless $hr1  =~ m/\b\d{2}\b/
        or                                          $hr1  =~ m/\b\d{6}\b/;
    die "Error. Invalid numhrs, $numhrs;"    unless $numhrs =~ m/\b\d+\b/;

    # check date/time values
    #-----------------------
    ($year, $month, $day) = ( $ymd1 =~ m/(\d{4})(\d{2})(\d{2})/ );
    $lastday = num_days_in_month($year, $month);
    $hr1 = substr($hr1, 0, 2) if $hr1 =~ m/\b\d{6}\b/;

    die "Error. Invalid month, $month;"       if $month < 1 or $month > 12;
    die "Error. Invalid day-of-month, $ymd1;" if $day   < 1 or $day   > $lastday;
    die "Error. Invalid hour, $hr1;"          if $hr1   < 0 or $hr1   > 24;

    # tick to end date/time
    #----------------------
    ($ymd2, $hhmmss2) = tick $ymd1, "${hr1}0000", 0, "${numhrs}0000";
    $hr2 = substr($hhmmss2, 0, 2);
}

#=======================================================================
# name - search_date
# purpose - parce line from obsys.rc to see if datetime is included
#=======================================================================
sub search_date {
    my ($line, $found, $begdatetime, $enddatetime);
    $line = shift @_;
    $found = 0;
    
    if ($line =~ m/\s+(\d{8})_(\d{2})z-(\d{8})_(\d{2})z/) {
        $begdatetime = "$1$2";
        $enddatetime = "$3$4";
        $found = 1 if "$ymd1$hr1" >= $begdatetime and "$ymd1$hr1" <= $enddatetime
            and       "$ymd2$hr2" >= $begdatetime and "$ymd2$hr2" <= $enddatetime;
    }
    return $found;
}

#=======================================================================
# name - usage
# purpose - print usage information
#=======================================================================
sub usage {
    my $script = basename($0);
    print << "EOF";

usage: $script obsclass yyyymmdd hh numhrs [options]

where
  obsclass   label(s) identifying obsclass to check (multiple data
             collection labels separated by commas with no spaces);
             if equal "all", then return all classes in the rc file which
             match the date criteria
  yyyymmdd   beginning year-month-day to check, e.g. 19980101
  hh         beginning hour to check, e.g. 12
  numhrs     number of hours to check

options
  -rc        obsys rc file name (default = obsys.rc)
  -sort      sort classes alphabetically in output string
  -h         print usage information

description:
This script takes a string of obsclass names (comma-separated, no spaces) and
returns a string with the subset of classes available for the time period
specified by the input arguments according to the obsys rcfile.

output printed to STDOUT:
=> obsStr: string containing comma-separated list of observation data sets from
           input obsclass which are available for the specified date/time range

EOF
;
    exit;
}
