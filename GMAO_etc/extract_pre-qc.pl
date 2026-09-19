#!/usr/bin/env perl
#--------------------------------------------------
#
# Purpose: create info files at given date/time
#          according to information in database. 
#
# Usage:
#
#  extract_pre-qc.pl [options] obsys.rc GSI_GridComp.rc
#
# !REVISION HISTORY:
#
#   28Sep2013 Todling  Initial code for append_gsigcrc.pl
#   05Mar2019 Sienkiewicz  modify code to print pre-qc file only
#
#--------------------------------------------------

use Env;                 # make env vars readily available
use FindBin;             # so we can find where this script resides
use File::Basename;      # for basename(), dirname()
use Getopt::Long;        # command line options

# look for perl packages in the following locations
#--------------------------------------------------
use lib ( "$FindBin::Bin", "$FVROOT/bin", "$ESMADIR/$ARCH/bin" );

GetOptions ( "debug",
             "h" );

# FVROOT is where the binaries have been installed
# ------------------------------------------------
$fvroot  = dirname($FindBin::Bin);
$fvroot  =~ s|/u/.realmounts/share|/share|;   # for portability across
                                              # NAS machines

# Initialize variables
# --------------------
  init();

# Execute this program
# --------------------
  run();

#....................................................................................
sub init {

   if ( $#ARGV  <  1 || $opt_h ) {
     print STDERR " Improper input parameters ; see usage:\n";
     usage();
   } else {              # required command line args
     $obsysrc = $ARGV[0];
     $gsigcrc = $ARGV[1];
   }

   die ">>>> ERROR <<< OBSCLASS env missing" unless $ENV{OBSCLASS};
   $obsclass = $ENV{OBSCLASS};

}

#....................................................................................
sub run {


# place name of table in file
open(MYJUNK,">>$gsigcrc") or
die ">>> ERROR <<< cannot write $gsigcrc";

@b = split(/,/, $obsclass);
foreach $oclass ( @b ) {
if ( "$oclass" =~ "pre-qc" ) {
 $row = ( `grep BEGIN $obsysrc | grep -v "#" | cut -c6- | grep $oclass` );
 $row =~ s/\s+$//; #remove trailing spaces
print  MYJUNK <<"EOF";
$row
EOF
}
}
 close(MYJUNK);

}

#....................................................................................
sub usage {

   print <<"EOF";

NAME
     extract_pre-qc.pl - create file with names of 'pre-qc' datasets

SYNOPSIS

     extract_pre-qc.pl [options] obsys.rc output-file

DESCRIPTION

     Based on the observation class and in the available classes
     in obsys.rc, this creates a file with templates for 'pre-qc'
     observation datasets for the experiment in question

    Optional Arguents:

    -h      help (echoes this usage)

    Required Env Vars:
      OBSCLASS  - set of obsclasses as defined in main ADAS, e.g.,
                  "merra2_cdas_pre-qc_bufr,merra2_ascat_pre-qc_bufr"

EOF

  exit(1);
}

