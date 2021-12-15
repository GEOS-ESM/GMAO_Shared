#!/usr/bin/env perl
#=======================================================================
# name - echorc.pl
# purpose - This is an alternative to the echorc.x program. 
#=======================================================================
use strict;
use warnings;

use File::Basename qw(basename);
use Getopt::Long qw(GetOptions);

# global variables
#-----------------
my (@tvals, $rcfile, $rcstring, $numcol, $fillFLG);
my (%month, %lastday);

%month = ("01" => "jan", "02" => "feb", "03" => "mar", "04" => "apr",
          "05" => "may", "06" => "jun", "07" => "jul", "08" => "aug",
          "09" => "sep", "10" => "oct", "11" => "nov", "12" => "dec" );

%lastday = ( 1 => 31,  2 => 0,   3 => 31,  4 => 30,
             5 => 31,  6 => 30,  7 => 31,  8 => 31,
             9 => 30, 10 => 31, 11 => 30, 12 => 31 );

#=======================================================================
# main program
#=======================================================================
{
    my ($value, $expid, $yyyy, $yy, $mm, $dd, $hh, $nn, $ss, $jjj);

    init();
    #----------------------
    #--echorc.x -rc test.rc -ncol 2 OBS_INPUT
    #--dtype:  ps  t  t  q
    #----------------------
    $value = getvalue();
    if (@tvals) {
        $expid = $tvals[0];
        ($yyyy, $yy, $mm, $dd) = ( $tvals[1] =~ m/^(\d{2}(\d{2}))(\d{2})(\d{2})$/ );
        ($hh, $nn, $ss) = ( $tvals[2] =~ m/^(\d{2})(\d{2})(\d{2})$/ );
        $jjj = julianday($yyyy, $mm, $dd);

        $value =~ s/%s/$expid/g;
        $value =~ s/%y4/$yyyy/g;
        $value =~ s/%y2/$yy/g;
        $value =~ s/%m2/$mm/g;
        $value =~ s/%m3/$month{$mm}/g;
        $value =~ s/%d2/$dd/g;
        $value =~ s/%h2/$hh/g;
        $value =~ s/%n2/$nn/g;
        $value =~ s/%j3/$jjj/g;
        $value =~ s/%c/?/g;
    }
    print "$value\n";
}

#=======================================================================
# name - julianday
# purpose -
#=======================================================================
sub julianday {
    my ($yyyy, $mm, $dd, $jjj, $month);

    $yyyy = shift @_;
    $mm = shift @_;
    $dd = shift @_;

    $lastday{2} = 28;
    $lastday{2} = 29 if $yyyy%4==0 and ($yyyy%100!=0 or $yyyy%400==0);

    $jjj = 0;
    foreach $month (1..$mm-1) { $jjj += $lastday{$month} }
    $jjj += $dd;
    $jjj = sprintf("%03d", $jjj);

    return $jjj;
}

#=======================================================================
# name - init
# purpose - get runtime parameters and options
#=======================================================================
sub init {
    my ($help, $debug);

    GetOptions( "template=s{3}" => \@tvals,
                "rc=s"          => \$rcfile,
                "ncol=i"        => \$numcol,
                "fill"          => \$fillFLG,
                "db|debug"      => \$debug,
                "help|h"        => \$help );
    usage() if $help;
    usage() unless scalar(@ARGV) == 1;
    $rcstring = shift @ARGV;
    $rcstring =~ s/\$PESTOROOT/$ENV{"PESTOROOT"}/;

    unless ($fillFLG) { $rcfile = "fvpsas.rc" unless $rcfile }
    $numcol = 1 unless $numcol;

    if (@tvals) {
        die "Error. Bad format for yyyymmdd, $tvals[1];"
            unless $tvals[1] =~ /^\d{8}$/;
        die "Error. Bad format for hhmmss, $tvals[2];"
            unless $tvals[2] =~ /^\d{6}$/;
    }
    checkinput() if $debug;
    usage() if $help;
}

#=======================================================================
# name - getvalue
# purpose - get a variable's value
#=======================================================================
sub getvalue {
    my ($variable, $line, $val, $value, @vals);

    return $rcstring if $fillFLG;

    die "Error. No rcfile specified;" unless $rcfile;
    die "Error. Cannot find rcfile, $rcfile;" unless -e $rcfile;

    open (RCF, "< $rcfile") or die "Error opening rcfile, $rcfile;";
    while (<RCF>) {

        # if namelist found, then strip out numcol
        #-----------------------------------------
        if (m/^\s*$rcstring\s*::/) {
            @vals = ();
            foreach $line (<RCF>) {
                $line =~ s/|^\s*|\s*$|//g; # remove leading and trailing blanks
                next if $line =~ m/^\!/;   # skip comment lines
                last if $line =~ m/^::$/;  # end of namelist

                $val = (split /\s+/, $line)[$numcol-1];
                push @vals, $val;
            }
            $value = "@vals";
        }

        # else, get variable's value from file
        #-------------------------------------
        else {
            $value = $1 if m/^\s*$rcstring\s*:\s*(.*)\s*$/;
        }
        last if $value;
    }
    close RCF;
    die "Error. Could not find $rcstring in $rcfile;" unless $value;
    return $value;
}

#=======================================================================
# name - checkinput
#=======================================================================
sub checkinput {
    if (@tvals) {
        print "check: \@tvals =";
        foreach (@tvals) { print " $_" }
        print "\n";
    }
    print "check: \$rcfile = $rcfile\n" if $rcfile;
    print "check: \$rcstring = $rcstring\n";
    print "check: \$numcol = $numcol\n" if $numcol;
    print "check: \$fillFLG = $fillFLG\n" if $fillFLG;
}

#=======================================================================
# name - usage
# purpose - print usage information
#=======================================================================
sub usage {
    my ($FH, $info);
    my $script = basename($0);

    $FH = select;
    open(INFO, ">", \$info);

    select INFO;
    print << "EOF";

Usage : $script [-template expid nymd nhms] [-rc rcfile|-fill] rcstring

where
 rcstring is string in rcfile whose value is to be echoed 
          or a GRADS template to be filled if the -fill flag is present
options:
 -template expid nymd nhms
                       values to fill into %s, the GRADS template;
                       where rcstring is the the GRADS template (with -fill option)
                       or rcstring\'s value extracted from rcfile is the GRADS template

 -rc rcfilename        rcfile name (default: fvpsas.rc)
 -ncol numcol          when reading table, allows getting given column (default: 1)
 -fill                 treat rcstring as GRADS template to use with -template option

Examples
--------
1. (using rc file)
    > $script -template myexp 20040101 060000 -rc myfile.rc templatename 

    result:  myexp.myfile.something.20040101_06z.bin 

2. (using default, fvpsas.rc)
    > $script variableName

    result: variable values

3. (direct template fill): 
    > $script -template myexp 20040101 060000 -fill %s.myfile.something.%y4%m2%d2_%h2z.bin 

    result:  myexp.myfile.something.20040101_06z.bin 

4. (-ncol option)
    > $script -rc gsi.rc -ncol 2 OBS_INPUT

    result: ps t t q

where

    myfile.rc contains the following line:
    -------------------------------------
templatename: %s.myfile.something.%y4%m2%d2_%h2z.bin 

    fvpsas.rc contains the following line:
    -------------------------------------
variableName: variable values

    gsi.rc contains the following lines:
    -----------------------------------
OBS_INPUT::
!  dfile          dtype       dplat       dsis                  dval    dthin  dsfcalc  obsclass
   prepbufr       ps          null        ps                    0.0     0      0        ncep_prep_bufr 
   prepbufr       t           null        t                     0.0     0      0        ncep_prep_bufr
   prepbufr_profl t           null        t                     0.0     0      0        ncep_acftpfl_bufr
   prepbufr       q           null        q                     0.0     0      0        ncep_prep_bufr
::
EOF
;
    close INFO;
    select $FH;

    system("echo \"$info\" | more");
    exit();
}
