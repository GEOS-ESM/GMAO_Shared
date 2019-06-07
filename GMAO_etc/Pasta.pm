package Pasta;
#=======================================================================
#  Revision History:
#  ----------------
#  26dec1999  da Silva  First version in c code
#  27dec1999  da Silva  Changed name from pesto_atmos to pasta
#  06jan2005  Owens     Added "%c" wildcard and some error trapping (PR 1260)
#  23jan2007  Owens     Added support for %y4 output with %y2 input
#  27jan2011  Stassi    Converted code to perl
#=======================================================================
use strict;
use File::Basename;

require 5.000;
require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw( pasta );

# global variables
#-----------------
my ($y4, $y2, $m2, $m3, $d2, $h2, $expid, $gradstr, $n2); # return values
my ($j3);
my $name = basename $0;

my %month = ("01" => "jan", "02" => "feb", "03" => "mar", "04" => "apr",
             "05" => "may", "06" => "jun", "07" => "jul", "08" => "aug",
             "09" => "sep", "10" => "oct", "11" => "nov", "12" => "dec" );
my %rmonth = reverse(%month);
my %lastday = ( 1 => 31,  2 => 28,  3 => 31,  4 => 30,
                5 => 31,  6 => 30,  7 => 31,  8 => 31,
                9 => 30, 10 => 31, 11 => 30, 12 => 31 );

#=======================================================================
# name - pasta
# purpose - use tokens in template to extract information from string
#=======================================================================
sub pasta {
    my ($string, $template, $numargs, $arg);

    # get calling parameters
    #-----------------------
    die ">> Error << wrong number of calling parameters;" if scalar(@_) != 2;
    ($string, $template) = @_;

    # extract information 
    #--------------------
    extract_vals($string, $template);
    fill_in_missing_vals();

    return ($y4, $y2, $m2, $m3, $d2, $h2, $expid, $gradstr, $n2);
}

#=======================================================================
# name - extract_vals
# purpose - use tokens in template to extract values from string
#
# input parameter
# => $string: input string
# => $template: template containing tokens
#=======================================================================
sub extract_vals {
    my ($string, $template);
    my (@specialchars, $sc, $tc);
    my ($fmt, $token, @list, @match, $label, $val);

    # input parameters
    #-----------------
    ($string, $template) = @_;

    # add slash in front of special matching characters in template
    #--------------------------------------------------------------
    $fmt = "";
    @specialchars = ("+", "?", ".", "*", "^", "(", ")", "[", "{", "|", "\\");
    outer: foreach (0..length($template)-1) {
        $tc = substr($template, $_, 1);

        foreach $sc (@specialchars) {
            if ($tc eq $sc) { $fmt .= "\\$tc"; next outer }
        }
        $fmt .= $tc;
    }
    $template = $fmt;

    # substitute perl character class patterns for tokens in the template
    #--------------------------------------------------------------------
    for (0..length($template)-2) {

        # 3-character tokens: (%y4, %y2, ... )
        #-------------------------------------
        $token = substr($template, $_, 3);

        if ($token eq "%y4") {
            push @list, "y4";
            $fmt =~ s/%y4/(\\d{4})/;
            next;
        }
        if ($token eq "%Y4") {
            push @list, "Y4";
            $fmt =~ s/%Y4/(\\d{4})/;
            next;
        }
        if ($token eq "%y2") {
            push @list, "y2";
            $fmt =~ s/%y2/(\\d{2})/;
            next;
        }
        if ($token eq "%Y2") {
            push @list, "Y2";
            $fmt =~ s/%Y2/(\\d{2})/;
            next;
        }
        if ($token eq "%m2") {
            push @list, "m2";
            $fmt =~ s/%m2/(\\d{2})/;
            next;
        }
        if ($token eq "%M2") {
            push @list, "M2";
            $fmt =~ s/%M2/(\\d{2})/;
            next;
        }
        if ($token eq "%d2") {
            push @list, "d2";
            $fmt =~ s/%d2/(\\d{2})/;
            next;
        }
        if ($token eq "%D2") {
            push @list, "D2";
            $fmt =~ s/%D2/(\\d{2})/;
            next;
        }
        if ($token eq "%j3") {
            push @list, "j3";
            $fmt =~ s/%j3/(\\d{3})/;
            next;
        }
        if ($token eq "%J3") {
            push @list, "J3";
            $fmt =~ s/%J3/(\\d{3})/;
            next;
        }
        if ($token eq "%h2") {
            push @list, "h2";
            $fmt =~ s/%h2/(\\d{2})/;
            next;
        }
        if ($token eq "%H2") {
            push @list, "H2";
            $fmt =~ s/%H2/(\\d{2})/;
            next;
        }
        if ($token eq "%n2") {
            push @list, "n2";
            $fmt =~ s/%n2/(\\d{2})/;
            next;
        }
        if ($token eq "%N2") {
            push @list, "N2";
            $fmt =~ s/%N2/(\\d{2})/;
            next;
        }
        if ($token eq "%s2") {
            push @list, "s2";
            $fmt =~ s/%s2/(\\d{2})/;
            next;
        }
        if ($token eq "%S2") {
            push @list, "S2";
            $fmt =~ s/%S2/(\\d{2})/;
            next;
        }
        if ($token eq "%mc") {
            push @list, "mc";
            $fmt =~ s/%mc/(\\w{3})/;
            next;
        }
        if ($token eq "%Mc") {
            push @list, "Mc";
            $fmt =~ s/%Mc/(\\w{3})/;
            next;
        }
        if ($token eq "%MC") {
            push @list, "MC";
            $fmt =~ s/%MC/(\\w{3})/;
            next;
        }

        # 2-character tokens: (%s, %c, and %n)
        #-------------------------------------
        $token = substr($template, $_, 2);

        if ($token eq "%s") {
            push @list, "expid";
            $fmt =~ s/%s/([^\\.]+)/; # accept any character except "."
            next;
        }
        if ($token eq "%c") {
            push @list, "char";
            $fmt =~ s/%c/(.)/; # accept any character
            next;
        }
        if ($token eq "%n") {
            push @list, "digit";
            $fmt =~ s/%n/(\\d)/; # accept any numeric
            next;
        }
    }

    # extract matched values; verify that expected number are extracted
    #------------------------------------------------------------------
    @match = ($string =~ /$fmt/);

    if (@list) {
        die "$name: error parsing token: $string $template;"
            unless scalar(@match) == scalar(@list);
    }

    # assign extracted values to variables
    #-------------------------------------
    init_vars_as_undef(\$expid, \$y4, \$y2, \$m2, \$m3, \$j3, \$d2, \$h2, \$n2);
    foreach $label (@list) {
        $val = shift @match;

        # if any variable assigned more than once, then take first
        #---------------------------------------------------------
        if ($label eq "expid") { $expid = $val unless $expid }
        elsif ($label eq "y4") {    $y4 = $val unless $y4    }
        elsif ($label eq "y2") {    $y2 = $val unless $y2    }
        elsif ($label eq "m2") {    $m2 = $val unless $m2    }
        elsif ($label eq "d2") {    $d2 = $val unless $d2    }
        elsif ($label eq "h2") {    $h2 = $val unless $h2    }
        elsif ($label eq "n2") {    $n2 = $val unless $n2    }
        elsif ($label eq "j3") {    $j3 = $val unless $j3    }

        if (lc($label) eq "mc") {
            die "Error: Unknown month label: $val;" unless $rmonth{lc($val)};
        }
        if ($label eq "mc") {
            unless ($m2 or $m3) {
                if ($rmonth{$val}) {
                    $m3 = $val;
                    $m2 = $rmonth{$val};
                }
                else { die "Error: month not lowercase: $val;" }
            }
        }
        elsif ($label eq "Mc") {
            unless ($m2 or $m3) {
                if ($val eq ucfirst(lc($val)) and $rmonth{lc($val)}) {
                    $m3 = $val;
                    $m2 = $rmonth{lc($val)};
                }
                else { die "Error: month not Capitalized: $val;" }
            }
        }
        elsif ($label eq "MC") {
            unless ($m2 or $m3) {
                if ($val eq uc($val) and $rmonth{lc($val)}) {
                    $m3 = $val;
                    $m2 = $rmonth{lc($val)};
                }
                else { die "Error: month not UPPERCASE: $val;" }
            }
        }
    }
}

#=======================================================================
# name - init_vars_as_undef
# purpose - initialize variables as undefined
#
# input: array of variable addresses
#=======================================================================
sub init_vars_as_undef { foreach (@_) { $$_ = undef }; return }

#=======================================================================
# name - fill_in_missing_vals
# purpose - fill in values that were not explicitly extracted from $string
#=======================================================================
sub fill_in_missing_vals {
    my ($cnt, $mm, $dd);

    # get $y2 from $y4
    #-----------------
    if ($y4 and ! $y2) { $y2 = substr($y4, 2, 2) }

    # get $y4 from $y2
    #--------------------------------------------------------------------------#
    # NOTE: After Y2K, nobody should be using 2-digit year representations !!! #
    #--------------------------------------------------------------------------#
    if ($y2 and ! $y4) {

        # use 1940 cutoff for consistency with c version of code
        #-------------------------------------------------------
        if ($y2 > 40) { $y4 = "19$y2" }
        else          { $y4 = "20$y2" }

        # uncomment for alternative algorithm
        #------------------------------------
        #--$y4 = y2_to_y4($y2);
    }

    # use $y4 and $j3 to determine $m2 and $d2
    #-----------------------------------------
    if ($y4 and $j3) {
        unless ($m2 and $d2) {
            if ($y4%4==0 and ($y4%100!=0 or $y4%400==0)) {
                $lastday{2} = 29;
            }
            $cnt = 0;
          mloop: for $mm (1..12) {
              for $dd (1..$lastday{$mm}) {
                  $m2 = sprintf("%02d", $mm);
                  $d2 = sprintf("%02d", $dd);
                  $cnt += 1;
                  last mloop if $cnt == $j3;
              }
          }
        }
    }

    # use $m2 to determine $m3
    #-------------------------
    unless ($m3) { $m3 = $month{$m2} if $m2 }

    # fill in defaults for missing values
    #------------------------------------
    $y4 = "yryr" unless $y4;
    $y2 = "yr"   unless $y2;
    $m2 = "mo"   unless $m2;
    $m3 = "mon"  unless $m3;
    $d2 = "dy"   unless $d2;
    $h2 = "hr"   unless $h2;
    $n2 = "00"   unless $n2;

    $expid = "expid" unless $expid;

    # assemble grads string from parts
    #---------------------------------
    $gradstr = "${h2}:${n2}Z${d2}${m3}${y4}";
}

#=======================================================================
# name - y2_to_y4
# purpose - convert 2-digit year representation to 4 digits
#
# input parameter
# => $yr2: 2-digit year
#
# return value
# => $yr4: 4-digit year
# => $delta: (optional) delta value to use in century determination
#
# Note:
# This routine uses the following algorithm
#   if $yr2 < (current year + $delta), then choose current century
#   if $yr2 > (current year + $delta), then choose previous century
#     where current year = realtime (i.e. time of running code)
#=======================================================================
sub y2_to_y4 {
    my ($yr2, $delta, $y4);
    my (@timeData, $currYR, $century, $currYR2);

    ($y2, $delta) = @_;
    $delta = 20 unless $delta;

    @timeData = localtime(time);
    $currYR  = 1900 + $timeData[5];
    $century = substr($currYR,0,2);
    $currYR2 = substr($currYR,2,2);

    if ($y2 > $currYR2 + $delta) { $y4 = ($century-1) .$y2 }
    else                         { $y4 = $century .$y2     }
}
1;
