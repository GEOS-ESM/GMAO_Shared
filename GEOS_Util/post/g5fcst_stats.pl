#!/usr/bin/env perl
#=======================================================================
# name - g5fcst_stats.pl 
# purpose - script to submit jobs to calculate forecast statistics
#
# key global variables -
# => $storedir: directory where output files get copied
#=======================================================================
use warnings;
use strict;

use FindBin qw($Bin);
use lib "$Bin";

use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Path qw(mkpath rmtree);
use Manipulate_time qw(tick);
use WriteLog qw(chdir_ query setprompt);

# global variables
#-----------------
my ($EXP_ARCHIVE);
my ($EXP_DATAMOVE_CONSTRAINT);
my ($anadir, $arcfile, $archiveFLG, $dasFLG, $nodes, $dryrun, $etcdir);
my ($expid, $fhours, $fhours_dflt, $fs_tag, $fv_etcdir);
my ($idate, $ihh, $jobdir, $localID, $ncsuffix, $ndays, $noprompt, $nver);
my ($offset_sec, $pesto, $progdir, $progtype, $secs_per_day, $secs_per_hour);
my ($statsX, $storedir, $tau_freq, $tau_fsec, $vanadir, $vanatype, $vexpid);
my (%opts, @acqIDs, @statsIDs);
my ($landonly, $landonly_dflt);
my ($nxonly, $nxonly_dflt);

$landonly_dflt = "no";
$nxonly_dflt = "no" ;
$fhours_dflt = 123;
$localID = $$;
$ncsuffix = ".nc4";
$opts{"verbose"} = 1;

$secs_per_hour = 60 * 60;
$secs_per_day  = 24 * $secs_per_hour;

# main program
#-------------
{
    my ($afile, $aname, $checkFLG, $climfilecnt, $cmd, $dd, $dd1, $dmls);
    my ($fdir, $file, $getlist, $ffile, $ffile0, $fname, $fstatswork);
    my ($mm, $mm1, $ndate, $ndd, $ndy, $nmm, $ntime, $nv, $nyyyy);
    my ($pid, $vdate, $vdate0, $vdate1, $vdd, $vhh, $vhh0, $vhh1);
    my ($vmm, $vmn, $vtime, $vyyyy, $yyyy, $yyyy1);
    my (%args, %fcsthash, %pesto_dirs, %vanahash);
    my (@alist, @ana_fnames, @climfiles, @fcst_fnames, @fcstlist, @fcstget);
    my (@pesto_dir_list, @rmlist, @vanalist, @vanaget);

    init();

    # key date and time variables
    #-----------------------------------------------
    # idate, ihh
    # ----------
    # > initial forecast date and time
    #
    # ndate, ntime
    # ------------
    # > current forecast date and time
    #   (when multiple days are being evaluated, i.e. ndays > 1)
    #
    # vdate, vtime, vhh
    # -----------------
    # > current verification date and time of current forecast
    #
    # vdate0, vhh0
    # ------------
    # > first verification date and time of first forecast;
    #   used for labeling the archive log file
    #
    # vdate1, vhh1
    # ------------
    # > first verification date and time of current forecast
    #   (when multiple days are being evaluated, i.e. ndays > 1)
    #   used for labeling outputs and log files
    #-----------------------------------------------

    # initial forecast date and time
    #-------------------------------
    $ndate = $idate;
    $ntime = "${ihh}0000";

    $vdate0 = 0;
    $vhh0 = 0;

    %pesto_dirs = ();

    write_g5fcst_stats_arc();
    verify_values();

    # loop over number of days to process
    #------------------------------------
    foreach $ndy (1..$ndays) {
        if ($dasFLG) {
            $checkFLG = das_check($ndate, $ntime);
            next unless $checkFLG;
        }

        # use a separate work directory for each day
        #-------------------------------------------
        $fstatswork = $ENV{"NOBACKUP"} ."/FSTATSWORK.$localID.$ndy";
        die "Work directory already exists;" if -e $fstatswork;
        mkpath($fstatswork, \%opts);

        # loop over verification hours
        #-----------------------------
        %fcsthash = ();
        %vanahash = ();

        # initial verification date and time for current forecast
        #--------------------------------------------------------
        ($vdate, $vtime) = tick($ndate, $ntime, $offset_sec);

        # save initial verification date/time from current forecast
        # for labeling output and log files
        #----------------------------------
        $vdate1 = $vdate;
        ($vhh1, $vmn) = extract_hh_mn($vtime);

        # save initial verification date/time from first forecast
        # for labeling archive log
        #-------------------------
        $vdate0 = $vdate1 unless $vdate0;
        $vhh0 = $vhh1 unless $vhh0;

        ($nyyyy, $nmm, $ndd) = extract_yyyy_mm_dd($ndate);

        foreach $nv (1..$nver) {
            ($vyyyy, $vmm, $vdd) = extract_yyyy_mm_dd($vdate);
            ($vhh, $vmn) = extract_hh_mn($vtime);

            if ($ihh eq "00" and $nv == 1) {
                $fdir = "$vanadir/Y$nyyyy/M$nmm";
                if ($nxonly eq "yes") {
                   $fname = "$vexpid.$vanatype.${vdate}_${vhh}${vmn}z";
                } else {
                   $fname = "$vexpid.$vanatype.inst3d_met_p.${vdate}_${vhh}${vmn}z";                 
                }     
                $ffile0 = "$fdir/$fname$ncsuffix";
            }
            else {

                # allows diffing with analysis (doing say EC-ana vs G5-ana)
                #----------------------------------------------------------
                if ($progtype eq "ana") {
                    $fdir  = "$anadir/Y$vyyyy/M$vmm";
                    $fname = "$expid.ana.inst3d_met_p.${vdate}_${vhh}${vmn}z";
                }
                else {
                    $fdir  = "$progdir/Y$nyyyy/M$nmm/D$ndd/H$ihh";
                    $fname = "$expid.prog.$progtype.${ndate}_${ihh}z+${vdate}_${vhh}${vmn}z";
                }
            }
            $ffile = "$fdir/$fname$ncsuffix";
            unless (-e $ffile) {
                rmtree($fstatswork) if -d $fstatswork;
                die "Error. Cannot find forecast : $ffile;";
            }
            $fcsthash{$ffile} = 1;

            # multiple naming convention options
            #-----------------------------------
            @alist = ();
            if ($vexpid eq "gfs" or $vexpid eq "ecmwf" or $vexpid eq "era5") {
                push @alist, "$vexpid.$vanatype.${vdate}_${vhh}00z";
                push @alist, "$vexpid.$vanatype.${vdate}_${vhh}z+${vdate}_${vhh}z";
                push @alist, "$vexpid.$vanatype.${vdate}_${vhh}z+${vdate}_${vhh}${vmn}z";
            }
            else {
                if ($nxonly eq "yes" ) {
                   push @alist, "$vexpid.$vanatype.${vdate}_${vhh}${vmn}z";
                } else {
                   push @alist, "$vexpid.$vanatype.inst3d_met_p.${vdate}_${vhh}${vmn}z"; 
                }
            }

            # loop through naming convention options
            #---------------------------------------
            {
                $aname = shift @alist;
                $afile = "$vanadir/Y$vyyyy/M$vmm/$aname$ncsuffix";
                if (-e $afile) {
                    print "Found ver analysis : $afile\n";
                }
                else {
                    print "\nWarning. Cannot find ver analysis : $afile\n";
                    redo if @alist;

                    rmtree($fstatswork) if -d $fstatswork;
                    die "Error. Cannot find ver analysis : $afile;";
                }
            }
            $vanahash{$afile} = 1;

            # increment lead-time of forecast (verification time)
            #----------------------------------------------------
            ($vdate, $vtime) = tick($vdate, $vtime, $tau_fsec);
        }
        if ($nxonly eq "yes") {
        @climfiles = (<$ENV{SHARE}/dao_ops/verification/stats/MERRA-2.inst1_2d_asm_Nx.198501_201412.clim_??z.576x361.data.nc4>);
        } else {
        @climfiles = (<$ENV{SHARE}/dao_ops/verification/stats/MERRA-2.inst3_3d_asm_Np.198501_201412.clim_??z.576x361.data.nc4>);
        #--@climfiles = (<$ENV{SHARE}/dao_ops/verification/stats/merrasc.197901-200812.clim_??z.144x91.data.nc>);
        }
        $climfilecnt = scalar(@climfiles);
        if ($climfilecnt < 4) {
            rmtree($fstatswork) if -d $fstatswork;
            die "Error. Only $climfilecnt out of 4 Climatology Files found - exiting;"
        }

        # get fcst and ana filenames
        #---------------------------
        if ($vexpid eq "g5ncep" or $vexpid eq "gfs" or $vexpid eq "ecmwf" or $vexpid eq "era5") {
            if ($progtype eq "ana") {
                push @fcst_fnames, "$expid.ana.inst3d_met_p.*$ncsuffix";
                push @ana_fnames, "$vexpid.$vanatype.*$ncsuffix";
            }
            else {
                push @fcst_fnames, "$expid.prog.$progtype.${ndate}_${ihh}z+*$ncsuffix";
                push @ana_fnames, "$vexpid.$vanatype.*$ncsuffix";
            }
        }
        elsif ($vanatype eq "asm") {
            if ($progtype eq "ana") {
                push @fcst_fnames, $ffile0 if $ffile0;
                push @fcst_fnames, "$expid.ana.inst3d_met_p.*$ncsuffix";
                push @ana_fnames, "$vexpid.asm.inst3d_met_p.*$ncsuffix";
            }
            else {
                push @fcst_fnames, $ffile0 if $ffile0;
                push @fcst_fnames, "$expid.prog.$progtype.${ndate}_${ihh}z+*$ncsuffix"; 
                if ($nxonly eq "yes" ) {
                     push @ana_fnames, "$vexpid.$vanatype.*$ncsuffix";
                } else {
                     push @ana_fnames, "$vexpid.asm.inst3d_met_p.*$ncsuffix";
                }
            }
        }
        else {
            push @fcst_fnames, $ffile0 if $ffile0;
            push @fcst_fnames, "$expid.prog.$progtype.${ndate}_${ihh}z+*$ncsuffix";
            if ($nxonly eq "yes" ) {
               push @ana_fnames, "$vexpid.$vanatype.*$ncsuffix" 
            } else {
               push @ana_fnames, "$vexpid.asm.inst3d_met_p.*$ncsuffix";
               push @ana_fnames, "$vexpid.ana.inst3d_met_p.*$ncsuffix";
            }
        }
        @fcstlist = sort keys %fcsthash;
        @vanalist = sort keys %vanahash;

        # dmget inputs
        #-------------
        @fcstget = ();
        foreach (@fcstlist) {
            $dmls = `dmls -lL $_`;
            push @fcstget, $_ if $dmls =~ / \(OFL\) /;
        }

        @vanaget = ();
        foreach (@vanalist) {
            $dmls = `dmls -lL $_`;
            push @vanaget, $_ if $dmls =~ / \(OFL\) /;
        }

        foreach $getlist (\@fcstget, \@vanaget) {
            next unless @$getlist;
            $cmd = "dmget @$getlist";
            print "\n$cmd\n";
            defined ($pid = fork) or die "Error. fork failed: $!";
            unless ($pid) {
                exec $cmd;
                die "Error. exec dmget failed: $!";
            }
        }

        # fetch inputs
        #-------------
        print "\nCopying files to directory: $fstatswork\n";
        foreach $file (@fcstlist, @vanalist) {
            print "=> " .basename($file) ."\n";
            copy $file, $fstatswork or die "Error. copy failed: $file; $!";
        }

        $nodes = "null" unless $nodes;

        # store calling arguments in a hash
        #----------------------------------
        $args{"ana_fnames_addr"}  = \@ana_fnames;
        $args{"climfiles_addr"}   = \@climfiles;
        $args{"fcst_fnames_addr"} = \@fcst_fnames;
        $args{"fcstlist_addr"}    = \@fcstlist;
        $args{"vanalist_addr"}    = \@vanalist;
        $args{"fstatswork"}       = $fstatswork;
        $args{"vdate0"}           = $vdate0;
        $args{"vdate1"}           = $vdate1;
        $args{"vhh0"}             = $vhh0;
        $args{"vhh1"}             = $vhh1;
        $args{"nodes"}            = $nodes;
        $args{"landonly"}         = $landonly;
        $args{"nxonly"}           = $nxonly;

        # submit job to calculate stats
        #------------------------------
        submit_calcjob(%args);

        # make list of search directories for archive job
        #------------------------------------------------
        ($yyyy1, $mm1, $dd1) = extract_yyyy_mm_dd($vdate1);
        $pesto_dirs{"prog/$fs_tag/Y$yyyy1/M$mm1"} = 1;
        $pesto_dirs{"etc/Y$yyyy1/M$mm1"} = 1;

        # increment to next day to be evaluated
        #--------------------------------------
        ($ndate, $ntime) = tick($ndate, $ntime, $secs_per_day);
    }

    # submit archive job
    #-------------------
    @pesto_dir_list = sort keys %pesto_dirs;
    $args{"pesto_dir_list_addr"} = \@pesto_dir_list;
    submit_archivejob(%args);
}

#=======================================================================
# name - init
# purpose - get runtime parameters and options
#=======================================================================
sub init {
    use Getopt::Long qw(GetOptions);
    my ($help, $offset, %opts);

    $ENV{"PATH"} = ".:$Bin:$ENV{PATH}";
    do "g5_modules_perl_wrapper";

    $EXP_ARCHIVE = $ENV{"ARCHIVE"};
    $EXP_ARCHIVE = $ENV{"FVARCH"} if $ENV{"FVARCH"};

    $EXP_DATAMOVE_CONSTRAINT = $ENV{"DATAMOVE_CONSTRAINT"};

    # flush buffer after each output operation
    #-----------------------------------------
    $| = 1;

    # get runtime options
    #--------------------
    GetOptions("ihh=i"     => \$ihh,
               "fhrs=i"    => \$fhours,

               "anadir=s"  => \$anadir,
               "progdir=s" => \$progdir,

               "vexpid=s"  => \$vexpid,
               "vanadir=s" => \$vanadir,

               "ptype=s"   => \$progtype,
               "vtype=s"   => \$vanatype,

               "storedir=s" => \$storedir,
               "archive!"   => \$archiveFLG,

               "nodes=s"  => \$nodes,

               "landonly=s" => \$landonly,

               "nxonly=s"   => \$nxonly,

               "das"      => \$dasFLG,
               "np"       => \$noprompt,

               "dryrun"   => \$dryrun,
               "h|help"   => \$help);
    usage() if $help;

    $noprompt = 1 if $dasFLG;
    setprompt(0) if $noprompt;

    # get runtime parameters
    #-----------------------
    $expid = shift @ARGV;
    $idate = shift @ARGV;
    $ndays = shift @ARGV;
    
    # landonly option
    #----------------
    $landonly = $landonly_dflt unless $landonly;
  
    #nxonly option
    #--------------
    $nxonly = $nxonly_dflt unless $nxonly;
           
    # initial fcst hour and offset
    #-----------------------------
    unless ($ihh) {
        if ($expid eq "a_flk_04") { $ihh =  0 }
        else                      { $ihh = 21 }
    }
    if ($ihh == 0 or $ihh == 6 or $ihh == 12 or $ihh == 18) { $offset = 0 }
    if ($ihh == 3 or $ihh == 9 or $ihh == 15 or $ihh == 21) { $offset = 3 }

    $ihh = sprintf("%02d", $ihh);
    $offset_sec = $offset * $secs_per_hour;

    # frequecy of verification (hours)
    #---------------------------------
    $tau_freq = 12;
    $tau_fsec = $tau_freq * $secs_per_hour;

    # length of fcst in hours
    #------------------------
    $fhours = $fhours_dflt unless $fhours;
    $fhours -= $offset;

    $nver = int($fhours/$tau_freq) + 1;

    # interactive runtime params
    #---------------------------
    until ($expid) {
        $expid = query("Forecast Experiment ID: ");
    }
    until ($idate and $idate =~ m/^\d{8}$/) {
        print "\nCannot decipher initial date: $idate\n" if $idate;
        $idate = query("initial date <yyyymmdd>: ");
    }
    until ($ndays and $ndays =~ m/^\d+$/) {
        print "\nCannot decipher ndays: $ndays\n" if $ndays;
        $ndays = query("number of days to process <n>: ");
    }
    usage() unless $idate and $ndays and $expid;

    # preliminary check for forecast and das hidden files
    #----------------------------------------------------
    das_check($idate, "${ihh}0000") if $dasFLG;

    # option defaults
    #----------------
    if ($dryrun) { $dryrun = "echo" }
    else         { $dryrun = "" }
    $archiveFLG = 1 unless defined($archiveFLG);

    # determine $storedir (the place where output gets copied prior to archival)
    #---------------------------------------------------------------------------
    $storedir = $ENV{"STOREDIR"} unless $storedir;
    unless ($storedir) { $storedir = dirname($ENV{"FVHOME"}) if $ENV{"FVHOME"} }
    unless ($storedir) { $storedir = $ENV{"NOBACKUP"} if $ENV{"NOBACKUP"} }
    die "Error. Storage directory not defined;" unless $storedir;

    mkpath($storedir, \%opts) unless -d $storedir;
    die "Error creating directory: $storedir;" unless -d $storedir;

    $jobdir = "$storedir/$expid/fstats";
    mkpath($jobdir, \%opts) unless -d $jobdir;

    $etcdir = "$storedir/$expid/etc";
    mkpath($etcdir, \%opts) unless -d $etcdir;

    # get $anadir and $progdir
    #-------------------------
    $anadir = "$EXP_ARCHIVE/$expid/ana" unless $anadir;
    until (-d $anadir) {
        unless ($noprompt) {
            print "ANA directory does not exist: $anadir\n\n";
            $anadir = query("verification ana directory:", $anadir);
        }
        #else { die "verification ana directory does not exist: $anadir;" }
    }
    $progdir = "$EXP_ARCHIVE/$expid/prog" unless $progdir;
    until (-d $progdir) {
        unless ($noprompt) {
            print "forecase prog directory does not exist: $progdir\n\n";
            $progdir = query("forecast prog directory:", $progdir);
        }
        else { die "forecast prog directory does not exist: $progdir;" }
    }

    # get vexpid and vanadir
    #-----------------------
    unless ($vexpid) {
        $vexpid = basename(basename($vanadir)) if $vanadir;
        $vexpid = $expid unless $vexpid;
        $vexpid = query("verifying experiment ID:", $vexpid);
    }
    unless ($vanadir) {
        $vanadir = dirname(dirname($anadir))."/$vexpid/ana";
        $vanadir = $anadir unless -d $vanadir;
        $vanadir = query("verifying experiment ana directory:", $vanadir);
    }
    until (-d $vanadir) {
        print "Cannot find VANA verifying directory: $vanadir\n\n";
        $vanadir = query("verifying experiment ana directory:");
        $vexpid = basename($vanadir);
    }

    # get progtype and vanatype
    #--------------------------
    unless ($progtype) {
        $progtype = "inst3d_met_p";
        $progtype = query("forecast prog type: ", $progtype);
    }
    unless ($vanatype) {
        $vanatype = "asm";
        $vanatype = query("verification data type: ", $vanatype);
    }

    # determine fs_tag
    #--------------------
    $fs_tag = "fstats_${vexpid}_$vanatype";

    if ($vanatype eq "inst3_3d_ana_Np") {
        $fs_tag = "fstats_ana"
    }

    if ($vexpid eq "era5") {
        $fs_tag = "fstats_era5";
        if ($progtype eq "ana") {
            $fs_tag = "astats_era5";
        }
    }
    elsif ($vexpid eq "ecmwf") {
        $fs_tag = "fstats_ecmwf";
        if ($progtype eq "ana") {
            $fs_tag = "astats_ecmwf";
        }
    }
    elsif ($vexpid eq "gfs") {
        $fs_tag = "fstats_gfs";
        if ($progtype eq "ana") {
            $fs_tag = "astats_gfs";
        }
    }
    elsif ($vexpid eq "rsens3") {
        $vanadir = "/archive/u/dndaescu/$vexpid/ana"
    }
    elsif ($vexpid eq "d5124_m2_jan00") {
        $vanadir = "/home/dao_ops/$vexpid/run/.../archive/ana";
    }
    elsif ($expid eq "a_flk_04") {
        $vanatype = "ana";
        $progdir = "$ENV{ARCHIVE}/geos4/$vexpid/prog" unless $progdir;
        $vanadir = "$ENV{ARCHIVE}/geos4/$vexpid/ana";
        $fs_tag = "fstats_$vanatype"
    }

    if ($progtype eq "ana" and $expid eq $vexpid) {
        $fs_tag = "astats_self";
    }

    if ($vexpid eq "d591_rpit3_jan11") {
        $vanadir = "/discover/nobackup/rtodling/d591_rpit3_jan11/ana";
    }

    if ($vexpid eq "g5ncep") {
        $vanadir = "/archive/u/rtodling/g5ncep/$vexpid";
        $fs_tag = "fstats_ncana";
    }

    # find stats.x and pesto
    #-----------------------
    chomp($statsX = `which stats.x`);
    die "Error. Cannot find stats.x program;" unless -e $statsX;

    chomp($pesto = `which pesto`);
    die "Error. Cannot find pesto;" unless -e $pesto;

    # find stats.rc
    #--------------
    $fv_etcdir = dirname($Bin) ."/etc";
    if ($nxonly eq "yes" ) {
        die "Error. Cannot find $fv_etcdir/statsNx.rc" unless -e "$fv_etcdir/statsNx.rc"; 
    } else {
        die "Error. Cannot find $fv_etcdir/stats.rc" unless -e "$fv_etcdir/stats.rc";
    }
}
#=======================================================================
# name - das_check
# purpose - check for hidden files in $FVHOME/fcst and $FVHOME directories
#           left by DAS jobs to determine whether to calculate statistics
#
# input parameters
# => $ndate: forecast date; format, yyyymmdd
# => $ntime: forecast time; format, hhmmss
#=======================================================================
sub das_check {
    my ($ndate, $ntime);
    my ($calculate_stats, $fdatetime, $FVHOME, $dotDONEFCST);
    my ($ddate, $dtime, $ddatetime, $dotDONE, %notfound);

    $ndate = shift @_;
    $ntime = shift @_;

    $calculate_stats = 1;

    # check that FVHOME variable is defined
    #--------------------------------------
    $FVHOME = $ENV{"FVHOME"};
    die "Error. FVHOME environment variable not defined;" unless $FVHOME;
    die "Error. Cannot find FVHOME directory: $FVHOME;" unless -d $FVHOME;

    # check for fcst DONE hidden file
    #-------------------------------------
    $fdatetime = "${ndate}_" .substr($ntime, 0, 2) ."z";
    $dotDONEFCST = "$ENV{FVHOME}/fcst/.DONE_FCST.$fdatetime";
    unless (-e $dotDONEFCST) {
        print "\nCannot find forecast hidden file:\n"
            . "- $dotDONEFCST\n\n";
        $calculate_stats = 0;
    }

    # check for das hidden DONE files
    #--------------------------------
    ($ddate, $dtime) = tick($ndate, $ntime, $offset_sec);
    for (1..$nver) {
        $ddatetime = $ddate .substr($dtime, 0, 2);
        $dotDONE = "$FVHOME/.DONE_CENTRAL_ADAS.$ddatetime";
        $notfound{$dotDONE} = 1 unless -e $dotDONE;
        ($ddate, $dtime) = tick($ddate, $dtime, $tau_fsec);
    }
    if (%notfound) {
        print "Cannot find ana hidden files:\n";
        foreach (sort keys %notfound) { print "- $_\n" }
        $calculate_stats = 0;
    }

    # quit if hidden files not found
    #-------------------------------
    unless ($calculate_stats) {
        print "\nForecast statistics will not be calculated for $fdatetime\n\n";
        exit(1) unless $ndays > 1;
    }
    return $calculate_stats;
}

#=======================================================================
# name - write_g5fcst_stats_arc
# purpose - write the g5fcst_stats.arc file
#=======================================================================
sub write_g5fcst_stats_arc {
    use File::Copy qw(move);
    my ($arcfile_tilde);

    $arcfile = "$jobdir/g5fcst.$fs_tag.arc";
    if (-e $arcfile) {
        $arcfile_tilde = "${arcfile}~";
        if (-e $arcfile_tilde) { unlink $arcfile }
        else                   { move $arcfile, $arcfile_tilde}
    }
    open ARC, "> $arcfile" or die "Error opening ARC file: $arcfile; $!";
    print ARC <<"EOF" or die "Error writing ARC file: $arcfile: $!";
\${PESTOROOT}%s/prog/$fs_tag/Y%y4/M%m2/%s.%H2z.globl.b%y4%m2%d2_%h2z.e%Y4%M2%D2_%H2z.ctl
\${PESTOROOT}%s/prog/$fs_tag/Y%y4/M%m2/%s.%H2z.globl.b%y4%m2%d2_%h2z.e%Y4%M2%D2_%H2z.data
\${PESTOROOT}%s/prog/$fs_tag/Y%y4/M%m2/%s.%H2z.stats.b%y4%m2%d2_%h2z.e%Y4%M2%D2_%H2z.ctl1
\${PESTOROOT}%s/prog/$fs_tag/Y%y4/M%m2/%s.%H2z.stats.b%y4%m2%d2_%h2z.e%Y4%M2%D2_%H2z.ctl2
\${PESTOROOT}%s/prog/$fs_tag/Y%y4/M%m2/%s.%H2z.stats.b%y4%m2%d2_%h2z.e%Y4%M2%D2_%H2z.data
\${PESTOROOT}%s/prog/$fs_tag/Y%y4/M%m2/%s.fstats.log.%y4%m2%d2_%h2z.txt
\${PESTOROOT}%s/etc/Y%y4/M%m2/%s.fstats_calc.log.%y4%m2%d2_%h2z.txt
EOF
;
    close ARC;
}

#=======================================================================
# name - submit_calcjob
# purpose - write and submit jobfile to run stats
#=======================================================================
sub submit_calcjob {
    my (%args, @ana_fnames, @climfiles, @fcst_fnames, @fcstlist, @vanalist);
    my ($fstatswork, $vdate1, $vhh1);
    my (@rmfilelist, $yyyy, $mm, $dd);
    my ($logdir, $logfile1, $logfile2, $jobname, $jobdate, $jobfile, $jobtype);
    my ($cmd, $jobID, $jobIDline);
    my (@levs, @levels_19, @levels_11, @levels_1);
    my ($mynodes,$usrnodes);
    my ($qos, $partition);
    my ($ntspn, $npn);
    my ($landonly, $landmaskdirfile); 
    my ($nxonly);
    my ($whichrc);
    
    @levels_19 = ( 1000.0, 975.0, 950.0, 925.0,
                    900.0, 850.0, 800.0, 750.0,
                    700.0, 600.0, 500.0, 400.0,
                    300.0, 250.0, 200.0, 150.0,
                    100.0,  70.0,  10.0 );

    @levels_11 = ( 1000.0, 925.0, 850.0, 700.0,
                    500.0, 400.0, 300.0, 250.0,
                    200.0, 150.0, 100.0 );

    @levels_1  =  ( 1000.0 ) ;

    if ($vexpid eq "ecmwf") { @levs = @levels_11 }
    else                    { @levs = @levels_19 }

    # input arguments
    #----------------
    %args = @_;

    @ana_fnames  = @{$args{"ana_fnames_addr"}};
    @climfiles   = @{$args{"climfiles_addr"}};
    @fcst_fnames = @{$args{"fcst_fnames_addr"}};
    @fcstlist    = @{$args{"fcstlist_addr"}};
    @vanalist    = @{$args{"vanalist_addr"}};

    $fstatswork  = $args{"fstatswork"};
    $vdate1      = $args{"vdate1"};
    $vhh1        = $args{"vhh1"};
    $usrnodes    = $args{"nodes"};
    $landonly    = $args{"landonly"};
    $nxonly      = $args{"nxonly"};
    print "nxonly: $nxonly\n";
    print "landonly: $landonly\n";

    if ($landonly eq "yes") {
     $landmaskdirfile = "$jobdir/landmaskfile.nc4";
     die "Error. Cannot find $landmaskdirfile" unless -e "$landmaskdirfile";
    }
    
    $whichrc = "stats.rc";
    if ($nxonly eq "yes") {
     @levs = @levels_1;
     $whichrc = "statsNx.rc";
     print "whichrc: $whichrc\n"; 
    }

    foreach (@fcstlist, @vanalist) { push @rmfilelist, basename($_) };

    ($yyyy, $mm, $dd) = extract_yyyy_mm_dd($vdate1);
    $logdir = "$etcdir/Y$yyyy/M$mm";
    mkpath($logdir, \%opts) unless -d $logdir;

    $jobtype = "fstats_calc";
    $jobdate = "${vdate1}_${vhh1}z";

    $jobname = "$jobtype.$jobdate";
    $jobfile = "$jobdir/$jobname.j";
    $logfile1 = "$jobdir/$jobtype.log.$jobdate.o%j.txt";
    $logfile2 = "$logdir/$expid.$jobtype.log.$jobdate.txt";
    $npn = `facter processorcount`; chomp($npn);
    if ( $npn == 40 ) {
      $mynodes = "sky";
      $ntspn   = 36;
      $qos     = "#SBATCH --qos=dastest";        # wired for now since only way to use SKY
      $partition = "#SBATCH --partition=preops"; # wired for now since only way to use SKY
    } elsif ( $npn == 48 ) {
      $mynodes = "cas";
      $ntspn   = 42;
      $qos     = "#SBATCH --qos=dastest";        # wired for now since only way to use CAS
      $partition = "#SBATCH --partition=preops"; # wired for now since only way to use CAS
    } else {
      $mynodes = "hasw";
      $ntspn   = 24;
      $qos     = "";
      $partition = "";
#     $qos     = "#SBATCH --qos=dastest";        # wired for now since only way to use HASW
#     $partition = "#SBATCH --partition=preops"; # wired for now since only way to use HASW
    }
    if ( $usrnodes ne "null" ) { $mynodes = $usrnodes }; # overwrite with specification from command line

    print "\nwriting jobfile: $jobfile\n";
    open FH, "> $jobfile" or die "Error opening $jobfile; $!";
    print FH <<"EOF" or die "Error writing to $jobfile: $!";
#!/usr/bin/csh
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --job-name=$jobname
#SBATCH --output=$logfile1
#SBATCH --export=NONE
#SBATCH --constraint=$mynodes
$qos
$partition

source $Bin/g5_modules
set echo
chdir $fstatswork

# slurm env is messed up (inconsistent)
#unsetenv SLURM_MEM_PER_CPU
unsetenv SLURM_MEM_PER_GPU
unsetenv SLURM_MEM_PER_NODE

if ( \$?I_MPI_ROOT ) then
  set mympi = "mpirun"
else
  set mympi = "mpiexec_mpt"
endif

$dryrun \$mympi -np 1 $statsX -fcst @fcst_fnames \\
                    -ana @ana_fnames \\
                    -cli @climfiles \\
                    -tag $expid.${ihh}z \\
                    -nfreq ${tau_freq}0000 \\
EOF
if ( $landonly eq "yes" ) {
print  FH <<"EOF";
                    -land ${landmaskdirfile} \\
EOF
}

print  FH <<"EOF";
                    -levs @levs \\
                    -o $expid.fstats.log.$jobdate.txt \\
                    -verif gmao \\
                    -fcsrc gmao \\
                    -fhour $fhours \\
                    -rc $fv_etcdir/$whichrc
EOF

print  FH <<"EOF";

@ calc_status = \$status

$pesto -arc $arcfile \\
       -expid $expid \\
       -d $fstatswork \\
       -r $storedir \\
       -l -v -clean
@ calc_status += \$status

if (! \$calc_status) then
   \\rm @rmfilelist
   ls -l
   chdir ~
   \\rm -r $fstatswork
   qalter -o $logfile2 \$SLURM_JOBID
endif

EOF
;
    close FH;

    print "submitting jobfile: $jobfile\n";
    $cmd = "sbatch $jobfile";
    print "> $cmd\n";
    chomp($jobIDline = `$cmd`);

    $jobID = (split /\s+/, $jobIDline)[-1];
    print "jobID = $jobID\n\n";
    push @statsIDs, $jobID;
}

#=======================================================================
# name - submit_archivejob
# purpose - write and submit jobfile for archiving stats output files
#=======================================================================
sub submit_archivejob {
    my (%args, @pesto_dir_list, $vdate0, $vhh0, $dtFLG);
    my ($yyyy, $mm, $dd, $logdir, $logfile1, $logfile2, $pdir);
    my ($jobname, $jobdate, $jobfile_0, $jobfile, $jobtype);
    my ($cmd, $deps, $dependFLG, $jobIDline, $jobID);
    my ($mynodes);

    # input arguments
    #----------------
    %args = @_;
    @pesto_dir_list = @{$args{"pesto_dir_list_addr"}};
    $vdate0 = $args{"vdate0"};
    $vhh0   = $args{"vhh0"};

    ($yyyy, $mm, $dd) = extract_yyyy_mm_dd($vdate0);
    $logdir = "$etcdir/Y$yyyy/M$mm";
    mkpath($logdir, \%opts) unless -d $logdir;

    $jobtype = "fstats_arch";
    $jobdate = "${vdate0}_${vhh0}z";

    $jobname = "$jobtype.$jobdate";
    $jobfile = "$jobdir/$jobname.j";
    $logfile1 = "$jobdir/$jobtype.log.$jobdate.o%j.txt";
    $logfile2 = "$jobdir/$jobtype.log.$jobdate.txt";

    $jobfile_0 = "${jobfile}_0";

    # if processing one day, then archive only one day
    #-------------------------------------------------
    if ($ndays == 1) { $dtFLG = "-date $vdate0 -syntime ${vhh0}0000" }
    else             { $dtFLG = "" }    

    # write archive jobfile
    #----------------------
    print "writing jobfile: $jobfile\n";
    open FH0, "> $jobfile_0" or die "Error opening $jobfile_0; $!";
    print FH0 <<"EOF" or die "Error writing to $jobfile_0: $!";
#!/usr/bin/csh
#SBATCH --time=1:00:00
#SBATCH --job-name=$jobname
#SBATCH --partition=datamove
#$EXP_DATAMOVE_CONSTRAINT
#SBATCH --output=$logfile1
#SBATCH --export=NONE

set echo
@ archive_status = 0
foreach dir ( \\
              )

    $pesto -arc $arcfile \\
           -expid $expid $dtFLG \\
           -d $storedir/$expid/\$dir \\
           -r $EXP_ARCHIVE \\
           -l -v -clean

    @ archive_status += \$status
end

if (! \$archive_status) then
   qalter -o $logfile2 \$SLURM_JOBID
endif
EOF
;
    close FH0;

    # add directories to search
    #--------------------------
    open FH0, "< $jobfile_0" or die "Error opening $jobfile_0: $!;";
    open FH,  "> $jobfile"   or die "Error opening $jobfile: $!;";
    while (<FH0>) {
        print FH $_;
        if (m/foreach dir/) {
            foreach $pdir (@pesto_dir_list) {
                printf FH " " x 13;
                print FH "$pdir \\\n";
            }
        }
    }
    return unless $archiveFLG;
    close FH0;
    close FH;
    unlink($jobfile_0);

    # submit archive job
    #-------------------
    $deps = "";
    foreach (@statsIDs) { $deps .= ":$_" if $_ }

    $dependFLG = "";
    $dependFLG = "--dependency=afterany$deps" if $deps;

    print "submitting jobfile: $jobfile\n";
    $cmd = "sbatch $dependFLG $jobfile";
    print "> $cmd\n";
    chomp($jobIDline = `$cmd`);

    $jobID = (split /\s+/, $jobIDline)[-1];
    print "jobID = $jobID\n\n";
}

#=======================================================================
# name: extract_hh_mn
# purpose: extract hh and mn values from hhmmss string
#
# input parameter:
# => $hhmmss: input time string
#
# return value:
# => $hh: extracted hour value
#=======================================================================
sub extract_hh_mn{
    my ($hhmmss, $hh, $mn);
    $hhmmss = shift @_;

    ($hh, $mn) = ($hhmmss =~ m/^(\d{2})(\d{2})\d{2}$/)
        or die "Error. Undecipherable time: $hhmmss;";
    return $hh, $mn;
}

#=======================================================================
# name: extract_yyyy_mm_dd
# purpose: extract yyyy, mm, and dd values from yyyymmdd string
#
# input parameter:
# => $yyyymmdd: input date string
#
# return values:
# => ($yyyy, $mm, $dd): extracted date values
#=======================================================================
sub extract_yyyy_mm_dd {
    my ($yyyymmdd, $yyyy, $mm, $dd);
    $yyyymmdd = shift @_;

    ($yyyy, $mm, $dd) = ($yyyymmdd =~ m/^(\d{4})(\d{2})(\d{2})$/)
        or die "Error. Undecipherable date: $yyyymmdd;";
    return ($yyyy, $mm, $dd);
}

#=======================================================================
# name - verify_values
# purpose - have user verify job values
#=======================================================================
sub verify_values {
    my ($ans);
    print "\nJob Inputs and Values\n";
    print   "---------------------\n";
    print   "localID: $localID\n";
    print   "---------------------\n";
    print "expid:       $expid\n";
    print "idate:       $idate\n";
    print "ndays:       $ndays\n";
    print "fhours:      $fhours\n\n";

    print "anadir:      $anadir\n";
    print "vanadir:     $vanadir\n";
    print "vanatype:    $vanatype\n";
    print "vexpid:      $vexpid\n\n";

    print "progdir:     $progdir\n";
    print "progtype:    $progtype\n";
    print "fs_tag:      $fs_tag\n\n";
    print "landonly:    $landonly\n\n";
    print "nxonly:      $nxonly\n\n";
    print "dryrun:      $dryrun\n\n" if $dryrun;

    $ans = query("Continue (y/n):", "y");
    if ($ans eq "n") {
        print "Exiting.\n";
        exit();
    }
}

#=======================================================================
# name - usage
# purpose - print usage information
#=======================================================================
sub usage {
    my $script = basename($0);
    print <<"EOF";

NAME
    $script

SYNOPSIS
    $script \$expid \$idate \$ndays [OPTIONS]

PARAMETERS
     expid             forecast experiment ID
     idate             initial date of forecast; format: yyyymmdd
     ndays             number of days to process

OPTIONS [defaults in brackets]
    -ihh ihh           initial hour of forecast [21]
    -fhrs fhours       forecast length in hrs [$fhours_dflt]

    -anadir anadir     forecast ana directory [\$ARCHIVE/\$expid/ana]
    -progdir progdir   forecast prog directory [\$ARCHIVE/\$expid/prog]

    -vexpid vexpid     verifying experiment ID [basename(basename(\$vanadir)) or \$expid]
    -vanadir vanadir   verifying experiment ana directory [basename(\$anadir)/\$vexpid]

    -ptype progtype    forecast prog type [inst3d_met_p]
    -vtype vanatype    verification data type [asm]

    -storedir storedir location to move outputs after processing, prior to archiving
                       [dirname(\$FVHOME) or \$NOBACKUP]
    -noarchive         do not archive outputs [archives by default]

    -landonly          calculate stats only over land surface (yes/no) 
                       [no by default]

    -nxonly            calculate stats of 2d Nx collection only (yes/no) 
                       [no by default] 

    -nodes nodesname   specify nodes (e.g., sky, hasw, or cas)
    -das               check for DAS hidden files before attempting to fetch files
                       and set no prompt; requires \$FVHOME environment variable;

    -np                no prompt; do not prompt for inputs
    -dryrun            dry run; show stats.x commands, but do not execute
    -h, -help          print usage information

EOF
exit();
}
