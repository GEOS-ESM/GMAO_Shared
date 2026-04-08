#!/usr/bin/env perl
#=======================================================================
# name - cmpdir
# purpose - this utility allows users to quickly and easily see
#           differences between files in two experiment directories.
#
# Note:
# 1. See usage subroutine for usage information
# 2. This script uses the diff command and xxdiff utility
#
# !Revision History
# 29Mar2010  Stassi  Initial version.
#=======================================================================
use strict;
use warnings;

use Cwd qw(abs_path cwd);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Path qw(mkpath rmtree);
use FindBin qw($Bin);
use Getopt::Long;

# global variables
#-----------------
my ($amperFLG, $bindiff, $BbwiFLG, $debug, $cdoX, $cmpX, $delim, $diffFLGs);
my ($dir1, $dir2, $dirA, $dirB, $dirL1, $dirL2, $dmgetX, $dn, $filemode);
my ($first, $follow, $ignoreFLG, $trivialFLG, $ignoreFile, $list, $listx);
my ($quiet, $recurse, $sortFLG, $subdir, $tarfile_dirA, $tarfile_dirB);
my ($tdirfile, $tmpdir, $verbose, $xxdiff_available, $xxFLG);
my (%ddelim, %different, %diffsBIN, %diffsTXT, %dir_display, %filesize);
my (%found, %identical, %ignore, %opts, %patterns, %vopts);
my (@exclude, @extEXCL, @extINC, @fileIDs, @files, @files1, @files2);
my (@p1, @p2, @unmatched1, @unmatched2);

$amperFLG = "&";
$delim = "="; # character to use when defining diffs to ignore
$tdirfile = "tarfile_directory_locationXXX"; # used when comparing tarfiles

# main program
#-------------
{
    init();
    cmp_files();
    show_results();

    if (-d $tmpdir) {
        underline("Removing temporary dir and files:");
        rmtree($tmpdir, \%vopts);
    }
}

#=======================================================================
# name - init
# purpose  - get input parameters and initialize global variables
#=======================================================================
sub init {
    my ($binX, $etcflg, $fcstflg, $help, $rsflg, $runflg, $xxdiff);
    my ($arrAddr, $modinit, $val, @values, @vlist);

    # use %vopts for verbose always on
    #---------------------------------
    %vopts = ("verbose" => 1);

    # create temporary directory
    #---------------------------
    $tmpdir = "$ENV{TMPDIR}/tmpdir.$$";
    rmtree($tmpdir, \%opts) if -d $tmpdir;

    # get runtime flags
    #------------------
    GetOptions("b"          => \$trivialFLG,
               "binX"       => \$binX,
               "d|subdir=s" => \$subdir,
               "db|debug"   => \$debug,
               "extX=s"     => \@extEXCL,
               "ext=s"      => \@extINC,
               "etc"        => \$etcflg,
               "fcst"       => \$fcstflg,
               "follow"     => \$follow,
               "h|help"     => \$help,
               "id=s"       => \@fileIDs,
               "list"       => \$list,
               "listx"      => \$listx,
               "p=s"        => \%patterns,
               "p1=s"       => \@p1,
               "p2=s"       => \@p2,
               "q"          => \$quiet,
               "r"          => \$recurse,
               "rs"         => \$rsflg,
               "run"        => \$runflg,
               "v"          => \$verbose,
               "X=s"        => \@exclude );
    usage() if $help;

    $follow = 0 unless $follow;
    $list = 1 if $listx;
    $recurse = 1 if $list;
    set_debug_vals();

    # exclude single- and double-dot directories
    #-------------------------------------------
    push @exclude, ".";
    push @exclude, "..";
    push @exclude, $tdirfile;

    # verbose option in %opts, if verbose mode
    #-----------------------------------------
    if ($verbose) { %opts = ("verbose" => 1) }
    else          { %opts = () }

    # extract comma-separated option values
    #--------------------------------------
    foreach $arrAddr (\@exclude, \@extEXCL, \@extINC, \@fileIDs, \@p1, \@p2) {
        @vlist = ();
        foreach (@$arrAddr) {
            @values = split ',', $_;
                
            foreach $val (@values) {

                # remove leading '.' from extension
                #----------------------------------
                if ($arrAddr == \@extEXCL or $arrAddr == \@extINC) {
                    while (index($val, ".") == 0) { $val = substr($val, 1) }
                }
                push @vlist, $val if $val;
            }
        }
        @$arrAddr = @vlist;
    }

    # extensions to exclude with -binX flag
    #--------------------------------------
    if ($binX) {
        foreach (qw(a d mod o pyc so x nc4 )) {
            push @extEXCL, $_;
        }
    }

    # shortcuts for checking etc, fcst, run, or rs directory 
    #-------------------------------------------------------
    unless ($subdir) {
        if ($etcflg)  { $subdir = "etc";  $recurse = 1 }
        if ($fcstflg) { $subdir = "fcst"; $recurse = 1 }
        if ($runflg)  { $subdir = "run";  $recurse = 1 }
        if ($rsflg)   { $subdir = "rs";   $recurse = 1 }
    }

    # runtime parameters
    #-------------------
    $dirA  = shift @ARGV;
    $dirB  = shift @ARGV;

    # check inputs
    #-------------
    usage() unless $dirA and $dirB;

    # initialize other variables
    #---------------------------
    $diffFLGs = "";
    $diffFLGs = "-Bbwi" if $trivialFLG;
    $xxFLG = 0;
    $first = 1;
    %ignore = ();
    %ddelim = ();

    # check for xxdiff availability
    #------------------------------
    chomp($xxdiff = `which xxdiff $dn`);
    unless (-x $xxdiff) {
        $modinit = "/usr/share/modules/init/perl";
        if (-e $modinit) {
            do $modinit;
            module("load xxdiff");
            chomp($xxdiff = `which xxdiff $dn`);
        }
    }

    if (-x $xxdiff) { $xxdiff_available = 1 }
    else {
        $xxdiff_available = 0;
        print "\nCannot find xxdiff utility.\n";
        print "All differences will be displayed with diff cmd.\n";
        pause();
    }

    # select program to use for $bindiff
    #-----------------------------------
    if (-x "/home/jstassi/bin/cdo") {
        $bindiff = "/home/jstassi/bin/cdo diffn";
        $cdoX = "/home/jstassi/bin/cdo";
    }
    else  {
        chomp($bindiff = `which h5diff`);
        $bindiff = "" unless -x $bindiff;
        $cdoX = 0;
    }

    # check for cmp program
    #----------------------
    chomp($cmpX = `which cmp`);
    $cmpX = 0 unless -x $cmpX;

    # look for dmget command
    #-----------------------
    $dmgetX = "/x/x/x/";
    chomp($dmgetX = `which dmget`);

    # default ignore list file
    #-------------------------
    $ignoreFile = "." .basename($0) .".ignore_list";
    if ($ENV{"HOME"}) { $ignoreFile = "$ENV{HOME}/$ignoreFile" }
    else              { $ignoreFile = cwd() ."/$ignoreFile" }

    # comparing two files or two directories?
    #----------------------------------------
    if (-f $dirA and -f $dirB) {
        $filemode = 1;
        init_filemode();
    }
    else {
        die ">> Error << Cannot find directory $dirA" unless -d $dirA;
        die ">> Error << Cannot find directory $dirB" unless -d $dirB;
        $filemode = 0;
        init_dirmode();
    }
}

#=======================================================================
# name - init_filemode
# purpose - initializations required when comparing two files
#=======================================================================
sub init_filemode {
    push @files1, $dirA;
    push @files2, $dirB;
    $found{$dirA} = $dirB;
}

#=======================================================================
# name - init_dirmode
# purpose - initializations required when comparing two directories
#=======================================================================
sub init_dirmode {
    my ($file1, $file2, $base1, $dirname1, $middir, $subdirname);

    # remove final slash from directory name
    #---------------------------------------
    $dirA = cleanDirName($dirA, 1);
    $dirB = cleanDirName($dirB, 1);

    # add user inputted patterns to pattern arrays
    #---------------------------------------------
    foreach (keys %patterns) {
        push @p1, $_;
        push @p2, $patterns{$_};
    }

    # add directory basenames to pattern arrays
    #------------------------------------------
    push @p1, basename $dirA;
    push @p2, basename $dirB;

    # check for @p1/@p2 correspondence
    #---------------------------------
    die ">> Error << unequal number of patterns (-p1/-p2)"
        unless scalar(@p1) == scalar(@p2);

    # check for requested subdirectory
    #---------------------------------
    if ($subdir) {
        die "not found: $dirA/$subdir\n" unless -d abs_path("$dirA/$subdir");
        die "not found: $dirB/$subdir\n" unless -d abs_path("$dirB/$subdir");
    }

    # get directories
    #----------------
    $dir1 = $dirA;
    $dir2 = $dirB;
    $dirL1 = basename $dir1;
    $dirL2 = basename $dir2;

    unless ($dirL1 eq $dirL2) {
        $ignore{$dirL1} = $dirL2;
        $ddelim{$dirL1} = "=";
    }

    if ($subdir) {
        $dir1  .= "/$subdir";
        $dir2  .= "/$subdir";
        $dirL1 .= "/$subdir";
        $dirL2 .= "/$subdir";
    }

    # get file lists
    #---------------
    @files1 = ();
    @files2 = ();

    $tarfile_dirA = get_filelist($dir1, \@files1);
    $tarfile_dirB = get_filelist($dir2, \@files2);
    die ">> Error << No files found in dir1: $dir1;" unless @files1;
    die ">> Error << No files found in dir2: $dir2;" unless @files2;
    if ($list) { cmp_lists(); return }

    show_file_counts(1);
    ignore_expid_diffs(\@files1, \@files2);

    # zero-out lists
    #---------------
    %diffsBIN = ();
    %diffsTXT = ();
    %different = ();
    %identical = ();
    @unmatched1 = ();
    @unmatched2 = ();
    %found = ();

    # which files are in both expdir1 and expdir2?
    #---------------------------------------------
    foreach $file1 (@files1) {

        $base1 = basename $file1;
        next if Xcluded($base1);

        $dirname1 = dirname $file1;
        $middir = "";
        $middir = $1 if ($dirname1 =~ m[^$dir1/(\S+)]);
        print "(1) checking $file1\n" if $verbose;

        if ($middir) { $file2 = namechange("$dir2/$middir/$base1") }
        else         { $file2 = namechange("$dir2/$base1")         }

        print "(2) checking $file2\n\n" if $verbose;

        if (-e $file2 or -l $file2) { $found{$file1} = $file2 }
        else                        { push @unmatched1, display($file1,1) }
    }
}

#=======================================================================
# name - cmp_files
# purpose - get list of files from both directories and compare files
#           to find those which are the same in both directories, those
#           which are not, and those which are in only one directory or
#           the other.
#=======================================================================
sub cmp_files {
    my ($file1, $file2, $base1, $base2);
    my ($max, $fmt1, $fmt2);
    my ($status, $indexB, $indexT, @tempArr, %Rfound);

    # send job to dmget files (just in case)
    #---------------------------------------
    dmget(1, keys %found);
    dmget(2, values %found);

    # compare files found in both directories
    #----------------------------------------
    @tempArr = (keys %found);
    $max = maxLength(\@tempArr, "basename");

    $fmt1 = "checking %s\n";
    $fmt2 = "checking %-${max}s <=> %s\n";

    foreach $file1 (sort keys %found) {
        $file2 = $found{$file1};
        $base1 = basename($file1);
        $base2 = basename($file2);
        unless ($quiet) {
            if ($base1 eq $base2) { printf $fmt1, $base1         }
            else                  { printf $fmt2, $base1, $base2 }
        }

        # get file sizes
        #---------------
        $filesize{$file1} = -s abs_path($file1);
        $filesize{$file2} = -s abs_path($file2);

        if ( text($file1, $file2) ) {
            $status = system_("diff $diffFLGs $file1 $file2 >& /dev/null");
        }
        elsif ($cdoX and ($file1 =~ /\.hdf/ or $file2 =~ /\.nc4/)) {
            $status = cdo_diff($file1, $file2);
            if ($status == 2 and $cmpX) {
                print "checking with cmp\n";
                $status = system_("$cmpX $file1 $file2 >& /dev/null");
            }
        }
        elsif ($cmpX) {
            $status = system_("$cmpX $file1 $file2 >& /dev/null");
        }
        else {
            $status = system_("diff $diffFLGs $file1 $file2 >& /dev/null");
        }
        unless ($status) {
            $identical{$file1} = $file2;
            next;
        }

        # files are different
        #--------------------
        $different{$file1} = $file2;
        if ( text($file1, $file2) ) { $diffsTXT{++$indexT} = $file1 }
        else                        { $diffsBIN{++$indexB} = $file1 }
    }

    # which $dir2 files are not in $dir1?
    #------------------------------------
    %Rfound = reverse %found;
    foreach $file2 (@files2) {
        push @unmatched2, display($file2,2) unless $Rfound{$file2};
    }

    @unmatched1 = sort @unmatched1;
    @unmatched2 = sort @unmatched2;
}

#=======================================================================
# name - cdo_diff
# purpose - use cdo utility to test equivalence of two files
#
# input parameters
# => file1: first file to compare
# => file2: second file to compare
#
# return values
# => 0: successful cdo comparison; zero diffs found
# => 1: successful cdo comparison; diffs found
# => 2: unsuccessful cdo comparison
#=======================================================================
sub cdo_diff {
    my ($file1, $file2) = @_;
    my ($line, $status);
    
    foreach $line (`$cdoX -s diffn $file1 $file2`) {
        $status = $?;
        return 2  if $status;
        return $1 if $line =~ m/(\d+)(\s+)of(\s+)(\d+) records differ/;
    }
    return 2;
}

#=======================================================================
# name - get_filelist
# purpose - get a list of files from a directory
#
# input parameters
# => $dirname: name of directory
# => $flAddr: address of array hold list of file names
#=======================================================================
sub get_filelist {
    my ($dirname, $flAddr);
    my (@dirs, $dir, @names, $name, $tardirfile, $tarfile_dir);

    $dirname = shift @_;
    $flAddr = shift @_;

    die "Error. $dirname is not a directory;" unless -d abs_path($dirname);
    die "Error. file list array address not given;" unless $flAddr;

    $tarfile_dir = "";
    $tardirfile = "$dirname/$tdirfile";
    chomp($tarfile_dir = `cat $tardirfile`) if -e $tardirfile;

    push @$flAddr, $dirname if $list;

    @dirs = ();
    @names = (<$dirname/*>, <$dirname/.*>);
    foreach $name (sort(@names)) {
        next if Xcluded($name);
        if    (dirOK($name))  { push @dirs,    $name }
        elsif (fileOK($name)) { push @$flAddr, $name }
    }
    if ($recurse) {
        foreach $dir (sort(@dirs)) { get_filelist($dir, $flAddr) }
    }
    return $tarfile_dir;
}

#=======================================================================
# name - cmp_lists
# purpose - compare lists of files in $dir1 and $dir2
#=======================================================================
sub cmp_lists {
    my ($lfile1, $lfile2);

    # dump lists to files
    #--------------------
    mkpath($tmpdir, \%vopts) unless -d $tmpdir;

    $lfile1 = "$tmpdir/cmpdir_list1_" .basename($dir1) .".$$";
    listdump(1, $dir1, $lfile1, @files1);
    print "\n(1) list written: $lfile1 (" .scalar(@files1) ." lines)\n";

    $lfile2 = "$tmpdir/cmpdir_list2_" .basename($dir2) .".$$";
    listdump(2, $dir2, $lfile2, @files2);
    print "(2) list written: $lfile2 (" .scalar(@files2) ." lines)\n\n";
        
    # set up global arrays and hashes to show differences in list files
    #------------------------------------------------------------------
    @files1 = ( $lfile1 );
    @files2 = ( $lfile2 );
    $found{$lfile1} = $lfile2;
    $different{$lfile1} = $lfile2;
    $diffsTXT{"1"} = $lfile1;
    $diffFLGs = "-Bbwi";
    return;
}

#=======================================================================
# name - ignore_expid_diffs
# purpose - check the file names to see if the files have a common expid
#           label. if so, then add expid labels to pattern arrays, if not
#           already present.
#
# input parameters
# => $arr1ADDR: address of 1st array of filenames
# => $arr1ADDR: address of 2nd array of filenames
#=======================================================================
sub ignore_expid_diffs {
    my ($arr1ADDR, $arr2ADDR, @arr1, @arr2);
    my ($base, @parts, %freq1, %freq2);
    my ($maxF, $expid1, $expid2);

    $arr1ADDR = shift @_;
    $arr2ADDR = shift @_;
    @arr1 = @$arr1ADDR;
    @arr2 = @$arr2ADDR;

    # look for common filename labels in each array
    #----------------------------------------------
    foreach (@arr1) {
        $base = basename $_;
        @parts = split /[.]/, $base;
        ++$freq1{$parts[0]} if scalar(@parts) > 1;
    }
    foreach (@arr2) {
        $base = basename $_;
        @parts = split /[.]/, $base;
        ++$freq2{$parts[0]} if scalar(@parts) > 1;
    }
    
    # assign expid to most common label
    #----------------------------------
    $maxF = 0; $expid1 = "";
    foreach (keys %freq1) {
        if ($freq1{$_} > $maxF) { $expid1 = $_; $maxF = $freq1{$_} }
    }
    #return unless $maxF*2 > scalar(@arr1);
    $expid1 = "" unless $maxF*2 > scalar(@arr1);

    $maxF = 0; $expid2 = "";
    foreach (keys %freq2) {
        if ($freq2{$_} > $maxF) { $expid2 = $_; $maxF = $freq2{$_} }
    }
    #return unless $maxF*2 > scalar(@arr2);
    $expid2 = "" unless $maxF*2 > scalar(@arr2);

    # add expids to pattern arrays if appropriate
    #--------------------------------------------
    if ($expid1 ne $expid2) {
        #unless ( in($expid1, @p1) and in($expid2, @p2) ) {
        unless (in($expid1, @p1)) {
            push @p1, $expid1;
            push @p2, $expid2;
            $ignore{$expid1} = $expid2;
            $ddelim{$expid1} = "=";
        }
    }
    if ($debug) {
        foreach (0..$#p1) { print "p1/p2 = $p1[$_]/$p2[$_]\n" };
        pause();
    }
    return;
}

#=======================================================================
# name - show_results
# purpose - give menu option for user to view results of comparison
#           between two directories
#=======================================================================
sub show_results {
    my ($opt, $dflt, $ans, $pid);

    $opt  = "";
    $dflt = 1;

    unless (@files1 or @files2) {
        print "No files were found\n";
        exit;
    }
    unless (%found) {
        print "No common files were found.\n"
            . "Continue? (y/n) [n] ";
        #====================
        chomp($ans = <STDIN>);
        #====================
        exit unless lc($ans) eq "y";
        $dflt = 3;
    }

    # display list differences and quit
    #----------------------------------
    if ($list) {
        $xxFLG = 1 if $xxdiff_available;
        #$amperFLG = "";
        show_text_diffs(%diffsTXT);
        return;
    }

    while (1) {

        # skip menu if comparing two files
        #---------------------------------
        if ($filemode) {
            if ($first) { $opt = 1; $first = 0 }
            else        { $opt = 0 }
        }

        # menu for comparing two directories
        #-----------------------------------
        else {
            underline("Make Selection",2);
            print " 1. differences\n"
                . " 2. identical\n"
                . " 3. not found\n\n"
                . " 4. file counts\n"
                . " 5. file lists\n\n"
                . " 0. quit\n\n"
                . "choose option: [$dflt] ";

            #====================
            chomp($opt = <STDIN>);
            #====================
            $opt =~ s/\s//g;
            $opt = $dflt if $opt eq "";

            # check for integer input
            #------------------------
            unless ($opt =~ m/\d+/) {
                print "\n$opt: Invalid option; Try again.\n";
                pause();
                next;
            }
        }
        return unless $opt;

        if (int($opt) > 2) { $dflt = 0 }
        else               { $dflt = $opt + 1 }

        if ($opt eq "1") { show_differences(); next }
        if ($opt eq "2") { show_identical();   next }
        if ($opt eq "3") { show_unmatched();   next }
        if ($opt eq "4") { show_file_counts(); next }
        if ($opt eq "5") { list_files();       next }

        print "\n$opt: Invalid option; Try again.\n";
        pause();
    }
}

#=======================================================================
# name - show differences
# purpose - print summary list of files which differ; give user option
#           to view differences in specific files using xxdiff utility
#=======================================================================
sub show_differences {
    unless (%diffsBIN or %diffsTXT) {
        print "\nNo differences found.\n\n";
        pause();
        return;
    }
    size_differences();
    show_binary_diffs() if %diffsBIN;
    show_text_diffs()   if %diffsTXT;
}

#=======================================================================
# name - size differences
# purpose - show size differences of common files
#=======================================================================
sub size_differences {
    my ($fmt, $file1, $file2, $label, $diffcnt);
    $fmt = "%s: %s (%d)\n";
    foreach $file1 (sort keys %filesize) {
        next unless $file2 = $found{$file1};
        if ($filesize{$file1} != $filesize{$file2}) {
            print "\n";
            print "These file sizes differ\n"
                . "-----------------------\n" unless $label;
            printf $fmt, "dir1", basename($file1), $filesize{$file1};
            printf $fmt, "dir2", basename($file2), $filesize{$file2};
            $diffcnt++;
            $label = 1;
        }
    }
    if ($diffcnt) { npause() unless $filemode }
}

#=======================================================================
# name - show_binary_diffs
# purpose - 
#=======================================================================
sub show_binary_diffs {
    my (%diffsB, @tempArr, $max, $align, $fmt1, $fmtB, $num, $show_menu, $ask);
    my ($file1, $file2, $base1, $base2, $ext1, $ext2, $dflt, $sel);

    if (@_) { %diffsB = @_ }
    else    { %diffsB = %diffsBIN }

    @tempArr = (values %diffsB);
    $max = maxLength(\@tempArr, "branch1");

    $fmt1 = "%3s. %s\n";
    $fmtB = "%3s. %-${max}s <=> %-s\n";

    $num = 1;
    $show_menu = 1;
    while (1) {

        if ($show_menu) {

            # select file to do binary diff
            #------------------------------
            $num--;
            underline("These binary files differ") if %diffsB;
            foreach (sort numeric keys %diffsB) {
                $file1 = $diffsB{$_};
                $base1 = branch($file1, "1");
                $base2 = branch($different{$file1}, "2");

                if ($base1 eq $base2) { printf $fmt1, $_, $base1 }
                else                  { printf $fmtB, $_, $base1, $base2 }
            }
        }

        unless ($bindiff) {
            pause();
            last;
        }

        print "\n";
        printf $fmt1, "0", "previous menu";
        printf $fmt1, "-1", "refresh menu\n";

        $show_menu = 0;
        $dflt = ++$num;

        while ($dflt) {
            unless ( $diffsB{$dflt} ) { $dflt = 0; last }
            last if $diffsB{$dflt} =~ /\.tar$/
                or  $diffsB{$dflt} =~ /\.png$/;

            if ($cdoX) {
                last if $diffsB{$dflt} =~ /\.hdf$/
                    or  $diffsB{$dflt} =~ /\.nc4$/
                    or  $diffsB{$dflt} =~ /\.ods$/;
            }
            else {
                last if $diffsB{$dflt} =~ /\.nc4$/;
            }
            $dflt++;
        }

        print "Compare files: [$dflt] ";
        #====================
        chomp($sel = <STDIN>);
        #====================
        $sel =~ s/\s//g;
        $sel = $dflt if $sel eq "";

        if ($sel ==  0) { return }
        if ($sel == -1) { $show_menu = 1; next }
            
        # show selected binary diff
        #--------------------------
        $num = $sel;
        unless ($diffsB{$num}) {
            print "Selection not found: $num; Try again.\n";
            $num = --$dflt;
            $sel = "";
            pause();
            next;
        }
        $file1 = $diffsB{$num};
        $file2 = $different{$file1};

        $ext1 = get_ext($file1);
        $ext2 = get_ext($file2);

        # compare tarfile contents
        #-------------------------
        if ($ext1 eq "tar" and $ext2 eq "tar") {
            $ask = query("Compare tarfile contents?", "y");
            if (yes($ask)) {
                cmp_tarfiles($file1, $file2);

                $base1 = basename($file1);
                $base2 = basename($file2);

                print "\n Tarfile comparison complete.";
                print "\n [$num] $base1";
                print " <=> $base2" if $base1 ne $base2;
                print "\n";
                npause(1)
            }
            $num++;
            $show_menu = 1;
        }

        # display Portable Network Graphic files
        #---------------------------------------
        if ($ext1 eq "png" and $ext2 eq "png") {
            system("display $file1 &");
            sleep(1);
            system("display $file2 &");
            $num++;
            $show_menu = 1;
        }

        # compare binary files
        #---------------------
        else { cmp_binary_files($num, $file1, $file2) }
    }
}

#=======================================================================
# name - cmp_binary_files
# purpose - compare two binary files
#
# input parameters
# => $num: index of %diff
# => $file1: 1st file to compare
# => $file2: 2nd file to compare
#=======================================================================
sub cmp_binary_files {
    my ($num, $file1, $file2, $base1, $base2, $status);

    $num   = shift @_;
    $file1 = shift @_;
    $file2 = shift @_;

    $base1 = branch($file1, "1");
    $base2 = branch($file2, "2");

    if ($base1 eq $base2) {
        printf "$bindiff (%d) %s\n", $num, $base1;
    } else {
        printf "$bindiff (%d) %s <=> %s\n", $num, $base1, $base2;
    }

    $status = system_("$bindiff $file1 $file2");
    unless ($bindiff =~ m/cdo/) {
        if ($status ) { print "FILES DIFFER\n" }
        else          { print "FILES MATCH\n"  }
    }
}

#=======================================================================
# name - cmp_tarfiles
# purpose - recursively call cmpdir.pl to compare contents of tarfiles
#=======================================================================
sub cmp_tarfiles {
    my  (@tarfile, @tmptardir, @tdir, $tardirfile);

    $tarfile[1] = shift @_;
    $tarfile[2] = shift @_;

    $tdir[1] = "$dirA\n";
    $tdir[2] = "$dirB\n";

    mkpath($tmpdir, \%vopts) unless -d $tmpdir;

    # untar into temporary directories
    #---------------------------------
    foreach (1..2) {
        $tmptardir[$_] = "$tmpdir/tdir$_.$$";
        rmtree($tmptardir[$_], \%opts) if -d $tmptardir[$_];
        mkpath($tmptardir[$_], \%vopts);

        system_("tar xf $tarfile[$_] -C $tmptardir[$_]")
            && die "Error untarring $tmptardir[$_];";

        # write tarfile directory location to file in tmpdir
        #---------------------------------------------------
        $tardirfile = "$tmptardir[$_]/$tdirfile";
        open tD, "> $tardirfile" or die "Error opening $tardirfile;";
        print tD $tdir[$_];
        close tD;
    }

    # use cmpdir.pl to compare temporary directories
    #-----------------------------------------------
    system("$Bin/cmpdir.pl -r $tmptardir[1] $tmptardir[2]");
}

#=======================================================================
# name - show_text_diffs
# purpose - show text differences
#=======================================================================
sub show_text_diffs {
    my (%diffsT, @tempArr, $max, $fmt0, $fmt1, $fmti, $fmtT, $num, $icnt);
    my ($file1, $firstpass, $base1, $base2, $dflt, $sel, $ncount);

    if (@_) { %diffsT = @_ }
    else    { %diffsT = %diffsTXT }

    @tempArr = (values %diffsT);
    $max = maxLength(\@tempArr, "branch1");

    $fmt0 = "%3s. %s";
    $fmt1 = "%3s. %s\n";
    $fmti = "%3s. %s (%d)\n";
    $fmtT = "%3s. %-${max}s <=> %-s\n";

    return unless %diffsT;
    $firstpass = 1;
    $sel = "";
    $num = 0;
    while (1) {

        # skip menu on first pass for list option
        #----------------------------------------
        if ($firstpass and $list) {
            $firstpass = 0;
            $num = 1;
            goto display;
        }

        # select which file to show differences
        #--------------------------------------
        underline("These text files differ");
        foreach (sort numeric keys %diffsT) {
            $file1 = $diffsT{$_};
            $base1 = branch($file1, "1");
            $base2 = branch($different{$file1}, "2");

            if ($base1 eq $base2) { printf $fmt1, $_, $base1 }
            else                  { printf $fmtT, $_, $base1, $base2 }
        }
        $dflt = ++$num;
        $dflt = 0 unless $diffsT{$dflt};
        print "\n";

        # display differences for multiple files
        #---------------------------------------
        printf $fmt1, "0", "previous menu\n";
        if (keys %diffsT > 1) {
            printf $fmt0, "a", "display all diffsT";
            print " (starting from $dflt)" if $dflt; print "\n";

            printf $fmt0, "aN", "display next N (int value) differences";
            print " (starting from $dflt)" if $dflt; print "\n";

            printf $fmt0, "-1", "reset default index to 1\n\n";
        }

        # toggle $trivialFLG, if trivial diffs not initially ignored
        #-----------------------------------------------------------
        if (! $trivialFLG) {
            if ($diffFLGs) { printf $fmt1, "b", "turn OFF -Bbwi diff flag" }
            else           { printf $fmt1, "b", "turn ON -Bbwi diff flag"  }
        }

        # edit ignore diffs
        #------------------
        if ($ignoreFLG) {
            $icnt = scalar(keys %ignore);
            if ($icnt) { printf $fmti, "e", "edit ignore diffs list", $icnt }
        }
        
        # toggle $ignoreFLG
        #------------------
        if ($ignoreFLG) { printf $fmt1, "i", "turn OFF ignore diffs"  }
        else            { printf $fmt1, "i", "turn ON ignore diffs"   }

        # toggle $sortFLG
        #----------------
        if ($sortFLG) { printf $fmt1, "s", "turn OFF sorted diff" }
        else          { printf $fmt1, "s", "turn ON sorted diff"  }

        # toggle $xxFLG
        #--------------
        if ($xxdiff_available) {
            if ($xxFLG) { printf $fmt1, "x", "turn OFF xxdiff, turn ON diff" }
            else        { printf $fmt1, "x", "turn ON xxdiff, turn OFF diff" }
        }
        print "\n";
        print "Make Selection: [$dflt] ";
        #====================
        chomp($sel = <STDIN>);
        #====================
        $sel =~ s/\s//g;
        $sel = $dflt if $sel eq "";
        $num = scalar(keys %diffsT) if $dflt == 0 and $sel !~ m/\d/;
        print "\n" unless $list;

        # return to previous menu
        #------------------------
        if ($sel eq "0") { system("clear") unless $list; return }

        # hidden option; toggle $debug
        #-----------------------------
        if ($sel eq "db") {
            $debug = ! $debug;
            set_debug_vals();

            print "debug = $debug\n";
            pause();
            $num -= 1; next;
        }

        # reset default to 1
        #-------------------
        if ($sel eq "-1") {
            $num = 0;
            next;
        }

        # show multiple diffs starting with current index
        #------------------------------------------------
        if (lc($sel) eq "a" or $sel =~ m/^a(\d*)$/) {
            $ncount = $1;
            $ncount = scalar(keys %diffsT) unless $ncount;

            $num = 0 unless $diffsT{$num};
            while ($diffsT{$num}) {
                display_text_diffs($num, %diffsT);
                last if $ncount and --$ncount <= 0;
                $num++;
            }
            next;
        }

        # toggle between xxdiff and diff comparisons
        #-------------------------------------------
        if ($sel eq "x") {
            if ($xxdiff_available) {
                $xxFLG = ! $xxFLG;
            }
            $num -= 1; next;
        }

        # toggle diff -Bbwi flag
        #-----------------------
        if ($sel eq "b") {
            $BbwiFLG = ! $BbwiFLG;
            if ($diffFLGs) { $diffFLGs = ""     }
            else           { $diffFLGs = "-Bbwi" }
            $num -= 1; next;
        }

        # edit ignore diffs list
        #-----------------------
        if ($ignoreFLG) {
            if ($sel eq "e" and %ignore) {
                diffs_to_ignore();
                $ignoreFLG = 0 unless %ignore;
                $num -= 1; next;
            }
        }
        
        # specify differences to ignore
        #------------------------------
        if ($sel eq "i") {
            $ignoreFLG = ! $ignoreFLG;
            if ($ignoreFLG and ! %ignore) {
                diffs_to_ignore();
                $ignoreFLG = 0 unless %ignore;
            }
            $num -= 1; next;
        }            

        # toggle $sortFLG
        #----------------
        if ($sel eq "s") {
            $sortFLG = ! $sortFLG;
            $num -= 1; next;
        }

        # show selected difference
        #-------------------------
        $num = $sel;
        unless ($diffsT{$num}) {
            print "Selection not found: $num; Try again.\n";
            $num = --$dflt;
            $sel = "";
            pause();
            next;
        }
      display:
        display_text_diffs($num, %diffsT);
    }
}

#=======================================================================
# name - display_text_diffs
# purpose - display the differences of two text files
#
# input parameters
# => $num: index number of difference to display (starting at 1)
#=======================================================================
sub display_text_diffs {
    my ($num, %diffsT);
    my ($file1, $file2, $base1, $base2);
    my ($srt, $f1, $f2, $cnt, $diff1, $diff2, $var, $var_);
    my ($status, $line);

    $num = shift @_;
    %diffsT = @_;

    $file1 = $diffsT{$num};
    $file2 = $different{$file1};
    $base1 = basename $file1;
    $base2 = basename $file2;

    if ($sortFLG) {
        mkpath($tmpdir, \%vopts) unless -d $tmpdir;
        $srt = "sorted ";

        $f1 = "$tmpdir/$base1.sorted.1";
        $f2 = "$tmpdir/$base2.sorted.2";
        unlink($f1) if -f $f1;
        unlink($f2) if -f $f2;

        print "Sorting files ... ";
        system_("sort $file1 > $f1") && die "Error: sort $file1 > $f1;";
        system_("sort $file2 > $f2") && die "Error: sort $file2 > $f2;";
        print "DONE\n";
    }
    else {
        $srt = "";

        $f1 = $file1;
        $f2 = $file2;
    }

    if ($ignoreFLG) {
        mkpath($tmpdir, \%vopts) unless -d $tmpdir;
        if (($f1 eq $file1) or ($f2 eq $file2)) {
            $f1 = "$tmpdir/$base1.ignore.1";
            $f2 = "$tmpdir/$base2.ignore.2";

            unlink($f1) if -f $f1;
            unlink($f2) if -f $f2;

            copy($file1, $f1) or die "Copy failed: $!;";
            copy($file2, $f2) or die "Copy failed: $!;";
        }

        $cnt = 0;
        print "Filtering out diffs to ignore ... " if %ignore;
        foreach $diff1 (sort keys %ignore) {
            $diff2 = $ignore{$diff1};
            $var = "var" .++$cnt;
            $var_ = "\\\$\{var" .$cnt ."\}";

            sed_subst($f1, $diff1, $var_);
            sed_subst($f2, $diff2, $var_);
            if ($xxFLG) {
                system_("sed -i \"${cnt}i$var = $diff1\" $f1")
                    && die "Error: sed cmd;";
                system_("sed -i \"${cnt}i$var = $diff2\" $f2")
                    && die "Error: sed cmd;";
            }
        }
        $cnt++;
        system_("sed -i \"${cnt}i==================\" $f1")
            && die "Error: sed cmd;";
        system_("sed -i \"${cnt}i==================\" $f2")
            && die "Error: sed cmd;";
        print "DONE\n";
    }

    if ($base1 eq $base2) {
        printf "showing ${srt}diffs for (%d) %s\n\n", $num, $base1;
    } else {
        printf "showing ${srt}diffs for (%d) %s <=> %s\n\n", $num, $base1, $base2;
    }

    # xxdiff text files
    #------------------
    if ($xxFLG) { system_("xxdiff --text $diffFLGs $f1 $f2 $dn $amperFLG") }
    else {
        system("clear");
        if ($base1 eq $base2) { underline("diffs for ($num) $base1", 3) }
        else       { underline("diffs for ($num) $base1 <=> $base2", 3) }

        # show ignore diffs
        #------------------
        if ($ignoreFLG) {
            $cnt = 0;
            $line = 0;
            foreach $diff1 (sort keys %ignore) {
                $diff2 = $ignore{$diff1};
                $line = "var" .++$cnt .": $diff1 <=> $diff2";
                print "$line\n";
            }
            print "-" x length($line) ."\n" if $line;
        }

        # show text diffs
        #----------------
        # it's possible that no differences found after
        # trivial diffs and ignored diffs filtered out)
        #----------------------------------------------
        $status = system_("diff --text $diffFLGs $f1 $f2");
        if ($status) {
            if ($base1 eq $base2) { underline("diffs for ($num) $base1", 3) }
            else       { underline("diffs for ($num) $base1 <=> $base2", 3) }
        }
        else { print "NO DIFFERENCES FOUND\n" unless $status }

        print "\n";
        pause();
    }
    
}

#=======================================================================
# name - diffs_to_ignore
# purpose - allow user to identify (or remove) differences to ignore
#=======================================================================
sub diffs_to_ignore {
    my ($addX, $ans, $dd, $delim1, $delim_dflt, $dflt);
    my ($diff1, $diff2, $diffs12, $diffstring, $fmt1, $fmt2);
    my ($index, $line, $max, $sel, $selA);
    my (%diff1Hash, @diff1List, %selOpts);

    $dflt = 0;
  outer: while (1) {

      @diff1List = (sort keys %ignore);
      $max = maxLength(\@diff1List);
      %selOpts = ();
      $fmt1 = "%3s. %-${max}s %s %s\n";
      $fmt2 = "%3s. %s\n";

      # print main menu for ignore diffs
      #---------------------------------
      $index = 0;

      underline("  Ignore diffs menu");
      for $diff1 (@diff1List) {
          $diff1Hash{++$index} = $diff1;
          printf $fmt1, $index, $diff1, $ddelim{$diff1}, $ignore{$diff1};
      }
      print "\n" if @diff1List;
      printf $fmt2, "a", "add new diffs to ignore" if @diff1List;
      printf $fmt2, "c", "change [$delim] delimiter";
      printf $fmt2, "r", "read ignore list from file";
      $selOpts{"c"} = 1;
      $selOpts{"r"} = 1;
      $selOpts{"db"} = 1;
      for (1..scalar(@diff1List)) { $selOpts{$_} = 1 }
      
      if (@diff1List) {
          printf $fmt2, "w", "write ignore list to file\n";
          printf $fmt2, "n", "remove n from ignore list, where 1<=n<=$index";
          printf $fmt2, "R", "remove all diffs from ignore list";
          $selOpts{"w"} = 1;
          $selOpts{"n"} = 1;
          $selOpts{"R"} = 1;
      }
      printf $fmt2, "0", "previous menu";

      # default to "add" option, unless diffs are defined;
      # addX=1 indicates ignore diff was added from main menu
      #------------------------------------------------------
      $selA = "";
      $addX = 0;
      if (@diff1List) {
          print "\n  Make Selection: [$dflt] ";
          #====================
          chomp($sel = <STDIN>);
          #====================
          $sel =~ s/\s//g;              # <==== $sel made here
          $sel = $dflt if $sel eq "";
          last unless $sel;

          # allow user to select main menu options
          #---------------------------------------
          $diffs12 = $sel;
          ($diff1, $diff2) = split /[$delim]/, $diffs12;
          if ($diff1 and $diff2) {
              $selA = $sel;
              $sel = "a";
              $addX = 1;
          }
      }
      else { print "\n"; $sel = "a" }

      # add diffs to ignore
      #--------------------
      if ($sel eq "a") {
          while (1) {
              unless ($addX) {
                  print "  Enter diff1${delim}diff2 [quit] ";
                  #=====================
                  chomp($selA = <STDIN>);
                  #=====================
                  $selA =~ s/\s//g;  # remove blank spaces
                  $sel = $selA;
              }
              if ($selA eq "") {
                  if (%ignore) { last }
                  else         { last outer }
              }

              # allow user to make menu choice from previous menu
              #--------------------------------------------------
              elsif ($selA eq "0")    { last outer }
              elsif ($selOpts{$selA}) { last }

              # extract ignore diff information
              #--------------------------------
              $diffs12 = $selA;
              ($diff1, $diff2) = ();
              ($diff1, $diff2) = split /[$delim]/, $diffs12;

              if ($diff1 and $diff2) {
                  $ignore{$diff1} = $diff2;
                  $ddelim{$diff1} = $delim;
              } else {
                  print "\n  Cannot decipher input: $selA; Try again.\n";
                  pause(2);
                  print "\n" unless $addX;
              }
              next outer if $addX;
          }
      }

      # write ignore list to file
      #--------------------------
      if ($sel eq "w") {
          if (%ignore) {
              print "  Write ignore list to file [$ignoreFile]: ";
              #====================
              chomp($ans = <STDIN>);
              #====================
              $ans =~ s/\s//g;  # remove blank spaces
              $ignoreFile = $ans if $ans;

              print "\n";
              open IGN, "> $ignoreFile"
                  or die "Error opening ignore file for write: $ignoreFile;";
              foreach $diff1 (sort keys %ignore) {
                  $dd = $ddelim{$diff1};
                  $diff2 = $ignore{$diff1};
                  $diffs12 = "$diff1$dd$diff2";
                  print "  WRITING: $dd,$diffs12\n";
                  print IGN "$dd,$diffs12\n";
              }
              close IGN;
              pause(2);
          }
          else { print "  The ignore list is empty. Nothing to write.\n" }
      }

      # delete all diffs from ignore list
      #----------------------------------
      elsif ($sel eq "R") {
          print "  Remove all diffs from ignore list. Are you sure (y/n)? [y] ";
          #====================
          chomp($ans = <STDIN>);
          #====================
          unless (lc($ans) eq "n") {
              %ignore = ();
              %ddelim = ();
              print "  All diffs removed!\n";
          }
          else { print "  No diffs removed!\n"; pause(2) }
      }

      # use different character for defining diffs to ignore
      #-----------------------------------------------------
      # this is necessary if "=" is in the
      # string to ignore or its substitute
      #-----------------------------------
      elsif ($sel eq "c") {
          if ($delim eq "=") { $delim_dflt = ":" }
          else               { $delim_dflt = "=" }
          print "  Character to use for ignore rules [$delim_dflt] ";
          #======================
          chomp($delim = <STDIN>);
          #======================
          $delim = $delim_dflt if $delim eq "";
      }

      # read ignore list from file
      #---------------------------
      elsif ($sel eq "r") {
          while (1) {
              print "  Read ignore list from file [$ignoreFile]: ";
              #====================
              chomp($ans = <STDIN>);
              #====================
              $ans =~ s/\s//g;  # remove blank spaces

              $ignoreFile = $ans if $ans;
              last if -e $ignoreFile;
              print "  File does not exist: $ignoreFile\n";
          }
          open IGN, "< $ignoreFile"
              or die "Error opening ignore file for read: $ignoreFile;";
          print "\n";
          foreach $line (<IGN>) {
              chomp($line);
              if ($line =~ m/^(.),(.*)\1(.*)$/) {
                  $delim1 = $1;
                  $diff1 = $2;
                  $diff2 = $3;
                  $ignore{$diff1} = $diff2;
                  $ddelim{$diff1} = $delim1;
                  print "  READING: $diff1 $ddelim{$diff1} $diff2\n";
              } else {
                  print "  ???????: $line\n"
              }
          }
          close IGN;
      }

      # hidden option; toggle $debug
      #-----------------------------
      elsif ($sel eq "db") {
          $debug = ! $debug;
          set_debug_vals();

          print "\n  debug = $debug\n";
          pause(2);
      }

      else {
          $diff1 = $diff1Hash{$sel} ? $diff1Hash{$sel} : "";
          $diff2 = $ignore{$diff1};

          # delete individual diffs from ignore list
          #-----------------------------------------
          if ($diff2) {
              $diffstring = "$diff1 $ddelim{$diff1} $diff2";
              print "  Remove diff: $diffstring. Are you sure (y/n)? [y] ";
              #====================
              chomp($ans = <STDIN>);
              #====================
              unless (lc($ans) eq "n") {
                  delete($ignore{$diff1});
                  delete($ddelim{$diff1});
                  print "  Diff removed!\n";

                  # update %diff1Hash
                  #------------------
                  while ($diff1Hash{++$sel}) {
                      $diff1Hash{$sel-1} = $diff1Hash{$sel};
                  }
                  delete($diff1Hash{--$sel});
              }
              else { print "  No diff removed!\n"; pause(2) }
          }

          # ??????
          #-------
          else {
              if ($sel) {
                  print "  Cannot decipher input: $sel; Try again.\n";
                  pause(2);
              }
          }
      }
  }
}

#=======================================================================
# name - sed_subst
# purpose - substitute $str2 for $str1 in $file
#
# input parameters
# => $file
# => $str1
# => $str2
#=======================================================================
sub sed_subst {
    my ($file, $str1, $str2);
    my ($D, @delimCharList, $char);

    $file = shift @_;
    $str1 = shift @_;
    $str2 = shift @_;

    # find delim character that is not in either string
    #--------------------------------------------------
    $D = "";
    @delimCharList =("/", "|", "^", "#");

    foreach $char (@delimCharList) {
        next if $str1 =~ m/$char/;
        next if $str2 =~ m/$char/;
        $D = $char;
        last;
    }
    return unless $D;

    # system sed substitute command
    #------------------------------
    system_("sed -i \"s${D}$str1${D}$str2${D}g\" $file")
        && die "Error: sed subst $str1 to $str2 in $file;";
}

#=======================================================================
# name - show_identical
# purpose - print summary list of files which are identical in both
#           directories
#=======================================================================
sub show_identical {
    my (@tempArr, $max, $num, $fmt1, $fmt2);
    my ($file1, $file2, $base1, $base2);

    @tempArr = (keys %identical);
    $max = maxLength(\@tempArr, "branch1");
    $fmt1 = "%2d. %s\n";
    $fmt2 = "%2d. %-${max}s <=> %-s\n";

    if (%identical) {
        $num = 0;
        underline("These files are identical in the two directories");
        foreach (sort keys %identical) {
            $file1 = $_;
            $file2 = $identical{$file1};
            $base1 = branch($file1, "1");
            $base2 = branch($file2, "2");
            if ($base1 eq $base2) { printf $fmt1, ++$num, $base1         }
            else                  { printf $fmt2, ++$num, $base1, $base2 }
        }
    } else {
        print "\nNo identical files were found in the two directories.\n";
    }
    npause();
}

#=======================================================================
# name - show_unmatched
# purpose - print summary lists of files which exist in one directory
#           but not in the other
#=======================================================================
sub show_unmatched {
    my ($num, $ddir1, $ddir2, $ddirL1, $ddirL2);
    my ($ask, $f1, $f2, $file1, $file2);
    my (%diffsB, %diffsT, $indexB, $indexT);
    my (%RdiffsT, %RdiffsB, $key);
    my ($menu, $len, $len3, $len7);
    my ($dflt1, $dflt2, $cnt1, $cnt2, $cnt, $delta);

    $ddir1  = display($dir1,1);
    $ddir2  = display($dir2,2);
    $ddirL1 = display($dirL1,1);
    $ddirL2 = display($dirL2,2);

    # list files found in dir1 but not in dir2
    #-----------------------------------------
    if (@unmatched1) {
        $num = 0;
        underline("FOUND in (1) $ddirL1 but NOT FOUND in (2) $ddirL2");
        foreach ( @unmatched1 ) { printf "%2d. %s\n", ++$num, $_ }
    } else {
        if (@files1) { print "\nAll files in $ddir1 are also in $ddir2\n"   }
        else         { print "\nNo files found in dir1: $ddir1\n" }
    }
    npause();

    # list files found in dir2 but not in dir1
    #-----------------------------------------
    if (@unmatched2) {
        $num = 0;
        underline("FOUND in (2) $ddirL2 but NOT FOUND in (1) $ddirL1");
        foreach ( @unmatched2 ) { printf "%2d. %s\n", ++$num, $_ }
    } else {
        if (@files2) { print "\nAll files in $ddir2 are also in $ddir1\n"   }
        else         { print "\nNo files found in dir2: $ddir1\n" }
    }
    npause() unless @unmatched1 and @unmatched2;

    # does user want to compare unmatched files?
    #-------------------------------------------
    if (@unmatched1 and @unmatched2) {
        $menu = 1;
        $cnt1 = scalar @unmatched1;
        $cnt2 = scalar @unmatched2;
        if ($cnt1 > $cnt2) { $cnt = $cnt1 }
        else               { $cnt = $cnt2 }

        $len = maxLength(\@unmatched1, "basename");
        $len = 15 if $len < 15;
        $len3 = $len + 3;
        $len7 = $len + 7;
        $dflt1 = 1;
        $delta = 0;

        %diffsB = (); $indexB = 0;
        %diffsT = (); $indexT = 0;

        $ask = query("\nCompare unmatched files from the two directories?", "n");

        if (yes($ask) and $#unmatched1 == 0 and $#unmatched2 == 0) {
            $file1 = $unmatched1[0];
            $file2 = $unmatched2[0];
            $different{$file1} = $file2;

            if (text($file1, $file2)) { $diffsT{++$indexT} = $file1 }
            else                      { $diffsB{++$indexB} = $file1 }
            $ask = "n";
        }
        while (yes($ask)) {
            dmget(1, @unmatched1);
            dmget(2, @unmatched2);

            if ($menu) {
                print "\n     dir1: $ddirL1\n";
                print   "     dir2: $ddirL2\n\n";

                printf  "     %-${len7}s %s\n", "dir1", "dir2";
                printf  "     %-${len7}s %s\n", "-"x3, "-"x3;
                foreach (1..$cnt) {
                    if ($_ > $cnt1) {
                        printf " "x$len7 ."  %2d. %s\n",
                        $_, basename($unmatched2[$_-1]);
                    }
                    elsif ($_ > $cnt2) {
                        printf " %2d. %-${len3}s\n",
                        $_, basename($unmatched1[$_-1]);
                    }
                    else {
                        printf " %2d. %-${len3}s %2d. %s\n",
                        $_, basename($unmatched1[$_-1]),
                        $_, basename($unmatched2[$_-1]);
                    }
                }
            }
            $menu = 0;

            print "\n -1. Reprint menu\n";
            print   "  0. Exit menu\n";

            if ($f1) {
                print "\nchoose file from dir1 (1-$cnt1): $f1\n";
            }
            else {
                {
                    $f1 = query("\nchoose file from dir1 (1-$cnt1):", $dflt1++);
                    if (($f1 !~ /^-?\d+$/) or ($f1 < -1) or ($f1 > $cnt1)) {
                        $dflt1--;
                        redo;
                    }
                }
                if ($f1 == -1) { $menu = 1; $f1 = ""; $dflt1--; redo }
                if ($f1 ==  0) { last }
            }
            $dflt2 = $f1 - $delta;
            $dflt2 = $cnt2 if $dflt2 > $cnt2;

            {
                $f2 = query("choose file from dir2 (1-$cnt2):", $dflt2);
                redo if ($f2 !~ /^-?\d+$/) or ($f2 < -1) or ($f2 > $cnt2);
            }
            if ($f2 == -1) { $menu = 1; $f1 = ""; $dflt1--; redo }
            if ($f2 ==  0) { last }
            $delta = $f1 - $f2;

            $file1 = $unmatched1[$f1-1];
            $file2 = $unmatched2[$f2-1];

            # cannot compare $file1 to more than one $file2;
            # remove duplicate $file1 entries from %diffs{T,B}
            #-------------------------------------------------
            %RdiffsT = reverse %diffsT;
            if ($key = $RdiffsT{$file1}) {
                delete $diffsT{$key};
                %diffsT = reindex(\%diffsT);
                $indexT--;
            }
            %RdiffsB = reverse %diffsB;
            if ($key = $RdiffsB{$file1}) {
                delete $diffsB{$key};
                %diffsB = reindex(\%diffsB);
                $indexB--;
            }

            # add new entries
            #----------------
            $different{$file1} = $file2;
            if (text($file1, $file2)) { $diffsT{++$indexT} = $file1 }
            else                      { $diffsB{++$indexB} = $file1 }

            $dflt1 = $f1 + 1;
            $dflt1 = 0 if $dflt1 > $cnt1;

            $f1 = "";
            $f2 = "";
        }
        show_binary_diffs(%diffsB) if %diffsB;
        show_text_diffs(%diffsT) if %diffsT;
    }
}

#=======================================================================
# name - show_file_counts
# purpose - show number of files in each of the two directories being compared
#
# input parameter
# => $flag: flag indicating whether to pause afterwards (=0 for pause)
#=======================================================================
sub show_file_counts {
    my ($flag, $len1, $len2, $max, $fmt);

    $flag = shift @_;

    $len1 = length($dir1);
    $len2 = length($dir2);

    if ($len1 > $len2) { $max = $len1 }
    else               { $max = $len2 }
    $fmt = "%s: %-${max}s (%d files)\n";

    underline("Directory file counts");
    if ($tarfile_dirA and $tarfile_dirB) {
        print "dirA: $tarfile_dirA\n";
        print "dirB: $tarfile_dirB\n\n";
    }
    printf $fmt, "dir1", $dir1, scalar(@files1);
    printf $fmt, "dir2", $dir2, scalar(@files2);
    print "\n";

    pause() unless $flag;
}

#=======================================================================
# name - list_files
# purpose - display list of files in the two directories being compared
#=======================================================================
sub list_files {
    my ($fmt, $cnt, $num, $base);

    $fmt = "%2d. %s\n";

    # print filenames in dir1
    #------------------------
    if (@files1) {
        underline("dir1: " .$dir1);
        print "$tarfile_dirA\n" if $tarfile_dirA;
        $cnt = 0;
        $num = 0;
        foreach (sort @files1) {
            printf $fmt, ++$num, branch($_, "1");
            if (++$cnt == 50) {
                npause();
                $cnt = 0;
            }
        }
    } else {
        print "\nNo files in dir1: $dir1\n";
    }
    npause();

    # print filenames in dir2
    #------------------------
    if (@files2) {
        underline("dir2: " .$dir2);
        print "$tarfile_dirB\n" if $tarfile_dirB;
        $cnt = 0;
        $num = 0;
        foreach (sort @files2) {
            printf $fmt, ++$num, branch($_, "2");
            if (++$cnt == 50) {
                npause();
                $cnt = 0;
            }
        }
    } else {
        print "\nNo files in dir2: $dir2\n";
    }
    npause();
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         UTILITY subroutines
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#=======================================================================
# name - cleanDirName
# purpose - remove the final slash from a directory name and find
#           the absolute path
#
# input parameters
# => $dir:  name of directory
# => $flag: return abbreviated name if ==1
#=======================================================================
sub cleanDirName {
    my ($dir, $flag, $abspath);

    $dir  = shift @_;
    $flag = shift @_;
    while ($dir =~ m[^(.*[^/])/+$]) { $dir = $1 };
    
    $abspath = abs_path($dir);
    $dir_display{$abspath} = $dir if $flag;
    return $abspath;
}

#=======================================================================
# name - dirOK
# purpose - return true (1) if dir is okay to include in our list
#
# input parameter
# => $name: name of dir to test for inclusion
#
# return value
# => 0 to not include dir in list
# => 1 to include dir in list
#=======================================================================
sub dirOK {
    my ($name);
    $name = shift @_;

    return 0 if -l $name and (!$follow);
    return 1 if -d abs_path($name);
    return 0
}

#=======================================================================
# name - display
# purpose - display path as given rather as absolute path
#
# input parameters
# => $name: filename to display
# => $flag: flag indicating whether files come from $dir1 or $dir2
#=======================================================================
sub display {
    my $name = shift @_;
    my $flag = shift @_;

    return $name if $filemode;
    if    ($flag == 1) { $name =~ s|$dirA|$dir_display{$dirA}| }
    elsif ($flag == 2) { $name =~ s|$dirB|$dir_display{$dirB}| }

    return $name;
}

#=======================================================================
# name - dmget
# purpose - send job to dmget files prior to comparing them
#
# input parameters
# => $flag: flag indicating whether files come from $dir1 or $dir2
# => @arr: array of files to dmget
#=======================================================================
sub dmget {
    my ($flag, @arr, $str, $name, $cnt, $max, $pid, $cmd);

    return unless -x $dmgetX;

    $flag = shift @_;
    @arr  = @_;
    
    $cnt = 0;
    $max = 150;

    # get list of files in archive directory
    #---------------------------------------
    $str = "";
    while (@arr) {
        $name = shift @arr;
        if (abs_path(dirname($name)) =~ /archive/) {
            $str .= " " .display($name, $flag);
            $cnt++;
        }

        # fork job to dmget files
        #------------------------
        if ($cnt > $max or (! @arr and $cnt)) {
            $cmd = "$dmgetX $str >& /dev/null";
            print "$cmd\n" if $verbose;
            defined($pid = fork) or die ">> Error << while forking: $!";
            unless ($pid) {
                exec $cmd;
                die ">> Error << $dmgetX command not executed: $!";
            }
            $str = "";
            $cnt = 0;
        }
    }
}

#=======================================================================
# name - fileOK
# purpose - return true (1) if file is okay to include in our list
#
# input parameter
# => $name: name of file to test for inclusion
#
# return value
# => 0 to not include file in list
# => 1 to include file in list
#=======================================================================
sub fileOK {
    my ($name, $okay);
    my ($ok_ext, $my_ext, $id_ok);

    $name = shift @_;

    return if -l $name and (!$follow);
    return unless -f abs_path($name);

    $okay = 1;

    # check file extension
    #---------------------
    if (@extINC) {
        foreach $ok_ext (@extINC) {
            $my_ext = get_ext($name);
            return 1 if $my_ext eq $ok_ext;
        }
        $okay = 0;
    }

    # check file ID
    #--------------
    if (@fileIDs) {
        foreach $id_ok (@fileIDs) {
            return 1 if basename($name) =~ m/$id_ok/;
        }
        $okay = 0;
    }
    return $okay;
}

#=======================================================================
# name - get_expid
# purpose - extract a guess at the expid from list of file names
#
# input parameter
# => @arr: list of file names
#=======================================================================
sub get_expid {
    my (@arr, $expid, $max);
    my ($fullname, $name, @dummy, %count);

    @arr = @_;
    $max = 1;
    $expid = "";

    foreach $fullname (@arr) {
        ($name, @dummy) = split /[.]/, basename($fullname);
        $count{$name}++;
    }
    foreach $name (sort keys %count) {
        if ($count{$name} == $max) {
            $expid = "";
            next;
        }
        if ($count{$name} > $max) {
            $expid = $name;
            $max = $count{$name};
        }
    }
    return $expid;
}

#=======================================================================
# name - get_ext
# purpose - return file name extension
#
# input parameter
# => $name: name of file, potentially including directory path
#=======================================================================
sub get_ext {
    my ($name, $ext);
    $name = shift @_;
    $ext = (split /[.]/, basename($name))[-1];
    $ext = "" if $ext eq basename($name);
    return $ext;
}

#=======================================================================
# name - in
# purpose - determine whether value is in an array
#
# input parameters
# => $value: value to search for in an array
# => @array: array to search
#=======================================================================
sub in {
    my ($value, @array, $flag);

    $value = shift @_;
    @array = @_;

    $flag = 0;
    foreach (@array) {
        if ($_ eq $value) { $flag = 1; last }
    }
    return $flag;
}

#=======================================================================
# name - listdump
# purpose - dump sorted contents of array to file
#
# input parameters
# => $num: number 1 or 2, indicating whether it is first or second list
# => $dir: name of top level directory for list
# => $lfile: name of file where the list will be written
# => $arr: array file names to dump
#=======================================================================
sub listdump {
    my ($num, $dir, $lfile, @arr);
    my ($ddir, $expid);

    $num = shift @_;
    $dir = shift @_;
    $lfile = shift @_;
    @arr = @_;

    $ddir = display($dir, $num);
    $expid = get_expid(@arr) if $listx;

    open LFILE, "> $lfile" or die "Error opening file, $lfile;";
    print LFILE "top directory: $ddir\n";
    print LFILE "EXPID: $expid\n" if $expid;

    foreach (sort @arr) {
        s/$dir//;
        s/$expid/\$EXPID/ if $expid;
        print LFILE "$_\n";
    }
    close LFILE;
}

#=======================================================================
# name - maxLength
# purpose - find maximum length of strings in an array
#
# input parameters
# => $arrAddr: address of array containing strings
# => $flag: (optional)
#    - if eq "basename" then find maximum length of string basenames
#    - if eq "branch1" then find max length of strings after removing $dir1
#    - if eq "branch2" then find max length of strings after removing $dir2
#=======================================================================
sub maxLength {
    my ($arrAddr, $flag, @arr, $string);
    my ($length, $max);

    $arrAddr = shift @_;
    $flag = shift @_;
    @arr = @$arrAddr;

    $max = 0;
    $flag = "" unless $flag;
    foreach (@arr) {
        if    ($flag eq "basename") { $string = basename($_) }
        elsif ($flag eq "branch1")  { $string = branch($_, 1) }
        elsif ($flag eq "branch2")  { $string = branch($_, 2) }
        else                        { $string = $_ }
        $length = length($string);
        $max = $length if $length > $max;
    }
    return $max;
}

#=======================================================================
# name - branch
# purpose - return filename with the top directory name removed
#
# note: This is not the same as taking the basename, since the basename
#       function will remove all the directories preceding the last name,
#       where this function will only remove the top directory, i.e. $dir1
#       or $dir2.
#
# input parameters
# => $name: full name of file, including path
# => $flag: flag indicating which directory path to remove from $name
#             if $flag eq "1", then remove $dir1
#             if $flag eq "2", then remove $dir2
#=======================================================================
sub branch {
    my ($name, $flag, $ddir, $base);

    $name = shift @_;
    $flag = shift @_;

    return $name if $filemode;

    if ($flag eq "1") {
        $ddir = display($dir1,1);
        $ddir = "\'.\'" if $ddir eq ".";
        $name =~ s/$dir1\///;
        $name =~ s/$ddir\///;
    }
    elsif ($flag eq "2") {
        $ddir = display($dir2,2);
        $ddir = "\'.\'" if $ddir eq ".";
        $name =~ s/$dir2\///;
        $name =~ s/$ddir\///;
    }
    $name =~ s/^$ENV{"TMPDIR"}\///;  # when comparing tarfiles
    return $name;
}

#=======================================================================
# name - namechange
# purpose - substitute patterns into a filename;
#           return the new filename, if the file exists
#
# input parameters
# => $name: name of file before name change
#=======================================================================
sub namechange {
    my ($name, $namesave, $flag);
    my ($dir, $base);

    $name = shift @_;
    $namesave = $name;

    # check each pattern substitution individually
    #---------------------------------------------
    foreach (0..$#p1) {
        last if -e $name;
        $dir = dirname $namesave;
        $base = basename $namesave;
        if ($p1[$_] eq "") {
            $base = "$p2[$_].$base";
        } elsif ($p2[$_] eq "") {
            $base =~ s/\b$p1[$_].//;
        } else {
            $base =~ s/\b$p1[$_]\b/$p2[$_]/g;
        }
        $name = "$dir/$base";
    }

    # try individual substitution in full pathname
    #---------------------------------------------
    unless (-e $name) {
        foreach (0..$#p1) {
            last if -e $name;
            $name = $namesave;
            $name =~ s/\b$p1[$_]\b/$p2[$_]/g;
        }
    }

    # check cumulative pattern substitution
    #--------------------------------------
    unless (-e $name) {
        $name = $namesave;
        foreach (0..$#p1) {
            last if -e $name;
            $dir = dirname $name;
            $base = basename $name;
            $base =~ s/\b$p1[$_]\b/$p2[$_]/g;
            $name = "$dir/$base";
        }
    }

    # try cumulative substitution in full pathname
    #---------------------------------------------
    unless (-e $name) {
        $name = $namesave;
        foreach (0..$#p1) {
            last if -e $name;
            $name =~ s/\b$p1[$_]\b/$p2[$_]/g;
        }
    }

    return $name;
}

#=======================================================================
# name - numeric
# purpose - used with perl sort command to do a numeric sort
#=======================================================================
sub numeric {
    return  1 if $a > $b;
    return -1 if $a < $b;
}

#=======================================================================
# name - pause
# purpose - pause processing until user input is detected
#
# input parameter
# => $s: number of spaces to put in front of print texte
#=======================================================================
sub pause {
    my ($s, $dummy);
    $s = shift @_;
    print " "x$s if $s;
    print "Hit <CR> to continue ... ";
    $dummy = <STDIN>;
}

#=======================================================================
# name - npause
# purpose - same as sub pause() except add "\n" to front of prompt
#
# input parameter
# => $s: number of spaces to put in front of print texte
#=======================================================================
sub npause {
    my ($s, $dummy);
    $s = shift @_;
    print "\n";
    print " "x$s if $s;
    print "Hit <CR> to continue ... ";
    $dummy = <STDIN>;
}

#=======================================================================
# name - reindex
# purpose - for a hash which uses an index number for the key, renumber
#           the keys (indices) if one has been removed
#
# input parameters
# => $hashAddr: the address of the hash to reindex
#
# return value
# => %new: reindexed hash
#=======================================================================
sub reindex {
    my ($hashAddr, %myHash, $index, %new);

    $hashAddr = shift @_;
    %myHash = %$hashAddr;

    %new = ();
    $index = 1;
    foreach (sort keys %myHash) { $new{$index++} = $myHash{$_}}

    return %new;
}

#=======================================================================
# name - system_
# purpose - call system command; print command to STDOUT if $debug
#
# input parameters
# => $cmd: command to send to system
# => $vflg: verbose flag; print $cmd to STDOUT if $vflg == 1
#=======================================================================
sub system_ {
    my ($cmd, $vflg, $status);
    $cmd = shift @_;
    $vflg = shift @_;
    print "$cmd\n" if $verbose or $vflg;
    ($status = system $cmd) /= 256;
    print "status = $status\n" if $debug;
    return $status;
}

#=======================================================================
# name - text
# purpose - determine whether both files are text (i.e. viewable)
#
# input parameters
# => $file1: 1st file
# => $file2: 2nd file
#=======================================================================
sub text {
    my ($file1, $file2, $type1, $type2, $txtFLG);
    $file1 = shift @_;
    $file2 = shift @_;

    $txtFLG = 0;
    $type1 = `file -L $file1`;
    $type2 = `file -L $file2`;

    $txtFLG = 1 if ($type1=~/ASCII/ or $type1=~/text/ or $type1=~/source/)
        and         ($type2=~/ASCII/ or $type2=~/text/ or $type2=~/source/);
    unless ($txtFLG) {
        $txtFLG = 1 if get_ext($file1) eq "txt"
            and         get_ext($file2) eq "txt";
    }
    return $txtFLG;
}

#=======================================================================
# name - set_debug_vals
# purpose - set debug values (or not)
#=======================================================================
sub set_debug_vals {
    if ($debug) { $verbose = 1; $dn = "" }
    else        { $verbose = 0; $dn = "2> /dev/null" }
}

#=======================================================================
# name - underline
# purpose - prints a string to stdout and underlines it
#
# input parameters
# => string: the string to underline
# => flag: (optional); defaults to =1
#           =1: underline only with '-'
#           =2: underline and overline with '='
#           =3: underline and overline with "-"
#=======================================================================
sub underline {
    my ($string, $flag, %pattern);
    my ($s1, $str, $len_s1, $len_str);

    $string = shift @_;
    $flag = shift @_;

    $pattern{1} = "-";
    $pattern{2} = "=";
    $pattern{3} = "-";

    $flag = 1 unless $flag;
    $flag = 1 unless $flag == 2 or $flag == 3;

    ($s1, $str) = ($string =~ m/^(\s*)(.+)(\s*)$/);
    $len_s1 = length($s1);
    $len_str = length($str);
    print "\n";
    print " "x$len_s1 .$pattern{$flag}x$len_str."\n" if $flag > 1;
    print " "x$len_s1 .$str."\n";
    print " "x$len_s1 .$pattern{$flag}x$len_str."\n";
}

#=======================================================================
# name - Xcluded
# purpose - return true (1) if file is on exclusion list, @exclude
#
# input parameter
# => $name: name of file to check for exclusion
#
# return value
# => 1 : do not include $name
#=======================================================================
sub Xcluded {
    my ($name, $ext, $xname);
    $name = shift @_;

    foreach $ext (@extEXCL) {
        return 1 if (split(/\./, $name))[-1] eq $ext;
    }
    foreach $xname (@exclude) {
        return 1 if basename($name) eq $xname;
    }
    return 0;
}

#=======================================================================
# name - query
# purpose - query user for a response and return the response
#
# input parameters
# => $str: use this line to prompt for a response
# => $dflt: (optional) default value to use for <cr> response
# 
#=======================================================================
sub query {
    my ($str, $dflt, $prompt, $ans);

    $str  = shift @_;
    $dflt = shift @_;

    # prepare prompt
    #---------------
    $prompt  = "$str ";
    $prompt .= "[$dflt] " unless blank($dflt);

    # get user response
    #------------------
    print $prompt;
    #====================
    chomp($ans = <STDIN>);
    #====================
    $ans =~ s/^\s+|\s+$//g;     # remove leading/trailing blanks from response
    $ans = expand_EnvVars($ans);
    if ( blank($ans) ) { $ans = $dflt unless blank($dflt) }

    return $ans;
}

#=======================================================================
# name - expand_EnvVars
# purpose - expand Environment Variables within a string
#=======================================================================
sub expand_EnvVars {
    my ($string, $cnt, $var, $name);

    $string = shift @_;
    $cnt = 0;

    # look for ${var} format
    #-----------------------
    while ($string =~ m/(\$\{(\w+)})/)   {
        $var = $1; $name = $2;
        $var =~ s/\$/\\\$/;
        $string =~ s/$var/$ENV{$name}/;
        die "Error. Infinite loop condition;" if ++$cnt > 100;
    }

    # look for $var format
    #---------------------
    while ($string =~ m/(\$\b(\w+)\b)/) {
        $var = $1; $name = $2;
        $var =~ s/\$/\\\$/;
        $string =~ s/$var/$ENV{$name}/;
        die "Error. Infinite loop condition;" if ++$cnt > 200;
    }
    return $string;
}

#=======================================================================
sub yes {
    my $str;
    $str = shift @_;
    $str = lc $str;           # make lowercase
    $str =~ s/^\s*|\s*$//g;   # remove leading/trailing blanks
    return 1 if $str eq "y" or $str eq "yes"
}

#=======================================================================
sub neg {
    my $str;
    $str = shift @_;
    $str = lc $str;           # make lowercase
    $str =~ s/^\s*|\s*$//g;   # remove leading/trailing blanks
    return 1 if $str eq "n" or $str eq "no"
}

#=======================================================================
sub blank {
    my $str;
    $str = shift @_;
    return 1 if $str =~ /^\s*$/;
}

#=======================================================================
# name - usage
# purpose - print script usage information
#=======================================================================
sub usage {    
    my $script = basename $0;
    print << "EOF";
Purpose: Compare two files or compare the contents of two directories

Compare two files
=================
Usage: $script file1 file2

where
  file1 = first file being compared
  file2 = second file being compared


Compare two directories
=======================
Usage: $script dir1 dir2 [options]

where
  dir1 = first directory being compared
  dir2 = second directory being compared

options
  -b                 ignore blanks and other trivial differences on
                     initial text file comparisons (using -Bbwi flags)
  -binX              exclude extensions: a, d, mod, nc4, o, pyc, so, x
                     i.e. same as -extX a,d,mod,nc4,o,pyc,so,x
  -d name            start comparison in specified subdirectory
  -extX extension    exclude all files with this extension (see Note 1)
  -ext extension     compare all files with this extension (see Note 1)
  -etc               shortcut for "-subdir etc -r"
  -fcst              shortcut for "-subdir fcst -r"
  -follow            follow symbollic links; default is to not follow
  -h(elp)            print usage information
  -id fileID         compare all files with \"fileID\" as part of its filename
                     (see Note 1)
  -list              compare list of file names in dir1 and dir2 (rather than
                     comparing the actual files)
  -listx             same as -list, except ignore expid name differences
  -q                 quiet mode
  -r                 recursively compare any subdirectories found
  -rs                shortcut for "-d rs -r"
  -run               shortcut for "-d run -r"
  -v                 verbose mode
  -X filename        exclude filename from list to compare (see Note 1)

pattern options
  -p1 pattern1       ignore these pattern differences in dir1/dir2 filenames
  -p2 pattern2       (see Notes 1-3)
or
  -p pattern1=pattern2

Notes:
1. Multiple values can be given for extension (-ext, -extX), fileID (-id),
   pattern1 (-p1), pattern2 (-p2), and str (-X) by separating values with
   a comma (no space) or by multiple use of the option flag.
2. Multiple pattern1=pattern2 values can be given by multiple uses of the -p flag
3. There must be a matching pattern2 for each pattern1, and vice versa.

EOF
exit;
}
