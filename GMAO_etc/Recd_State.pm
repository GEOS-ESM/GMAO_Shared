package Recd_State;
use Exporter;
@ISA = ('Exporter');
@EXPORT = (recd_state);

# ***************** Record job status to file **************
#
#     The included modules handle past data from job, format data
#  and call scheduler utility - status to write status record to
#  file task_status.
#
#  Inputs:
#     @pramt - array of strings. Some of them will be written in file:
#       [0] -- job name;
#       [1] -- job status (FAILED or COMPLETE);
#       [2] -- part of record in status file that are:
#               sched_cnfg, sched_id, sched_synp, sched_c_dt
#       [3] -- scheduler's top dir;
#       [4] -- schedule file local or temprary dir. It would be used when
#              job status was COMPLETE and write to status file false;
#
#  Output:
#    Call scheduler utility "status" with these variables:
#       $pramt[0] $pramt[1] $pramt[2] $pramt[3]
#
#  Example:
#
#    recd_state( "get_ssmi.pl",
#                "FAILED",
#                "FLK, GET-SSMI-01, 06:00, 2004-12-07",
#                "/home/dao_ops/Scheduler",
#                "/home/dao_ops" )
#
#  Author:  Xiaochuan Huang
#
#  History: Initial code 12/14/04
#
#  *****************************************************

sub recd_state {
 
    my ($jb, $state, $keys, $sch_dir, $local_dir) = @_;
 
    #  $keys carry 4 data of strings. It will be write in task_status
    #  file that request data should be seperated by comma and white
    #  space like "FLK, TOP-FLK-DAS, 12:00, 2004-12-10". Count comma
    #  and space for checking. If there was no space between strings,
    #  add it to.
    ##
    $cnt_space = $keys =~ tr/ //;
    $cnt_cama = $keys =~ tr/,//;
    if( $cnt_space < $cnt_cama )
    {
        print "\n!!! add space in.";
 
        #  Add white space after camma.
        ##
        $keys =~ s/\,/, /g;
    }
   
    $comd_wrt = "$sch_dir/utility/status";
    $args = "$jb $state $keys $sch_dir";
    print "\n args ==$args";
 
    $signal = system "$comd_wrt $args";
    print "\nwrite status to schedule file signal = $signal";
 
    # If it was false in writing table information to the status file,
    # write the table information to a local or temporary file.
 
    if ( $signal != 0 && $state =~ "COMP" )
    {
        my $stat_dt_tm = `date -u +%Y-%m-%d-%H:%M`;
        chomp( $stat_dt_tm );
        $tmp_stat = "$keys, $state, $stat_dt_tm, task $jb";
 
        #  if local file dir is not existe, create it.
        ##
        my $lis_dir = `ls -l $local_dir`;
        if ( $lis_dir eq "" || $! =~ "No" )
        {
            print "\n>>>> make local dir. tmp_stat =$tmp_stat";
            `mkdir $local_dir`;
            `chmod 777 $local_dir`;
        }
        $stat_fl = "$local_dir/task_status";
        `echo "$tmp_stat" >>$stat_fl`;
        `chmod -R 666 $stat_fl`;
    }
 
}
# end recd_state

