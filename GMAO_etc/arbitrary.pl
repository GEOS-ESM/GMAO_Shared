#!/usr/bin/env perl
# Lay out task in given node
# e.g., arbitray.pl 4,1 
# lays out 4 tasks in first node and 1 task in second.
# From: https://slurm.schedmd.com/faq.html
# By: R. Todling, 02Jun2020
my @tasks = split(',', $ARGV[0]);
my ($layout);

$layout = arbitrary(@tasks);
print join("\n", split(/,/,$layout));print("\n");

sub arbitrary {

my @tasks = @_;
my @nodes = `scontrol show hostnames $SLURM_JOB_NODELIST`;
my $node_cnt = $#nodes + 1;
my $task_cnt = $#tasks + 1;

if ($node_cnt < $task_cnt) {
  print STDERR "ERROR: You only have $node_cnt nodes, but requested layout on $task_cnt nodes.\n";
  $task_cnt = $node_cnt;
  exit(1);
}

my $cnt = 0;
my $layout;
foreach my $task (@tasks) {
  my $node = $nodes[$cnt];
  last if !$node;
  chomp($node);
  for(my $i=0; $i < $task; $i++) {
      $layout .= "," if $layout;
      $layout .= "$node";
  }
  $cnt++;
}
return $layout;
}
