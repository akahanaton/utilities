use strict;
use warnings;
use POSIX ":sys_wait_h";
use Data::Dumper;

my $max_children = 8;
my %work = map { $_ => 1 } 1 .. 1000;
print scalar keys %work, "\n";
my @works = keys %work;

my $ppid = getppid();
print "Parent Process ID $ppid\n";

my %pids;
while (%work) {
    #while there are still empty slots
    while (@works and keys %pids < $max_children) {
#--------------------------------------------------
#     while (@work) {
#-------------------------------------------------- 
	  #get some work for the child to do
	  my $cur_work = shift @works;
	  die "could not fork" unless defined(my $pid = fork);

	  #parent
	  if ($pid) {
		$pids{$pid} = 1;
		next;
	  }

	  #child
	#--------------------------------------------------
	#   print "$$ doing work $cur_work\n";
	#-------------------------------------------------- 
	#--------------------------------------------------
	#   system("blat query.fa.cut/query.fa.$work tigr7.0.rmsk/all.chrs.con query.fa.cut/query.fa.$work.blast9 -out=blast9");
	#-------------------------------------------------- 
	  sleep .1;
	#--------------------------------------------------
	#   print "$$ done doing work $cur_work\n";
	#-------------------------------------------------- 
	  print "$cur_work\t$pid\t", scalar keys %pids,"\t",scalar keys %work,"\t", scalar @works,"\n";
	  exit $cur_work;
    }

    my $pid = waitpid -1, WNOHANG;
#--------------------------------------------------
#     print "$ppid\t$pid\n";
#-------------------------------------------------- 
    if ($pid > 0) {
	  delete $pids{$pid};
	#--------------------------------------------------
	#   print "\$\?:\t $?\t";
	#-------------------------------------------------- 
	  my $rc = $? >> 8; #get the exit status
	#--------------------------------------------------
	#   print "\$rc:\t $rc\n";
	#-------------------------------------------------- 
	#--------------------------------------------------
	#   print "saw $pid was done with $rc\n";
	#-------------------------------------------------- 
	  delete $work{$rc};
	#--------------------------------------------------
	#   print "work left: ", join(", ", sort keys %work), "\n";
	#-------------------------------------------------- 
    }elsif ($pid  == -1){
	  last;
    }
}
