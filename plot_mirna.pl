use MyMod::Bio::Seq::mirna;

my $seq=$ARGV[0];
my $mature=$ARGV[1];

my $mirna=MyMod::Bio::Seq::mirna->new(
	-display_id=>"test",
	-seq=>$seq,
	-mature=>$mature
);

print $mirna->plot_txt_simple()."\n";
my $loc = $mirna->structure_loc();
print $loc->start();
