use 5.010;
use strict;
use warnings;

use Bio::SeqIO;

my $fcds="$ARGV[1]/$ARGV[0]";
my $outfile="$ARGV[1]/megahitCoassembly_$ARGV[2].fasta";
my $sampleID="$ARGV[2]";

my $in  = Bio::SeqIO-> new ( -file   => $fcds,
                             -format => 'fasta' );
my $out = Bio::SeqIO-> new ( -file   => ">>  $outfile",
                              -format => 'fasta' );

while (my $record = $in->next_seq()) {
          my $contig_id = $record->display_id();
          my $replace = "$sampleID"."-ContigCoA-$contig_id";
          $record->display_id($replace);
          $record->desc("");
          $out->write_seq($record);
}

