#!/usr/bin/perl -w

# Use this to parse the output from the genotyping tool MIRNA_SEARCH

use strict;
use Getopt::Long; 


my ($genotype_in);
my $read=10;# need at least 10 reads
my $perc=2;# the number of reads need to be 2% of what is mapped to that gene
my $nb_genes=2; # at least 2 genes need to multiple strains for it to count

&GetOptions(
	    'in:s'      => \$genotype_in,#the genotype file from MIRNA_SEARCH
	    'read:i'   => \$read,#the read threshold (default = 10)
	    'perc:i'  => \$perc, #perc of reads mapped (default =2)
	    'genes:i'  => \$nb_genes, # nunmber of genes where the multiple strains have been observed
           );

open(IN,"<$genotype_in")||die "Can't open $genotype_in\n";
my (%mapped,%gene_genotype);

while(<IN>){
  chomp($_);
  my @element=split(/\t/,$_);
  my $gene=(split(/-/,$element[2]))[0];
  $gene=~s/>//g;
  if ($element[2] !~/\^/){
    if ($mapped{$gene}){
      $mapped{$gene}=$mapped{$gene}+$element[0];
    }else{
      $mapped{$gene}=$element[0];
    }
  }
  if ($element[0]>$read && $element[2] !~/\^/){# the genotypes with ^ correspond to recombinants that are not taken into account
    #print "$element[0] $element[2]\n";
    $gene_genotype{$gene}{$element[2]}=$element[0];
  }
}

# calculate the percentage
my (%genotype_cnt);
for my $gene (keys %gene_genotype){
  for my $genotype (keys %{$gene_genotype{$gene}}){
    my $perc_genotype=($gene_genotype{$gene}{$genotype}/$mapped{$gene})*100;
    #print "$genotype $perc_genotype\n";
    if ($perc_genotype>2){
      $genotype_cnt{$gene}++;    
    }
  }
}

# count the number of times you see genes with 1 or more genotype
# see genes with 2 or more, 3 or more etc...
my @cnt=(1..50); # the max number of genotypes is 14 at the moment
my $max_strain=0;
my %strain_cnt;
foreach my $cnt (@cnt){
  foreach my $gene (keys %genotype_cnt){
    #print "$gene $genotype_cnt{$gene}\n";
    if ($genotype_cnt{$gene}>=$cnt){
      $strain_cnt{$cnt}++;
    }
  }
}

foreach my $cnt (keys %strain_cnt){
  if ($strain_cnt{$cnt}>=$nb_genes && $cnt>$max_strain){
    #print "$cnt $strain_cnt{$cnt}\n";
    $max_strain=$cnt;
  }  
}
print "$max_strain\n";