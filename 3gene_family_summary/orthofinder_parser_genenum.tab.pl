#Ortholog        Lnix    Ttub    Hhol    Tcor    Juli    Gmae    Plat    Smar    Total
#OG_10000        426     0       0       0       0       0       0       0       426

open GENENUM,"Orthogroups.GeneCount.tsv";
my @sp;
my $Single_copy_in_all_species = 0;
my $Present_in_all_species = 0;
my $total_sp = 0;
my %sum;
my %number;
while(<GENENUM>){
        chomp;
        if (/^Orthogroup/){
                @sp = split;
                shift @sp;
                pop @sp;
        #       print join ("\t",@sp);
                $total_sp = $#sp+1;
                #print "$total_sp\n";
        }else{
                my @OG_number = split;
                my $OG_id = shift @OG_number;
                pop @OG_number;
                #print join ("\t",@OG_number),"\n";
                my $Single_copy_marker = 0;
                my $nonPresent_marker = 0;
                foreach (0..$#OG_number){
                        $Single_copy_marker++ if ($OG_number[$_] ==1);
                        $nonPresent_marker++ if ($OG_number[$_] ==0);
                        $number{$OG_id}{$sp[$_]} = $OG_number[$_];
                }
                foreach (0..$#sp){
                        if ($Single_copy_marker == $total_sp){
                                $sum{$sp[$_]}{'Single_copy_in_all_species'}++;
                        }elsif(!$nonPresent_marker){
                                $sum{$sp[$_]}{'Present_in_all_species'} += $number{$OG_id}{$sp[$_]};
                        }elsif($nonPresent_marker == $#sp){
                                $sum{$sp[$_]}{'No_significant_homology'}++ if ($number{$OG_id}{$sp[$_]} == 1);
                                $sum{$sp[$_]}{'Self_homology_only'} += $number{$OG_id}{$sp[$_]} if ($number{$OG_id}{$sp[$_]} > 1);
                        }else{
                                $sum{$sp[$_]}{'With_ortholog_in_at_least_1_other_species'} += $number{$OG_id}{$sp[$_]} if ($number{$OG_id}{$sp[$_]});
                        }

                }
        }
}

print "Species\tSingle_copy_in_all_species\tPresent_in_all_species\tWith_ortholog_in_at_least_1_other_species\tSelf_homology_only\tNo_significant_homology\n";

foreach (0..$#sp){
        my $sn = $_;
        open OUT,">Ortholog.$sp[$sn].xls";
        foreach my $OG_id(keys %number){
                print OUT "$OG_id\n" if ($number{$OG_id}{$sp[$sn]} > 0);
        }
        close OUT;
        print "$sp[$sn]\t";
        print "$sum{$sp[$sn]}{'Single_copy_in_all_species'}\t";
        print "$sum{$sp[$sn]}{'Present_in_all_species'}\t";
        print "$sum{$sp[$sn]}{'With_ortholog_in_at_least_1_other_species'}\t";
        print "$sum{$sp[$sn]}{'Self_homology_only'}\t";
        print "$sum{$sp[$sn]}{'No_significant_homology'}\n";
}
