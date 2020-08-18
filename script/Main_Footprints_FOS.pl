####nohup perl /users/jwang3/RetReg/Programs/Main_Footprints_FOS.pl &
#/!perl-w
@Dir1='/dcl02/qian/data/jwang3/RetReg/ATAC'; 
@Dir2='/users/jwang3/RetReg';

@Spec1='Zf'; @ATACcond1='nmdahPldtP'; 
@ExpTech0=('scRNA', 'RNA', 'ScBk'); @ExpTech01=$ExpTech0[2]; 
#@Spec1='Mm';@ATACcond1='nmdaPldP'; 
@Spec2=lc($Spec1[0]); @Tech1='ATAC'; @GFP1='P'; 

@ATACFile1='ATAC_Files_RefinedF';
@File1=join('', ($Spec2[0],$Tech1[0],$ATACcond1[0])); @Fdr1=('05');
@PeaksFps1="$File1[0]_PFF_$ExpTech01[0]_Cand";
@ModuleF1="zfRNAldtPnmdaP_zfAdzfNMDAzfLDzfTR_ScBk_G40ClusterMT";
@Wid1[0]=2;

if($Spec1[0] eq 'Zf'){ $Dir1[0]=join('', ($Dir1[0],'/Zebrafish'));
	@PeaksDir1=join('', ($Dir1[0],'/Peaks'));
	@FootprintsD1=join('', ($Dir1[0],'/Footprints'));
	@Genome1='GRCz10Chr'; @RNAcond1='ldPnmdahP'; @Cond3=@RNAcond1;
}elsif($Spec1[0] eq 'Mm'){ $Dir1[0]=join('', ($Dir1[0],'/Mouse')); 
	@Genome1='GRCm38Chr'; @RNAcond1='ldPnmdaP'; @Cond3='nmdaP';
}else{ print('Unknown species') }
chdir $FootprintsD1[0];

=pod
print "Get partial ATACseq samples\n";
open(INFILE1,"$Dir2[0]/ATAC/$ATACFile1[0].txt")||die;
open(OUTFILE1,">$Dir2[0]/ATAC/$ATACFile1[0]_$Spec1[0].txt")||die;
@con1=<INFILE1>;

for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]);
  if($acc1[1]=~/$Spec2[0]$Tech1[0]/){ print "$acc1[1]\n";
  	print OUTFILE1 "$con1[$i]\n";
} }
close INFILE1; close OUTFILE1;
=cut

=pod
print "Get insertions for footprint regions\n";
open(INFILE1,"$Dir2[0]/ATAC/$ATACFile1[0]_$Spec1[0].txt")||die;
open(OUTFILE1,">$Dir2[0]/Programs/Tool_comm.txt")||die;
@con1=<INFILE1>;

for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]); print "$acc1[1]\n";
  print OUTFILE1 "nohup /users/jwang3/RetReg/Programs/dnase_wig_tracks_both2.py $FootprintsD1[0]/$PeaksFps1[0].bed $PeaksDir1[0]/$acc1[1]_CuFiQ10No_sorted.bam $FootprintsD1[0]/$File1[0]_Fimo2/$acc1[1]_$PeaksFps1[0]_cuts.wig &\n";
}
close INFILE1; close OUTFILE1;
system "sh $Dir2[0]/Programs/Tool_comm.txt";
=cut

=pod
print "Generate commands to transfer wig to multiple rows\n";
open(INFILE1,"$Dir2[0]/ATAC/$ATACFile1[0]_$Spec1[0].txt")||die;
open(OUTFILE1,">$Dir2[0]/Programs/Tool_comm.txt")||die;
open(OUTFILE2,">$FootprintsD1[0]/$File1[0]_Fimo2/Temp.txt")||die;
@con1=<INFILE1>;

for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]);
  print OUTFILE1 "nohup perl /users/jwang3/RetReg/Programs/Trans_WigToMultirows.pl $FootprintsD1[0]/$File1[0]_Fimo2 $acc1[1]_$PeaksFps1[0]_cuts wig &\n";
}
close INFILE1;
system "sh $Dir2[0]/Programs/Tool_comm.txt";
=cut

=pod
print "Add the size of motifs\n";
open(INFILE1,"$Dir2[0]/ATAC/$ATACFile1[0]_$Spec1[0].txt")||die;
open(INFILE2,"$PeaksFps1[0].txt")||die;
open(OUTFILE2,">Errors.txt")||die;
@con1=<INFILE1>; @con2=<INFILE2>;

for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]); print "$acc1[1]\n";
	open(INFILE3,"$File1[0]_Fimo2/$acc1[1]_$PeaksFps1[0]_cutsp.bed")||die;
	open(OUTFILE1,">$File1[0]_Fimo2/$acc1[1]_$PeaksFps1[0]_cutsp2.bed")||die;
	%hash3=();
	while(<INFILE3>){	$var3=$_; chomp $var3; @acc3=split(/	/,$var3); $var31=join('	',@acc3[0..2]); 
		$var32=join(',',@acc3[3..$#acc3]); $hash3{$var31}=$var32;
	}
	
	for($j=0;$j<=$#con2;$j++){ chomp $con2[$j]; @acc2=split(/	/,$con2[$j]); $var2=join('	',@acc2[0,4,5]); $size1=$acc2[2]-$acc2[1]+1;
		if(exists $hash3{$var2}){
			print OUTFILE1 "$acc2[$#acc2]	$acc2[3]	$size1	$acc2[0]	$acc2[1]	$acc2[2]	$hash3{$var2}\n";
		}else{ print OUTFILE2 "$acc1[1]	$j\n"; }
	}
	close INFILE3; close OUTFILE1;
}
close INFILE1; close INFILE2; close OUTFILE2;
=cut

=pod
print "Calculate footprint occupancy score (FOS)\n";
open(INFILE1,"$Dir2[0]/ATAC/$ATACFile1[0]_$Spec1[0].txt")||die;
open(OUTFILE1,">$Dir2[0]/Programs/Tool_comm.txt")||die;
open(OUTFILE2,">$FootprintsD1[0]/$File1[0]_Fimo2/Temp.txt")||die;
@con1=<INFILE1>;

for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]); print "$acc1[1]\n";
	print OUTFILE1 "nohup Rscript $Dir2[0]/Programs/Cal_Footprints_FOS2.R $FootprintsD1[0]/$File1[0]_Fimo2/$acc1[1]_$PeaksFps1[0]_cutsp2.bed &\n";
}
close INFILE1; close OUTFILE1;
system "sh $Dir2[0]/Programs/Tool_comm.txt";
=cut

###Go to Main_Footprints_FOS.R
##system "Rscript $Dir2[0]/Programs/Main_Footprints_FOS.R Cal_Footprints_FOS Spec1";
##system "Rscript $Dir2[0]/Programs/Main_Footprints_FOS.R Combine_Footprints_FOS Spec1";
##system "Rscript $Dir2[0]/Programs/Main_Footprints_FOS.R Filter_Footprints Spec1";

=pod
print "Get the potential regulation between candidate TFs and genes\n";
open(INFILE1,"$PeaksFps1[0]_FOSF.txt")||die;
open(OUTFILE1,">$PeaksFps1[0]_FOSF_Reg.txt")||die;
@con1=<INFILE1>;

print OUTFILE1 "TF	Motif	Target	TypeChrStartEndFOS\n";
for($i=1;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/, $con1[$i]); 
	@acc11=split(/\|/, $acc1[1]); @acc12=split(/;/, $acc1[2]); @acc121=();
	for($j=0;$j<=$#acc12;$j++){ @acc122=$acc12[$j]=~/(.*?)\|(.*)$/; @acc1222=split(/,/,@acc122[0]);
		$var122=join(',', ($acc122[1],@acc1[4..6,$#acc1]));
		for($k=0;$k<=$#acc1222;$k++){
		  @acc121=(@acc121, join('	',($acc1222[$k],$var122)));
	}	}
	
	for($j=0;$j<=$#acc11;$j++){ @acc112=split(/;/, $acc11[$j]);
		for($k=1;$k<=$#acc112;$k++){
		  for($n=0;$n<=$#acc121;$n++){
		  	print OUTFILE1 "$acc112[$k]	$acc112[0]	$acc121[$n]\n";
		} }		
}	}
close INFILE1; close OUTFILE1; 
=cut

=pod
print "Merge the same pairs of TFs and genes\n";
open(INFILE1,"$PeaksFps1[0]_FOSF_Reg.txt")||die;
open(OUTFILE1,">$PeaksFps1[0]_FOSF_RegM.txt")||die;
@con1=<INFILE1>;

%hash1=%hash12=();
for($i=1;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/, $con1[$i]); $var1=join('	', @acc1[0,2]); $var2=join(',', @acc1[1,3]);
	if(exists $hash1{$var1}){ $hash1{$var1}=$hash1{$var1}+1; $hash12{$var1}{$hash1{$var1}}=$var2;
  }else{ $hash1{$var1}=0; $hash12{$var1}{$hash1{$var1}}=$var2; }
}

print OUTFILE1 "TF	Target	MotifTypeChrStartEndFOS\n";
@key1=keys %hash1; @key2=sort{$a cmp $b}@key1;
for($i=0;$i<=$#key2;$i++){ @acc1=();
	for($j=0;$j<=$hash1{$key2[$i]};$j++){ @acc1=(@acc1, $hash12{$key2[$i]}{$j}); }
	$var1=join(';', @acc1); print OUTFILE1 "$key2[$i]	$var1\n";
}
close INFILE1; close OUTFILE1; 
=cut

=pod
print "Merge TFs and genes\n";
open(INFILE1,"$PeaksFps1[0]_FOSF_RegM.txt")||die;
open(OUTFILE11,">$PeaksFps1[0]_FOSF_RegM_TF.txt")||die;
open(OUTFILE12,">$PeaksFps1[0]_FOSF_RegM_Targ.txt")||die;
@con1=<INFILE1>;

%hash1=%hash12=%hash2=%hash22=();
for($i=1;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/, $con1[$i]); $var1=join('	',@acc1[1..2]); $var2=join('	',@acc1[0,2]);
	if(exists $hash1{$acc1[0]}){ $hash1{$acc1[0]}=$hash1{$acc1[0]}+1; $hash12{$acc1[0]}{$hash1{$acc1[0]}}=$var1;
  }else{ $hash1{$acc1[0]}=0; $hash12{$acc1[0]}{$hash1{$acc1[0]}}=$var1; }
  	
  if(exists $hash2{$acc1[1]}){ $hash2{$acc1[1]}=$hash2{$acc1[1]}+1; $hash22{$acc1[1]}{$hash2{$acc1[1]}}=$var2;
  }else{ $hash2{$acc1[1]}=0; $hash22{$acc1[1]}{$hash2{$acc1[1]}}=$var2; }
}

print OUTFILE11 "TF	Target	MotifTypeChrStartEndFOS\n";
print OUTFILE12 "Target	TF	MotifTypeChrStartEndFOS\n";
for($i1=0;$i1<=1;$i1++){
	if($i1==0){ %hash3=%hash1; %hash32=%hash12; 
	}else{ %hash3=%hash2; %hash32=%hash22; }
	
  @key1=keys %hash3; @key2=sort{$a cmp $b}@key1; 
  print "$#key2\n";
  for($i=0;$i<=$#key2;$i++){ @acc11=@acc12=();
  	for($j=0;$j<=$hash3{$key2[$i]};$j++){ @acc1=split(/	/, $hash32{$key2[$i]}{$j});
  		@acc11=(@acc11, $acc1[0]); @acc12=(@acc12, $acc1[1]);
  	}
  	$var11=join('|', @acc11); $var12=join('|', @acc12); 
  	if($i1==0){ print OUTFILE11 "$key2[$i]	$var11	$var12\n";
    }else{ print OUTFILE12 "$key2[$i]	$var11	$var12\n"; }
} }  
close INFILE1; close OUTFILE11; close OUTFILE12; 
=cut

=pod
print "Get regulation of TFs to modules\n";
open(INFILE1,"$ModuleF1[0]_EnTFs5.txt")||die;
open(OUTFILE1,">$ModuleF1[0]_EnTFs5_Reg.txt")||die;
@con1=<INFILE1>;

$Thr1=5;
chomp $con1[0]; @acc01=split(/	/, $con1[0]);
chomp $con1[1]; @acc012=split(/	/, $con1[1]);
@Ind1=@Name1=();
if($#acc01 != $#acc012){ @acc01=('ID', @acc01); }
for($i=0;$i<=$#acc01;$i++){
	if($acc01[$i]=~/^([PN])fdr(\d+)$/){
		@Name1=(@Name1, join('', ($1,$2))); @Ind1=(@Ind1, $i);
}	}

print OUTFILE1 "TF	TFSymbol	TFGroup	TFScBk	TFPeakExpTime	TargetModule	TargetGroup	Regulation	Number	Nlogfdr\n";
for($i=1;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/, $con1[$i]);
	for($j=0;$j<=$#Ind1;$j++){	$Fdr1=-log($acc1[$Ind1[$j]])/log(10);
		if($Fdr1 > $Thr1){
			$var1=join('	', @acc1[0..4]); @Name2=$Name1[$j]=~/^([PN])(\d+)$/;
			if($Name2[0]eq'P'){ $Name3='Positive'; }elsif($Name2[0]eq'N'){ $Name3='Negative'; }else{ print "Errors: $Name1[$j]\n"; }
			print OUTFILE1 "$var1	$Name2[1]	$Name2[1]	$Name3	$acc1[$Ind1[$j]-2]	$Fdr1\n";
}	}	}
=cut

close INFILE1;close INFILE2;close OUTFILE1;close OUTFILE2;


sub HashC{
	my($ind00,$ind01,$ind02,$ind03,@con0)=@_;my %hash01;
	for(my $i0=0;$i0<=$#con0;$i0++){chomp $con0[$i0];my @acc01=split(/	/,$con0[$i0]);
	  my $var00=join('	',@acc01[$ind00..$ind01]);my $var01=join('	',@acc01[$ind02..$ind03]);$hash01{$var00}=$var01;	  
	}
	return (%hash01);
}

sub HashR{
	my($ind00,$ind01,$ind02,@con0)=@_;my %hash01;
	chomp $con0[0] ;my @acc01=split(/	/,$con0[$ind00]);
	for(my $i01=0;$i01<=$#acc01;$i01++){my @acc03=();
	  for(my $i02=$ind01;$i02<=$ind02;$i02++){chomp $con0[$i02];my @acc02=split(/	/,$con0[$i02]);
	    @acc03=(@acc03,$acc02[$i01]);
	  }
	  $hash01{$acc01[$i01]}=join('	',@acc03);
	}
	return (%hash01);
}

sub PosFea{
	my($fea0,@con0)=@_;my @acc00=();my $no00=0;
	for(my $i0=0;$i0<=$#con0;$i0++){if($con0[$i0]=~/$fea0/){$acc00[$no00]=$i0;$no00++;}}
	$acc00[$no00]=$#con0+1;return ($no00,@acc00);
}

sub Nonred{
	my(@con0)=@_;my @acc00=();$acc00[$no00]=0;my $no00=0;
	for(my $i0=1;$i0<=$#con0;$i0++){my $flag0=0;
		for(my $j0=0;$j0<=$no00;$j0++){if($con0[$i0]eq$con0[$acc00[$j0]]){$flag0=1;last;}}
		if(!$flag0){$acc00[$no00+1]=$i0;$no00++;}
	}
	return (@acc00);
}

sub ChrFeaH{###$num0:the column of chr;$num01:the number of splitting
	my($num0,$num01,@con0)=@_;my @pos0=();my %hash0=();
	if($con0[0]=~/	[Cc]hrom.*?	/){$pos0=1;}else{$pos0=0;}
	chomp $con0[$pos00[0]];my @acc00=split(/	/,$con0[$pos00[0]]);my @chr00=();
	my $no00=1;my @temp1=split(//,$acc00[$num0+1]);my $ID0;
	if($num01==-1){$ID0=$acc00[$num0];}else{$ID0=join('',($acc00[$num0],'_',@temp1[0..$num01],(0)x($#temp1-$num01)));}
	for(my $i0=$pos0+1;$i0<=$#con0;$i0++){chomp $con0[$i0];my @acc01=split(/	/,$con0[$i0]);
		my @temp1=split(//,$acc01[$num0+1]);my $ID01;
		if($num01==-1){$ID01=$acc01[$num0];}else{$ID01=join('',($acc01[$num0],'_',@temp1[0..$num01],(0)x($#temp1-$num01)));}
		if($ID0 ne $ID01){$hash0{$ID0}=join('	',($pos0,$i0));$pos0=$i0;$ID0=$ID01;}
	}
	$hash0{$ID0}=join('	',($pos0,$#con0+1));
	return (%hash0);
}

sub ChrFea{
	my($fea0,@con0)=@_;my $pos00=0;my %hash0;
	my @acc00=$con0[$pos00]=~/$fea0/;
	for(my $i0=$pos00+1;$i0<=$#con0;$i0++){my @acc01=$con0[$i0]=~/$fea0/;
		if($acc01[0]=~/^\d+$/&$acc01[0]!=$acc00[0]|$acc01[0]ne$acc00[0])
		{$hash0{$acc00[0]}=join('	',($pos00,$i0));@acc00=@acc01;$pos00=$i0;}
	}
	$hash0{$acc00[0]}=join('	',($pos00,$#con0+1));
	return %hash0;
}

sub Search{
	my($var0,$step0,@con0)=@_;my $ind0=-1;
	for(my $i0=0;$i0<=$#con0;$i0=$i0+$step0){
		if($con0[$i0]=~/^$var0\W/i|$con0[$i0]=~/\W$var0\W/i|$con0[$i0]=~/\W$var0$/i|$con0[$i0]=~/^$var0$/i)
		{$ind0=$i0;last;}
	}
	return ($ind0);
}

sub RevSeq{
  my ($seq00)=@_;my @seq0=split(//,$seq00);my $leng0=$#seq0;my @seq02=();my $no0=0;
  for(my $i0=$leng0;$i0>=0;$i0=$i0-1){
  	if($seq0[$i0]eq'A'){$seq02[$no0]='T';}elsif($seq0[$i0]eq'a'){$seq02[$no0]='t';}
  	elsif($seq0[$i0]eq'T'){$seq02[$no0]='A';}elsif($seq0[$i0]eq't'){$seq02[$no0]='a';}
  	elsif($seq0[$i0]eq'C'){$seq02[$no0]='G';}elsif($seq0[$i0]eq'c'){$seq02[$no0]='g';}
  	elsif($seq0[$i0]eq'G'){$seq02[$no0]='C';}elsif($seq0[$i0]eq'g'){$seq02[$no0]='c';}
  	elsif($seq0[$i0]eq'-'){$seq02[$no0]='-';}elsif($seq0[$i0]eq'N'|$seq0[$i0]eq'n'){$seq02[$no0]=$seq0[$i0];}
  	else{print $seq0[$i0];}$no0++;
  }
  my $seq03=join('',@seq02);return $seq03;
}

sub Min1{
  my (@con0)=@_;$var0=$con0[0];
  for(my $i0=1;$i0<=$#con0;$i0++){
    if($var0>$con0[$i0]){$var0=$con0[$i0];}
  }
  return $var0;
}

sub Max1{
  my (@con0)=@_; $var0=$con0[0];
  for(my $i0=1;$i0<=$#con0;$i0++){
    if($var0<$con0[$i0]){$var0=$con0[$i0];}
  }
  return $var0;
}

sub RandN{
  my ($rang0,$num0)=@_;my %has0=();my @seq0=();
  if($rang0<$num0){print "Errors\n";return -1;}
  else{
    for(my $i0=0;$i0<$num0;$i0++){
     my $range=100;my $var0=int(rand($rang0))+1;
     if(!exists $has0{$var0}){@seq0=(@seq0,$var0);$has0{$var0}=1;}
     else{$num0++;}
    }
    return @seq0;
} }
