####nohup perl /users/jwang3/RetReg/Programs/Main_Footprints.pl &
#/!perl-w
$Dir1='/dcl02/qian/data/jwang3/RetReg/ATAC'; chdir $Dir1;
$Dir2='/users/jwang3/RetReg'; $ATACDir1="$Dir2/ATAC"; chdir $ATACDir1;

@ExpTech0=('scRNA', 'RNA', 'ScBk'); @ExpTech01=$ExpTech0[2]; 
@Spec1='Zf'; @ATACcond1='nmdahPldtP';
#@Spec1='Mm'; @ATACcond1='nmdaPldP';
@Spec2=lc($Spec1[0]); @Tech1='ATAC'; @GFP1='P'; 

@ATACFile1='ATAC_Files_RefinedF';
@File1=join('', ($Spec2[0],$Tech1[0],$ATACcond1[0])); @Fdr1='05';
@SourceFile1="$File1[0]_CuFiQ10No_sorted_fdr$Fdr1[0]"; 
@Wid1[0]=2; @Diff1='DiffCorMG';

if($Spec1[0] eq 'Zf'){ $Dir1=join('',($Dir1,'/Zebrafish/Footprints'));
	if($ExpTech01[0] eq 'RNA'){ @RNAcond1='nmdahPldP'; @RNAcond2=join('', ($RNAcond1[0],'5')); 
		@TargDir1="$Dir2/RNA";
	}elsif($ExpTech01[0] eq 'scRNA'){ @RNAcond1='scRNANMDALDTR'; @RNAcond2=join('', ($RNAcond1[0],'01')); 
		@TargDir1='/users/jwang3/RetReg/scRNA/Zebrafish';
	  @TargGene1='zfAdzfNMDA_zfAdzfLD_zfAdzfTR_PseudotimeDiffGenes01_CommonGenes';
	}else{ @RNAcond2='ScBk'; @TargDir1="$Dir2/RNA";
		@TargGene1='zfRNAldtPnmdaP_zfAdzfNMDAzfLDzfTR_ScBk_DiffCorTFs_Genes';
	}
	@CandTF1="Tranfac201803_$Spec1[0]_MotifTFsF_$RNAcond2[0]_$Diff1[0]"; @Genome1='GRCz10Chr'; 
}elsif($Spec1[0] eq 'Mm'){ $Dir1=join('',($Dir1,'/Mouse/Footprints')); 
	@RNAcond1='nmdaPldP'; @Genome1='GRCm38Chr'; 
}else{ print('Unknown species') }

=pod
print "generate command to merge bam files for calling footprints\n";
open(INFILE1,"$ATACFile1[0].txt")||die;
open(OUTFILE1,">$Dir2/Programs/Tool_comm.txt")||die;
@con1=<INFILE1>;

$var1='';
for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]);
	if($acc1[1]=~/$Spec1[0].*ld/i){
		$var11=join('', ($acc1[1],'_CuFiQ10No_nsorted.bam'));
	  $var1=join(' ', ($var1, $var11));
	}
}
print OUTFILE1 "cd $Dir1\n";
print OUTFILE1 "samtools cat $var1 >$Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No.bam\n";
print OUTFILE1 "samtools sort -@ 20 $Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No.bam -o $Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted.bam\n";
print OUTFILE1 "samtools index $Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted.bam\n";
print OUTFILE1 "/users/jwang3/Software/bam_compact_split_util/sam2compact.pl $Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted.sam\n";
print OUTFILE1 "perl /users/jwang3/Software/bam_compact_split_util/split.pl $Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted_compact.txt\n";
print OUTFILE1 "calcDFT /dcl02/qian/data/jwang3/Public/$Genome1[0]/ $Dir1/$Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted/$Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted >$Dir1/$Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted/$Spec2[0]ATAC$RNAcond1[0]_CuFiQ10No_sorted_dft.txt\n";
=cut

###run Main_Footprints.R to get footprints

=pod
print "merge footprints whose distance is less than 4\n";
open(INFILE1,"$File1[0]_CuFiQ10No_sorted_fdr0.$Fdr1[0]0000.bed")||die;
open(OUTFILE1,">$File1[0]_CuFiQ10No_sorted_fdr$Fdr1[0]0.txt")||die;
@con1=<INFILE1>;

$Wid2=$Wid1[0]*2; chomp $con1[1];@acc1=split(/	/,$con1[1]);
$str1=$acc1[0]; $start1=$acc1[1];$end1=$acc1[2];$pva1=$acc1[3];$no1=1;
for($i=2;$i<=$#con1;$i++){chomp $con1[$i]; @acc2=split(/	/,$con1[$i]);
  if($str1 eq $acc2[0]&($start1>=$acc2[1]&$start1<=$acc2[2]|$end1>=$acc2[1]&$end1<=$acc2[2]|$start1<=$acc2[1]&$end1>=$acc2[2]|$acc2[1]>$end1&$acc2[1]-$end1<=$Wid2|$start1>$acc2[2]&$start1-$acc2[2]<=$Wid2)){
    @acc22=sort{$a<=>$b}($start1,$end1,@acc2[1,2]);$start1=$acc22[0];$end1=$acc22[3];$pva1=$pva1+$acc2[3];$no1++;
  }
  else{$pva2=$pva1/$no1;print OUTFILE1 "$str1	$start1	$end1	$pva2\n";$str1=$acc2[0];$start1=$acc2[1];$end1=$acc2[2];@acc1=@acc2;$pva1=$acc1[3];$no1=1;}
}
$pva2=$pva1/$no1;print OUTFILE1 "$str1	$start1	$end1	$pva2\n";
=cut

=pod
print "revise results from dnase2tf\n";
open(INFILE1,"$File1[0]_CuFiQ10No_sorted_fdr$Fdr1[0]0.txt")||die;
open(OUTFILE1,">$File1[0]_CuFiQ10No_sorted_fdr$Fdr1[0].txt")||die;
@con1=<INFILE1>;

for($i=1;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]);
	$start1=$acc1[1]-$Wid1[0]; $end1=$acc1[2]+$Wid1[0];
	if($start1<1){ $start1=1; }
	print OUTFILE1 "$acc1[0]	$start1	$end1	$acc1[3]\n";
}
=cut

=pod
print "get sequences of the footprints\n";
system "bedtools getfasta -fo $Dir1/$File1[0]_CuFiQ10No_sorted_fdr$Fdr1[0].fa -fi /dcl02/qian/data/jwang3/Public/$Genome1[0]/Genome/$Genome1[0].fa -bed $Dir1/$File1[0]_CuFiQ10No_sorted_fdr$Fdr1[0].txt";
=cut

=pod
print "find motifs in the footprints\n";
open(INFILE1,"Errors.txt")||die; #$CandTF1[0]
open(OUTFILE2,">$Dir2/Programs/Fimo_All.txt")||die;
@con1=<INFILE1>;
mkdir "$File1[0]_Fimo"; $no1=1; $Step1=20;

for($i=1;$i<=$#con1;$i=$i+$Step1){ print "$i\n"; 
	print OUTFILE2 "nohup sh $Dir2/Programs/Fimo$no1.txt &\n";
	open(OUTFILE1,">$Dir2/Programs/Fimo$no1.txt")||die; $no1++;
	for($j=$i;$j<=$i+$Step1-1;$j++){ chomp $con1[$j]; @acc1=split(/	/,$con1[$j]); 
	  if($acc1[0]ne''){ print OUTFILE1 "fimo --parse-genomic-coord --max-stored-scores 2000000  --text >$Dir1/$File1[0]_Fimo/$acc1[0].txt /users/jwang3/Public/Mememotif/$acc1[0].txt $Dir1/$File1[0]_CuFiQ10No_sorted_fdr$Fdr1[0].fa\n"; }
	}
	close OUTFILE1;
}
system "sh /users/jwang3/RetReg/Programs/Fimo_All.txt";
=cut

=pod
print "check whether all motifs are included from $CandTF1[0]\n";
open(INFILE1,"$CandTF1[0].txt")||die;
open(OUTFILE1,">Errors.txt")||die;
@con1=<INFILE1>;

opendir(FOLDER1, "$Dir1/$File1[0]_Fimo")||die;
@DirFile1=grep(!/^\.\.?$/, readdir FOLDER1); 

for($i=0;$i<=$#DirFile1;$i++){ $hash1{$DirFile1[$i]}=$DirFile1[$i]; }

print OUTFILE1 $con1[0];
for($i=1;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]); $var1=join('',($acc1[0],'.txt'));
	if(!exists $hash1{$var1}){ print "$acc1[1] no\n"; print OUTFILE1 "$con1[$i]\n"; }
}
=cut

chdir $Dir1;
=pod
print "Combine all footprints of motifs from $CandTF1[0]\n";
print "Get size of $Spec1[0] $ATACcond1[0] motifs\n";
open(INFILE1,"$ATACDir1/$CandTF1[0].txt")||die; 
open(INFILE2,"$Dir2/Public/Tranfac201803_MotifPWM.txt")||die;
open(OUTFILE1,">$File1[0]_Fimo_$ExpTech01[0]_$Diff1[0].txt")||die;
open(OUTFILE2,">Errors.txt")||die;
@con1=<INFILE1>; @con2=<INFILE2>; ($no02,@acc02)=PosFea('>',@con2);

for($i=0;$i<$no02;$i++){
	if($con2[$acc02[$i]]=~/>(.*?) /){
	  $hash2{$1}=$acc02[$i+1]-$acc02[$i]-1;
	}else{ print OUTFILE2 $con2[$acc02[$i]]; }
}

chomp $con1[0]; $con1[0]=~s/\r//; print OUTFILE12 "$con1[0]	MotifSize\n";
for($i=1;$i<=$#con1;$i++){ chomp $con1[$i]; $con1[$i]=~s/\r//; @acc1=split(/	/,$con1[$i]); print "$i $acc1[0]\n"; 
	open(INFILE12,"$File1[0]_Fimo/$acc1[0].txt")||die; $var02=<INFILE12>;
  while(<INFILE12>){ chomp $_; @acc2=split(/	/,$_); $var2=join('	',@acc2[2..5,7,9,0]);
  	print OUTFILE1 "$var2\n";
  }
  close INFILE12;
  
  if(!exists $hash2{$acc1[0]}){ print OUTFILE2 "$acc1[0]\n"; }
}
=cut

=pod
print "Get bed file for differential or correlated peaks";
open(INFILE1,"$ATACDir1/$File1[0]_Diff_Fc1Fdr1_DiffCor.txt")||die; 
open(OUTFILE1,">$ATACDir1/$File1[0]_Diff_Fc1Fdr1_DiffCor.bed")||die;

$var01=<INFILE1>; @acc01=split(/	/,$var01);
for($i=0;$i<=$#acc01;$i++){
	if($acc01[$i]eq'chr'){ $Ind1=$i; last; }
}

while(<INFILE1>){ $var1=$_; chomp $var1; @acc1=split(/	/,$var1);
	$var12=join('	',@acc1[$Ind1+1..$Ind1+3]); print OUTFILE1 "$var12\n";
}
=cut

=pod
print "Overlap differential peaks and motif footprints\n";
system "bedtools intersect -a $Dir1/$File1[0]_Fimo_$ExpTech01[0]_$Diff1[0].txt -b $ATACDir1/$File1[0]_Diff_Fc1Fdr1_DiffCor.bed -wa -wb >$Dir1/$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0].txt";

print "Merge and extend footprint regions from $File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]\n";
open(INFILE1,"$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0].txt")||die; 
open(INFILE2,"$ATACDir1/$CandTF1[0].txt")||die; 
open(OUTFILE1,">$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Merged.bed")||die;
open(OUTFILE2,">$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Merged.txt")||die;
open(OUTFILE3,">Errors.txt")||die;
@con1=<INFILE1>; @con2=<INFILE2>; %hash2=HashC(0,0,4,4,@con2);

%hash1=%hash12=();
for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]); 
	$var1=join('	',@acc1[0..3]); 
	if(exists $hash2{$acc1[6]}){ $var12=join(';', (@acc1[6],$hash2{$acc1[6]})); 
	}else{ print OUTFILE3 "$acc1[6]\n"; }
	if(exists $hash1{$var1}){ $hash1{$var1}=$hash1{$var1}+1; $hash12{$var1}{$hash1{$var1}}=$var12;
  }else{ $hash1{$var1}=0; $hash12{$var1}{$hash1{$var1}}=$var12; }
}

@key1=keys %hash1; @key2=sort{$a cmp $b}@key1;
for($i=0;$i<=$#key2;$i++){ @key21=split(/	/,$key2[$i]); @acc2=();
	$num1=int(($key21[1]+$key21[2])/2); $size1=$key21[2]-$key21[1]+1; 
  $num11=$num1-$size1*5; $num12=$num1+$size1*5;
	for($j=0;$j<=$hash1{$key2[$i]};$j++){
		@acc2=(@acc2, $hash12{$key2[$i]}{$j});
  }
	$var2=join('|',@acc2); print OUTFILE1 "$key21[0]	$num11	$num12\n";
	print OUTFILE2 "$key21[0]	$key21[1]	$key21[2]	$key21[3]	$num11	$num12	$var2\n";
}
close INFILE1; close INFILE2; close OUTFILE1; close OUTFILE2; close OUTFILE3;
=cut

=pod
print "Get footprint-related genes\n";
system "perl /users/jwang3/RetReg/Programs/Deal_annovar.pl $Dir1 $File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Merged_NoStr txt 0 $Spec2[0] F T";
=cut

=pod
print "Get strand-specific footprints\n";
open(INFILE1,"$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Merged.txt")||die; 
open(OUTFILE1,">$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Merged_StrG.txt")||die;
open(OUTFILE2,">$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Merged_Str.bed")||die;

while(<INFILE1>){ chomp $_; @acc1=split(/	/, $_); @acc10=split(/;/, $acc1[0]); 
	@acc11=split(/;/, $acc1[1]); @acc12=split(/;/, $acc1[2]); @acc20=@acc21=();
	for($i=0;$i<=$#acc11;$i++){ @acc111=split(/,/, $acc11[$i]); @acc121=split(/,/, $acc12[$i]); @acc211=();
		for($j=0;$j<=$#acc121;$j++){
			if($acc121[$j] eq $acc1[6]){ @acc211=(@acc211, @acc111[$j]); }
		}
		if($#acc211>=0){ @acc20=(@acc20, @acc10[$i]); @acc21=(@acc21, join(',', @acc211)); }
	}
	if($#acc20>=0){ $var20=join(';', @acc20); $var21=join(';', @acc21); 
		$var1=join('	', ($var20, $var21, @acc1[3..5,7..8,6,9]));
		print OUTFILE1 "$var1\n"; print OUTFILE2 "$acc1[3]	$acc1[7]	$acc1[8]\n"; 
}	}
close INFILE1; close OUTFILE1; close OUTFILE2;
=cut

#=pod
print "Get candidate genes/TFs-related peaks\n";
open(INFILE1,"$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_MergedG.txt")||die; 
open(INFILE12,"$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Merged.bed")||die; 
open(INFILE21,"$ATACDir1/$CandTF1[0]_TFList.txt")||die;
open(INFILE22,"$TargDir1[0]/$TargGene1[0].txt")||die;
open(OUTFILE1,">$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Candidates.txt")||die;
open(OUTFILE12,">$File1[0]_PeaksFootprints_Fimo_$ExpTech01[0]_Candidates.bed")||die;
@con1=<INFILE1>; @con12=<INFILE12>; @con21=<INFILE21>;  @con22=<INFILE22>;
@con2=(@con21, @con22); %hash2=HashC(0,0,1,1,@con2);

for($i=0;$i<=$#con1;$i++){ chomp $con1[$i]; @acc1=split(/	/,$con1[$i]); @acc11=split(/;/,@acc1[0]); @acc12=split(/;/,@acc1[1]); @acc131=@acc132=();
	for($j=0;$j<=$#acc11;$j++){ @acc121=split(/,/,$acc12[$j]); @acc122=();
		for($k=0;$k<=$#acc121;$k++){
			if(exists $hash2{$acc121[$k]}){
				@acc122=(@acc122, $acc121[$k]);
		}	}
		if($#acc122>=0){ $var11=join('|', (join(',',@acc122), $acc11[$j]));
			@acc131=(@acc131, $var11);
	}	}
	if($#acc131>=0){ $var131=join(';', @acc131);
		$var1=join('	', (@acc1[2..4],$var131,@acc1[5..$#acc1]));
		print OUTFILE1 "$var1\n"; print OUTFILE12 $con12[$i]; 
}	}
#=cut

###Go to Main_Footprints_FOS.pl

close INFILE1; close INFILE12; close INFILE2; close INFILE21; close INFILE22; close OUTFILE1; close OUTFILE12; close OUTFILE2;


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
