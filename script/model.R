setwd("/users/fwang/work/mmRetNet/model/")
setwd("E://work//RetReg//model//")
fimo=read.table("zfATACnmdahPldtP_Fimo_ScBk_DiffCorMG_CandidateTFs.txt",sep = '\t')
fos2=read.table("zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp2_TFsTargets.bed",sep = '\t')
fos=read.table("zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp_FootprintInserstions.bed",sep = '\t')
#get factors and target genes binding information
fos2_2=fos2[1,]
mark=1
for(i in 1:nrow(fos2))
{
	factors=unlist(strsplit(as.character(fos2[i,1]),"[|]"))
	if(length(factors)>1)
	{
		fos2[i,1]=factors[1]
		for(j in 2:length(factors))
		{
			fos2_2[mark,]=fos2[i,]
			fos2_2[mark,1]=factors[j]
			mark=mark+1
		}
	}
	#print(i)
}

fos2=rbind(fos2,fos2_2)
save(fos2,file="zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp2_TFsTargets_splited_motifs.RData")

###Calculate footprint occupancy score (FOS)

Cal_Footprints_FOS <- function(FP1, FlankFold1=3){
  FP2 <- apply(FP1, 1, function(X1){ MotifSize1 <- as.numeric(X1[3])
     WidthL1 <- floor(MotifSize1/2); WidthR1 <- MotifSize1-WidthL1; 
     WidthL2 <- floor(MotifSize1*(2*FlankFold1+1)/2); WidthR2 <- MotifSize1*(2*FlankFold1+1)-WidthL2;
     Insertion1 <- strsplit(X1[7], ',')[[1]]; Midpoint1 <- floor(length(Insertion1)/2)
     Insertion2 <- as.numeric(Insertion1[c((Midpoint1-WidthL2+1):(Midpoint1+WidthR2))])
     
     L1 <- sum(Insertion2[1:(WidthL2-WidthL1)])/FlankFold1
     M1 <- sum(Insertion2[(WidthL2-WidthL1+1):(WidthL2+WidthR1)])
     R1 <- sum(Insertion2[(WidthL2+WidthR1+1):length(Insertion2)])/FlankFold1
     FOS1 <- min(1-(M1+1)/(L1+1), 1-(M1+1)/(R1+1));
     if(FOS1<0)
     {
     	FOS1=0
     }
     return(c(X1[1:6], as.numeric(FOS1)))
     } )
     
  return(t(FP2))  
}

FOS=Cal_Footprints_FOS(fos2)

FOS2=as.data.frame(array(0,dim=c(nrow(FOS),9)))
mark=1
for(i in 426272:nrow(FOS))
{
	factors=unlist(strsplit(FOS[i,1],";"))
	genes=unlist(strsplit(FOS[i,2],"[|]"))
	if(length(factors)==2)
	{
		FOS2[mark,1:2]=factors
		FOS2[mark,3:4]=genes
		FOS2[mark,5:9]=FOS[i,3:7]
		mark=mark+1
	}
	if(length(factors)>2)
	{	
		for(j in 2:length(factors))
		{
			FOS2[mark,1]=factors[1]
			FOS2[mark,2]=factors[j]
			FOS2[mark,3:4]=genes
			FOS2[mark,5:9]=FOS[i,3:7]
			mark=mark+1
		}
	}
}
save(FOS2,file="zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp2_TFsTargets_splited_motifs2.RData")

FOS=FOS2
FOS=unique(FOS)
save(FOS,file="zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp2_TFsTargets_splited_motifs3.RData")

#get PB matrix, transcription factors binding strength in each gene, row is gene, col is factor
load("~/work/mmRetNet/model/zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp2_TFsTargets_splited_motifs3.RData")
factors=names(table(FOS[,2]))
genes=names(table(FOS[,3]))
#remove multi genes
genes_new=c()
mark=1
for(i in 1:length(genes))
{
	gene_name=unlist(strsplit(genes[i],","))
	if(length(gene_name)==1)
	{
		genes_new[mark]=genes[i]
		mark=mark+1
	}
}

PB=array(0,dim=c(length(genes),length(factors)))
for(i in 1:length(genes))
{
	for(j in 1:length(factors))
	{
		FOS2=FOS[FOS[,2]==factors[j] & FOS[,3]==genes[i],]
		if(length(FOS2)==9)
		{
			PB[i,j]=as.numeric(FOS2[9])
		}
		if(length(FOS2)>9)
		{
			FOS_value=c(as.numeric(FOS2[,9]))
			FOS_value=1-prod(1-FOS_value)
			PB[i,j]=FOS_value
		}
	}
	print(i)
}
save(PB,file="PB.RData")

#get factor expression matrix
load("~/work/mmRetNet/model/zfAdzfNMDA_zfLD_zfTR_pbmcSubC_CondSub_MG_umap_seurat.RData")
#order cells by pseudotime
cells_inf=MG@meta.data[order(MG@meta.data[,"Pseudotime"]),]
cells_pse=cells_inf[,"Pseudotime"]
expression_matrix=MG@data[,rownames(cells_inf)]
factor_exp=expression_matrix[factors,]
save(factor_exp,file="zf_factor_exp.RData")

#get gene expression matrix, this is for each single cell. Most neighbor cells have no expression difference
gene_exp=expression_matrix[genes_new,]
gene_exp_dt=gene_exp[,1:2]
gene_exp_dt[,2]=gene_exp[,2]-gene_exp[,1]
for(i in 3:ncol(gene_exp))
{
	gene_exp_dt=cbind(gene_exp_dt,gene_exp[,i]-gene_exp[,(i-1)])
	print(i)
}
save(gene_exp_dt,file="zf_gene_exp_dt.RData")

#cut pseudotime-ordered cells into bins, calculate average expression levels of cells within each bin
gene_exp=as.matrix(gene_exp)
gene_exp_bin=array(0,dim=c(nrow(gene_exp),500))

for(i in 1:nrow(gene_exp))
{
	bin=56
	end=0
	for(j in 1:206)
	{
		start=end+1
		end=start+bin-1
		gene_exp_bin[i,j]=sum(gene_exp[i,start:end])/bin
	}
	bin=57
	for(j in 207:500)
	{
		start=end+1
		end=start+bin-1
		gene_exp_bin[i,j]=sum(gene_exp[i,start:end])/bin
	}
	print(i)
}
save(gene_exp_bin,file="zf_gene_exp_500_bin.RData")

#get nerghbor bins expression difference
gene_exp_bin_dt=gene_exp_bin[,1:2]
gene_exp_bin_dt[,2]=gene_exp_bin[,2]-gene_exp_bin[,1]
for(i in 3:ncol(gene_exp_bin))
{
	gene_exp_bin_dt=cbind(gene_exp_bin_dt,gene_exp_bin[,i]-gene_exp_bin[,(i-1)])
	print(i)
}
save(gene_exp_bin_dt,file="zf_gene_exp_500_bin_dt.RData")
gene_exp_bin_dt=gene_exp_bin_dt[,2:ncol(gene_exp_bin_dt)]

#get factors expression level within each pseudotime bin
factor_exp=as.matrix(factor_exp)
factor_exp_bin=array(0,dim=c(nrow(factor_exp),500))

for(i in 1:nrow(factor_exp))
{
	bin=56
	end=0
	for(j in 1:206)
	{
		start=end+1
		end=start+bin-1
		factor_exp_bin[i,j]=sum(factor_exp[i,start:end])/bin
	}
	bin=57
	for(j in 207:500)
	{
		start=end+1
		end=start+bin-1
		factor_exp_bin[i,j]=sum(factor_exp[i,start:end])/bin
	}
	print(i)
}
save(factor_exp_bin,file="zf_factor_exp_500_bin.RData")
factor_exp_bin=factor_exp_bin[,1:(ncol(factor_exp_bin)-1)]

#construct model for all genes, add loop
load("~/PB.RData")
load("~/work/mmRetNet/model/zf_gene_exp_500_bin_dt.RData")
load("~/work/mmRetNet/model/zf_gene_exp_500_bin.RData")
load("~/zf_factor_exp_500_bin.RData")
setwd("/users/fwang/work/mmRetNet/model/")
cor=read.table("zfATACnmdahPldtP_PFF_ScBk_Cand_FOSF_RegMT_Cor.txt",sep='\t')
rownames(PB)=genes
colnames(PB)=factors
PB2=round(PB[genes_new,],2)

weight_parm_all=c()
factor_lable_all=c()
for(i in 1:nrow(factor_exp_bin))
{
	weight_parm_all[i]=paste("w",i,sep="")
	factor_lable_all[i]=paste("factor",i,sep="")
}
nls_mod=list()
for(i in 1:10)
{
	dG=gene_exp_bin[i,2:500]
	gene_exp_pre=gene_exp_bin[i,1:499]
	factor_pos=which(PB2[i,]!=0)
	TF_binding_exp=factor_exp_bin[factor_pos,1:499]
	factor_lable=factor_lable_all[factor_pos]
	if(length(factor_lable)==1)
	{
		data=as.data.frame(cbind(dG,TF_binding_exp,gene_exp_pre))
	}
	if(length(factor_lable)>1)
	{
		data=as.data.frame(cbind(dG,t(TF_binding_exp),gene_exp_pre))
	}
	colnames(data)=c("gene",factor_lable,"gene_pre")
	weight_parm=weight_parm_all[factor_pos]
	PB_level=PB2[i,factor_pos]
	#construct regulatory formula for each gene
	if(length(factor_lable)>0)
	{
		formula_left="gene"
		formula_right1=c()
		for(j in 1:length(factor_pos))
		{
			formula_right1[j]=paste(PB_level[j],weight_parm[j],factor_lable[j],sep="*")
		}
		formula_right2=paste(formula_right1,collapse="+")
		formula_right=paste("(",formula_right2,")/(1+",formula_right2,")+w0*gene_pre+e",sep="")
		regulatory_formula=as.formula(paste(formula_left,formula_right,sep=" ~ "))
		start=c(rep(1,length(factor_lable)),1,1)
		start=as.list(start)
		names(start)=c(weight_parm,"w0","e")
		nls_mod[[i]]=nls(regulatory_formula,data=data,control=nls.control(warnOnly = T),algorithm="port",start=start,lower=c(rep(-1,length(factor_lable)),-1,-1),upper=c(rep(1,length(factor_lable)),1,1))
		#optim(c(rep(0,14),-1,0),fn=regulatory_formula,lower=rep(-1,16),upper=rep(1,16),data=data, method = "L-BFGS-B") 
	}
	print(i)
}



#construct model example, for gene ENSDARG00000000002
weight_parm=c()
for(i in 1:nrow(factor_exp_bin))
{
	weight_parm[i]=paste("w",i,sep="")
}

dG=gene_exp_bin[4325,2:500]
gene_exp_pre=gene_exp_bin[4325,1:499]
factor_pos=which(PB2[i,]!=0)
TF_binding_exp=factor_exp_bin[factor_pos,]
factor_lable=factor_lable_all[factor_pos]
data=as.data.frame(cbind(dG,t(TF_binding_exp),gene_exp_pre))
colnames(data)=c("gene",factor_lable,"gene_pre")

mod2 <- nls(gene1 ~ (w1*factor1+w2*factor2+w3*factor3+w4*factor4+w5*factor5+w6*factor6+w7*factor7+w8*factor8+w9*factor9+w10*factor10+w11*factor11)/(1+w1*factor1+w2*factor2+w3*factor3+w4*factor4+w5*factor5+w6*factor6+w7*factor7+w8*factor8+w9*factor9+w10*factor10+w11*factor11)-w0*gene1_pre+e,
	data=data)

start=c(w0=1,e=0,w1=1,w2=1,w3=1,w4=1,w5=1,w6=1,w7=1,w8=1,w9=1,w10=1,w11=1)

gene1=(w1*factor1+w2*factor2+w3*factor3+w4*factor4+w5*factor5+w6*factor6+w7*factor7+w8*factor8+w9*factor9+w10*factor10+w11*factor11+w12*factor12+w13*factor13+w14*factor14+w15*factor15+w16*factor16+w17*factor17+w18*factor18+w19*factor19+w20*factor20+w21*factor21+w22*factor22+w23*factor23+w24*factor24+w25*factor25+w26*factor26+w27*factor27+w28*factor28+w29*factor29+w30*factor30+w31*factor31+w32*factor32+w33*factor33+w34*factor34+w35*factor35+w36*factor36+w37*factor37+w38*factor38+w39*factor39+w40*factor40+w41*factor41+w42*factor42+w43*factor43+w44*factor44)/(1+w1*factor1+w2*factor2+w3*factor3+w4*factor4+w5*factor5+w6*factor6+w7*factor7+w8*factor8+w9*factor9+w10*factor10+w11*factor11+w12*factor12+w13*factor13+w14*factor14+w15*factor15+w16*factor16+w17*factor17+w18*factor18+w19*factor19+w20*factor20+w21*factor21+w22*factor22+w23*factor23+w24*factor24+w25*factor25+w26*factor26+w27*factor27+w28*factor28+w29*factor29+w30*factor30+w31*factor31+w32*factor32+w33*factor33+w34*factor34+w35*factor35+w36*factor36+w37*factor37+w38*factor38+w39*factor39+w40*factor40+w41*factor41+w42*factor42+w43*factor43+w44*factor44)-w0*gene1_pre+e
mod2 <- nls(gene1 ~ (w1*factor1+w2*factor2+w3*factor3+w4*factor4+w5*factor5+w6*factor6+w7*factor7+w8*factor8+w9*factor9+w10*factor10+w11*factor11+w12*factor12+w13*factor13+w14*factor14+w15*factor15+w16*factor16+w17*factor17+w18*factor18+w19*factor19+w20*factor20+w21*factor21+w22*factor22+w23*factor23+w24*factor24+w25*factor25+w26*factor26+w27*factor27+w28*factor28+w29*factor29+w30*factor30+w31*factor31+w32*factor32+w33*factor33+w34*factor34+w35*factor35+w36*factor36+w37*factor37+w38*factor38+w39*factor39+w40*factor40+w41*factor41+w42*factor42+w43*factor43+w44*factor44)/(1+w1*factor1+w2*factor2+w3*factor3+w4*factor4+w5*factor5+w6*factor6+w7*factor7+w8*factor8+w9*factor9+w10*factor10+w11*factor11+w12*factor12+w13*factor13+w14*factor14+w15*factor15+w16*factor16+w17*factor17+w18*factor18+w19*factor19+w20*factor20+w21*factor21+w22*factor22+w23*factor23+w24*factor24+w25*factor25+w26*factor26+w27*factor27+w28*factor28+w29*factor29+w30*factor30+w31*factor31+w32*factor32+w33*factor33+w34*factor34+w35*factor35+w36*factor36+w37*factor37+w38*factor38+w39*factor39+w40*factor40+w41*factor41+w42*factor42+w43*factor43+w44*factor44)-w0*gene1_pre+e,
	data=data,start=c(w0=1,e=0,w1=1,w2=1,w3=1,w4=1,w5=1,w6=1,w7=1,w8=1,w9=1,w10=1,w11=1,w12=1,w13=1,w14=1,w15=1,w16=1,w17=1,w18=1,w19=1,w20=1,w21=1,w22=1,w23=1,w24=1,w25=1,w26=1,w27=1,w28=1,w29=1,w30=1,w31=1,w32=1,w33=1,w34=1,w35=1,w36=1,w37=1,w38=1,w39=1,w40=1,w41=1,w42=1,w43=1,w44=1))


#plot gene and factor expression distribution
data2=as.data.frame(cbind(unlist(smooth_exp[genes_new[i],]),unlist(smooth_exp[factors[217],])))
ggplot(data = data2, mapping = aes(x = V1,
                                  y = V2)) + 
  geom_point(colour = "#426671", size = 2) +  geom_smooth(method = loess,colour='#764C29',fill='#E7E1D7')




cor(data[,1],data[,4])


