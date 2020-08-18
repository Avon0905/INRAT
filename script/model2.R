#construct model for all genes, add loop
setwd("E://work//RetReg//model")
load("PB.RData")
load("zf_gene_exp_500_bin_dt.RData")
load("zf_gene_exp_500_bin.RData")
load("zf_factor_exp_500_bin.RData")
smooth_exp=read.table("zfAdzfNMDA_zfLD_zfTR_pbmcSubC_CondSub_MG_Bin50_SmoothExp.txt",sep='\t')
load("zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp2_TFsTargets_splited_motifs3.RData")
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

i=4325
for(i in 1:10)
{
	dG=gene_exp_bin_dt[i,2:500]
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
	#get correlation in scRNA-seq, use correlation as start value
	cor=c()
	if(length(factor_pos)>1)
	{
		for(j in 1:length(factor_pos))
		{
			factor_smooth_exp=smooth_exp[factors[factor_pos[j]],]
			target_smooth_exp=smooth_exp[genes_new[i],]
			cor[j]=cor(unlist(factor_smooth_exp),unlist(target_smooth_exp))
		}
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
		formula_right3=c()
		for(j in 1:length(factor_pos))
		{
			formula_right3[j]=paste("abs(",formula_right1[j],")",sep="")
		}
		formula_right4=paste(formula_right3,collapse="+")
		formula_right=paste("alpha*(",formula_right2,")/(1+",formula_right4,")+w0*gene_pre+e",sep="")
		regulatory_formula=as.formula(paste(formula_left,formula_right,sep=" ~ "))
		start=c(1,cor,1,0)
		start=as.list(start)
		names(start)=c("alpha",weight_parm,"w0","e")
		nls_mod[[i]]=nls(regulatory_formula,data=data,control=nls.control(warnOnly = T),algorithm="port",start=start,lower=c(-Inf,rep(-1,length(factor_lable)),-1,-1),upper=c(Inf,rep(1,length(factor_lable)),1,1))
		#optim(c(rep(0,14),-1,0),fn=regulatory_formula,lower=rep(-1,16),upper=rep(1,16),data=data, method = "L-BFGS-B") 
	}
	print(i)
}



