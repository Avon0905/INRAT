##### Main process #####
setwd("G://data//model//zebrafish")
#setwd("/zp1/data/fwang/model")
source("G://data//model//script//main.R")

#load("zfAdzfNMDA_zfLD_zfTR_pbmcSubC_CondSub_MG_umap_seurat.RData")
# Generate PB matrix, transcription factors binding strength in each gene, row is gene, col is factor
load("zfATACldtP00R1_zfATACnmdahPldtP_PeaksFootprints_Fimo_ScBk_Candidates_cutsp2_TFsTargets_splited_motifs3.RData")
factors_ID=names(table(FOS[,2]))
genes_ID=names(table(FOS[,3]))

load("PB.RData")

rownames(PB)=genes_ID
colnames(PB)=factors_ID

#remove multi genes
genes_ID_new=c()
mark=1
for(i in 1:length(genes_ID))
{
	gene_name=unlist(strsplit(genes_ID[i],","))
	if(length(gene_name)==1)
	{
		genes_ID_new[mark]=genes_ID[i]
		mark=mark+1
	}
}
genes_ID=genes_ID_new
#PB=Generate_PB_matrix(FOS, factors_ID, genes_ID)
PB=PB[genes_ID,]

# Generate factors' and target genes' fitted smooth expression data matrix
Pse_Bin_num=100
ByBin='Equal.Pseudotime'
condition="zfNMDA"
#Pseudotime1=MG@meta.data
load("Pseudotime1.RData")
Pseudotime1=Pseudotime1[Pseudotime1$Treatment %in% c("zfAd",condition),]
Pseudotime=Pseudotime1[order(Pseudotime1[,"Pseudotime"]),"Pseudotime"]

#factors_exp=Generate_exp_matrix(MG, factors_ID, condition)
load("factors_exp.RData")
#diff_factors=Differential_genes(factors_exp, Pseudotime1)
#factors_exp=factors_exp[diff_factors$gene_ID,]
#factors_smooth_exp=Smooth_exp_by_bin(factors_exp, Pseudotime1, NULL, Pse_Bin_num, ByBin)
factors_smooth_fitted_exp=Extract_fitting_value(factors_exp, Pseudotime, fitting_method = "spline", Pse_Bin_num)
#factors_smooth_fitted_exp=factors_smooth_fitted_exp[,1:(ncol(factors_smooth_fitted_exp)-1)]

#genes_exp=Generate_exp_matrix(MG, genes_ID, condition)
load("genes_exp.RData")
#diff_genes=Differential_genes(genes_exp, Pseudotime1)
#genes_smooth_exp=Smooth_exp_by_bin(genes_exp, Pseudotime1, NULL, Pse_Bin_num, ByBin)
genes_smooth_fitted_exp=Extract_fitting_value(genes_exp, Pseudotime, fitting_method = "spline", Pse_Bin_num)
genes_smooth_fitted_ratio=Change_ration_pse(genes_exp, Pseudotime, fitting_method = "spline", Pse_Bin_num)
#gene_exp_pre=genes_smooth_fitted_exp[,1:(ncol(genes_smooth_fitted_exp)-1)]

### Generate binding matrix between target genes and factors ###
PB=PB[genes_ID, rownames(factors_exp)]

# Generate all parameters
marker=data.frame(index1=c(1:nrow(factors_smooth_fitted_exp)))
weight_parm_all=apply(marker,2,function(x1){paste("w",x1,sep="")})[,1] ## Generate all parameters and symbols
factor_label_all=apply(marker,2,function(x1){paste("factor",x1,sep="")})[,1] ## Generate all parameters and symbols


cor_method="spearman"
cor_threshold=0.5
direction_mark=matrix(0,nrow(genes_smooth_fitted_exp),2)
final_regulatory_formula=list()
for(index in 1:nrow(genes_smooth_fitted_exp))
{
	##### Generate fitting data for estimate parameters #####
	factor_ind=which(PB[index,]!=0) ## index of candidate binding factors for target gene

	target_exp=genes_smooth_fitted_exp[index,]
	factor_exp=factors_smooth_fitted_exp[factor_ind,]

	dG=genes_smooth_fitted_ratio[index,] ## expression ratio of target gene
	dG_pre=genes_smooth_fitted_exp[index,] ## previous time expression of target gene
	PB_level=PB[index,factor_ind] ## binding FOS value between factors and target genes
	TF_binding_exp=factors_smooth_fitted_exp[factor_ind,] ## expression matrix for corresponding factors
	factor_label=factor_label_all[factor_ind] ## label for corresponding factors
	weight_parm=weight_parm_all[factor_ind] ## parameters for corresponding factors

	if(length(factor_label)==1)
	{
		data=as.data.frame(cbind(dG,TF_binding_exp,dG_pre))
	}
	if(length(factor_label)>1)
	{
		data=as.data.frame(cbind(dG,t(TF_binding_exp),dG_pre))
	}
	colnames(data)=c("gene",factor_label,"gene_pre")

	#initial_factor_direction=Initial_factor_direction(target_exp, factor_exp, cor_method="spearman", cor_threshold=0.5)
	#initial_factor_direction[,1]=rep(1,length(factor_label))
	#factors_direction=Estimate_regulation_diretion(factor_label, initial_factor_direction, data)

	factors_direction = Binary_GA_direction(factor_label, data)

	direction_mark[index,1]=round(factors_direction[length(factors_direction)-1])
	direction_mark[index,2]=factors_direction[length(factors_direction)]
	factors_direction=round(factors_direction[1:(length(factors_direction)-2)])
	
	final_regulatory_formula=Regulatory_formula(factor_label, weight_parm, PB_level, factors_direction)

	final_regulatory_parameters=Fitting_parameters(data, final_regulatory_formula)
}






index=4325
ENSDARG00000087554
plot(genes_smooth_fitted_exp[index,1:99],factors_smooth_fitted_exp[242,])






