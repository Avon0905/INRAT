##### Calculate new FOS values of footprints #####
Cal_Footprints_FOS <- function(FP1, FlankFold1=3)
{
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

##### Generate factors and target genes binding information #####
Generate_PB_matrix = function(FOS, factors, genes)
{
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
	rownames(PB)=genes
	colnames(PB)=factors
	return(PB)
	#save(PB,file="PB.RData")
}

##### Generate genes expression matrix #####
Generate_exp_matrix = function(MG, gene_list)
{
	cells_inf=MG@meta.data[order(MG@meta.data[,"Pseudotime"]),]
	cells_pse=cells_inf[,"Pseudotime"]
	expression_matrix=MG@data[,rownames(cells_inf)]
	genes_exp=expression_matrix[gene_list,]
	genes_exp=as.matrix(genes_exp)
	colnames(genes_exp)=rownames(cells_inf)
	rownames(genes_exp)=gene_list
	return(genes_exp)
	#save(genes_exp,file="zf_genes_exp.RData")
}

##### Cut expression value into pseudotime bins #####
Smooth_exp_by_bin = function(genes_exp, Pseudotime1, PseudotimeRange1=NULL, SmoothLength1=100, ByBin1=c('Equal.Pseudotime','Equal.Cells'))
{
	if(ncol(genes_exp)!=nrow(Pseudotime1)){
       stop('The length of pseudotime is not equal to the number of cells')
    }else if(!is.element('Pseudotime', colnames(Pseudotime1))){
       stop('No pseudotime inforamtion in variable Pseudotime1')
    }
	if(ByBin1[1]=='Equal.Pseudotime'){
      if(is.null(PseudotimeRange1)){
        PseudotimeBin1 <- seq(min(Pseudotime1$Pseudotime), max(Pseudotime1$Pseudotime), length.out = SmoothLength1+1)
      }else{ PseudotimeBin1 <- seq(PseudotimeRange1[1], PseudotimeRange1[2], length.out = SmoothLength1+1) }
    }else{
      Pseudotime1 <- Pseudotime1[order(Pseudotime1$Pseudotime), ]
      Bin1 <- ceiling(nrow(Pseudotime1)/SmoothLength1)
      PseudotimeBin1 <- seq(1, nrow(Pseudotime1), by = Bin1)
      PseudotimeBin1[length(PseudotimeBin1)+1] <- nrow(Pseudotime1)
    }
    Smooth_Exp <- array(0, dim=c(nrow(genes_exp), SmoothLength1))
    for(i in 1:(length(PseudotimeBin1)-1)){
      if(ByBin1[1]=='Equal.Pseudotime'){
        if(i==length(PseudotimeBin1)-1){
          Cells1 <- Pseudotime1[Pseudotime1$Pseudotime>=PseudotimeBin1[i] & Pseudotime1$Pseudotime<=PseudotimeBin1[i+1], ]
        }else{ Cells1 <- Pseudotime1[Pseudotime1$Pseudotime>=PseudotimeBin1[i] & Pseudotime1$Pseudotime<PseudotimeBin1[i+1], ] }
      }else{
        if(i==length(PseudotimeBin1)-1){
          Cells1 <- Pseudotime1[PseudotimeBin1[i]:PseudotimeBin1[i+1], ]
        }else{ Cells1 <- Pseudotime1[PseudotimeBin1[i]:(PseudotimeBin1[i+1]-1), ] }
      }
      
      if(nrow(Cells1)>1){
        Smooth_Exp[, i] <- rowMeans(genes_exp[, match(rownames(Cells1), colnames(genes_exp))])
      }else if(nrow(Cells1)==1){
        Smooth_Exp[, i] <- genes_exp[, match(rownames(Cells1), colnames(genes_exp))]
      }
    }
    rownames(Smooth_Exp) <- rownames(genes_exp)
    return(Smooth_Exp)
}

##### Extract values of fitting line of genes' smooth bin expression #####
Extract_fitting_value = function(Smooth_Exp, fitting_method="loess", Pse_Bin)
{
	if(fitting_method=="loess")
	{
		Smooth_Exp_Fitting=t(apply(Smooth_Exp, 1, function(x1){
			data1=as.data.frame(cbind(c(1:Pse_Bin),x1))
			colnames(data1)=c("Pse_Bin","Smooth_Exp")
			smooth_vals = predict(loess(Smooth_Exp~Pse_Bin,data1), data1$Pse_Bin)
			return(smooth_vals)
			}))
		return(Smooth_Exp_Fitting)
	}
}

##### Get change ratio for target genes in each psedotime_bin point #####
Change_ration_pse = function(Exp)
{
	change_ration_exp=t(apply(Exp, 1, function(x1){
		return(x1[2:length(x1)]-x1[1:(length(x1)-1)])	
		})) #ignore the first psedotime_bin point
	return(change_ration_exp)
}

##### Construct regulatory parameters and formula #####
index=data.frame(index1=c(1:nrow(factor_exp_bin)))
weight_parm_all=apply(index,2,function(x1){paste("w",x1,sep="")})[,1]
factor_lable_all=apply(index,2,function(x1){paste("factor",x1,sep="")})[,1]
Regulatory_formula = function(PB, pos, factor_lable_all, weight_parm_all, factors_direction)
{
	factor_pos=which(PB[pos,]!=0)
	factor_lable=factor_lable_all[factor_pos]	
	weight_parm=weight_parm_all[factor_pos]
	PB_level=PB[pos,factor_pos]
	#construct regulatory formula for each gene
	if(length(factor_lable)>0)
	{
		formula_left="gene"
		formula_right1=c()
		positive_index=which(factors_direction==1)
		for(j in positive_index)
		{
			formula_right1[j]=paste(PB_level[j],weight_parm[j],factor_lable[j],sep="*")
		}
		formula_right2=paste(formula_right1[positive_index],collapse="+")
		formula_right3=c()
		for(j in 1:length(factor_pos))
		{
			formula_right3[j]=paste(PB_level[j],weight_parm[j],factor_lable[j],sep="*")
		}
		formula_right4=paste(formula_right3,collapse="+")
		formula_right=paste("alpha*(",formula_right2,")/(1+",formula_right4,")+belta*gene_pre+e",sep="")
		regulatory_formula=as.formula(paste(formula_left,formula_right,sep=" ~ "))
		return(regulatory_formula)
	}else{
		return(NA)
	}
}

##### Fitting parameters for each model #####
nls_mod=list()
Fitting_parameters = function(factor_exp_bin, gene_exp_bin, gene_exp_pre, PB, pos, factor_lable_all, regulatory_formula, start)
{
	dG=gene_exp_bin[pos,]
	dG_pre=gene_exp_pre[pos,]
	factor_pos=which(PB[pos,]!=0)
	TF_binding_exp=factor_exp_bin[factor_pos,]
	factor_lable=factor_lable_all[factor_pos]	
	if(length(factor_lable)==1)
	{
		data=as.data.frame(cbind(dG,TF_binding_exp,dG_pre))
	}
	if(length(factor_lable)>1)
	{
		data=as.data.frame(cbind(dG,t(TF_binding_exp),dG_pre))
	}
	colnames(data)=c("gene",factor_lable,"gene_pre")
	start=c(1,rep(1,length(factor_lable)),1,0)
	start=as.list(start)
	names(start)=c("alpha",weight_parm,"belta","e")
	nls_mod=nls(regulatory_formula,data=data,control=nls.control(warnOnly = T),algorithm="port",start=start,lower=c(-Inf,rep(0,length(factor_lable)),-Inf,-Inf),upper=c(Inf,rep(1,length(factor_lable)),Inf,Inf))
	return(nls_mod)
	#optim(c(rep(0,14),-1,0),fn=regulatory_formula,lower=rep(-1,16),upper=rep(1,16),data=data, method = "L-BFGS-B") 
}



#source("E://model//regulatory_model.R")
#用旧的模型筛选factor








