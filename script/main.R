library(monocle)
library(Seurat)
##### Calculate new FOS values of footprints #####
Cal_Footprints_FOS = function(FP1, FlankFold1=3)
{
  	FP2 = apply(FP1, 1, function(X1){ MotifSize1 = as.numeric(X1[3])
	WidthL1 = floor(MotifSize1/2); WidthR1 = MotifSize1-WidthL1; 
	WidthL2 = floor(MotifSize1*(2*FlankFold1+1)/2); WidthR2 = MotifSize1*(2*FlankFold1+1)-WidthL2;
	Insertion1 = strsplit(X1[7], ',')[[1]]; Midpoint1 = floor(length(Insertion1)/2)
	Insertion2 = as.numeric(Insertion1[c((Midpoint1-WidthL2+1):(Midpoint1+WidthR2))])

	L1 = sum(Insertion2[1:(WidthL2-WidthL1)])/FlankFold1
	M1 = sum(Insertion2[(WidthL2-WidthL1+1):(WidthL2+WidthR1)])
	R1 = sum(Insertion2[(WidthL2+WidthR1+1):length(Insertion2)])/FlankFold1
	FOS1 = min(1-(M1+1)/(L1+1), 1-(M1+1)/(R1+1));
	if(FOS1<0)
	{
		FOS1=0
	}
	return(c(X1[1:6], as.numeric(FOS1)))
	} )
	return(t(FP2))  
}

##### Generate factors and target genes binding information #####
Generate_PB_matrix = function(FOS, factors_ID, genes_ID)
{
	PB=array(0,dim=c(length(genes_ID),length(factors_ID)))
	for(i in 1:length(genes_ID))
	{
		for(j in 1:length(factors_ID))
		{
			FOS2=FOS[FOS[,2]==factors_ID[j] & FOS[,3]==genes_ID[i],]
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
		#print(i)
	}
	rownames(PB)=genes_ID
	colnames(PB)=factors_ID
	PB=round(PB,2)
	return(PB)
	#save(PB,file="PB.RData")
}

##### Generate genes expression matrix #####
Generate_exp_matrix = function(MG, gene_ID, condition="zfNMDA")
{
	cells_inf=MG@meta.data[order(MG@meta.data[,"Pseudotime"]),]
	cells_inf=cells_inf[cells_inf$Treatment %in% c("zfAd",condition),]
	cells_pse=cells_inf[,"Pseudotime"]
	genes_exp=MG@raw.data[gene_ID,rownames(cells_inf)]
	genes_exp=as.matrix(genes_exp)
	colnames(genes_exp)=rownames(cells_inf)
	rownames(genes_exp)=gene_ID
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
        PseudotimeBin1 = seq(min(Pseudotime1$Pseudotime), max(Pseudotime1$Pseudotime), length.out = SmoothLength1+1)
      }else{ PseudotimeBin1 = seq(PseudotimeRange1[1], PseudotimeRange1[2], length.out = SmoothLength1+1) }
    }else{
      Pseudotime1 = Pseudotime1[order(Pseudotime1$Pseudotime), ]
      Bin1 = ceiling(nrow(Pseudotime1)/SmoothLength1)
      PseudotimeBin1 = seq(1, nrow(Pseudotime1), by = Bin1)
      PseudotimeBin1[length(PseudotimeBin1)+1] = nrow(Pseudotime1)
    }
    Smooth_Exp = array(0, dim=c(nrow(genes_exp), SmoothLength1))
    for(i in 1:(length(PseudotimeBin1)-1)){
      if(ByBin1[1]=='Equal.Pseudotime'){
        if(i==length(PseudotimeBin1)-1){
          Cells1 = Pseudotime1[Pseudotime1$Pseudotime>=PseudotimeBin1[i] & Pseudotime1$Pseudotime<=PseudotimeBin1[i+1], ]
        }else{ Cells1 = Pseudotime1[Pseudotime1$Pseudotime>=PseudotimeBin1[i] & Pseudotime1$Pseudotime<PseudotimeBin1[i+1], ] }
      }else{
        if(i==length(PseudotimeBin1)-1){
          Cells1 = Pseudotime1[PseudotimeBin1[i]:PseudotimeBin1[i+1], ]
        }else{ Cells1 = Pseudotime1[PseudotimeBin1[i]:(PseudotimeBin1[i+1]-1), ] }
      }
      
      if(nrow(Cells1)>1){
        Smooth_Exp[, i] = rowMeans(genes_exp[, match(rownames(Cells1), colnames(genes_exp))])
      }else if(nrow(Cells1)==1){
        Smooth_Exp[, i] = genes_exp[, match(rownames(Cells1), colnames(genes_exp))]
      }
    }
    rownames(Smooth_Exp) = rownames(genes_exp)
    return(Smooth_Exp)
}

##### Extract values of fitting line of genes' smooth bin expression #####
Extract_fitting_value = function(exp_data, Pseudotime, fitting_method=c("loess", "spline"), Pse_Bin_num)
{
	if(fitting_method == "spline")
	{
		##### Use spline to simulate single cell gene expression data along pseudotime #####
		Smooth_Exp_Fitting = t(apply(exp_data, 1, function(x1){
			ss10 = smooth.spline(Pseudotime, log(x1+1), df = 10)
			return(predict(ss10, x = Pseudotime)$y)
		}))
		return(Smooth_Exp_Fitting)
		
		#plot(Pseudotime, log(genes_exp[index,]+1), pch = 20, col = "steelblue")
		#lines(ss10, lty = 2, col = "red")

		#t = 0.0001
		#t0 = 3
		#velocity = (predict(ss10,x=(t0+t))$y-predict(ss10,x=t0)$y)/t
	}

	if(fitting_method=="loess")
	{
		Smooth_Exp_Fitting=t(apply(exp_data, 1, function(x1){
			data1=as.data.frame(cbind(c(1:Pse_Bin_num),x1))
			colnames(data1)=c("Pse_Bin_num","Smooth_Exp")
			smooth_vals = predict(loess(exp_data~Pse_Bin_num,data1), data1$Pse_Bin_num)
			return(smooth_vals)
			}))
		return(Smooth_Exp_Fitting)
	}
}

##### Get change ratio for target genes in each psedotime_bin point #####
Change_ration_pse = function(exp_data, Pseudotime, fitting_method=c("loess", "spline"), Pse_Bin_num)
{
	if(fitting_method == "spline")
	{
		##### Use spline to simulate single cell gene expression data along pseudotime #####
		change_ratio = t(apply(exp_data, 1, function(x1){
			ss10 = smooth.spline(Pseudotime, log(x1+1), df = 10)
			delt_t = 0.0001
			velocity = (predict(ss10,x=(Pseudotime+delt_t))$y-predict(ss10,x=Pseudotime)$y)/delt_t
			return(velocity)
		}))
		return(change_ratio)
		
		#plot(Pseudotime, log(genes_exp[index,]+1), pch = 20, col = "steelblue")
		#lines(ss10, lty = 2, col = "red")

		#t = 0.0001
		#t0 = 3
		#velocity = (predict(ss10,x=(t0+t))$y-predict(ss10,x=t0)$y)/t
	}
}
# Change_ration_pse = function(Exp)
# {
# 	change_ration_exp=t(apply(Exp, 1, function(x1){
# 		return(x1[2:length(x1)]-x1[1:(length(x1)-1)])	
# 		})) #ignore the first psedotime_bin point
# 	return(change_ration_exp)
# }

##### Generate initial regulation directions of factors according to correlation
Initial_factor_direction = function(target_exp, factor_exp, cor_method=c("pearson","spearman"), cor_threshold=0.5)
{
	direction_parm=rep(0,nrow(factor_exp))
	direction_mark=rep(0,nrow(factor_exp))
	cor_val=apply(factor_exp,1,function(x1){return(cor(x1,target_exp,method=cor_method))})
	for(j in 1:length(cor_val))
	{
		if(cor_val[j]>0)
		{
			direction_parm[j]=1
		}
		if(cor_val[j]<0)
		{
			direction_parm[j]=0
		}
		if(abs(cor_val[j])>cor_threshold)
		{
			direction_mark[j]=1
		}
	}
	factor_direction=cbind(direction_parm,direction_mark)
	return(factor_direction)
}

##### Construct random regulatory parameters and formula for predicting regulation direction #####
Regulatory_formula2 = function(factor_label, factor_direction)
{
	#construct regulatory formula for each gene
	if(length(factor_label)>0)
	{
		formula_left="gene"
		formula_right1=c()
		for(j in 1:length(factor_label))
		{
			formula_right1[j]=paste(factor_label[j],factor_direction[j],sep="^")
		}
		formula_right2=paste(formula_right1,collapse="+")
		formula_right4=paste(factor_label,collapse="+")
		formula_right=paste("alpha*(",formula_right2,")/(1+",formula_right4,")+belta*gene_pre+e",sep="")
		regulatory_formula=as.formula(paste(formula_left,formula_right,sep=" ~ "))
		return(regulatory_formula)
	}else{
		return(NA)
	}
}

##### Nonlinear fitting to estimate regulation direction of factors #####
Fitting_parameters2 = function(fitting_data, regulatory_formula)
{
	start=c(1,1,0)
	start=as.list(start)
	names(start)=c("alpha","belta","e")
	nls_mod=nls(regulatory_formula,data=fitting_data,control=nls.control(warnOnly = T),start=start)
	return(nls_mod)
}

##### Generate new leaf node, change only one (won't change nodes with high correlation) ##### 
Leaf_node = function(factor_direction)
{
	parent_node=factor_direction[factor_direction[,2]==0,1]
	leaf_direction=array(factor_direction[,1],dim=c(nrow(factor_direction),length(parent_node)))
	for(i in 1:length(parent_node))
	{
		leaf_node=parent_node
		leaf_node[i]=abs(parent_node[i]-1)
		leaf_direction[factor_direction[,2]==0,i]=leaf_node
	}
	return(leaf_direction)
}

##### Generate new leaf node, change only one (will change all nodes) ##### 
Leaf_node2 = function(factor_direction)
{
	parent_node=factor_direction[,1]
	leaf_direction=array(factor_direction[,1],dim=c(nrow(factor_direction),length(parent_node)))
	for(i in 1:length(parent_node))
	{
		leaf_node=parent_node
		leaf_node[i]=abs(parent_node[i]-1)
		leaf_direction[,i]=leaf_node
	}
	return(leaf_direction)
}

##### Calculate coefficient of determination to evaluate the model #####
Cal_R_square = function(y_true, y_predict, n, p)
{
	y_mean=mean(y_true)
	SS_tot=sum((y_true-y_mean)^2)
	SS_res=sum((y_true-y_predict)^2)
	R2=(1-SS_res/SS_tot)
	R2_adjusted=1-(1-R2)*(n-1)/(n-p-1)
	return(R2_adjusted)
}

##### Iterated fitting parameters and select the best regulation direction combination of factors #####
Estimate_regulation_diretion = function(factor_label, initial_factor_direction, fitting_data)
{
	initial_regulatory_formula=Regulatory_formula2(factor_label, initial_factor_direction)
	initial_nls_mod=Fitting_parameters2(fitting_data, initial_regulatory_formula)
	n=nrow(fitting_data)
	p=length(factor_label)
	mark=Cal_R_square(fitting_data$gene,predict(initial_nls_mod),n,p)
	leaf_direction=Leaf_node2(initial_factor_direction)
	final_factor_direction=initial_factor_direction[,1]
	record=0
	while(TRUE)
	{
		record=record+1
		leaf_R2_adj=apply(leaf_direction,2,function(x){
			leaf_factor_direction=cbind(x,initial_factor_direction[,2])
			leaf_regulatory_formula=Regulatory_formula2(factor_label, leaf_factor_direction)
			leaf_nls_mod=Fitting_parameters2(fitting_data, leaf_regulatory_formula)
			return(Cal_R_square(fitting_data$gene,predict(leaf_nls_mod),n,p))
		})
		max_index=which.max(leaf_R2_adj)
		if(leaf_R2_adj[max_index]>mark)
		{
			mark=leaf_R2_adj[max_index]
			final_factor_direction=leaf_direction[,max_index]
			leaf_direction=Leaf_node2(cbind(leaf_direction[,max_index],initial_factor_direction[,2]))
		}else{
			break
		}
	}
	print(paste0("Leaf loop number: ",record))
	print(paste0("Adjusted R2: ",round(mark,3)))
	return(c(final_factor_direction,record,round(mark,3)))
}

##### Binary GA algorithm to estimate regulation direction #####
Binary_GA_direction = function(factor_label, fitting_data)
{
	library(genalg)
	### Construct evaluate function for GA algorithm
	Eval_func = function(x)
	{
		n=nrow(fitting_data)
		p=length(factor_label)
		regulatory_formula=Regulatory_formula2(factor_label, x)
		nls_mod=Fitting_parameters2(fitting_data, regulatory_formula)
		return(Cal_R_square(fitting_data$gene,predict(nls_mod),n,p))
	}
	Size = length(factor_label) #the number of genes in the chromosome (The number of candidate factors)

	GAmodel = rbga.bin(size = Size, mutationChance = 0.01, evalFunc = Eval_func)
	result = summary(GAmodel)
	result = unlist(strsplit(unlist(strsplit(result, split="Best Solution : "))[2], split=" \n"))
	result = as.numeric(unlist(strsplit(result, split=" ")))
	return(result)
}

##### Construct regulatory parameters and formula #####
Regulatory_formula = function(factor_label, weight_parm, PB_level, factors_direction)
{
	#construct regulatory formula for each gene
	if(length(factor_label)>0)
	{
		formula_left="dG"
		formula_right1=c()
		positive_index=which(factors_direction==1)
		for(j in positive_index)
		{
			formula_right1[j]=paste(PB_level[j],weight_parm[j],factor_label[j],sep="*")
		}
		formula_right2=paste(formula_right1[positive_index],collapse="+")
		formula_right3=c()
		for(j in 1:length(factor_label))
		{
			formula_right3[j]=paste(PB_level[j],weight_parm[j],factor_label[j],sep="*")
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
Fitting_parameters = function(fitting_data, regulatory_formula)
{
	parameter_num=ncol(fitting_data)-2
	start=c(1,rep(1,parameter_num),1,0)
	start=as.list(start)
	names(start)=c("alpha",weight_parm,"belta","e")
	nls_mod=nls(regulatory_formula,data=fitting_data,control=nls.control(warnOnly = T),algorithm="port",start=start,lower=c(-Inf,rep(0,parameter_num),-Inf,-Inf),upper=c(Inf,rep(1,parameter_num),Inf,Inf))
	return(nls_mod)
	#optim(c(rep(0,14),-1,0),fn=regulatory_formula,lower=rep(-1,16),upper=rep(1,16),data=data, method = "L-BFGS-B") 
}

##### Identify differential factors or genes #####
Differential_genes = function(exp_matrix, Pseudotime1)
{
	library(monocle)

	pheno.data.df <- Pseudotime1[colnames(exp_matrix),] # Must be data frame object
	pd <- new('AnnotatedDataFrame', data = pheno.data.df) 
	rownames(pd) <- colnames(exp_matrix)

	gene_annotation <- array(0,dim = c(nrow(exp_matrix),2))
	gene_annotation[,1] <- rownames(exp_matrix)
	gene_annotation[,2] <- rownames(exp_matrix)
	gene_annotation <- as.data.frame(gene_annotation)
	colnames(gene_annotation) <- c("gene_ID", "gene_short_name")
	fd <- new("AnnotatedDataFrame", data = gene_annotation)
	rownames(fd) <- rownames(exp_matrix)

	mef <- newCellDataSet(as.matrix(exp_matrix), phenoData = pd, featureData = fd, expressionFamily = VGAM::negbinomial.size())
	mef <- estimateSizeFactors(mef)
	mef <- estimateDispersions(mef)
	mef <- detectGenes(mef, min_expr = 0.1)

	diff_test_res <- differentialGeneTest(mef[rownames(exp_matrix),], fullModelFormulaStr = "~sm.ns(Pseudotime)")
	diff_genes <- subset(diff_test_res,num_cells_expressed >= round(0.01*nrow(pData(mef))))
	diff_genes <- diff_genes[which(diff_genes$qval<0.001),]
}