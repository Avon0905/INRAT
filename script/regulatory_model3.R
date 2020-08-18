##### Use correlation to predict the regulation direction of factors #####
Initial_factor_direction = function(genes_smooth_exp, factors_smooth_exp, PB, pos, cor_method=c("pearson","spearman"), cor_threshold=0.5)
{
	target_exp=genes_smooth_exp[pos,]
	factor_pos=which(PB[pos,]!=0)
	factor_exp=factors_smooth_exp[factor_pos,]

	direction_parm=rep(0,length(factor_pos))
	direction_mark=rep(0,length(factor_pos))

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
Regulatory_formula2 = function(PB, pos, factor_lable_all, factor_direction)
{
	factor_pos=which(PB[pos,]!=0)
	factor_lable=factor_lable_all[factor_pos]
	#construct regulatory formula for each gene
	if(length(factor_lable)>0)
	{
		formula_left="gene"
		formula_right1=c()
		for(j in 1:length(factor_pos))
		{
			formula_right1[j]=paste(factor_lable[j],factor_direction[j,1],sep="*")
		}
		formula_right2=paste(formula_right1,collapse="+")
		formula_right4=paste(factor_lable,collapse="+")
		formula_right=paste("alpha*(",formula_right2,")/(1+",formula_right4,")+belta*gene_pre+e",sep="")
		regulatory_formula=as.formula(paste(formula_left,formula_right,sep=" ~ "))
		return(regulatory_formula)
	}else{
		return(NA)
	}
}

##### Nonlinear fitting to estimate regulation direction of factors #####
Fitting_parameters2 = function(factor_exp_bin, gene_exp_bin_dt, gene_exp_pre, PB, pos, factor_lable_all, regulatory_formula)
{
	dG=gene_exp_bin_dt[pos,]
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
	start=c(1,1,0)
	start=as.list(start)
	names(start)=c("alpha","belta","e")

	nls_mod=nls(regulatory_formula,data=data,control=nls.control(warnOnly = T),start=start)
	return(nls_mod)
}

##### Generate new leaf node, change only one ##### 
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

##### Calculate coefficient of determination to evaluate the model #####
Cal_R_square=function(y_true,y_predict,n,p)
{
	y_mean=mean(y_true)
	SS_tot=sum((y_true-y_mean)^2)
	SS_res=sum((y_true-y_predict)^2)
	R2=(1-SS_res/SS_tot)
	R2_adjusted=1-(1-R2)*(n-1)/(n-p-1)
	return(R2_adjusted)
}

##### Iterated fitting parameters and select the best regulation direction combination of factors #####
Estimate_regulation_diretion = function(factor_exp_bin, gene_exp_bin_dt, gene_exp_pre, PB, pos, factor_lable_all, initial_factor_direction)
{
	regulatory_formula=Regulatory_formula2(PB, pos, factor_lable_all, initial_factor_direction)
	nls_mod=Fitting_parameters2(factor_exp_bin, gene_exp_bin_dt, gene_exp_pre, PB, pos, factor_lable_all, regulatory_formula)
	n=ncol(gene_exp_bin_dt)
	p=length(which(PB[pos,]!=0))
	mark=Cal_R_square(gene_exp_bin_dt[pos,],predict(nls_mod),n,p)
	leaf_direction=Leaf_node(initial_factor_direction)
	final_factor_direction=initial_factor_direction[,1]
	record=0
	while(TRUE)
	{
		record=record+1
		leaf_R2_adj=c(1:ncol(leaf_direction))
		for(i in 1:ncol(leaf_direction))
		{
			leaf_factor_direction=cbind(leaf_direction[,i],initial_factor_direction[,2])
			leaf_regulatory_formula=Regulatory_formula2(PB, pos, factor_lable_all, leaf_factor_direction)
			leaf_nls_mod=Fitting_parameters2(factor_exp_bin, gene_exp_bin_dt, gene_exp_pre, PB, pos, factor_lable_all, leaf_regulatory_formula)
			leaf_R2_adj[i]=Cal_R_square(gene_exp_bin_dt[pos,],predict(leaf_nls_mod),n,p)
		}
		index=which.max(leaf_R2_adj)
		if(leaf_R2_adj[index]>mark)
		{
			mark=leaf_R2_adj[index]
			final_factor_direction=leaf_direction[,index]
			leaf_direction=Leaf_node(cbind(leaf_direction[,index],initial_factor_direction[,2]))
		}else{
			break
		}
	}
	print(paste0(rownames(gene_exp_bin_dt)[pos]," Done!"))
	print(paste0("Leaf loop number: ",record))
	print(paste0("Adjusted R2: ",round(mark,3)))
	return(c(final_factor_direction,record,round(mark,3)))
}




