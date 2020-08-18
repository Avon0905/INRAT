##### Get gene regulation pairs from simulated networks #####
Get_gene_pair = function(interaction_CID)
{
	grn = data.frame(target = 0, factor = 0, parm = 0, type = 0)
	index = 1
	for(i in 1:nrow(interaction_CID))
	{
		interaction_temp = as.numeric(unlist(strsplit(as.character(interaction_CID[i,1]), ",")))
		for(j in 1:interaction_temp[2])
		{
			grn[index,"target"] = round(interaction_temp[1])
			grn[index,"factor"] = round(interaction_temp[j+2])
			grn[index,"parm"] = interaction_temp[j+2+interaction_temp[2]]
			grn[index,"type"] = round(interaction_temp[j+2+interaction_temp[2]]/abs(interaction_temp[j+2+interaction_temp[2]]))
			index = index+1
		}
	}
	return(grn)
}

##### Construct monocle2 trajectory #####
Construct_monocle2_trajectory = function(exp_matrix, cell_type = c(1:7), cells_per_type = 300)
{
	library(monocle)

	pd = as.matrix(cbind(colnames(exp_matrix),colnames(exp_matrix)))
	cell_type_all = c()
	for(i in 1:length(cell_type))
	{
		cell_type_all = c(cell_type_all,rep(cell_type[i],cells_per_type))
	}
	pd = cbind(pd,cell_type_all)
	colnames(pd) = c("cell_id","cell_name","cell_type")
	pheno.data.df = as.data.frame(pd)
	pd <- new('AnnotatedDataFrame', data = pheno.data.df) 
	rownames(pd)=colnames(exp_matrix)


	fd = as.matrix(cbind(rownames(exp_matrix),rownames(exp_matrix)))
	colnames(fd) = c("gene_ID","gene_short_name")
	gene_annotation = as.data.frame(fd)
	fd <- new("AnnotatedDataFrame", data = gene_annotation)
	rownames(fd)=rownames(exp_matrix)

	mef <- newCellDataSet(as.matrix(exp_matrix), phenoData = pd, featureData = fd)
	mef <- estimateSizeFactors(mef)
	#mef <- estimateDispersions(mef)
	mef <- detectGenes(mef, min_expr = 0.1)

	diff_test_res <- differentialGeneTest(mef, fullModelFormulaStr = "~cell_type")
	ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

	mef <- setOrderingFilter(mef, ordering_genes)
	mef <- reduceDimension(mef, max_components = 2, reduction_method = "DDRTree")
	mef <- orderCells(mef)
	#save(mef,file=paste0(dir,species,"/",condition,"/",cell_type,"/",species,"_",condition,"_",cell_type,"_",find_var_genes_method,"_monocle2.RData"))

	return(mef)
}

##### Change single cell data to bulk data #####
SC_to_bulk = function(exp_matrix, number_bins, number_sc)
{
	bulk_exp_matrix = matrix(0, nrow=nrow(exp_matrix), ncol=number_bins)
	for(i in 1:number_bins)
	{
		bulk_exp_matrix[,i] = apply(exp_matrix, 1, function(x)
		{
			mean(x[(i*100-99):(i*100)])
		})
	}
	rownames(bulk_exp_matrix) = rownames(exp_matrix)
	colnames(bulk_exp_matrix) = c(1:number_bins)
	return(bulk_exp_matrix)
}

##### Use correlation to identify 
Regulation_prediction = function(exp_matrix, targets, factors, method = c("cor", "GENIE3", "SCODE"), cor_thre = 0.6, genie3_thre = 200, dir = sim_dir)
{
	if(method == "cor")
	{
		grn_cor = data.frame(target = 0, factor = 0, type = 0)
		temp = 1
		for(i in 1:length(targets))
		{
			for(j in 1:length(factors))
			{
				correlation = cor(exp_matrix[targets[i],],exp_matrix[factors[j],])
				if(correlation > thre)
				{
					grn_cor[temp,1] = targets[i]
					grn_cor[temp,2] = factors[j]
					grn_cor[temp,3] = "1"
					temp = temp+1
				}
				if(correlation < (-thre))
				{
					grn_cor[temp,1] = targets[i]
					grn_cor[temp,2] = factors[j]
					grn_cor[temp,3] = "-1"
					temp = temp+1
				}
			}
		}
		grn_predict = grn_cor
	}
	if(method == "GENIE3")
	{
		library(GENIE3)
		set.seed(123)
		exprMatr = SC_to_bulk(exp_matrix, number_bins, number_sc)
		weightMat = GENIE3(exprMatr, regulators=factors, targets=targets)
		linkList = getLinkList(weightMat, reportMax=200)
		grn_predict = cbind(linkList$targetGene, linkList$regulatoryGene, linkList$weight)
		return(grn_predict)
	}
	if(method == "SCODE")
	{
		source("G://data//model//SCODE//SCODE-master//SCODE.R")
		out_dir = paste0(dir, "SCODE//")
		SCODE(exp_matrix, pData(sim_mef)[,"Pseudotime"], out_dir, tfnum = nrow(exp_matrix), pnum = 4, cnum = ncol(exp_matrix), maxite = 100)
		A = read.table(paste0(out_dir,"A.txt"), sep='\t')
		rownames(A) = rownames(exp_matrix)
		colnames(A) = rownames(exp_matrix)
		grn_predict = data.frame(target = 0, factor = 0, parm = 0, type = 0)
		mark = 1
		for(i in 1:length(targets))
		{
			for(j in 1:length(factors))
			{
				if(A[targets[i],factors[j]] > 0.5)
				{
					grn_predict[mark,"target"] = targets[i]
					grn_predict[mark,"factor"] = factors[j]
					grn_predict[mark,"parm"] = A[targets[i],factors[j]]
					grn_predict[mark,"type"] = "1"
					mark = mark+1
				}
				if(A[targets[i],factors[j]] < (-0.5))
				{
					grn_predict[mark,"target"] = targets[i]
					grn_predict[mark,"factor"] = factors[j]
					grn_predict[mark,"parm"] = A[targets[i],factors[j]]
					grn_predict[mark,"type"] = "-1"
					mark = mark+1
				}
			}
		}
		return(grn_predict)
	}
}

##### Calculate accuracy of methods #####
Accuracy_of_method = function(grn1, grn2, targets, factors, method = c("cor", "GENIE3", "SCODE"))
{
	if(method == "cor" | method == "SCODE")
	{
		### real interactions
		grn_index1 = apply(grn1, 1, function(x)
		{
			index = paste0(x, collapse = ";")
		})
		### predicted interactions
		grn_index2 = apply(grn2, 1, function(x)
		{
			index = paste0(x, collapse = ";")
		})
		### all possible interactions
		interaction_type = c("1","-1")
		all_index = unlist(lapply(targets, function(x1)
		{
			lapply(factors, function(x2)
			{
				lapply(interaction_type, function(x3)
				{
					return(paste0(x1,";",x2,";",x3))
				})
			})
		}))

		TP = length(intersect(grn_index1,grn_index2)) ### True positive (a predicted edge which is actually present)
		FP = length(grn_index2) - TP ### False positive (a predicted edge which is not actually present)
		TN = length(setdiff(all_index,union(grn_index1,grn_index2))) ### True Negative (an unpredicted edge which is actually not present)
		FN = length(grn_index1) - TP ### False Negative (an unpredicted edge which is actually present)
		TPR = TP/(TP+FN) ### True Positive Rate: the proportion of predicted edges that are actually present
		FPR = FP/(FP+TN) ### False Positive Rate: the proportion of predicted edges that are not actually present
		F_score = 2*TP/(2*TP+FP+FN)
		return(F_score)
	}
	if(method == "GENIE3")
	{
		### real interactions
		grn_index1 = apply(grn1, 1, function(x)
		{
			index = paste0(x[1:2], collapse = ";")
		})
		### predicted interactions
		grn_index2 = apply(grn2, 1, function(x)
		{
			index = paste0(x[1:2], collapse = ";")
		})
		### all possible interactions
		all_index = unlist(lapply(targets, function(x1)
		{
			lapply(factors, function(x2)
			{
				return(paste0(x1,";",x2))
			})
		}))

		TP = length(intersect(grn_index1, grn_index2)) ### True positive (a predicted edge which is actually present)
		FP = length(grn_index2) - TP ### False positive (a predicted edge which is not actually present)
		TN = length(setdiff(all_index,union(grn_index1,grn_index2))) ### True Negative (an unpredicted edge which is actually not present)
		FN = length(grn_index1) - TP ### False Negative (an unpredicted edge which is actually present)
		TPR = TP/(TP+FN) ### True Positive Rate: the proportion of predicted edges that are actually present
		FPR = FP/(FP+TN) ### False Positive Rate: the proportion of predicted edges that are not actually present
		F_score = 2*TP/(2*TP+FP+FN)
		return(F_score)
	}
}

#predict_method = c("cor","GENIE3","SCODE")
predict_method = "SCODE"
number_genes = 100
number_bins = 7
number_sc = 100
data_type = "clean"
dir = "G://data//model//SERGIO//sim_datasets//"
sim_dir = paste0(dir,as.character(number_genes),'G_',as.character(number_bins),'T_',as.character(number_sc),'cPerT_dynamics//')

##### Get gene regulation pairs from simulated networks #####
interaction_CID = read.table(paste0(sim_dir,"Interaction_cID.txt"))
grn = Get_gene_pair(interaction_CID)
factors = names(table(grn$factor))
targets = names(table(grn$target))
#write.csv(grn,paste0(sim_dir,"grn.csv"))

##### Get pseudotime of single cells (use monocle2)#####
if(data_type == "clean")
{
	filename = "exprS_clean.csv"
}
if(data_type == "noise")
{
	filename = "count_matrix_S_O_L_D.csv"
}
sim_exp = read.csv(paste0(sim_dir,filename),header = FALSE)
rownames(sim_exp) = c(0:(number_genes-1))
cell_type = c(1:number_bins)
cells_per_type = number_sc

#sim_mef = Construct_monocle2_trajectory(sim_exp, cell_type, cells_per_type)

##### Use correlation to identify gene regulation relationship #####
sim_exp_mat = as.matrix(sim_exp)
rownames(sim_exp_mat) = rownames(sim_exp)
colnames(sim_exp_mat) = colnames(sim_exp)
grn_predict = Regulation_prediction(sim_exp_mat, targets, factors, method = predict_method, cor_thre = 0.6, genie3_thre = nrow(grn), dir = sim_dir)

##### Check the accuracy of method #####
acc = Accuracy_of_method(grn1 = grn[,c(1,2,4)], grn2 = grn_predict[,c(1,2,4)], targets, factors, method = predict_method)


index = 99
ss10 = smooth.spline(sim_exp_mat[1,], sim_exp_mat[index,], df = 10)
plot(sim_exp_mat[1,], sim_exp_mat[index,], pch = 20, col = "steelblue")
lines(ss10, lty = 2, col = "red", type="b")

