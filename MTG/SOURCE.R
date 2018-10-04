# ----------
# Paul's R Source File
# ----------

# ----------
# Libraries
# ----------
if(FALSE){
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(mclust)) # for naive clustering
suppressPackageStartupMessages(library(matrixcalc)) # to calculate rcond??
suppressPackageStartupMessages(library(numDeriv)) # to approx gradient and hessian
suppressPackageStartupMessages(library(survival)) # to run survival regressions
suppressPackageStartupMessages(library(splines)) # for possible cubic spline fitting
}

# ----------
# General Functions
# ----------
smart_date = function(){
	paste(strsplit(date()," ")[[1]][c(2,3,5)],collapse="-")
}
smart_hist = function(x,...){
	hist(x,col="gray",freq=FALSE,...)
	lines(density(x,na.rm=TRUE),lwd=1.5,lty=2,col="blue")
}
smart_merge = function(x,y,mess=NULL,...){
	if( !is.null(mess) ){
		intersect_vars = paste(intersect(names(x),names(y)),collapse=", ")
		cat(paste0("Merging dataframes on variables = { ",intersect_vars," }\n"))
	}
	
	merge(x,y,by=intersect(names(x),names(y)),...)
}
smart_table = function(...){
	table(...,useNA='ifany')
}
smart_RT = function(...){
	read.table(...,stringsAsFactors=FALSE)
}
smart_WT = function(...){
	write.table(...,row.names=FALSE,quote=FALSE)
}
smart_df = function(...){
	data.frame(...,stringsAsFactors=FALSE)
}
smart_sapply = function(...){
	sapply(...,USE.NAMES=FALSE,simplify=FALSE)
}
convert_obj = function(x){
	# convert x (vector,matrix) of character elements to numbers
	# x = char_mat
	# x = c("a","b","c","a")
	conv_df = data.frame(xx = unique(sort(x)),stringsAsFactors = FALSE)
	conv_df$num = seq(nrow(conv_df))
	
	if( length(dim(x)) == 0 ){ # aka a vector
		orig_df = data.frame(xx = x,stringsAsFactors = FALSE)
		orig_df$order = seq(nrow(orig_df))
		orig_df = smart_merge(orig_df,conv_df)
		orig_df = orig_df[order(orig_df$order),]
		orig_df$num
	} else if( length(dim(x)) == 2 ){ # aka a matrix
		orig_df = data.frame(xx = as.vector(x),stringsAsFactors = FALSE)
		orig_df$order = seq(nrow(orig_df))
		orig_df = smart_merge(orig_df,conv_df)
		orig_df = orig_df[order(orig_df$order),]
		matrix(orig_df$num,nrow(x),ncol(x),byrow=FALSE)
	}
}
convert_color = function(x,y){
	# x = demos$sex; y = c("red","green")
	x2 = unique(sort(x))
	if( length(x2) > length(y) ){
		stop("Mismatch in input colors")
	}
	
	conv_df = data.frame(xx = x2, col = y[seq(length(x2))], stringsAsFactors = FALSE)
	orig_df = data.frame(xx = x, stringsAsFactors = FALSE)
	orig_df$order = seq(nrow(orig_df))
	orig_df = smart_merge(orig_df,conv_df)
	orig_df = orig_df[order(orig_df$order),]
	orig_df$col
}

smart_image = function(orig_mat,border=FALSE,axes=FALSE,...){
	if(FALSE){
		# orig_mat = c(1,1,2,2,2,3)
		# orig_mat = tmp_mat
		orig_mat = blah$MAT_COLORS
		border = FALSE; axes = FALSE
	}
	
	# Code for creating the heatmap
	if( !is.character(orig_mat) && is.numeric(orig_mat) ){
		# No border
		if(border == FALSE){
			par(mar = rep(0,4))
		}
		# Orienting the matrix to how its organized
		if(length(dim(orig_mat)) == 2 && nrow(orig_mat) > 1){
			image_mat = t(orig_mat[rev(seq(nrow(orig_mat))),])
		} else if(length(dim(orig_mat)) == 0 && is.null(nrow(orig_mat)) ){
			image_mat = as.matrix(orig_mat)
		}
		# Plot
		if(length(dim(image_mat)) == 2){
			num_rows = nrow(image_mat); num_cols = ncol(image_mat)
			x_coor = seq(0,1,1/(num_rows+1))[-c(1,num_rows+2)] - 0.5 /(num_rows+1)
			xlims = as.numeric(summary(x_coor)[c(1,6)] + c(-1,1) * (0.5 /(num_rows+1)))
			y_coor = seq(0,1,1/(num_cols+1))[-c(1,num_cols+2)] - 0.5 /(num_cols+1)
			ylims = as.numeric(summary(y_coor)[c(1,6)] + c(-1,1) * (0.5 /(num_cols+1)))
			
			image(x = x_coor, y = y_coor, z = image_mat, 
				axes = axes,
				xaxs = "i",yaxs = "i",
				xlim = xlims, ylim = ylims,
				...)

		} else if( length(dim(image_mat)) == 0 ){
			
			stop("No code for vector input yet")
		}
		par(mar=c(5,4,4,2)+0.1)
	}
	
	# If orig_mat is a character color matrix, convert then plot
	if( is.character(orig_mat) && !is.numeric(orig_mat) ){
		# orig_mat = bb_col
		numeric_mat = orig_mat
		numeric_mat = matrix(NA,nrow(orig_mat),ncol(orig_mat))
		mat_names = names(table(orig_mat))
		mat_colors = rep(NA,length(mat_names))
		count = 1
		for(mat_name in mat_names){
			mat_colors[count] = mat_name
			numeric_mat[orig_mat == mat_name] = count
			count = count + 1
		}
		# numeric_mat = matrix(as.numeric(numeric_mat), nrow = dim(numeric_mat)[1])
		smart_image(numeric_mat,col=mat_colors)
	}
	
}
convert_MAT_COLORS = function(MAT){
	# MAT = as.matrix(bb)
	MAT_COLORS = matrix(NA,nrow(MAT),ncol(MAT))
	colnames(MAT_COLORS) = colnames(MAT)
	rownames(MAT_COLORS) = rownames(MAT)
	LEGENDS = list()
	
	for(ii in colnames(MAT)){
		# ii = names(bb)[1]
		VARS = MAT[,ii]
		uVARS = sort(unique(VARS))
		nVARS = length(uVARS)
		if( nVARS == 2 ){
			COLORS = c("black","white")
		} else {
			COLORS = rainbow(nVARS)
		}
		map_VARS_COLORS = smart_df(VARS = uVARS,COLORS)
		LEGENDS[[ii]] = map_VARS_COLORS
		dfVARS = smart_df(VARS,order=seq(length(VARS)))
		dfVARS = smart_merge(dfVARS,map_VARS_COLORS)
		dfVARS = dfVARS[order(dfVARS$order),]
		MAT_COLORS[,ii] = dfVARS$COLORS
	}
	
	list(MAT_COLORS=MAT_COLORS,LEGENDS=LEGENDS)
}

make_dummy = function(x){
	# x = gsub(" ","_",x)
	len_x = length(x)
	all_factors = sort(unique(x))
	num_dummy = length(all_factors) - 1
	fact_matrix = data.frame(matrix(0,nrow=len_x,ncol=max(1,num_dummy)))
	
	if(num_dummy == 0){
		stop("Variable vector is constant")
	} else {
		names(fact_matrix) = paste0(all_factors[-1],"_vs_",all_factors[1])
	}
	
	if(num_dummy == 0){
		for(i in 1:len_x){
			fact_matrix[i,1] = 1
		}
	} else {
		for(i in 1:len_x){
			pos = which(x[i]==all_factors) - 1
			if(is.na(x[i])) fact_matrix[i,] = NA
			else if(pos>0) fact_matrix[i,pos] = 1
		}
	}
	
	fact_matrix
}
smart_uniq_df = function(input_df,vars){
	input_df[!duplicated(input_df[,vars]),]
}
smart_perm = function(x){
	len_x = length(x)
	matrix(unlist(combinat::permn(x)),ncol=len_x,byrow=TRUE)
}
collapse_var = function(ORIG_VAR,ORIG_VALUES,NEW_VALUE){
	ORIG_VAR[which(ORIG_VAR %in% ORIG_VALUES)] = NEW_VALUE
	ORIG_VAR
}
name_change = function(DATA,ORIG_NAME,NEW_NAME){
	index = which(names(DATA) == ORIG_NAME)
	if( length(index) > 0 ){
		names(DATA)[index] = NEW_NAME
	}
	DATA
}
complete_cases = function(DATA,NAMES){
	for(ii in NAMES){
		DATA = DATA[which(!is.na(DATA[,ii])),]
	}
	DATA
}
smart_SN = function(x,digits=2){
	# For scientific notation
	formatC(x,format = "e",digits = digits)
}
smart_digits = function(x,digits=2){
	sprintf(paste0("%.",digits,"f"),round(x,digits))
}
smart_label = function(DATA,VAR,VD){
	if(FALSE){
		DATA = DATA_cc
		VAR = model_vars[6]; VAR
		VD = VD
	}
	
	# Determine if each variable in coefficients table is continuous, categorical or interaction
	INT = ifelse(length(grep(":",VAR)) > 0,TRUE,FALSE)
	if( !INT ){
		NUM = ifelse(class(DATA[,VAR]) == "numeric",TRUE,FALSE)
		if( NUM ){
			DISP = VD$D[which(VD$VAR == VAR)]
		} else {
			tmp_vec = sort(unique(DATA[,VAR]))
			DISP = paste0(VD$D[which(VD$VAR == VAR)],":",paste0(tmp_vec[-1]," vs ",tmp_vec[1]))
		}
	} else {
		VARS = strsplit(VAR,":")[[1]]
		DISPS = sapply(VARS,function(x) smart_label(DATA=DATA,VAR=x,VD=VD))
		if(class(DISPS) == "list"){
			DISP = paste0(DISPS[[1]]," by ",DISPS[[2]])
		} else {
			DISP = paste0(DISPS[1]," by ",DISPS[2])
		}
	}
	DISP
}
rdsn = function(x,rdnum){
	# round quantities, if smaller than epsilon
	ifelse(abs(x) <= 10^(-rdnum),
		smart_SN(x,digits=1),
		format(round(x,rdnum),nsmall=rdnum))
}


# ----------
# File manipulation
# ----------
smart_append = function(new_fn,vec_fn){
	# new_fn = name of new file to create
	# vec_fn = vector of files to append together in the order they're stored
	aa = file.create(new_fn)
	bb = sapply(vec_fn,function(x) file.append(new_fn,x))
}
smart_remove = function(FN){
	# remove file
	if( file.exists(FN) ){
		blah = file.remove(FN)
	} else {
		warning(paste0(FN," doesn't exist"))
	}
}
smart_Rcpp = function(new_Rcpp_fn,vec_Rcpp_fn){
	smart_remove(new_Rcpp_fn)
	smart_append(new_Rcpp_fn,vec_Rcpp_fn)
		# Create merged Rcpp file
	cat(paste0("Create ",new_Rcpp_fn,"\n"))
	sourceCpp(file = new_Rcpp_fn,showOutput = TRUE)
}
smart_mkdir = function(input_dir){
	if( !file.exists(input_dir) || !dir.exists(input_dir) ){
		dir.create(input_dir)
	}
}
smart_sprintf = function(...){
	orig = sprintf(...)
	orig = gsub("\n","",orig)
	orig = gsub("\t","",orig)
	orig
}


# ----------
# Simulation Summary Functions
# ----------
calc_BIAS = function(input_df,expected_MEAN){
	# Assume input_df is two columns, each row is a replicate with an estimate and ASE
	mean(input_df[,1] - expected_MEAN)
}
calc_SE = function(input_df){
	# SE = standard error of the parameter estimator
	# Assume input_df is two columns, each replicates estimate and ASE estimate
	sd(input_df[,1])
}
calc_SEE = function(input_df){
	# SEE = empirical average of the standard error estimator
	# Assume input_df is two columns, each replicates estimate and ASE estimate
	mean(input_df[,2])
}
calc_CP = function(input_df,prop=0.95,expected_MEAN){
	# CP = empirical coverage percentage of the 95% confidence interval
	# Assume input_df is two columns, each replicates estimate and ASE estimate
	prop_quantile = 1 - (1 - prop)/2
	all_conf_int = matrix(0,nrow=nrow(input_df),ncol=2)
	all_conf_int[,1] = input_df[,1] - qnorm(prop_quantile)*input_df[,2]
	all_conf_int[,2] = input_df[,1] + qnorm(prop_quantile)*input_df[,2]
	mean(ifelse(expected_MEAN >= all_conf_int[,1] 
		& expected_MEAN <= all_conf_int[,2],1,0))
}
calc_POWER = function(input_df,H0,prop=0.95){
	# input_df = blah2; prop = 0.95
	# Check proportion of times we reject null H0
	vec_z = abs(input_df[,1] - H0) / input_df[,2]
	vec_p = 2*(1-pnorm(vec_z))
	mean(ifelse(vec_p < 1-prop,1,0))
}


# ----------
# Plotting Functions
# ----------
show_png = function(){
	cat("# ----------\n")
	cat("# PNG Template\n")
	cat("# ----------\n")
	cat("png(file.path(),units='px',height=,width=,res=250,type='cairo')\n")
	cat("...\n")
	cat("dev.off()\n")
}
plot_KM = function(DATA,VARS,OUTCOME="Time to ...",my_lwd=3,...){
	if(FALSE){
		DATA = blah
		VARS = "IDH"
		VARS = "histological_type"
		VARS = c("stage")
		VARS = c("new_stage","binE")
	}
	
	if( !all(c("Time","Delta") %in% names(DATA)) ){
		stop("Missing Time and Delta")
	}
	
	# Construct formula
	control_vars = paste(VARS,collapse=" + ")
	my_formula = formula(paste0("Surv(Time,Delta) ~ ",control_vars))
	
	# Get complete cases
	DATA = complete_cases(DATA,VARS)

	# Perform logrank test
	logrank_out = survdiff(my_formula,data = DATA)
	pval = 1 - pchisq(logrank_out$chisq,df=length(logrank_out$obs)-1)
	
	# Get uniq values/colors for legend
	if( length(VARS) == 1 ){
		uniq_df = smart_df(sort(unique(DATA[,VARS])))
	} else {
		uniq_df = smart_df(unique(DATA[,VARS]))
	}
	uniq_df = smart_df(uniq_df[do.call(order,uniq_df),])
	uniq_df$values = as.character(apply(uniq_df,1,function(x) paste(x,collapse=',')))
	uniq_df$col = seq(nrow(uniq_df))
	# uniq_df = uniq_df[order(uniq_df$values),]
	num_values = nrow(as.matrix(unique(DATA[,VARS])))
	
	# Plot
	plot(survfit(my_formula,data=DATA),mark.time=TRUE,
		# main=paste0(OUTCOME," ~ ",control_vars,"; P = ",round(pval,3)),
		#xlab="Days",ylab="P(Survival)",
		# cex.main=0.8,
		lwd=my_lwd,
		col=uniq_df$col,...)
	legend(x="bottomleft",col=uniq_df$col,lwd=my_lwd,legend=uniq_df$values,cex=0.75)

}
my_ggplot_KM = function(fit,nrow=2,fs=20,fs_leg=12){
	if(FALSE){
		fit = survfit(Surv(Time,Delta) ~ bin_rMB,data = aa)
	}
	my_fontsize = fs
	# ggsurvplot(fit)
	ggsurvplot(fit,risk.table=FALSE,conf.int=FALSE,
		xlab = "Time (days)",ylab = "Survival Probability",
		legend = "bottom",size = 2,
		#ggtheme = theme_grey() + 
		#	theme(plot.title = element_text(hjust = 0.5,face = "bold",size=15)),
		font.x=my_fontsize,font.tickslab = my_fontsize,font.y=my_fontsize,
		font.legend = list(size=fs_leg,face="bold")) +
		guides(colour = guide_legend(nrow=nrow))
}
make_yx = function(...){
	abline(a=0,b=1,lty=2,...)
}
get_colors = function(color){
	colors()[grep(color,colors())]
}
show_colors = function(color){
	# color = "blue"
	tmp_df = smart_df(col=get_colors(color))
	tmp_df$num = seq(nrow(tmp_df))
	tmp_df
	smart_image(as.matrix(tmp_df$num),col=tmp_df$col)
	num_rows = nrow(tmp_df)
	text(smart_df(x=0.25,
		y=rev(seq(0,1,1/(num_rows+1))[-c(1,num_rows+2)] - 0.5 /(num_rows+1))),
		labels=tmp_df$col,cex=0.75)
}


# ----------
# Statistical Functions
# ----------
smart_solve = function(MATRIX){
	# Invert MATRIX
	mat_rcond = rcond(MATRIX)
	if( mat_rcond == 0 ){
		NA
	} else {
		solve(MATRIX,tol=0.1*mat_rcond)
	}
}
smart_optim = function(...){
	my_control = list(fnscale = -1,maxit = 2e3,abstol = 1e-40,reltol = 1e-40)
	
	tryCatch({
		optim(...,control = my_control)
		}, error = function(err){
			list(par=NA)
		}
	)
}
grad_descent = function(params,LL_func,gr_func=NULL,max_iter=2e3,my_epsilon=1e-5){
	if(FALSE){
		params = init_params
		LL_func = joint_wLL_1
		gr_func = joint_wGRAD_1
		# gr_func = NULL
		max_iter = 2e3; my_epsilon = 1e-6
	}

	iter = 0; curr_LL = NA; curr_params = NA

	# Newer code
	if( is.null(gr_func) ){
		my_grad = function(input_params){
			grad(LL_func,input_params)
		}
		my_hess = function(input_params){
			hessian(LL_func,input_params)
		}
	} else {
		my_grad = function(input_params){
			gr_func(input_params)
		}
		my_hess = function(input_params){
			jacobian(gr_func,input_params)
		}
	}

	while(iter < max_iter){
		old_LL = LL_func(params)
		old_grad = my_grad(params)
		norm_old_grad = Rcpp_norm(old_grad)
		
		if( any(is.na(old_grad)) || is.na(old_LL) ){
			iter = max_iter
			break
		}
		
		check_update = 0
		for(ss in seq(0,40)){
			# ss = 0
			new_params = params + 1/4^ss * old_grad / max(c(norm_old_grad,1))
			new_LL = LL_func(new_params)
			if( !is.na(new_LL) && new_LL > old_LL ){
				old_LL = new_LL
				params = new_params
				check_update = 1
				break
			}
		}
		
		# After 1 iteration, 
		# if norm(old_LL,curr_LL) < epsilon & norm(curr_params,params) < epsilon,
		# if yes check norm_eff_score = norm(solve(-hessian) %*% score)
		if(iter > 0){
			
			if( !is.na(old_LL) && abs(curr_LL - old_LL) < my_epsilon 
				&& Rcpp_norm(curr_params - params) < my_epsilon ){
				
				tmp_score = my_grad(params)
				tmp_hess = my_hess(params)
				tmp_eff_score = as.numeric(solve(tmp_hess,
					tol=0.1*rcond(tmp_hess)) %*% tmp_score)
				
				if( Rcpp_norm(tmp_eff_score) < my_epsilon ){
					break
				}
				
			}
		
		}
		
		curr_LL = old_LL; curr_params = params
		iter = iter + 1
	}

	final_LL = LL_func(params)
	final_score = my_grad(params)
	final_hess = my_hess(params)
	final_eff_score = as.numeric(solve(final_hess,
		tol=0.1*rcond(final_hess)) %*% final_score)
	
	list(conv = ifelse(iter < max_iter,"yes","no"), 
		iter = iter, LL = final_LL,
		norm_score = Rcpp_norm(final_score),
		norm_eff_score = Rcpp_norm(final_eff_score),
		params = params)
}
bin_cont_var = function(VAR,NUM_GROUPS,ROUND=3,binNUM=FALSE){
	if(FALSE){
		# VAR = sort(runif(50))
		# VAR = aa$AscatM_E
		VAR = aa$LOCI
		NUM_GROUPS = 3
		ROUND = 3
	}
	
	my_quantiles = as.numeric(quantile(x = VAR,
		probs = seq(NUM_GROUPS-1)/NUM_GROUPS,
		na.rm = TRUE))
	
	out_VAR = rep(NA,length(VAR))
	for(ii in seq(NUM_GROUPS)){
		if(ii == 1){
			if(binNUM){
				out_VAR[which(VAR <= my_quantiles[ii])] = ii
			} else {
				out_VAR[which(VAR <= my_quantiles[ii])] = paste0(ii,") ",round(min(VAR,na.rm=TRUE),ROUND),
					"-",round(my_quantiles[ii],ROUND))
			}
		} else if(ii == NUM_GROUPS){
			if(binNUM){
				out_VAR[which(VAR > my_quantiles[ii-1])] = ii
			} else {
				out_VAR[which(VAR > my_quantiles[ii-1])] = paste0(ii,") ",round(my_quantiles[ii-1],ROUND),
					"-",round(max(VAR,na.rm=TRUE),ROUND))
			}
		} else {
			if(binNUM){
				out_VAR[which(VAR > my_quantiles[ii-1] & VAR <= my_quantiles[ii])] = ii
			} else {
				out_VAR[which(VAR > my_quantiles[ii-1] & VAR <= my_quantiles[ii])] = paste0(ii,") ",round(my_quantiles[ii-1],ROUND),
					"-",round(my_quantiles[ii],ROUND))
			}
		}
	}
	
	# smart_df(VAR,out_VAR)
	if(binNUM) out_VAR = as.character(out_VAR)
	
	out_VAR
}
count_freq = function(VAR){
	# VAR = clin_list$clin_overall$stage
	tmp_tab = smart_table(VAR)
	names(tmp_tab)[which(is.na(names(tmp_tab)))] = "N/A"
	smart_df(value=names(tmp_tab),freq=as.numeric(tmp_tab),prop=round(as.numeric(tmp_tab)/length(VAR),3))
}
smart_trans = function(x,trans){
	if(trans == "log"){
		log(x)
	} else if(trans == "inv_log"){
		exp(x)
	} else if(trans == "logit"){
		log(x/(1-sum(x)))
	} else if(trans == "inv_logit"){
		exp(x)/(1+sum(exp(x)))
	} else if(trans == "q"){
		# assumes input is unc_q
		x2 = smart_trans(x,"inv_logit")
		c(x2,1-sum(x2))
	} else if(trans == "unc_q"){
		x2 = x[-length(x)]
		smart_trans(x2,"logit")
	}
}
smart_max_matrix = function(x){
	row_max = which.max(apply(x,1,max))
	col_max = which.max(apply(x,2,max))
	c(row_max,col_max)
}


# ----------
# Latex Output Functions
# ----------
format_latex = function(INPUT){
	# INPUT = "optE_AIC%"
	INPUT2 = gsub("%","\\\\%",INPUT)
	INPUT2 = gsub("_","\\\\_",INPUT2)
	INPUT2
}
clean_repeats = function(VEC){
	if(FALSE){
		VEC = c(rep("a",2),rep("b",2),"a","c")
		VEC
	}
	
	curr_string = NA
	for(ii in seq(length(VEC))){
		# ii = 1
		if(ii == 1){
			curr_string = VEC[ii]
		} else {
			if(VEC[ii] == curr_string){
				VEC[ii] = ""
			} else {
				curr_string = VEC[ii]
			}
		}
	}
	
	VEC
}
print_latex_table = function(DATA,repeat_VARS=NULL,my_align=NULL,
	add_table=FALSE,fontsize=NULL,caption=NULL,label=NULL,midrule1=NULL,...){
	
	if(FALSE){
		# DATA = res_cont_nothing; repeat_VARS = c("N","entropy"); add_table = TRUE; fontsize = "tiny"
		DATA = res_cont_knowE; repeat_VARS = "entropy"; add_table = TRUE; fontsize = ""
		DATA = est; repeat_VARS = NULL; add_table = !TRUE; fontsize = ""; file = "C:/Users/Admin/Desktop/hola.txt"
	}

	# New code: convert every column to character
	if(class(DATA) == "data.frame"){
		orig_names = names(DATA) 
	} else if(class(DATA) == "matrix"){
		orig_names = colnames(DATA)
	}
	
	DATA = smart_df(apply(DATA,2,as.character))
	
	if( !is.null(repeat_VARS) && length(repeat_VARS) > 0 ){ # loop thru vector(column) to find repeats and replace with ""
		# print("blah")
		tmp_index = which(orig_names %in% repeat_VARS)
		DATA[,tmp_index] = apply(DATA[,tmp_index,drop=FALSE],2,clean_repeats)
	}
	
	prep_DATA = DATA
	
	cat("\n",...)
	if(add_table){
		cat(paste0("\\begin{table}[!htbp] \n\\centering\n"),...)
		if( !is.null(fontsize) ) cat(paste0("\\",fontsize,"\n"),...)
		if( !is.null(caption) ) cat(paste0("\\caption{",caption,"}\n"),...)
		if( !is.null(label) ) cat(paste0("\\label{tab:",label,"}\n"),...)
	}
	
	if( is.null(my_align) ){
		cat(paste0("\\begin{tabular}{l",paste(rep("c",ncol(prep_DATA)-1),collapse=""),"}\n"),...)
	} else {
		cat(paste0("\\begin{tabular}{",my_align,"}\n"),...)
	}
	cat("\\toprule\n",...)
	# cat(paste0(paste(sapply(names(prep_DATA),format_latex),collapse=" & ")," \\\\\n"))
	cat(paste0(paste(sapply(orig_names,format_latex),collapse=" & ")," \\\\\n"),...)
	if( is.null(midrule1) ){
		cat("\\midrule\n",...)
	} else {
		cat(paste0(midrule1,"\n"),...)
	}
	apply(prep_DATA,1,function(x) cat(paste0(paste(sapply(x,format_latex),collapse=" & ")," \\\\\n"),...))
	cat("\\bottomrule\n\\end{tabular}\n",...)

	if(add_table){
		# if( !is.null(caption) ) cat(paste0("\\caption{",caption,"}\n"),...)
		# if( !is.null(label) ) cat(paste0("\\label{tab:",label,"}\n"),...)
		cat(paste0("\\end{table}\n"),...)
	}

	cat("\n",...)
}

proto_LATEX_table = function(orig_TABLE,center=FALSE,fontsize=NULL,table=FALSE,...){
	if(FALSE){
		# orig_TABLE = smart_table(my_ith[,c("myP1_SCP2","myP2_SCP2")])
		orig_TABLE = prop_repro
		center = FALSE
		fontsize = NULL
		table = FALSE
	}
	table_names = format_latex(names(dimnames(orig_TABLE))) # (row,col)
	TABLE = orig_TABLE
	TABLE = rbind(colnames(TABLE),TABLE); colnames(TABLE) = NULL
	TABLE = cbind("",rownames(TABLE),TABLE); rownames(TABLE) = NULL
	num_cols = ncol(TABLE)
	num_rows = nrow(TABLE)
	if(table) cat("\\begin{table}[!htbp]\n",...)
	if(table && center) cat("\\centering\n",...)
	if( table && !is.null(fontsize) ) cat(paste0("\\",fontsize,"\n"),...)
	cat(paste0("\\begin{tabular}{",paste(c(rep("c",2),"|",rep("c",num_cols-2)),collapse=""),"}\n"),...)
	cat("\\toprule\n",...)
	cat(paste0(" && \\multicolumn{",num_cols-2,"}{c}{",table_names[2],"} \\\\\n"),...)
	cat(paste0(paste(TABLE[1,],collapse=" & ")," \\\\\n"),...)
	cat("\\midrule\n",...)
	cat(paste0("\\multirow{",num_rows-1,"}{*}{\\rotatebox[origin=c]{90}{",table_names[1],"}}\n"),...)
	apply(TABLE[-1,],1,function(x) cat(paste0(paste(format_latex(x),collapse=" & ")," \\\\\n"),...))
	cat("\\bottomrule\n\\end{tabular}\n",...)
	if(table) cat("\\end{table}\n",...)
}
# proto_LATEX_table(smart_table(my_ith[,c("myP1_SCP2","myP2_SCP2")]),center=TRUE)


# ----------
# Experimental
# ----------
if(FALSE){ # To run negative binomial fit

library(MASS)
xx = rnbinom(1e3,mu=1)

}
if(FALSE){ # Simulate coxph

my_hazard = function(input,func,choice){
	if(choice == 1){
		cc = exp(2)
		dd = 2
		if(func == "cumhaz"){
			cumhaz = cc * input^dd
			cumhaz
		} else if(func == "inv_cumhaz"){
			inv_cumhaz = (input / cc)^(1/dd)
			inv_cumhaz
		}
	}
}

set.seed(100)
nn = 2e3
true_beta = c(3,2)
true_choice = 1
obs_data = smart_df(V = runif(nn), V2 = runif(nn),x1 = rnorm(nn,0,1),x2 = rbinom(nn,1,0.3))
obs_data$T = my_hazard(input = -log(obs_data$V) / exp(as.numeric(as.matrix(obs_data[,c("x1","x2")]) %*% true_beta)),
	func = "inv_cumhaz",choice = true_choice)
obs_data$C = -log(obs_data$V2) / 0.28
obs_data$Delta = 1*(obs_data$T <= obs_data$C)
smart_table(obs_data$Delta)
obs_data$X = pmin(obs_data$T,obs_data$C)
# out_surv = survfit(Surv(X,Delta) ~ x2,data = obs_data); plot(out_surv)
coxph_out = coxph(Surv(X,Delta) ~ x1 + x2,data = obs_data)
summary(coxph_out)
basehaz_df = basehaz(coxph_out,centered=FALSE)
plot(basehaz_df[,c("time","hazard")],type="s",col="red",lwd=2)
lines(basehaz_df$time,my_hazard(basehaz_df$time,func = "cumhaz",choice = true_choice),type="l",lty=2,col="blue",lwd=2)
lines(basehaz_df$time,my_hazard(basehaz_df$time,func = "cumhaz",choice = true_choice),type="s",lty=2,col="green",lwd=2)

exp_out = summary(survreg(Surv(X,Delta) ~ x1 + x2,data = obs_data,dist = "exponential"))
exp_out$table[,c(1,3)] = -1 * exp_out$table[,c(1,3)]
exp_out

}
if(FALSE){
set.seed(1)
uniq_mm = matrix(runif(6),3,2,byrow=TRUE)
mm = uniq_mm[sample(seq(nrow(uniq_mm)),5,replace=TRUE),]

get_uniq = function(mm){
	uniq_mat = matrix(NA,nrow(mm),ncol(mm))
	num_uniq = 1
	uniq_mat[num_uniq,] = mm[num_uniq,]
	for(rr in seq(2,nrow(mm))){
		# rr = 2
		cc = 0
		for(uu in seq(1,num_uniq)){
			# uu = 1
			if( all(mm[rr,] == uniq_mat[uu,]) ){
				cc = 1
			}
		}
		
		if(cc == 0){
			num_uniq = num_uniq + 1
			uniq_mat[num_uniq,] = mm[rr,]
		}
	}

	uniq_mat[seq(num_uniq),]
}

}






####

