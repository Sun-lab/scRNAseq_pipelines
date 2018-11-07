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


# ----------
# scRNA Template Functions
# ----------
run_barcodeRanks_emptyDrops = function(sce){
	# ref = "https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-3-tenx.html#calling-cells-from-empty-droplets"
	# cat(paste0(ref,"\n"))

	# Calling cells from empty droplets
	bcrank = DropletUtils::barcodeRanks(counts(sce))

	# Only show unique points for plotting speed.
	uniq = !duplicated(bcrank$rank)

	par(mfrow=c(1,1),mar=c(5,4,2,1),bty="n")
	plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", 
		xlab="Rank", ylab="Total UMI count", cex=0.5, cex.lab=1.2)
	abline(h=bcrank$inflection, col="darkgreen", lty=2,lwd=2)
	abline(h=bcrank$knee, col="dodgerblue", lty=2,lwd=2)
	legend("left",legend=c("Inflection","Knee"), bty="n", 
		col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2,lwd=2)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
	
	cat(paste0("Inflection = ",bcrank$inflection,"\n"))
	cat(paste0("Knee = ",bcrank$knee,"\n"))

	# summary(bcrank$total)
	# table(bcrank$total >= bcrank$knee)
	# table(bcrank$total >= bcrank$inflection)

	set.seed(100)
	print(date())
	e_out = DropletUtils::emptyDrops(counts(sce))
	print(date())
	# length(unique(e_out$FDR))
	# table(e_out$FDR)

	# tapply(e_out$Total, e_out$FDR, summary)
	list(bcrank = bcrank,e_out = e_out)
}
run_QC_filter = function(work_dir,sce){
	# Checks
	check_chr_name = "chromosome_name" %in% names(rowData(sce))
	if( !check_chr_name ) stop("'chromosome_name' is missing from rowData(sce)")
	check_mito = length(which(rowData(sce)$chromosome_name) == "MT") > 0
	if( !check_mito ) stop("Recode mitochondria contig as 'MT'")
	
	file_link = "https://www.genenames.org/cgi-bin/genefamilies/set/1054/download/branch"
	file_name = strsplit(file_link,"/")[[1]]
	file_name = file_name[length(file_name)]
	ribo_fn = file.path(work_dir,file_name)
	if( !file.exists(ribo_fn) ){
		system(sprintf("cd %s; wget %s",work_dir,file_link))
	}

	ribo = smart_RT(ribo_fn,sep='\t',header=TRUE)
	# ribo[1:2,]
	
	is_mito = which(rowData(sce)$chromosome_name == "MT")
	is_ribo = which(rowData(sce)$gene %in% ribo$Approved.Symbol)
	# length(is_mito)
	# length(is_ribo)
	
	sce = calculateQCMetrics(sce,feature_controls=list(Mt=is_mito, Ri=is_ribo))
	
	par(mfrow=c(2,2), mar=c(5, 4, 1, 1), bty="n")
	smart_hist(log10(sce$total_counts),xlab="log10(Library sizes)",main="", 
		breaks=20,ylab="Number of cells")
	smart_hist(log10(sce$total_features),xlab="log10(# of expressed genes)", 
		main="",breaks=20,ylab="Number of cells")
	smart_hist(sce$pct_counts_Ri,xlab="Ribosome prop. (%)",
		ylab="Number of cells",breaks=40,main="")
	smart_hist(sce$pct_counts_Mt,xlab="Mitochondrial prop. (%)", 
		ylab="Number of cells",breaks=80,main="")
	smoothScatter(log10(sce$total_counts),log10(sce$total_features), 
		xlab="log10(Library sizes)",ylab="log10(# of expressed genes)", 
		nrpoints=500,cex=0.5)
	smoothScatter(log10(sce$total_counts),sce$pct_counts_Ri,
		xlab="log10(Library sizes)", ylab="Ribosome prop. (%)",
		nrpoints=500,cex=0.5)
	smoothScatter(log10(sce$total_counts),sce$pct_counts_Mt,
		xlab="log10(Library sizes)", ylab="Mitochondrial prop. (%)",
		nrpoints=500,cex=0.5)
	smoothScatter(x=sce$pct_counts_Ri,y=sce$pct_counts_Mt,
		xlab="Ribosome prop. (%)", ylab="Mitochondrial prop. (%)",
		nrpoints=500,cex=0.5)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
	libsize_drop = isOutlier(sce$total_counts,nmads=3,type="lower",log=TRUE)
	feature_drop = isOutlier(sce$total_features_by_counts,nmads=3,type="lower",log=TRUE)
	mito_drop = isOutlier(sce$pct_counts_Mt,nmads=3,type="higher")
	ribo_drop = isOutlier(sce$pct_counts_Ri,nmads=3,type="higher")

	keep = !(libsize_drop | feature_drop | mito_drop | ribo_drop)
	filter_summary = smart_df(ByLibSize=sum(libsize_drop),
		ByFeature=sum(feature_drop),ByMito=sum(mito_drop),
		ByRibo=sum(ribo_drop),Remaining=sum(keep))
	
	list(sce=sce,keep=keep,filter_summary=filter_summary)
}
run_normalization = function(sce,min_mean){
	print(date())
	clusters = quickCluster(sce,min.mean=min_mean,method="igraph")
	print(table(clusters))
	print(date())
	sce = computeSumFactors(sce,cluster=clusters,min.mean=min_mean)
	print(date())
	print(summary(sizeFactors(sce)))

	# Remove cells with negative or very small size factors
	sce = sce[,which(sizeFactors(sce) > 0)]

	par(mfrow=c(1,2),mar=c(5,4,2,1),bty="n")
	smoothScatter(sce$total_counts, sizeFactors(sce), log="xy", 
		xlab="total counts", ylab="size factors")
	plot(sce$total_counts, sizeFactors(sce), log="xy", 
		xlab="total counts", ylab="size factors", 
		cex=0.3, pch=20, col=rgb(0.1,0.2,0.7,0.3))
	par(mfrow=c(1,1))

	sce = normalize(sce)
	
	sce
}
run_dimension_reduction = function(sce,span=0.05){
	print(date())
	new_trend = makeTechTrend(x=sce)
	print(date())
	fit = trendVar(sce,use.spikes=FALSE,loess.args=list(span=span))

	par(mfrow=c(1,1), mar=c(5,4,2,1), bty="n")
	plot(fit$mean, fit$vars, pch=20, col=rgb(0.1,0.2,0.7,0.6), 
		xlab="log(mean)", ylab="var")
	curve(new_trend(x), col="red", lwd=2, add=TRUE)
	curve(fit$trend(x), col="orange", lwd=2, add=TRUE)
	legend("topright", legend=c("Poisson noise", "observed trend"), 
		lty=1, lwd=2, col=c("red", "orange"), bty="n")
	
	# For Poisson trend
	fit_pois = fit
	fit_pois$trend = new_trend
	dec_pois = decomposeVar(fit=fit_pois)
	top_dec_pois = dec_pois[order(dec_pois$bio,decreasing=TRUE),]
	par(mfrow=c(2,2),oma=c(0,0,2,0))
	smart_hist(dec_pois$bio,breaks=30,xlab="Biological Variance")
	smart_hist(dec_pois$FDR,breaks=30,xlab="FDR")
	smart_hist(log10(dec_pois$FDR + 1e-6),breaks=30,xlab="log10(FDR + 1e-6)")
	mtext("Poisson Trend",outer=TRUE,cex=1.2)
	par(mfrow=c(1,1),oma=rep(0,4))
	plotExpression(sce, features=rownames(top_dec_pois)[1:10])
	
	# For Observed trend
	dec_obs = decomposeVar(fit=fit)
	top_dec_obs = dec_obs[order(dec_obs$bio,decreasing=TRUE),]
	par(mfrow=c(2,2),oma=c(0,0,2,0))
	smart_hist(dec_obs$bio,breaks=30,xlab="Biological Variance")
	smart_hist(dec_obs$FDR,breaks=30,xlab="FDR")
	smart_hist(log10(dec_obs$FDR + 1e-6),breaks=30,xlab="log10(FDR + 1e-6)")
	mtext("Observed Trend",outer=TRUE,cex=1.2)
	par(mfrow=c(1,1),oma=rep(0,4))
	plotExpression(sce, features=rownames(top_dec_obs)[1:10])
	
	list(fit=fit,new_trend=new_trend,dec_pois=dec_pois,dec_obs=dec_obs)
}


# ----------
# ggplot Functions
# ----------
gg.heatmap <- function(t2, ylab){
  #create a new variable from incidence
  if(max(t2$value) > 1000){
    t2$valueFactor = cut(t2$value, 
                         breaks = c(-1,0,5,10,50,100,1000,max(t2$value)),
                         labels=c("0","1-5","6-10","11-50","51-100","101-1000",">1000"))
  }else{
    t2$valueFactor = cut(t2$value, 
                         breaks = c(-1,0,5,10,50,100,max(t2$value)),
                         labels=c("0","1-5","6-10","11-50","51-100","101-1000"))
  }
  
  #change level order
  t2$valueFactor = factor(as.character(t2$valueFactor),
                          levels=rev(levels(t2$valueFactor)))
  
  y.axis.size = 8
  if(nrow(t2) > 700)   y.axis.size = 7
  if(nrow(t2) > 800)   y.axis.size = 6
  
  col1 = c("#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da")
  col2 = c("#d53e4f", col1)
  
  if(max(t2$value) > 1000){
    col2use = col2
  }else{
    col2use = col1
  }
  
  g1 = ggplot(t2, aes(Cell_Type, cluster)) + 
    theme_bw() +
    geom_tile(aes(fill = valueFactor), color = "white", size=0.25) +
    scale_fill_manual(values=col2use) +  
    ylab(ylab) + xlab("cell types") +
    theme(legend.title=element_blank(),
          legend.text=element_text(colour="grey40",size=9,face="bold"),
          legend.key.height=grid::unit(0.8,"cm"),
          legend.key.width=grid::unit(0.4,"cm"),
          plot.title = element_text(size=16),
          axis.title = element_text(size=14, face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1), 
          axis.text.y = element_text(size=y.axis.size), 
          panel.border=element_blank())
  g1
}
ggplot_custom = function(DATA,X,Y,COL){
	ggplot(DATA,aes_string(x = X,y = Y,col = COL)) + 
		geom_point(size=0.2,alpha=0.6) + theme_classic() + 
		guides(color = guide_legend(override.aes = list(size=3),
			ncol = 10,byrow = TRUE)) +
		theme(legend.position="bottom")
}

