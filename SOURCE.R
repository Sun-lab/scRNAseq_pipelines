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

