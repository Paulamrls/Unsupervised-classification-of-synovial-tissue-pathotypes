# Local functional enrichment of positive markers for numeric clusters 1..5.
# Gene lists derived from samples are never transmitted to a web service.
project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(project_root, ".r_libs", "4.5"), .libPaths()))
suppressPackageStartupMessages({library(msigdbr); library(ggplot2); library(patchwork)})
old_results <- file.path(project_root, "analysis_output", "results")
out <- file.path(project_root, "results")
dir.create(file.path(out, "enrichment"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "figures"), recursive=TRUE, showWarnings=FALSE)
old_labels <- setNames(c("Endothelial","Fibroid","Myeloid","Lymphoid","IFN_high"), 1:5)
deg <- lapply(old_labels, function(x) read.csv(file.path(old_results,"DEGs",paste0("DEGs_ALL_",x,".csv")), check.names=FALSE))
universe <- Reduce(intersect, lapply(deg, function(x) unique(na.omit(x$gene))))
reactome <- msigdbr(db_species="HS", species="human", collection="C2", subcollection="CP:REACTOME")
db_version <- unique(reactome$db_version); stopifnot(length(db_version)==1)
sets <- lapply(split(reactome$gene_symbol, reactome$gs_name), function(x) intersect(unique(x),universe))
sets <- sets[lengths(sets)>=10 & lengths(sets)<=500]

ora <- function(id) {
  d <- deg[[id]]
  genes <- intersect(unique(d$gene[!is.na(d$padj) & d$padj<0.05 & d$log2FoldChange>=1]), universe)
  ovgenes <- lapply(sets, intersect, y=genes); ov <- lengths(ovgenes); keep <- ov>=3
  s <- sets[keep]; ovgenes <- ovgenes[keep]; ov <- ov[keep]; bg <- lengths(s)
  pv <- phyper(ov-1, bg, length(universe)-bg, length(genes), lower.tail=FALSE)
  z <- data.frame(cluster_id=as.integer(id), pathway=names(s), overlap=ov,
    input_genes=length(genes), pathway_genes_in_universe=bg,
    universe_genes=length(universe), gene_ratio=ov/length(genes),
    background_ratio=bg/length(universe),
    fold_enrichment=(ov/length(genes))/(bg/length(universe)), p_value=pv)
  z$fdr <- p.adjust(z$p_value,"BH")
  z$overlap_genes <- vapply(ovgenes,paste,collapse=";",FUN.VALUE=character(1))
  z <- z[order(z$fdr,-z$fold_enrichment,-z$overlap),]
  write.csv(z,file.path(out,"enrichment",paste0("Reactome_cluster_",id,".csv")),row.names=FALSE)
  z
}
all_results <- lapply(names(old_labels), ora)
combined <- do.call(rbind,all_results)
write.csv(combined,file.path(out,"enrichment","Reactome_all_clusters.csv"),row.names=FALSE)
top <- do.call(rbind,lapply(split(combined,combined$cluster_id),function(x) head(x[x$fdr<0.05,],5)))
stopifnot(nrow(top)>0,all(1:5 %in% top$cluster_id))
top$cluster_label <- factor(paste("Cluster",top$cluster_id),levels=paste("Cluster",1:5))
top$minus_log10_fdr <- -log10(pmax(top$fdr,.Machine$double.xmin))
top$term_short <- gsub("_"," ",sub("^REACTOME_","",top$pathway))
top$term_short <- vapply(top$term_short,
                         function(x) paste(strwrap(x, width=34), collapse="\n"),
                         FUN.VALUE=character(1))
fold_limits <- range(top$fold_enrichment,na.rm=TRUE)
fold_breaks <- pretty(fold_limits,n=5)
fold_breaks <- fold_breaks[fold_breaks>=fold_limits[1] & fold_breaks<=fold_limits[2]]
overlap_limits <- range(top$overlap,na.rm=TRUE)
overlap_breaks <- pretty(overlap_limits,n=5)
overlap_breaks <- overlap_breaks[overlap_breaks>=overlap_limits[1] & overlap_breaks<=overlap_limits[2]]
make_panel <- function(cluster_id) {
  panel_df <- top[top$cluster_id==cluster_id,]
  panel_df$term_short <- factor(panel_df$term_short,levels=rev(panel_df$term_short))
  ggplot(panel_df,aes(minus_log10_fdr,term_short,size=overlap,colour=fold_enrichment))+
    geom_point(alpha=.92)+
    scale_colour_gradient(low="#56B1F7",high="#CA0020",limits=fold_limits,breaks=fold_breaks)+
    scale_size_continuous(range=c(4.5,10),limits=overlap_limits,breaks=overlap_breaks)+
    guides(
      colour=guide_colourbar(title="Fold enrichment",title.position="top",
        barwidth=grid::unit(6.5,"cm"),barheight=grid::unit(.65,"cm")),
      size=guide_legend(title="Gene overlap",title.position="top",override.aes=list(alpha=1))
    )+
    labs(title=paste("Cluster",cluster_id),x="-log10(FDR)",y=NULL)+
    theme_bw(base_size=14)+
    theme(plot.title=element_text(face="bold",size=15,hjust=.5,colour="black"),
      axis.text=element_text(size=12.5,colour="black"),
      axis.text.y=element_text(size=12.5,lineheight=.95,colour="black"),
      axis.title.x=element_text(size=13.5,colour="black"),panel.grid.minor=element_blank(),
      legend.title=element_text(size=14,face="bold",colour="black"),
      legend.text=element_text(size=13,colour="black"))
}
plots <- lapply(1:5,make_panel)
layout <- c(area(t=1,l=1,b=1,r=2),area(t=1,l=3,b=1,r=4),
  area(t=2,l=1,b=2,r=2),area(t=2,l=3,b=2,r=4),area(t=3,l=2,b=3,r=3))
p <- wrap_plots(plots,design=layout,guides="collect")+
  plot_annotation(
    title="Reactome pathway enrichment of upregulated genes",
    subtitle=paste0("MSigDB Reactome ",db_version,
      "; local hypergeometric test; BH correction; input FDR < 0.05 and log2FC >= 1"),
    theme=theme(plot.title=element_text(size=21,face="bold",hjust=.5,colour="black"),
      plot.subtitle=element_text(size=14,hjust=.5,colour="black"))) &
  theme(legend.position="bottom",legend.box="horizontal",
    legend.title=element_text(size=14,face="bold",colour="black"),
    legend.text=element_text(size=13,colour="black"))
ggsave(file.path(out,"figures","reactome_enrichment_all_clusters.png"),p,width=16,height=16,dpi=450,bg="white")
write.csv(top[,c("cluster_id","pathway","overlap","input_genes","gene_ratio","fold_enrichment","p_value","fdr","overlap_genes")],
          file.path(out,"enrichment","Reactome_top5_by_cluster.csv"),row.names=FALSE)
metadata <- data.frame(item=c("method","multiple_testing","gene_set_source","msigdb_version","collection","input_definition","universe_definition","min_pathway_size","max_pathway_size","min_overlap"),
 value=c("Local over-representation analysis (hypergeometric test)","Benjamini-Hochberg","MSigDB Reactome via msigdbr",db_version,"C2:CP:REACTOME","DEG FDR < 0.05 and log2FC >= 1","Genes present in all five one-vs-rest DEG result tables",10,500,3))
write.csv(metadata,file.path(out,"enrichment","enrichment_metadata.csv"),row.names=FALSE)
cat("Universe:",length(universe),"genes; MSigDB version:",db_version,"\n")
cat("Significant pathways per cluster:\n"); print(table(combined$cluster_id[combined$fdr<0.05]))
