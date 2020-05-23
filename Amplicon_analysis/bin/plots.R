library("phyloseq"); #packageVersion("phyloseq")
library("ggplot2")    # for high quality graphics
library(argparser, quietly=TRUE)

p <- arg_parser("Plot in R")

p = add_argument(p, "-l", help="Enter list file [MiSeq_16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list]", type='character')
p = add_argument(p, "-t", help="Enter phylip tree file [MiSeq_16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.otu_modified.phylip.tre]", type='character')
p = add_argument(p, "-s", help="Enter shared file file [MiSeq_16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared]", type='character')
p = add_argument(p, "-tx", help="Enter taxonomy file [MiSeq_16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy]", type='character')
p = add_argument(p, "-m",  help="Enter taxonomy file [mouse.dpw.metadata]", type='character')
p = add_argument(p, "-o",  help="Enter path for output folder", type='character')


# Parse the command line arguments
argv <- parse_args(p)
ul_list = argv$l
tre = argv$t
shared_file = argv$s
taxonomy_file = argv$tx
g_metadata_file = argv$m

g_mothur_ps <- import_mothur(mothur_list_file =  ul_list,
                           mothur_group_file = NULL,
                           mothur_tree_file = tre,
                           cutoff = NULL,
                           mothur_shared_file = shared_file,
                           mothur_constaxonomy_file = taxonomy_file,
                           parseFunction = parse_taxonomy_default)


obj <- g_mothur_ps
g_metadata_file = argv$m
metadata=read.delim(g_metadata_file, header=T)
rownames(metadata) = metadata[,1]
sample_data(obj) <- metadata

col_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(obj)
rank_colnames <- c(Rank1 = col_names[1], Rank2 = col_names[2], Rank3 = col_names[3], Rank4 = col_names[4], Rank5 = col_names[5], Rank6 = col_names[6], Rank7 = col_names[7])
colnames(tax_table(obj)) <- rank_colnames[1:length(colnames(tax_table(obj)))]
any(taxa_sums(obj) == 0)
obj1 <- prune_taxa(taxa_sums(obj) > 0, obj)

#rarefaction
source("https://raw.githubusercontent.com/gauravsk/ranacapa/master/R/ggrare.R")
rare_plot = ggrare(obj1,
  step=20,
  se=TRUE,
  plot=TRUE,
  label="dpw",
  color="group")
ggsave(paste(argv$o,"/RarefactionCurve.png",sep=""), plot = rare_plot, width = 30, height = 20, units = "cm")

#boxplot
bp=plot_richness(obj1, x="group", color="dpw", measures="Chao1") +
  geom_boxplot(size=0.5,
  alpha=0.1,
  fill = "white",
  notch = FALSE,
  outlier.colour = "red",
  outlier.fill= "red",
  outlier.shape = 2,
  outlier.size=5,
  outlier.stroke=1,
  outlier.alpha=0.8)
ggsave(paste(argv$o,"/ChaosRichness.png",sep=""), plot = bp, width = 30, height = 20, units = "cm")

#tree
 tree=plot_tree(obj1,
   color="Phylum",
   label.tips="taxa_names",
   ladderize="left",
   plot.margin=0.6,
   size="abundance",
   min.abundance=20) +
   coord_polar(theta="y")
ggsave(paste(argv$o,"/PhyloTree.png",sep=""), plot = tree, width = 30, height = 20, units = "cm")


## Calculating realtive abundance of the data
## In the following steps data is transformed to relative abundance
   GPr = transform_sample_counts(obj1, function(x) x / sum(x) )
   GPfr = filter_taxa(GPr, function(x) var(x) > 1e-5, TRUE)


#creating stacked bar graphs without border
   stacked_bar=plot_bar(GPfr, fill = "Genus") +
   geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
ggsave(paste(argv$o,"/StackedBar.png",sep=""), plot = stacked_bar, width = 30, height = 20, units = "cm")
