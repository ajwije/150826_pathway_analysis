library(ggplot2)
library(reshape2)
first_up <- unique(unlist(kegg.gs[rownames(native_kegg_fc$greater)[1]]))
aminoacyl_tRNA <- norm_counts[row.names(norm_counts) %in% first_upset, ]


aminoacyl_format <- melt(Aminoacyl_tRNA)
head(aminoacyl_format)
colnames(aminoacyl_format) <- c("Gene_names",
                                "Genotype", 
                                "norm_counts" )

ggplot(data = aminoacyl_format, aes(x=Genotype, y=Gene_names, fill=norm_counts)) + 
        geom_tile()
