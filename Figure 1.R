'R version 4.1.3'
library(ggplot2)
library(cowplot)
library(ggpubr)

# read ratios exported from FlowJo X
ratios <- read.table("data/Ratios of T cell subpopulations.txt", header = T, row.names = 1)
rownames(ratios)
ratios <- ratios[1:30, ]
ratios$group <- unlist(lapply(rownames(ratios), function(x) strsplit(x, "_")[[1]][1]))

# melt data
data <- data.table::melt(ratios, id.vars = "group")

head(data)

# plot list
plist <- lapply(unique(data[, "variable"]), function(i) {
    subdat = data[data[, "variable"] == i,]
    #
    p = ggplot(subdat, aes(group, value, color = group)) +
        geom_violin(show.legend = "none", trim = F, draw_quantiles = .5)
    p = p +
        geom_jitter(width = 0.2, show.legend = 'none', alpha = .5) +
        xlab('') + ylab('') + ggtitle(i) + theme_classic() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 15),
            axis.text.y = element_text(size = 15)
        )
    #
    sig = compare_means(value ~ group, subdat)
    sig = sig[sig$p.signif != 'ns',]
    if (nrow(sig) > 0) {
        compare_list = lapply(1:nrow(sig), function(x) c(sig$group1[x], sig$group2[x]))
        p = p + stat_compare_means(comparisons = compare_list, method = 'wilcox.test',
                                    aes(label = ..p.signif..), hide.ns = T,
                                    vjust = 0.5, step.increase = 0.1)
    } else print(paste0('No sig result found in ', i))
    p
})

#
plot_grid(plotlist = plist, ncol = 5, nrow = 2)
