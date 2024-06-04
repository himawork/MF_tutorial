'R version 4.1.3'

library(Seurat)
library(dynoPlot)
library(dynwrap)

source('FUN.R', encoding = 'UTF-8')

# Figure 2B
# T-cells substates
{
  load('data/ddd.Skin.mat.rda')
  ddd$Group = ddd$Sample
  ddd1 = ddd
  load('data/ddd.Blood.mat.rda')
  ddd$Group = ddd$Sample
  ddd2 = ddd
}

# cobmine
{
  # ddd = rbind(ddd1, ddd2)
  ddd = ddd1
  ddd$Group[which(grepl('CTCL', ddd$Group))] = 'CTCL'
  ddd$Group[which(grepl('HC1_B|HB', ddd$Group))] = 'HB'
  ddd$Group[which(grepl('Healthy|HC', ddd$Group))] = 'HC'
  ddd$Group[which(grepl('Blood|PBMC', ddd$Group))] = 'MS_PBMC'
  ddd$Group[which(grepl('MF|lesion|patch', ddd$Group))] = 'MF'
  ddd$Group[which(grepl('SS|SZ', ddd$Group))] = 'SS'

  table(ddd$Group)
}

head(ddd)

# plot list
plist <- lapply(unique(ddd[, "CellTypes"]), function(i) {
    subdat = ddd[ddd[, "CellTypes"] == i,]
    #
    p = ggplot(subdat, aes(Group, Count, color = Group)) +
        geom_violin(show.legend = "none", trim = F, draw_quantiles = .5)
    p = p +
        geom_jitter(width = 0.2, show.legend = 'none', alpha = .5) +
        xlab('') + ylab('') + ggtitle(i) + theme_classic() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 15),
            axis.text.y = element_text(size = 15)
        )
    #
    sig = compare_means(Count ~ Group, subdat)
    sig = sig[sig$p.signif != 'ns',]
    if (nrow(sig) > 0) {
        compare_list = lapply(1:nrow(sig), function(x) c(sig$group1[x], sig$group2[x]))
        p = p + stat_compare_means(comparisons = compare_list, method = 't.test',
                                    aes(label = ..p.signif..), hide.ns = T,
                                    vjust = 0.5, step.increase = 0.1)
    } else print(paste0('No sig result found in ', i))
    p
})

#
plot_grid(plotlist = plist, ncol = 3, nrow = 3)


# Figure 2C
load('data/cluster_slingshot_SEOB_model.rda')
levels(Idents(SEOB))
colnames(SEOB@meta.data)

#
plot_dimred(model, expression_source = dataset$expression,
            grouping = dataset$grouping, color_density = "grouping",
            label_milestones = T, plot_milestone_network = T)  +
  guides(color = 'none', fill = guide_legend(title = 'Clusters'))


# Figure 2D
load('data/Tcells_functional.cluster_markers.rda')
markerVocalno(markers = Tmarkers, topn = 5)


# Figure 2E
load('data/Biopsy_slingshot_SEOB_model.rda')
levels(Idents(SEOB))
colnames(SEOB@meta.data)

#
plot_dimred(model, expression_source = dataset$expression,
            grouping = dataset$grouping, color_density = "grouping",
            label_milestones = T, plot_milestone_network = T)  +
  guides(color = 'none', fill = guide_legend(title = 'Biopsies'))


# Figure 2F
# ... branch markers
model$milestone_network
branching_milestone <- model$milestone_network %>%
  group_by(from) %>% filter(n() > 1) %>% pull(from) %>% unique()

# top 30 features
branch_feature_importance <-
  calculate_branching_point_feature_importance(model, expression_source = dataset$expression,
                                               milestones_oi = branching_milestone[1])

branching_point_features <- branch_feature_importance %>%
  top_n(30, importance) %>% pull(feature_id)

# heatmap
plot_heatmap(
  model,
  expression_source = dataset$expression,
  features_oi = branching_point_features
) |> print()
