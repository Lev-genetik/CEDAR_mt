# Compare enterotype predictions
# 2 enterotypes
# IL
comp.12.34.IL <- CompEntPred.pair(ent.Raes12.2e.sep.by.loc.IL[[2]],ent.Raes34.2e.sep.by.loc.IL[[2]])
comp.12.56.IL <- CompEntPred.pair(ent.Raes12.2e.sep.by.loc.IL[[2]],ent.Raes56.2e.sep.by.loc.IL[[2]])
comp.34.56.IL <- CompEntPred.pair(ent.Raes34.2e.sep.by.loc.IL[[2]],ent.Raes56.2e.sep.by.loc.IL[[2]])
# TR
comp.12.34.TR <- CompEntPred.pair(ent.Raes12.2e.sep.by.loc.TR[[2]],ent.Raes34.2e.sep.by.loc.TR[[2]])
comp.12.56.TR <- CompEntPred.pair(ent.Raes12.2e.sep.by.loc.TR[[2]],ent.Raes56.2e.sep.by.loc.TR[[2]])
comp.34.56.TR <- CompEntPred.pair(ent.Raes34.2e.sep.by.loc.TR[[2]],ent.Raes56.2e.sep.by.loc.TR[[2]])
# RE
comp.12.34.RE <- CompEntPred.pair(ent.Raes12.2e.sep.by.loc.RE[[2]],ent.Raes34.2e.sep.by.loc.RE[[2]])
comp.12.56.RE <- CompEntPred.pair(ent.Raes12.2e.sep.by.loc.RE[[2]],ent.Raes56.2e.sep.by.loc.RE[[2]])
comp.34.56.RE <- CompEntPred.pair(ent.Raes34.2e.sep.by.loc.RE[[2]],ent.Raes56.2e.sep.by.loc.RE[[2]])

comp.12.34.IL$Maximum_number_samples_clustering_together[,2]
comp.12.56.IL$Maximum_number_samples_clustering_together[,2]
comp.34.56.IL$Maximum_number_samples_clustering_together[,2]

comp.12.34.TR$Maximum_number_samples_clustering_together[,2]
comp.12.56.TR$Maximum_number_samples_clustering_together[,2]
comp.34.56.TR$Maximum_number_samples_clustering_together[,2]

comp.12.34.RE$Maximum_number_samples_clustering_together[,2]
comp.12.56.RE$Maximum_number_samples_clustering_together[,2]
comp.34.56.RE$Maximum_number_samples_clustering_together[,2]

# Number of samples
comp.12.34.IL$Maximum_number_samples_clustering_together[,3]
comp.12.56.IL$Maximum_number_samples_clustering_together[,3]
comp.34.56.IL$Maximum_number_samples_clustering_together[,3]

comp.12.34.TR$Maximum_number_samples_clustering_together[,3]
comp.12.56.TR$Maximum_number_samples_clustering_together[,3]
comp.34.56.TR$Maximum_number_samples_clustering_together[,3]

comp.12.34.RE$Maximum_number_samples_clustering_together[,3]
comp.12.56.RE$Maximum_number_samples_clustering_together[,3]
comp.34.56.RE$Maximum_number_samples_clustering_together[,3]


# cluster correspondences
comp.12.34.IL$Optimal_cluster_correspondences
comp.12.56.IL$Optimal_cluster_correspondences
comp.34.56.IL$Optimal_cluster_correspondences

comp.12.34.TR$Optimal_cluster_correspondences
comp.12.56.TR$Optimal_cluster_correspondences
comp.34.56.TR$Optimal_cluster_correspondences

comp.12.34.RE$Optimal_cluster_correspondences
comp.12.56.RE$Optimal_cluster_correspondences
comp.34.56.RE$Optimal_cluster_correspondences