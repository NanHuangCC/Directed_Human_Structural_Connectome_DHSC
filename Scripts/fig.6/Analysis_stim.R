library(pheatmap)
library(ggplot2)

# load basic imformation
# load coor.file #
conn_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Data"
conn_coor <- read.csv(paste0(conn_path, "/Nodes_informations.csv"), header = T)
colnames(conn_coor)[c(3,4)] <- c("region_names", "Network")
region_name <- conn_coor$region_names
# define color list for brain networks
network_colors <- c("#CC3333","#FFFF99", "#339933",
                    "#FFFF00" ,"#66CC33","#E2593D",
                    "#6633CC", "#BC249F", "#DDA036",
                    "gray", "#993399","#6EB4E1","#3366CC")
names(network_colors) <- names(table(conn_coor$NETWORK))

#load trophic levels 
data_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Data"
trophic <- read.csv(paste0(data_path,"/Nodes_informations.csv"))
colnames(trophic)[c(3,4)] <- c("region_names", "Network")
trophic$Network <- factor(trophic$Network, levels = c("Visual1", "Visual2", "Orbito-Affective",
                                                      "Posterior-Multimodal", "Ventral-Multimodal", "Auditory","Language",
                                                      "Dorsal-Attention", "Cingulo-Opercular",
                                                      "subcortex", "Default", "Frontoparietal", "Somatomotor"))


# read time sires
ts_list <- readRDS("D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Scripts/fig.6/average_ts.RDS")
conn <- "directedsymm"
stim <- "V1"
ts <- ts_list[[conn]][[stim]]$bold$mean
dim(ts)

library(pheatmap)
library(viridis)
bk <- seq(0.0,0.1,0.001)

P <- ggplot() +
  #geom_density(aes(x=trophic$Trophic_level[ts[30,] > 0.02])) + 
  geom_density(aes(x=trophic$Trophic_level[ts[32,] > 0.02]), color = mako(10)[9]) +
  geom_density(aes(x=trophic$Trophic_level[ts[34,] > 0.02]), color = mako(10)[8]) +
  geom_density(aes(x=trophic$Trophic_level[ts[36,] > 0.02]), color = mako(10)[7]) +
  geom_density(aes(x=trophic$Trophic_level[ts[38,] > 0.02]), color = mako(10)[6]) +
  geom_density(aes(x=trophic$Trophic_level[ts[40,] > 0.02]), color = mako(10)[5]) +
  geom_density(aes(x=trophic$Trophic_level[ts[42,] > 0.02]), color = mako(10)[4]) +
  geom_segment(data=trophic, aes(x=Trophic_level, xend = Trophic_level, y = -0.2,yend = -0.1, color=as.character(Network)), size =1) +
  scale_color_manual(values = network_colors) +
  theme_classic() +
  scale_x_continuous(limits = c(0,2)) + scale_y_continuous(limits = c(-0.2,4),expand = c(0,0))
P
save_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/CCEP_stim/synthetic/20260330/results"
file <- paste0(save_path,"/","active_distribution",conn,"_",stim,".pdf")
pdf(file,width=8,height=2)
print(P)
dev.off()

# classify time-sires pattern
# library(uwot)
# library(dbscan)
# start_time <- 28
# end_time <- 80
# conn_list <- names(ts_list)
# stim <- "V1"
# conn <- conn_list[1]
# data_tmp <- ts_list[[conn]][[stim]]$bold$mean
# data_tmp <- t(data_tmp)
# data_tmp [is.na(data_tmp )] <- 0
# data_tmp <- data_tmp[,start_time:end_time]
# data <- data_tmp[1:2,]
# for (conn in conn_list) {
#   for (stim in c("V1","A1")) {
#    data_tmp <- ts_list[[conn]][[stim]]$bold$mean
#   data_tmp <- t(data_tmp)
#    data_tmp [is.na(data_tmp )] <- 0
#    data_tmp <- data_tmp[,start_time:end_time]
#    rownames(data_tmp) <- paste0(conn,"_",stim,"_",rownames(data_tmp))
#    data <- rbind(data, data_tmp)
#  }
# }
#data <- data[c(-1,-2),]
#
# conn_id <- c()
# stim_id <- c()
# for(i in 1:nrow(data)){
#   conn_id <- c(conn_id, strsplit(rownames(data)[i],"_")[[1]][1])
#  stim_id <- c(stim_id, strsplit(rownames(data)[i],"_")[[1]][2])
# }
#
# ncomponents <- 3
# umap_data <- umap(data, 
#                   n_neighbors = 100, # 5 ~ 15 for local structure / 50~200 for global sturcture
#                   min_dist = 0.05,
#                  n_components = ncomponents,
#                   metric = "cosin")
#
# colnames(umap_data) <- paste0("components", 1:ncomponents)
# umap_full <- as.data.frame(umap_data)
# Full_umap_clus <- hdbscan(umap_data, minPts=100, gen_simplified_tree = F, gen_hdbscan_tree = F)
# labels_cosin <- as.character(Full_umap_clus$cluster)
#
# ncomponents <- 3
# umap_data <- umap(data, 
#                   n_neighbors = 50, # 5 ~ 15 for local structure / 50~200 for global sturcture
#                   min_dist = 0.05,
#                   n_components = ncomponents,
#                   metric = "euclidean")
#
# colnames(umap_data) <- paste0("components", 1:ncomponents)
# umap_full <- as.data.frame(umap_data)
# Full_umap_clus <- hdbscan(umap_data, minPts=10, gen_simplified_tree = F, gen_hdbscan_tree = F)
# labels_euc <- as.character(Full_umap_clus$cluster)
#
# umap_full$labels_euc <- labels_euc
# umap_full$labels_cosin <- labels_cosin
#
# umap_full$labels_intergrate <- 0 
# umap_full$labels_intergrate[umap_full$labels_cosin==1] <- 1
# umap_full$labels_intergrate[umap_full$labels_cosin==2] <- 1
# umap_full$labels_intergrate[umap_full$labels_cosin==3] <- 2
# umap_full$labels_intergrate[umap_full$labels_euc==6] <- 3
# umap_full$labels_intergrate[umap_full$labels_euc==9] <- 4
# umap_full$labels_intergrate[umap_full$labels_euc==1] <- 4
# umap_full$labels_intergrate <- as.character(umap_full$labels_intergrate)
# labels_intergrate <- umap_full$labels_intergrate
# 
# umap_full$labels_intergrate <- labels_intergrate
# umap_full$labels_intergrate[umap_full$components1 < -3] <- 1 
# umap_full$conn <- conn_id
# umap_full$stim  <- stim_id 
# nrow(umap_full)

# for (conn in conn_list) {
#   for (stim in c("V1", "A1")) {
#     label <- as.numeric(umap_full[umap_full$stim == stim & umap_full$conn == conn,"labels_intergrate"])
#     ts <- ts_list[[conn]][[stim]]$bold$mean[start_time:end_time,]
#     if(stim == "V1"){label[c(1,213 +1)] <- 1}
#     if(stim == "A1"){label[c(24,213 + 24)] <- 1}
#    amp_init <- c()
#     amp_mant <- c()
#   for (i in 1:426) {
#      region_label <- label[i]
#      if(region_label == 1){
#        ts_init <- ts[1:14, i]
#        amp_init_r <- max(ts_init )
#      } else {
#        amp_init_r <- 0
#      }
#      ts_mant_r <- ts[28:42, i]
#      amp_mant_r <- mean(ts_mant_r)
#      amp_init <- c(amp_init,amp_init_r)
#      amp_mant <- c(amp_mant,amp_mant_r)
#   }
#    feature_df[,paste0(conn, "_",stim, "_Info")] <- label
#    feature_df[,paste0(conn, "_",stim, "_initial_amp")] <- amp_init
#    feature_df[,paste0(conn, "_",stim, "_maintain_amp")] <- amp_mant
#  }
#}
#feature_df <- as.data.frame(feature_df)
#feature_df$network <- conn_coor$NETWORK
#feature_df$trophic <- trophic$Trophic.level

# input coordinate from 
#1. Rottschy, C. et al. Modelling neural correlates of working memory: A coordinate-based meta-analysis. NeuroImage 60, 830–846 (2012).
#WM_df <- data.frame(
# Macroanatomical_location = c(
#    "Left anterior insula",
#    "Left inferior frontal gyrus pars opercularis",
#    "Left caudal lateral prefrontal cortex",
#    "Left rostral lateral prefrontal cortex",
#    "Right anterior insula",
#    "Right inferior frontal gyrus pars triangularis",
#    "Right caudal lateral prefrontal cortex",
#    "Right rostral lateral prefrontal cortex",
#    "Posterior medial frontal cortex",
#    "Left posterior superior frontal gyrus",
#    "Right posterior superior frontal gyrus",
#    "Left intraparietal sulcus",
#    "Left superior parietal lobule/intraparietal sulcus",
#    "Left posterior superior parietal lobule",
#    "Right intraparietal sulcus",
#    "Right intraparietal sulcus",
#    "Right posterior superior parietal lobule",
#    "Left thalamus (prefrontal/temporal)",
#    "Left basal ganglia (caudate)",
#    "Left basal ganglia (putamen)",
#    "Left basal ganglia (pallidum)",
#    "Right thalamus (prefrontal/temporal)",
#    "Right thalamus (prefrontal)",
#    "Left cerebellum / left fusiform gyrus",
#    "Right cerebellum / right fusiform gyrus"
#  ),
#  
#  Cytoarchitectonic_location = c(
#   NA,
#    "Area 44/45",
#    NA,
#    NA,
#    NA,
#    "Area 44/45",
#    NA,
#   NA,
#    NA,
#    NA,
#    NA,
#    "hIP3, hIP2, hIP1",
#    "7PC, 7A, hIP3",
#    "7A",
#    "hIP2/hIP3/hIP1",
#    "hIP3/hIP1",
#    "7P/7A",
#    NA,
#    NA,
#   NA,
#    NA,
#   NA,
#    NA,
#    "Lobules VI / VIIa Crus I",
#    "Lobule VI / VIIa Crus I"
#  ),
#  
#  x = c(
#    -32, -48, -46, -38, 36, 50, 44, 38, 2, -28, 30,
#    -42, -34, -24, 42, 32, 16, -12, -16, -18, -16,
#    12, 8, -34, 32
#  ),
#  
#  y = c(
#    22, 10, 26, 50, 22, 14, 34, 54, 18, 0, 2,
#    -42, -52, -66, -44, -58, -66, -12, 2, 4, 0,
#    -10, -18, -66, -64
#  ),
#  
#  z = c(
#    -2, 26, 24, 10, -6, 24, 32, 6, 48, 56, 56,
#    46, 48, 54, 44, 48, 56, 12, 14, 6, 2,
#    10, 4, -20, -18
#  ),
#  
#  z_score = c(
#    8.30, 8.25, 8.23, 6.33, 8.26, 8.22, 7.80, 4.28,
#    8.29, 7.60, 7.06, 8.25, 8.25, 6.64, 8.22, 8.22,
#    5.17, 5.47, 5.30, 6.33, 5.40, 4.08, 3.64, 4.94, 5.44
#  ),
#  
#  stringsAsFactors = FALSE
#)

# mapping WM score to HCP_MMP_nodes
# effect_score <- matrix(0, nrow = nrow(conn_coor), ncol = nrow(WM_df), dimnames = list(conn_coor$HCP_MMP_name, WM_df$Macroanatomical_location))
# for (WM_region in 1:nrow(WM_df)) {
#  for (HCP_region in 1:nrow(conn_coor)) {
#    coor_WM <- WM_df[WM_region,c("x", "y", "z")]
#    coor_HCP <- conn_coor[HCP_region,c("x", "y", "z")]
#    coor_diff <-  coor_WM - coor_HCP
#    distance <-  sum(coor_diff^2)^0.5
#    distance <- max(c(1,distance)) # values not larger than orignal value
#    effect <- WM_df[WM_region,c("z_score")]/(distance)
#    effect_score[HCP_region,WM_region] <- effect
#  }
#}
#effect_score <-  effect_score / apply(effect_score, 2, FUN = sum) 
#effect_score_sum <-  apply(effect_score, 1, FUN = sum) 
#effect_score_sum <-scale(effect_score_sum)
#feature_df$effect_score_sum <- effect_score_sum

feature_df <- read.csv(paste0("D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Scripts/fig.6/signal_features.csv"), row.names = 1)

hist(feature_df$directed_V1_maintain_amp, breaks = 50)
feature_df$network <- factor(trophic$Network, levels = c("Visual1", "Visual2", "Orbito-Affective",
                                                      "Posterior-Multimodal", "Ventral-Multimodal", "Auditory","Language",
                                                      "Dorsal-Attention", "Cingulo-Opercular",
                                                      "subcortex", "Default", "Frontoparietal", "Somatomotor"))


library(ggplot2)
library(viridis)

# Fig.6f -------------------------------------------------------------------------------------------------------------- 
data_ini_V1 <- feature_df[-c(1,214),]
data_ini_V1 <- data_ini_V1[, c(grep("initial_amp",colnames(feature_df)), grep("maintain_amp",colnames(feature_df)))]
data_mean <- aggregate(data_ini_V1, by = list(feature_df$trophic_group[-c(1,214)]), mean)
data_sd <- aggregate(data_ini_V1, by = list(feature_df$trophic_group[-c(1,214)]), sd)
data_mean <- data.frame(data_mean, data_sd)
P <- ggplot() + 
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directed_V1_initial_amp + directed_V1_initial_amp.1/(42^0.5), ymin = directed_V1_initial_amp - directed_V1_initial_amp.1/(42^0.5)),fill = "#E39745", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = directed_V1_initial_amp), color = "#E39745") + 
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directed_V1_maintain_amp + directed_V1_maintain_amp.1/(42^0.5), ymin = directed_V1_maintain_amp - directed_V1_maintain_amp.1/(42^0.5)), fill = "#E39745", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = directed_V1_maintain_amp), color = "#E39745") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directedsymm_V1_initial_amp + directedsymm_V1_initial_amp.1/(42^0.5), ymin = directedsymm_V1_initial_amp - directedsymm_V1_initial_amp.1/(42^0.5)), fill  = "#407277", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = directedsymm_V1_initial_amp), color = "#407277") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directedsymm_V1_maintain_amp + directedsymm_V1_maintain_amp.1/(42^0.5), ymin = directedsymm_V1_maintain_amp - directedsymm_V1_maintain_amp.1/(42^0.5)), fill  = "#407277", alpha=0.2) +
  geom_line(data = data_mean, aes(x=Group.1, y = directedsymm_V1_maintain_amp), color = "#407277") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = threshold_V1_initial_amp + threshold_V1_initial_amp.1/(42^0.5), ymin = threshold_V1_initial_amp - threshold_V1_initial_amp.1/(42^0.5)), fill  = "gray", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = threshold_V1_initial_amp), color = "gray") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = threshold_V1_maintain_amp + threshold_V1_maintain_amp.1/(42^0.5), ymin = threshold_V1_maintain_amp - threshold_V1_maintain_amp.1/(42^0.5)), fill  = "gray", alpha=0.2)  +
  geom_line(data = data_mean, aes(x=Group.1, y = threshold_V1_maintain_amp), color = "gray") +
  theme_classic() + scale_x_continuous(limits = c(1,10), expand = c(0,0), breaks = c(1:10)) + scale_y_continuous(limits = c(-0.001,0.08), expand = c(0,0))
P

data_ini_A1 <- feature_df[-c(24,24+213),]
data_ini_A1 <- data_ini_A1[, c(grep("initial_amp",colnames(feature_df)), grep("maintain_amp",colnames(feature_df)))]
data_mean <- aggregate(data_ini_A1, by = list(feature_df$trophic_group[-c(24,24+213)]), mean)
data_sd <- aggregate(data_ini_A1, by = list(feature_df$trophic_group[-c(24,24+213)]), sd)
data_mean <- data.frame(data_mean, data_sd)
P <- ggplot() + 
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directed_A1_initial_amp + directed_A1_initial_amp.1/(42^0.5), ymin = directed_A1_initial_amp - directed_A1_initial_amp.1/(42^0.5)),fill = "#E39745", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = directed_A1_initial_amp), color = "#E39745") + 
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directed_A1_maintain_amp + directed_A1_maintain_amp.1/(42^0.5), ymin = directed_A1_maintain_amp - directed_A1_maintain_amp.1/(42^0.5)), fill = "#E39745", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = directed_A1_maintain_amp), color = "#E39745") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directedsymm_A1_initial_amp + directedsymm_A1_initial_amp.1/(42^0.5), ymin = directedsymm_A1_initial_amp - directedsymm_A1_initial_amp.1/(42^0.5)), fill  = "#407277", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = directedsymm_A1_initial_amp), color = "#407277") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = directedsymm_A1_maintain_amp + directedsymm_A1_maintain_amp.1/(42^0.5), ymin = directedsymm_A1_maintain_amp - directedsymm_A1_maintain_amp.1/(42^0.5)), fill  = "#407277", alpha=0.2) +
  geom_line(data = data_mean, aes(x=Group.1, y = directedsymm_A1_maintain_amp), color = "#407277") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = threshold_A1_initial_amp + threshold_A1_initial_amp.1/(42^0.5), ymin = threshold_A1_initial_amp - threshold_A1_initial_amp.1/(42^0.5)), fill  = "gray", alpha=0.2) + 
  geom_line(data = data_mean, aes(x=Group.1, y = threshold_A1_initial_amp), color = "gray") +
  geom_ribbon(data = data_mean, aes(x=Group.1, ymax = threshold_A1_maintain_amp + threshold_A1_maintain_amp.1/(42^0.5), ymin = threshold_A1_maintain_amp - threshold_A1_maintain_amp.1/(42^0.5)), fill  = "gray", alpha=0.2)  +
  geom_line(data = data_mean, aes(x=Group.1, y = threshold_A1_maintain_amp), color = "gray") +
  theme_classic() + scale_x_continuous(limits = c(1,10), expand = c(0,0), breaks = c(1:10)) + scale_y_continuous(limits = c(-0.001,0.08), expand = c(0,0))
P
#------------------------------------------------------------------------------------------

# fig.6h ----------------------------------------------------------------------------------
df <- feature_df[c(,)]
#df <- df[df$effect_score_sum <3,]
P <- ggplot() +
  geom_point(data =df, aes(x=effect_score_sum, y=(directed_V1_maintain_amp  + directed_A1_maintain_amp)/2, color = trophic)) +
  scale_color_gradientn(colors =  c(viridis(10)[1],viridis(10), viridis(10)[10])) +
  geom_smooth(data = df, aes(x=effect_score_sum, y=(directed_V1_maintain_amp  + directed_A1_maintain_amp)/2), method = "lm", color = "black") + 
  theme_bw() + scale_x_continuous(limits = c(-2.0, 3.0), expand = c(0.01,0.01)) + 
  scale_y_continuous(limits = c(0, 0.09), expand = c(0.002,0.002), breaks = c(0,0.03,0.06,0.09)) +
  coord_fixed(5/(0.09*2))
P
cor.test(df$effect_score_sum, (df$directed_V1_maintain_amp  + df$directed_A1_maintain_amp)/2)

P <- ggplot() +
  geom_point(data =df, aes(x=effect_score_sum, y=(threshold_V1_maintain_amp  + threshold_A1_maintain_amp)/2, color = trophic)) +
  scale_color_gradientn(colors =  c(viridis(10)[1],viridis(10), viridis(10)[10])) +
  geom_smooth(data = df, aes(x=effect_score_sum, y=(threshold_V1_maintain_amp  + threshold_A1_maintain_amp)/2), method = "lm", color = "black") + 
  theme_bw() + scale_x_continuous(limits = c(-2.0, 3.0), expand = c(0.01,0.01)) + 
  scale_y_continuous(limits = c(0, 0.03), expand = c(0.001,0.001), breaks = c(0,0.01,0.02,0.03)) +
  coord_fixed(5/(0.03*2))
P
cor.test(df$effect_score_sum, (df$threshold_V1_maintain_amp  + df$threshold_A1_maintain_amp)/2)

feature_df$region_name <- row.names(feature_df)
# plot different path
conn_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Data"
directional_conn <- read.csv(paste0(conn_path, "/connections/DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv"), header = F, check.names = F)
directional_conn <- as.matrix(directional_conn)
dimnames(directional_conn) <- list(region_name, region_name)
feature_df$V1_conn <- as.numeric(directional_conn["L_V1",]) + as.numeric(directional_conn["R_V1",])
df <- feature_df[-c(214,1),]

# fig.6e ----------------------------------------------------------------------------------------------
feature_df$V1_conn <- as.numeric(directional_conn["L_V1",]) + as.numeric(directional_conn["R_V1",])
df <- feature_df[-c(1,214),]
df <- df[df$V1_conn > 0,]
P <- ggplot() +
  geom_smooth(data = df, aes(x=V1_conn, y =directed_V1_initial_amp), method = "lm", color = "#103380") + 
  geom_point(data = df[,], aes(x=V1_conn, y =directed_V1_initial_amp), color = "#103380") +
  geom_smooth(data = df, aes(x=V1_conn, y =directed_V1_maintain_amp), method = "lm", color = "#E3C42B") +
  geom_point(data = df[,], aes(x=V1_conn, y =directed_V1_maintain_amp), color = "#E3C42B") + 
  scale_y_continuous(limits = c(-0.02, 0.12), expand = c(0.0,0.0), breaks = c(0,0.04,0.08,0.12)) +
  scale_x_continuous(limits = c(0, 8), expand = c(0.0,0.0), breaks = c(0,2,4,6,8)) +
  coord_fixed(4/0.14) +
  theme_bw()
P
cor.test(df$V1_conn, df$directed_V1_maintain_amp)
feature_df$A1_conn <- as.numeric(directional_conn["L_A1",]) + as.numeric(directional_conn["R_A1",])
df <- feature_df[-c(24,24 +213),]
df <- df[df$A1_conn > 0,]
P <- ggplot() +
  geom_smooth(data = df, aes(x=A1_conn, y =directed_A1_initial_amp), method = "lm", color = "#E7211A") + 
  geom_point(data = df[,], aes(x=A1_conn, y =directed_A1_initial_amp), color = "#E7211A") +
  geom_smooth(data = df, aes(x=A1_conn, y =directed_A1_maintain_amp), method = "lm", color = "#E3C42B") +
  geom_point(data = df[,], aes(x=A1_conn, y =directed_A1_maintain_amp), color = "#E3C42B") + 
  scale_y_continuous(limits = c(-0.02, 0.08), expand = c(0.0,0.0), breaks = c(0,0.04,0.08,0.12)) +
  scale_x_continuous(limits = c(0, 3), expand = c(0.0,0.0), breaks = c(0,1,2,3)) +
  coord_fixed(1.5/0.10) +
  theme_bw()
P
cor.test(df$A1_conn, df$directed_A1_maintain_amp)
# ------------------------------------------------------------------------------------------------------

# fig.6d -----------------------------------------------------------------------------------------------
# venn plot foor initual 
library(UpSetR)
# same list as above
my_list <- list(
  D_V1 = feature_df[feature_df$directed_V1_Info == 1,"idx"],
  D_A1 = feature_df[feature_df$directed_A1_Info == 1,"idx"],
  DS_V1 = feature_df[feature_df$directedsymm_V1_Info == 1,"idx"],
  DS_A1 = feature_df[feature_df$directedsymm_A1_Info == 1,"idx"],
  T_V1 = feature_df[feature_df$threshold_V1_Info == 1,"idx"],
  T_A1 = feature_df[feature_df$threshold_A1_Info == 1,"idx"]
)

# Convert list to binary matrix
data <- fromList(my_list)
# Plot
upset(data, nsets = 6, keep.order = T, order.by = c("freq", "degree"))



# fig.6j----------------------------------------------------------------------------------------------
# node plot
library(igraph)
library(viridis)
df <- feature_df[,] # left hemisphere
# select regions response to stimulus
region_select <- (df$directed_V1_Info != 2 & df$directed_V1_Info != 0) | (df$directed_A1_Info != 2 & df$directed_A1_Info != 0)
directed_matrix <- directional_conn[,]
directed_matrix <- directed_matrix[region_select,region_select]
directed_matrix <- as.matrix(directed_matrix)
directed_matrix[directed_matrix < 2.0] <- 0 # view strong connections make plot readable
region_select <- (apply(directed_matrix, 1, sum) > 0) | (apply(directed_matrix, 1, sum) > 0)
directed_matrix <- directed_matrix[region_select,region_select]

# 创建有向图
df <- df[rownames(directed_matrix),]
digraph <- graph_from_adjacency_matrix(directed_matrix, mode = "directed")

# re-scale colors 
df$directed_V1_initial_amp[df$directed_V1_initial_amp > 0.1] <- 0.1
df$directed_V1_initial_amp <- df$directed_V1_initial_amp/0.1
df$directed_A1_initial_amp[df$directed_A1_initial_amp > 0.02] <- 0.03
df$directed_A1_initial_amp <- df$directed_A1_initial_amp/0.03
df$directed_A1_maintain_amp[df$directed_A1_maintain_amp > 0.1] <- 0.1
df$directed_A1_maintain_amp <- df$directed_A1_maintain_amp/0.1
df$directed_V1_maintain_amp[df$directed_V1_maintain_amp > 0.1] <- 0.1
df$directed_V1_maintain_amp <- df$directed_V1_maintain_amp/0.1

# amp based color
# 3 type infomation
colors <- rep("gray", nrow(directed_matrix))
for (i in 1:nrow(df)) {
  if(df[i,"directed_V1_Info"] == 1){
    colors[i] <- rgb(16/255,51/255,128/255,df[i,"directed_V1_initial_amp"])
  }
  if(df[i,"directed_A1_Info"] == 1){
    colors[i] <- rgb(184/255,28/255,37/255,df[i,"directed_A1_initial_amp"])
  }
  if(df[i,"directed_V1_Info"] == 1 & df[i,"directed_A1_Info"] == 1){
    colors[i] <- rgb(102/255,51/255,153/255,(df[i,"directed_A1_initial_amp"] + df[i,"directed_V1_initial_amp"])/2)
  }
  if(df[i,"directed_V1_Info"] == 4 & df[i,"directed_A1_Info"] == 4){
    colors[i] <- rgb(227/255,196/255,43/255,(df[i,"directed_A1_maintain_amp"] + df[i,"directed_V1_maintain_amp"])/2)
  }
}



# trophic color
colors <- rep("gray", nrow(directed_matrix))
for (i in 1:nrow(df)) {
  colors[i] <- viridis(10)[df[i,"trophic_group"]]
}



# network color
colors <- rep("gray", nrow(directed_matrix))
for (i in 1:nrow(df)) {
  colors[i] <- network_colors[as.character(df[i,"network"])]
}


colors <- rep("gray", nrow(directed_matrix))
for (i in 1:nrow(df)) {
  if(((df[i,"directed_A1_maintain_amp"] + df[i,"directed_V1_maintain_amp"])/2) > 0){
    colors[i] <- rgb(227/255,196/255,43/255,(df[i,"directed_A1_maintain_amp"] + df[i,"directed_V1_maintain_amp"])/2)
  }
  if(((df[i,"directed_A1_maintain_amp"] + df[i,"directed_V1_maintain_amp"])/2) <= 0){
    colors[i] <- "gray"
  }
}


set.seed(5)
plot(digraph,
     vertex.size = 5,
     vertex.color = colors,
     edge.arrow.size = 0,
     vertex.label = "",
     layout = layout_with_fr,
     vertex.frame.color ="black")

plot(digraph,
     vertex.size = 5,
     vertex.color = colors,
     edge.arrow.size = 0,
     vertex.label.cex = 0.2,
     layout = layout_with_fr,
     vertex.frame.color ="black")


# fig.6g-----------------------------------------------------------------
# plot correlations
df <- feature_df[,]
P <- ggplot() +
  geom_point(data =df, aes(x=directed_V1_maintain_amp, y=directed_A1_maintain_amp, color = trophic)) +
  scale_color_gradientn(colors =  c(viridis(10)[1],viridis(10), viridis(10)[10])) +
  geom_smooth(data = df, aes(x=directed_V1_maintain_amp, y=directed_A1_maintain_amp), method = "lm", color = "black") + 
  theme_bw() + scale_x_continuous(limits = c(0, 0.08), expand = c(0.002,0.002)) + 
  scale_y_continuous(limits = c(0, 0.08), expand = c(0.002,0.002), breaks = c(0,0.02,0.04,0.06,0.08)) +
  geom_text(data = df[df$directed_V1_maintain_amp > 0.075 | df$directed_A1_maintain_amp > 0.075 ,], aes(x=directed_V1_maintain_amp, y =directed_A1_maintain_amp, label = region_name)) +
  coord_fixed(1/2)
P

P <- ggplot() +
  geom_point(data =df, aes(x=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, y=threshold_V1_maintain_amp, color = trophic)) +
  scale_color_gradientn(colors =  c(viridis(10)[1],viridis(10), viridis(10)[10])) +
  geom_smooth(data = df, aes(x=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, y=threshold_V1_maintain_amp), method = "lm", color = "black") + 
  theme_bw() + scale_x_continuous(limits = c(0, 0.08), expand = c(0.002,0.002)) + 
  scale_y_continuous(limits = c(0, 0.08), expand = c(0.002,0.002), breaks = c(0,0.02,0.04,0.06,0.08)) +
  coord_fixed(1/2)
P

# fig.5i -------------------------------------------------------------
df$in_degree <- apply(directional_conn, 2, FUN = sum)
P <- ggplot() +
  geom_smooth(data = df, aes(x=in_degree, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2), method = "loess", color = "black") + 
  geom_point(data =df, aes(x=in_degree, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, color = trophic)) +
  scale_color_gradientn(colors =  c(viridis(10)[1],viridis(10), viridis(10)[10])) +
  theme_bw() + scale_x_continuous(limits = c(0, 450), expand = c(0.0,0.0)) + 
  scale_y_continuous(limits = c(0, 0.08), expand = c(0.002,0.002), breaks = c(0,0.02,0.04,0.06,0.08)) +
  geom_text(data = df[df$in_degree >300,], aes(x=in_degree, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, label = region_name)) +
  coord_fixed(450/0.08) 
P

df$out_degree <- apply(directional_conn, 1, FUN = sum)
P <- ggplot() +
  geom_smooth(data = df, aes(x=out_degree, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2), method = "loess", color = "black") + 
  geom_point(data =df, aes(x=out_degree, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, color = trophic)) +
  scale_color_gradientn(colors =  c(viridis(10)[1],viridis(10), viridis(10)[10])) +
  theme_bw() + scale_x_continuous(limits = c(0, 350), expand = c(0.0,0.0)) + 
  scale_y_continuous(limits = c(0, 0.08), expand = c(0.002,0.002), breaks = c(0,0.02,0.04,0.06,0.08)) +
  geom_text(data = df[df$out_degree >250,], aes(x=out_degree, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, label = region_name)) +
  coord_fixed(350/0.08) 
P

P <- ggplot() +
  geom_smooth(data = df, aes(x=trophic, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2), method = "loess", color = "black") + 
  geom_point(data =df, aes(x=trophic, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, color = trophic)) +
  scale_color_gradientn(colors =  c(viridis(10)[1],viridis(10), viridis(10)[10])) +
  theme_bw() + scale_x_continuous(limits = c(0, 1.8), expand = c(0.0,0.0)) + 
  scale_y_continuous(limits = c(0, 0.08), expand = c(0.002,0.002), breaks = c(0,0.02,0.04,0.06,0.08)) +
  geom_text(data = df[df$trophic >1.5,], aes(x=trophic, y=(directed_A1_maintain_amp + directed_V1_maintain_amp)/2, label = region_name)) +
  coord_fixed(1.8/0.08) 
P
save_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/CCEP_stim/synthetic/20260330/results"
file <- paste0(save_path,"/VA_trophic.pdf")
pdf(file,width=8,height=8)
print(P)
dev.off()

feature_df$V1_conn <- as.numeric(directional_conn["L_V1",]) + as.numeric(directional_conn["R_V1",])
df <- feature_df[-c(214,1),]
df <- df[df$V1_conn > 0,]
P <- ggplot() +
  geom_smooth(data = df, aes(x=V1_conn, y =directed_V1_initial_amp), method = "lm", color = "#103380") + 
  geom_point(data = df[,], aes(x=V1_conn, y =directed_V1_initial_amp), color = "#103380") +
  geom_smooth(data = df, aes(x=V1_conn, y =directed_V1_maintain_amp), method = "lm", color = "#E3C42B") +
  geom_point(data = df[,], aes(x=V1_conn, y =directed_V1_maintain_amp), color = "#E3C42B") + 
  scale_y_continuous(limits = c(-0.02, 0.12), expand = c(0.0,0.0), breaks = c(0,0.04,0.08,0.12)) +
  scale_x_continuous(limits = c(0, 8), expand = c(0.0,0.0), breaks = c(0,2,4,6,8)) +
  coord_fixed(4/0.14) +
  theme_bw()
P

#------------------------------------------------------------------------------------------------


# anova trophic group & conn
stim <- c("V1", "A1")[1]
condition <- c("initial", "maintain")[1]
model_list <- list()
for (stim in c("V1", "A1")) {
  for (condition in c("initial", "maintain")) {
    conn <- conn_list[1]
    data_tmp1 <- feature_df[,c(paste0(conn,"_",stim,"_",condition,"_amp") ,"trophic_group")]
    data_tmp1$conn <- conn
    colnames(data_tmp1) <- c("values", "trophic", "conn")
    
    conn <- conn_list[2]
    data_tmp2 <- feature_df[,c(paste0(conn,"_",stim,"_",condition,"_amp") ,"trophic_group")]
    data_tmp2$conn <- conn
    colnames(data_tmp2) <- c("values", "trophic", "conn")
    
    conn <- conn_list[3]
    data_tmp3 <- feature_df[,c(paste0(conn,"_",stim,"_",condition,"_amp") ,"trophic_group")]
    data_tmp3$conn <- conn
    colnames(data_tmp3) <- c("values", "trophic", "conn")
    
    data_anova <- rbind(data_tmp1,rbind(data_tmp2,data_tmp3))
    
    model_full <- aov(values ~ trophic * conn, data = data_anova)
    model_list[[stim]][[condition]] <- model_full
  }
}


# same list as above
my_list <- list(
  D_V1 = row.names(feature_df)[feature_df$directed_V1_Info == 1],
  D_A1 = row.names(feature_df)[feature_df$directed_A1_Info == 1],
  DS_V1 = row.names(feature_df)[feature_df$directedsymm_V1_Info == 1],
  DS_A1 = row.names(feature_df)[feature_df$directedsymm_A1_Info == 1],
  T_V1 = row.names(feature_df)[feature_df$threshold_V1_Info == 1],
  T_A1 = row.names(feature_df)[feature_df$threshold_A1_Info == 1]
)






# fig.6d----------------------------------------------------------------------
# Convert list to binary matrix
data <- fromList(my_list)


list_v <- list(
  D_V1 = row.names(feature_df)[feature_df$directed_V1_Info == 1],
  DS_V1 = row.names(feature_df)[feature_df$directedsymm_V1_Info == 1],
  T_V1 = row.names(feature_df)[feature_df$threshold_V1_Info == 1]
)

list_a <- list(
  D_A1 = row.names(feature_df)[feature_df$directed_A1_Info == 1],
  DS_A1 = row.names(feature_df)[feature_df$directedsymm_A1_Info == 1],
  T_A1 = row.names(feature_df)[feature_df$threshold_A1_Info == 1]
)

# Function to compute all intersections
get_intersections <- function(sets) {
  n <- length(sets)
  set_names <- names(sets)
  results <- list()
  # loop over all combinations (excluding empty set)
  for (k in 1:n) {
    combs <- combn(set_names, k, simplify = FALSE)
    
    for (cmb in combs) {
      # intersection
      inter <- Reduce(intersect, sets[cmb])
      
      # subtract elements belonging to other sets (exclusive region)
      others <- setdiff(set_names, cmb)
      if (length(others) > 0) {
        other_union <- Reduce(union, sets[others])
        inter <- setdiff(inter, other_union)
      }
      
      name <- paste(cmb, collapse = "&")
      results[[name]] <- length(inter)
    }
  }
  unlist(results)
}

fit_data <- get_intersections (list_v)
fit_data <- fit_data[fit_data > 0]
fit <- euler(fit_data)
plot(fit)
get_region_elements <- function(sets) {
  n <- length(sets)
  set_names <- names(sets)
  
  region_list <- list()
  
  for (k in 1:n) {
    combs <- combn(set_names, k, simplify = FALSE)
    
    for (cmb in combs) {
      # elements in all selected sets
      inter <- Reduce(intersect, sets[cmb])
      
      # remove elements that also appear in other sets
      others <- setdiff(set_names, cmb)
      if (length(others) > 0) {
        other_union <- Reduce(union, sets[others])
        inter <- setdiff(inter, other_union)
      }
      
      name <- paste(cmb, collapse = "&")
      region_list[[name]] <- inter
    }
  }
  
  return(region_list)
}

regions <- get_region_elements(my_list)
#--------------------------------------------------------------------------------
