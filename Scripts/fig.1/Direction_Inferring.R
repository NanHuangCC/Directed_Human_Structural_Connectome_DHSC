# date 20250220
# annalys for blueprints # 
# ref: 
#1. Lu, Y. et al. Macaque Brainnetome Atlas: A multifaceted brain map with parcellation, connection, and histology. Science Bulletin 69, 2241–2259 (2024).
#2. Mars, R. B. et al. Whole brain comparative anatomy using connectivity blueprints. eLife 7, e35237 (2018).
# read data blue prints results #

#############################################################
# coding which include raw data treating were marked as Notes
#############################################################
blueprint_path <- "../Data/blueprints"
tract_order <- read.table(paste0(blueprint_path, "/tract_order.txt"))
##############################################################
blueprint_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Data/blueprints"
# read index 
tract_order <- read.table(paste0(blueprint_path, "/tract_order.txt"))
tracts <- c(grep("_l", tract_order$V1), grep("_r", tract_order$V1))
tracts <- c(tracts,setdiff(1:42, tracts))
# read region level blueprints
blueprints.list <- list()
# read blueprint human
c_bp <- read.csv(paste0(blueprint_path, "/blueprint_human_HCPex_cortex_HCP-YA100.csv"))
rownames(c_bp) <-1:360 
c_bp <- c_bp[,-1]
c_bp <- c_bp[,tract_order$V1[tracts]]
s_bp <- read.csv(paste0(blueprint_path, "/blueprint_human_HCPex_subcortex_HCP-YA100.csv"))
rownames(s_bp) <- s_bp[,1]
s_bp <- s_bp[,-1]
s_bp <- s_bp[,tract_order$V1[tracts]]
blueprints.list[["HCPex"]][["cortex"]] <- c_bp
blueprints.list[["HCPex"]][["subcortex"]] <- s_bp

# read blueprint macaque
c_bp <- read.csv(paste0(blueprint_path, "/blueprint_macaque_D99_cortex_BNA.csv"))
rownames(c_bp) <- c_bp[,1]
c_bp <- c_bp[,-1]
c_bp <- c_bp[,tract_order$V1[tracts]]
s_bp <- read.csv(paste0(blueprint_path, "/blueprint_macaque_D99_subcortex_BNA.csv"))
rownames(s_bp) <- s_bp[,1]
s_bp <- s_bp[,-1]
s_bp <- s_bp[,tract_order$V1[tracts]]
blueprints.list[["D99"]][["cortex"]] <- c_bp
blueprints.list[["D99"]][["subcortex"]] <- s_bp

# Fig. 1b ----------------------------------------------------
# plot bp 
library(pheatmap)
bk = seq(0,0.4,0.02)
nrow(blueprints.list[["HCPex"]][["subcortex"]])
rownames(blueprints.list[["HCPex"]][["subcortex"]])
pheatmap(blueprints.list[["D99"]][["cortex"]],
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         breaks = bk,
         color = colorRampPalette(c("white","#94BCC6","#2D6C7C"))(length(bk)),
         cellheight = 0.5,
         cellwidth =5,
         show_rownames = F,
         show_colnames = F,
         border_color = F,
         clustering_method = "ward.D2")
#------------------------------------------------------------


# Fig. 1c ---------------------------------------------------
# combine all blueprint
bp_all <- blueprints.list[[1]][[1]][1,]
struc_list <- c()
subj_list <- c()
for(subj in names(blueprints.list)){
  for(struc in names(blueprints.list[[1]])){
    t_bp <- blueprints.list[[subj]][[struc]]
    rownames(t_bp) <- paste0(subj,"_", rownames(t_bp))
    struc_list <- c(struc_list, rep(struc, nrow(t_bp)))
    subj_list <- c(subj_list, rep(subj, nrow(t_bp)))
    bp_all <- rbind(bp_all, t_bp)
  }
}
nrow(bp_all)
bp_all <- bp_all[-1,]
rownames(bp_all)

library(ggplot2)
# cortex region: V1(M01_322, H_181) M1 8 
332 +522
H_region <- "HCPex_181"
M_region <- "D99_854"
df <- data.frame(tracts = colnames(bp_all),
                 H=as.numeric(bp_all[grep(H_region,row.names(bp_all)),]), 
                 M=as.numeric(bp_all[grep(M_region,row.names(bp_all)),])) 
df$tracts <- factor(df$tracts, levels = df$tracts)

P <- ggplot() + 
  geom_point(data=df, aes(x=tracts, y=H), color="#F2B038") +
  geom_line(data=df, aes(x=tracts, y=H), group=1, color="#F2B038") +
  geom_point(data=df, aes(x=tracts, y=M), color ="#541B13") +
  geom_line(data=df, aes(x=tracts, y=M), group=1, color="#541B13") +
  theme_classic() + scale_y_continuous(expand = c(0,0) ,limits = c(0, 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
P
#---------------------------------------------------------------

# calculate dierverce matrix between human and macaque
# 2. Mars, R. B. et al. Whole brain comparative anatomy using connectivity blueprints. eLife 7, e35237 (2018).
# Dij = sum(Mik * log(Mik/Hjk)) + sum(Hjk * log(Hjk/Mik))
mapping_mat <- function(bp_h, bp_m) {
  c_h <- bp_h
  c_m <- bp_m
  map_mat <- matrix(NA, nrow = nrow(c_h), ncol = nrow(c_m))
  dimnames(map_mat) <- list(rownames(c_h), row.names(c_m))
  # calculate the 1st value for mapping
  # intergate blueprints
  for(m_reg in colnames(map_mat)){
    print(m_reg)
    for(h_reg in row.names(map_mat)){
      tmp_m <- c_m[m_reg,]
      tmp_h <- c_h[h_reg,]
      log_r1 <- (log2(tmp_m/tmp_h))
      log_r1[log_r1 < -50] <- -50 # avoid Inf
      log_r1[log_r1 > 50] <- 50
      sum1 <- 0
      for(i in 1:length(tmp_h)){
        sum1 <- sum1 + as.numeric(tmp_m[i])*as.numeric(log_r1[i])
      }
      sum1
      log_r2 <- (log2(tmp_h/tmp_m))
      log_r2[log_r2 < -50] <- -50
      log_r2[log_r2 > 50] <- 50
      sum2 <- 0
      for(i in 1:length(tmp_h)){
        sum2 <- sum2 + as.numeric(tmp_h[i])*as.numeric(log_r2[i])
      }
      sum = sum1+sum2
      map_mat[h_reg,m_reg] <- as.numeric(sum)
    }
  }
  return(map_mat)
} 

# cortex structures
c_h <- bp_all[subj_list == "HCPex" & struc_list == "cortex",]
c_m <- bp_all[subj_list == "D99" & struc_list == "cortex",]
map_mat <- mapping_mat(c_h, c_m)

# subcortex structures 
s_h <- bp_all[subj_list == "HCPex" & struc_list == "subcortex",]
s_m <- bp_all[subj_list == "D99" & struc_list == "subcortex",]
map_mat_s <- mapping_mat(s_h, s_m)

#reload mapping matrix and directional connection in macaque######################################
# map tracer based connectome in macaque th huamn 
# ref: Modha, D. S. & Singh, R. Network architecture of the long-distance pathways in the macaque brain.
# Proceedings of the National Academy of Sciences 107, 13485–13490 (2010). URL https://pnas.org/doi/full/10.1073/pnas.1008054107.
# read index
conn_path <- "../Data/connections"
D99_conn <- read.csv(paste0(conn_path ,"/ConnectionMatrix_Mac_D99.csv"),header = T, row.names = 1)
# read D99 to HCPex mapping 
mapping_c <- map_mat
mapping_s <- map_mat_s
# left hemisphere
idx_h <- as.data.frame(strsplit(row.names(mapping_c),"_"))[2,]
idx_h <- as.numeric(idx_h)
idx_m <- as.data.frame(strsplit(colnames(mapping_c),"_"))[2,]
idx_m <- as.numeric(idx_m)
mapping_c_l <- mapping_c[idx_h <= 180, idx_m <=522]
# left subcortical structure
idx_hs <- data.frame(idx=0, hemi="A", name="A")
tmp <- strsplit(row.names(mapping_s),"_")
for(i in 1:length(tmp)){
  tmp_r <- data.frame(idx=tmp[[i]][2], hemi=tmp[[i]][3], name=tmp[[i]][4])
  idx_hs <- rbind(idx_hs,tmp_r)
}
idx_hs<- idx_hs[-1,]

mapping_s <- mapping_s[order(idx_hs$idx),]
idx_hs <- idx_hs[order(idx_hs$idx),]
idx_ms <- as.data.frame(strsplit(colnames(mapping_s),"_"))[2,]
idx_ms <- as.numeric(idx_ms)
mapping_s_l <- mapping_s[idx_hs$hemi == "L", idx_ms <=522]
colnames(mapping_s_l) <- paste0("D99_",idx_ms[idx_ms <= 522])

#############################################
# mapping process ###########################
#############################################
# left hemisphere
proj_mat <- matrix(0, nrow = nrow(mapping_c_l) + nrow(mapping_s_l), ncol = ncol(D99_conn))
dimnames(proj_mat) <- list(c(rownames(mapping_c_l), rownames(mapping_s_l)), colnames(D99_conn))
para_y <- -4
for (reg_h in rownames(mapping_c_l)) {
  D_list <- mapping_c_l[reg_h,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,ncol(D99_conn)) 
  for(reg in names(w_list)){
    reg_conn <- reg 
    tmp <- as.numeric(w_list[reg]) * as.numeric(D99_conn[reg_conn,]) / sum(w_list)
    if(is.na(tmp[1])){print(paste0(reg_conn, ":",sum(w_list)))}
    proj_tmp <- proj_tmp + tmp
  }
  proj_mat[reg_h,] <-proj_tmp 
}

for (reg_h in rownames(mapping_s_l)) {
  D_list <- mapping_s_l[reg_h ,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,ncol(D99_conn)) 
  for(reg in names(w_list)){
    reg_conn <- reg 
    tmp <- as.numeric(w_list[reg]) * as.numeric(D99_conn[reg_conn,]) / sum(w_list)
    proj_tmp <- proj_tmp + tmp
  }
  proj_mat[reg_h,] <-proj_tmp 
}

# calculate direct mat
direct_mat <- matrix(0, nrow = nrow(mapping_c_l) + nrow(mapping_s_l), ncol = nrow(mapping_c_l) + nrow(mapping_s_l))
dimnames(direct_mat) <- list(c(rownames(mapping_c_l), rownames(mapping_s_l)), c(rownames(mapping_c_l), rownames(mapping_s_l)))
dim(direct_mat)
para_y <- -4
for (reg_h in rownames(mapping_c_l)) {
  print(reg_h)
  D_list <- mapping_c_l[reg_h,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,nrow(proj_mat)) 
  for(reg in names(w_list)){
    reg_conn <- reg
    tmp <- as.numeric(w_list[reg]) * as.numeric(proj_mat[,reg_conn]) / sum(w_list)
    proj_tmp <- proj_tmp + tmp
  }
  direct_mat[,reg_h] <- proj_tmp 
}

for (reg_h in rownames(mapping_s_l)) {
  print(reg_h)
  D_list <- mapping_s_l[reg_h,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,nrow(proj_mat)) 
  for(reg in names(w_list)){
    reg_conn <- reg
    tmp <- as.numeric(w_list[reg]) * as.numeric(proj_mat[,reg_conn]) / sum(w_list)
    proj_tmp <- proj_tmp + tmp
  }
  direct_mat[,reg_h] <- proj_tmp 
}

############################################################################################
# right hemisphere
# left hemisphere
mapping_c_r <- mapping_c[idx_h > 180, idx_m > 522]
# left subcortical structure
mapping_s_r <- mapping_s[idx_hs$hemi == "R", idx_ms > 522]
colnames(mapping_s_r) <- paste0("M01_",idx_ms[idx_ms > 522])
# project matrix generate
proj_mat <- matrix(0, nrow = nrow(mapping_c_r) + nrow(mapping_s_r), ncol = ncol(D99_conn))
dimnames(proj_mat) <- list(c(rownames(mapping_c_r), rownames(mapping_s_r)), colnames(D99_conn))
para_y <- -4
for (reg_h in rownames(mapping_c_r)) {
  print(reg_h)
  D_list <- mapping_c_r[reg_h ,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,ncol(D99_conn)) 
  for(reg in names(w_list)){
    id <- as.numeric(strsplit(reg,"_")[[1]][2]) - 522
    reg_conn <- reg
    tmp <- as.numeric(w_list[reg]) * as.numeric(D99_conn[reg_conn,]) / sum(w_list)
    proj_tmp <- proj_tmp + tmp
  }
  proj_mat[reg_h,] <-proj_tmp 
}


for (reg_h in rownames(mapping_s_r)) {
  print(reg_h)
  D_list <- mapping_s_r[reg_h ,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,ncol(D99_conn)) 
  for(reg in names(w_list)){
    id <- as.numeric(strsplit(reg,"_")[[1]][2]) - 522
    reg_conn <- reg
    tmp <- as.numeric(w_list[reg]) * as.numeric(D99_conn[reg_conn,]) / sum(w_list)
    proj_tmp <- proj_tmp + tmp
  }
  proj_mat[reg_h,] <-proj_tmp 
}


# calculate direct mat
direct_mat <- matrix(0, nrow = nrow(mapping_c_r) + nrow(mapping_s_r), ncol = nrow(mapping_c_r) + nrow(mapping_s_r))
dimnames(direct_mat) <- list(c(rownames(mapping_c_r), rownames(mapping_s_r)), c(rownames(mapping_c_r), rownames(mapping_s_r)))
dim(direct_mat)
para_y <- -4
for (reg_h in rownames(mapping_c_r)) {
  print(reg_h)
  D_list <- mapping_c_r[reg_h,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,nrow(proj_mat)) 
  for(reg in names(w_list)){
    id <- as.numeric(strsplit(reg,"_")[[1]][2]) - 522
    reg_conn <- reg
    tmp <- as.numeric(w_list[reg]) * as.numeric(proj_mat[,reg_conn]) / sum(w_list)
    proj_tmp <- proj_tmp + tmp
  }
  direct_mat[,reg_h] <- proj_tmp 
}

for (reg_h in rownames(mapping_s_r)) {
  print(reg_h)
  D_list <- mapping_s_r[reg_h,]
  w_list <- D_list^para_y
  proj_tmp <- rep(0,nrow(proj_mat)) 
  for(reg in names(w_list)){
    id <- as.numeric(strsplit(reg,"_")[[1]][2]) - 522
    reg_conn <- reg
    tmp <- as.numeric(w_list[reg]) * as.numeric(proj_mat[,reg_conn]) / sum(w_list)
    proj_tmp <- proj_tmp + tmp
  }
  direct_mat[,reg_h] <- proj_tmp 
}

# result provided in Data folder
#########################################################################
# Binarization
mapping_path <- "../Data/connections"
mapping_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Data/connections"
# calculate 20 density of right and left hemisphere connectivity
# left
direct_mat <- read.csv(paste0(mapping_path, "/Lh_connectionP.csv"), header = T, check.names = F, row.names = )
rownames(direct_mat) <- direct_mat[,1]
direct_mat <- direct_mat[,-1] 
direct_mat <- direct_mat/max(direct_mat) # scale
conn_values <- as.numeric(as.matrix(direct_mat))
conn_values <- conn_values[order(-conn_values)]
thr <- conn_values[as.integer(length(conn_values) * 0.20)]

direct_mat_20 <- direct_mat 
direct_mat_20[direct_mat_20 < thr] <- 0
direct_mat_20[direct_mat_20 > 0] <- 1

# right
direct_mat <- read.csv(paste0(mapping_path, "/Rh_connectionP.csv"), header = T, check.names = F)
rownames(direct_mat) <- direct_mat[,1]
direct_mat <- direct_mat[,-1]

direct_mat <- direct_mat/max(direct_mat)
hist(as.matrix(direct_mat), breaks = 200)

conn_values <- as.numeric(as.matrix(direct_mat))
conn_values <- conn_values[order(-conn_values)]
thr <- conn_values[as.integer(length(conn_values) * 0.20)]

direct_mat_20 <- direct_mat 
direct_mat_20[direct_mat_20 < thr] <- 0
direct_mat_20[direct_mat_20 > 0] <- 1

#write.csv(direct_mat_20, paste0(mapping_path, "/directional_conn_M01H01_20dens.csv"))
anno <- data.frame(row.names = row.names(direct_mat), structure = c(rep("cortex",180),rep("subcortex",212-180)))  
bk = seq(0,0.5,0.01)
library(pheatmap)
pheatmap(direct_mat_20,
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         breaks = bk,
         color = colorRampPalette(c("white","#94BCC6","#2D6C7C"))(length(bk)),
         cellheight = 1,
         cellwidth = 1,
         show_rownames = F,
         show_colnames = F,
         clustering_method = "complete")





