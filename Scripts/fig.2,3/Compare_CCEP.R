# load basic imformation
flatten_mat <- function(mat, k = 0) {
  result <- mat
  result <- result[col(mat) != row(mat)]
  return(result)
}
# load coor.filmat# load coor.file #
conn_path <- "../Data"
conn_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Data"
conn_coor <- read.csv(paste0(conn_path, "/Nodes_informations.csv"), header = T)
region_name <- conn_coor[,3]

# define color list for brain networks
network_colors <- c("#CC3333","#FFFF99", "#339933",
                    "#FFFF00" ,"#66CC33","#E2593D",
                    "#6633CC", "#BC249F", "#DDA036",
                    "gray", "#993399","#6EB4E1","#3366CC")



# load all connectomes
# read posibility directed conn L hemi-sphere
directional_conn <- read.csv(paste0(conn_path, "/connections/Lh_connectionP.csv"), header = T, check.names = F, row.names = 1)
directional_conn_LC <- directional_conn[1:180,1:180]
dimnames(directional_conn_LC) <- list(region_name[1:180], region_name[1:180]) 

# CCEP 
# Read the file, skipping comment lines
CCEP_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Dataset/Human/CCEP_f_tract/"
data <- read.table(paste0(conn_path,"/connections/ccep_probability.txt") , 
                   header = FALSE,      # no column headers
                   comment.char = "#",  # ignore lines starting with #
                   na.strings = "nan")  # treat 'nan' as NA
CCEP_region_name <- read.table(paste0(conn_path,"/connections/ccep_MNI-HCP-MMP1.txt"))  # treat 'nan' as NA
CCEP_region_name <- CCEP_region_name$V1
# Convert to a numeric matrix
mat <- as.matrix(data)
dimnames(mat) <- list(CCEP_region_name, CCEP_region_name)
# select cortex region
CCEP_mat_LC <- mat[row.names(directional_conn_LC),row.names(directional_conn_LC)]
directional_conn_LC <- directional_conn_LC/max(directional_conn_LC)

# DWI based
DTI <- read.csv(paste0(conn_path, "/connections/DWI_connectivity_HCP-YA100_streamline.csv"), header = T, row.names = 1)
DTI_c <- DTI[1:180,1:180]
DTI_c <- log10(DTI_c+1)
DTI_c <- DTI_c/max(DTI_c)
DTI_LC <- as.matrix(DTI_c)
dimnames(DTI_LC) <- list(region_name[1:180], region_name[1:180])

diag(DTI_LC) <- 0
diag(CCEP_mat_LC) <- 0
diag(directional_conn_LC) <- 0

library(viridis)
library(pheatmap)
bk <- seq(0,3,0.1)
CCEP_mat_LC[is.na(CCEP_mat_LC)] <- 0
anno_row <- data.frame(row.names = conn_coor[,3], network = conn_coor[,4])
names(network_colors) <- names(table(conn_coor[,4]))
ann_colors=list(networks=network_colors)
#"magma" (or "A") "inferno" (or "B") "plasma" (or "C") "viridis" (or "D") "cividis" (or "E") "rocket" (or "F") "mako" (or "G") "turbo" (or "H")
bk <- seq(0,0.6,0.01)
conn_coor[,4] <- factor(conn_coor[,4], levels = c("Visual1", "Visual2", "Orbito-Affective",
                                                      "Posterior-Multimodal", "Ventral-Multimodal", "Auditory","Language",
                                                      "Dorsal-Attention", "Cingulo-Opercular",
                                                      "subcortex", "Default", "Frontoparietal", "Somatomotor"))


# Fig. 2a ----------------------------------------------------------------------------------------
bk <- seq(0,0.6,0.01)
P <- pheatmap(CCEP_mat_LC[order(conn_coor[,4][1:180]),order(conn_coor[,4][1:180])],
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         breaks = bk,
         color = colorRampPalette(c("white","white","#9CCCD6","#307284","black"))(length(bk)),
         cellheight = 1,
         cellwidth = 1,
         show_rownames = F,
         show_colnames = F,
         clustering_method = "complete",
         annotation_row = anno_row,
         annotation_col = anno_row,
         annotation_colors = ann_colors)

bk <- seq(0,0.4,0.01)
P <- pheatmap(directional_conn_LC[order(conn_coor[,4][1:180]),order(conn_coor[,4][1:180])],
              scale = "none",
              cluster_rows = F,
              cluster_cols = F,
              breaks = bk,
              color = colorRampPalette(c("white","white","#9CCCD6","#307284","black"))(length(bk)),
              cellheight = 1,
              cellwidth = 1,
              show_rownames = F,
              show_colnames = F,
              clustering_method = "complete",
              annotation_row = anno_row,
              annotation_col = anno_row,
              annotation_colors = ann_colors)
# -------------------------------------------------------------------------------------------------

################################################################################
# define function for p calculation
top_mask <- function(matrix, ratio = 0.2){
  rank <- rev(sort(matrix))
  thr_pos <- round((length(rank) * ratio), 0) 
  thr <- rank[thr_pos]
  matrix_new <- matrix
  matrix_new[matrix_new >= thr] <- 1
  matrix_new[matrix_new < thr] <- 0
  return(matrix_new)
}

# use different matrices to evaluate overlap bettween connectomes
overlap_evalue <- function(from_mat, check_mat){
  matrix <- from_mat + 2 * check_mat
  counts <- c(sum(as.numeric(matrix == 1)),
              sum(as.numeric(matrix == 2)),
              sum(as.numeric(matrix == 3)))
  I <- counts[3]
  nA <- counts[1] + I
  nB <- counts[2] + I
  N <- length(matrix)
  K <- (nA + nB)/2
  mean_hyp <- K * (K / N)
  var_hyp <- K*(K/N)*(1-K/N)*((N-K)/(N-1))  # 等价公式，或使用 phyper 的参数
  sd_hyp <- sqrt(var_hyp)
  jaccard <- I / (nA + nB - I)
  dice <- 2*I / (nA + nB)
  fold_enrichment <- (I * N) / (nA * nB)
  z_score <- (I - ((nA*nB)/N))/sd_hyp 
  p_log <- phyper(I, m=nA, n=N-nA, k=nB, lower.tail=FALSE, log.p=TRUE)
  if("matrix" %in% is(from_mat)){
    cor <- cor(flatten_mat(from_mat), flatten_mat(check_mat))
  }else{
    cor <- cor(from_mat, check_mat)
  }
  features <- c(I, jaccard, dice, fold_enrichment, z_score, p_log, cor)
  return(features)
}

###############################################################################
# calculate the posibility 
directional_conn_LC <- as.matrix(directional_conn_LC)
sweep_ratio <- seq(0,1,0.002)
features_mat <- matrix(NA, nrow = length(sweep_ratio), ncol = 3*7)
colnames(features_mat) <- paste0(rep(c("CE_D", "DTI_D", "DTI_CE"),each=7), "_",
                                 c("I", "jaccard", "dice", "foldEnrichment", "zScore", "pLog", "cor"))


for(i in 1:length(sweep_ratio)){
  ratio = sweep_ratio[i]
  CCEP_tmp <- top_mask(CCEP_mat_LC, ratio = ratio)
  D_tmp <- top_mask(directional_conn_LC, ratio = ratio)
  DTI_tmp <- top_mask(DTI_LC, ratio = ratio)
  
  # calculate features 
  overlap_features_CE_D <- overlap_evalue(CCEP_tmp, D_tmp)
  overlap_features_DTI_D <- overlap_evalue(DTI_tmp, D_tmp)
  overlap_features_DTI_CE <- overlap_evalue(DTI_tmp, CCEP_tmp)
  features_mat[i,1:7] <- overlap_features_CE_D
  features_mat[i,8:14] <- overlap_features_DTI_D
  features_mat[i,15:21] <- overlap_features_DTI_CE
}
features_mat <- as.data.frame(features_mat)
features_mat$ratio <- sweep_ratio 
features_mat

cor_D_fract_max <- c()
for(i in 1:length(sweep_ratio)){
  ratio = sweep_ratio[i]
  CCEP_tmp <- top_mask(CCEP_mat_LC, ratio = 0.45)
  D_tmp <- top_mask(directional_conn_LC, ratio = ratio)
  cor_D_fract_max <-c(cor_D_fract_max, cor(flatten_mat(CCEP_tmp),flatten_mat(D_tmp)))  
}
features_mat$CE_D_cor_max <- cor_D_fract_max 
ratio = 0.2


# fig. 2b -----------------------------------------------------------------------------------
library(ggplot2)
ggplot() + 
  geom_line(data=features_mat[2:500,], aes(x=ratio, y=CE_D_cor), color = "#307284") +
  geom_line(data=features_mat[225:500,], aes(x=ratio, y=CE_D_cor_max), color = "#607284") +
  geom_line(data=features_mat[2:500,], aes(x=ratio, y=DTI_D_cor), color ="gray") +
  #geom_line(data=features_mat[2:500,], aes(x=ratio, y=DTI_CE_cor), color = "blue") +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0)) +
  theme_classic()

################################################################################

# fig. 2c -----------------------------------------------------------------------------------
# conpare top 20% connection matrix
D_20 <- top_mask(directional_conn_LC, ratio = 0.2)
C_20 <- top_mask(CCEP_mat_LC, ratio = 0.2)
T_20 <- top_mask(DTI_LC, ratio = 0.2)

merge_m <- D_20 + 2*C_20
bk <- seq(0,3,0.1)
P <- pheatmap(merge_m[order(conn_coor[,4][1:180]),order(conn_coor[,4][1:180])],
              scale = "none",
              cluster_rows = F,
              cluster_cols = F,
              breaks = bk,
              color = colorRampPalette(c("white","#ea902f","#307284","black"))(length(bk)),
              cellheight = 1,
              cellwidth = 1,
              show_rownames = F,
              show_colnames = F,
              clustering_method = "complete",
              annotation_row = anno_row,
              annotation_col = anno_row,
              border_color = F,
              annotation_colors = ann_colors)
#----------------------------------------------------------------------------------------------

# fig. 2e -------------------------------------------------------------------------------------
brain_region_info <- conn_coor[1:180, c(1:4)]
brain_region_info$D_out <- apply(D_20, 1, sum)
brain_region_info$D_in <- apply(D_20, 2, sum)
brain_region_info$C_out <- apply(C_20, 1, sum)
brain_region_info$C_in <- apply(C_20, 2, sum)

P <- ggplot() + 
  geom_point(data=brain_region_info, aes(x=D_in, y=D_out, colour = NETWORK..Ji..J..L..et.al..2019.),size =2) + 
  scale_color_manual(values = network_colors) +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,100)) + 
  theme_classic() +
  coord_fixed()
P <- ggplot() + 
  geom_point(data=brain_region_info, aes(x=C_in, y=C_out, colour = NETWORK..Ji..J..L..et.al..2019.),size =2) + 
  scale_color_manual(values = network_colors) +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,100)) + 
  theme_classic() +
  coord_fixed()
#-------------------------------------------------------------------------------------

# fig. 2d ----------------------------------------------------------------------------
DTI_length <- read.csv(paste0(conn_path, "/connections/DWI_connectivity_HCP-YA100_streamline_length.csv"), header = F)
DTI_length <- as.matrix(DTI_length)
diag(DTI_length) <- 0
DTI_length_c <- DTI_length[1:180,1:180]
ggplot() +
  geom_histogram(aes(DTI_length_c[D_20>0]), bins = 50, fill = "red", alpha = 0.4) +
  geom_histogram(aes(DTI_length_c[T_20>0]), bins = 50, fill = "gray", alpha = 0.4) +
  geom_histogram(aes(DTI_length_c[C_20>0]), bins = 50, fill = "blue", alpha = 0.4) +
  theme_bw() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
#-------------------------------------------------------------------------------------

# fig. 2f ----------------------------------------------------------------------------
# graph_analysis
brain_region_info$C_in_out_ratio <- brain_region_info$C_in/brain_region_info$C_out
brain_region_info$D_in_out_ratio <- brain_region_info$D_in/brain_region_info$D_out

P <- ggplot() +
  geom_histogram(data=brain_region_info, aes(x=log2(D_in_out_ratio)), fill = "#ea902f", bins = 40) +
  geom_histogram(data=brain_region_info, aes(x=log2(C_in_out_ratio)), fill = "#307284", bins = 40) +
  theme_bw() + scale_y_continuous(limits =  c(0,60), expand = c(0,0)) + scale_x_continuous(limits =  c(-5,5), expand = c(0,0))
P
#-------------------------------------------------------------------------------------------

# fig. 3 data generate code ----------------------------------------------------------------
# sensitive analysis for: small world / asymm ratio / Henricl index / schur index
anno_row <- data.frame(row.names = trophic$Region, networks =  trophic$Network)
names(network_colors) <- names(table(conn_coor$NETWORK))
ann_colors=list(networks=network_colors)
# import calculation methods
library(igraph)
library(pheatmap)
library(pracma)
install.packages("expm")
library(expm)
#select top matrix
top_mask <- function(matrix, ratio = 0.2){
  rank <- rev(sort(matrix))
  thr_pos <- round((length(rank) * ratio), 0) 
  thr <- rank[thr_pos]
  matrix_new <- matrix
  matrix_new[matrix_new >= thr] <- 1
  matrix_new[matrix_new < thr] <- 0
  return(matrix_new)
}

calculate_directed_lattice_cc <- function(n, k) {
  return((3 * (k - 2)) / (4 * (k - 1)) * 1.0) 
}

create_directed_lattice <- function(n, k) {
  grid_size <- round(sqrt(n))
  g <- make_lattice(
    dimvector = c(grid_size, grid_size),
    nei = 1,
    directed = TRUE,
    mutual = FALSE,
    circular = TRUE
  )
  if (vcount(g) > n) {
    g <- delete_vertices(g, (n+1):vcount(g))
  }
  out_degrees <- degree(g, mode = "out")
  current_avg_out <- mean(out_degrees)
  if (current_avg_out < k) {
    needed_edges <- round((k - current_avg_out) * n)
    for(i in 1:needed_edges) {
      from <- sample(1:vcount(g), 1)
      to <- sample(setdiff(1:vcount(g), from), 1)
      if (!are_adjacent(g, from, to)) {
        g <- add_edges(g, c(from, to))
      }
    }
  }
  return(g)
}

calculate_directed_omega <- function(graph, num_random = 10) {
  n <- vcount(graph)  
  avg_out_degree <- mean(degree(graph, mode = "out")) 
  k <- round(avg_out_degree)
  C_actual <- transitivity(graph, type = "average") 
  L_actual <- average.path.length(graph, directed = TRUE) 
  L_random_sum <- 0
  valid_random_count <- 0
  for (i in 1:num_random) {
    # Erdős–Rényi model
    random_g <- erdos.renyi.game(n, p.or.m = (avg_out_degree/(n-1)), 
                                 type = "gnp", directed = TRUE)
    if (is.connected(random_g, mode = "strong") || 
        is.connected(random_g, mode = "weak")) {
      if (!is.connected(random_g, mode = "strong")) {
        random_g <- giant.component(random_g, mode = "weak")
      }
      L_random_i <- average.path.length(random_g, directed = TRUE)
      L_random_sum <- L_random_sum + L_random_i
      valid_random_count <- valid_random_count + 1
    }
  }
  L_random <- L_random_sum / valid_random_count
  g <- make_ring(n)
  g <- connect.neighborhood(g, order = k/2)
  C_lattice <- transitivity(g, type = "average")

  # Omega (ω)
  omega <- (L_random / L_actual) - (C_actual / C_lattice)
  result <- list(
    omega = omega,
    metrics = list(
      C_actual = C_actual,
      L_actual = L_actual,
      L_random = L_random,
      C_lattice = C_lattice,
      n = n,
      avg_out_degree = avg_out_degree,
      avg_in_degree = mean(degree(graph, mode = "in")),
      valid_random_networks = valid_random_count
    )
  )
  return(result)
}
# calculate asymm connection rate
calculate_uniconn_rato <- function(matrix) {
  matrix = matrix + t(matrix)
  res <- table(matrix)
  ratio <- (res["1"]/2)/(res["1"]/2 + res["2"])
  return(ratio)
}
# calculate henrici
henrici_index <- function(A) {
  # Ensure matrix
  A <- as.matrix(A)
  At <- Conj(t(A))
  comm <- At %*% A - A %*% At
  num <- norm(comm, type = "F")
  den <- norm(A, type = "F")^2
  return(num / den)
}
# calculate schur
schur_upper_index <- function(A) {
  A <- as.matrix(A)
  sch <- Schur(A)   # complex Schur
  T <- sch$T
  # strictly upper triangular part
  U <- T
  U[lower.tri(U, diag = TRUE)] <- 0
  num <- norm(U, type = "F")
  den <- norm(T, type = "F")
  return(num / den)
}

step <- 0.02
sweep_ratio <- seq(step,1,step)
features_mat <- matrix(NA, nrow = length(sweep_ratio), ncol = 8)
colnames(features_mat) <- c("omega", "L_actual", "L_random", "C_actual", "C_lattice",
                            "uni_ratio", "henrici","schur")

for(i in 1:length(sweep_ratio)){
  ratio = sweep_ratio[i]
  D_tmp <- top_mask(directional_conn, ratio =ratio)
  # calculate features 
  
  tryCatch({g <- graph_from_adjacency_matrix(adjmatrix = D_tmp,
                                   mode = "directed", 
                                   weighted = NULL) 
  result <- calculate_directed_omega(g, num_random = 10)
  features_mat[i,"omega"] <- result$omega
  features_mat[i,"L_actual"] <- result$metrics$L_actual
  features_mat[i,"L_random"] <- result$metrics$L_random
  features_mat[i,"C_actual"] <- result$metrics$C_actual
  features_mat[i,"C_lattice"] <- result$metrics$C_lattice}, error = function(e) {
    # 错误处理
    message(paste("loop", i, "error:", e$message))
  })
  
  features_mat[i,"uni_ratio"] <- calculate_uniconn_rato(D_tmp)
  features_mat[i,"henrici"] <- henrici_index(D_tmp)
  features_mat[i,"schur"] <- schur_upper_index(D_tmp)
}
features_mat <- as.data.frame(features_mat)
features_mat$ratio <- sweep_ratio 





