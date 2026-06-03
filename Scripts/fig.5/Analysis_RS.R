library(pracma)
library(gsignal)
library(pheatmap)
library(ggplot2)
library(uwot)
library(fitdistrplus)
library(dbscan)
library (MASS)
library (convey)
###########################################
# define fcuntions
# resting state annalysis
# load brain region information
scale01 <- function(ts){
  ts_s <- (ts - min(ts))/(max(ts) - min(ts))
}

extract_upper_triangle <- function(mat, k = 0) {
  result <- mat
  result <- result[col(mat) - row(mat) >= k]
  return(result)
}

flatten_mat <- function(mat, k = 0) {
  result <- mat
  result <- result[col(mat) != row(mat)]
  return(result)
}

read_BOLD <- function(file_path, region_name){
  cortex_L <- 181:360
  cortex_R <- 1:180
  subcor_L <- 361:(361+32)
  subcor_R <- 394:(394+32)
  region_rank <- c(cortex_L, subcor_L, cortex_R, subcor_R)
  # load BOLD singnal
  BOLD_HCP <- read.csv(file_path, header = F)
  BOLD_HCP <- BOLD_HCP[,region_rank]
  colnames(BOLD_HCP) <- conn_coor$HCP_MMP_name
  return(BOLD_HCP)
}

# GSR ananlysis for noise flite
perform_gsr <- function(data_matrix) {
  global_signal <- colMeans(data_matrix, na.rm = TRUE)
  global_signal_deriv <- c(0, diff(global_signal))
  design_matrix <- cbind(global_signal, global_signal_deriv)
  gsr_data <- matrix(NA, nrow(data_matrix), ncol(data_matrix))
  for (i in 1:nrow(data_matrix)) {
    if (all(is.na(data_matrix[i, ]))) {
      gsr_data[i, ] <- NA
      next
    }
    fit <- lm(data_matrix[i, ] ~ design_matrix)
    gsr_data[i, ] <- residuals(fit)
  }
  return(gsr_data)
}


# 1. preprocessing of time sires
# freq fliter
preprocess_fmri <- function(ts, tr = 0.72, high_freq = 0.08, low_freq = 0.008) {
  detrended <- apply(ts, 2, function(x) {
    lm(x ~ seq_along(x))$residuals
  })
  nyq_freq <- 1/(2*tr)
  bf <- butter(2, c(low_freq/nyq_freq, high_freq/nyq_freq), type = "pass")
  # filtered
  filtered <- apply(detrended, 2, function(x) {
    filtfilt(bf, x)
  })
  
  z_scores <- apply(filtered, 2, scale)
  return(t(z_scores)) # shape : brain region * time points
} 

preprocess_fmri_nopass <- function(ts) {
  detrended <- apply(ts, 2, function(x) {
    lm(x ~ seq_along(x))$residuals
  })
  #nyq_freq <- 1/(2*tr)
  #bf <- butter(2, c(low_freq/nyq_freq, high_freq/nyq_freq), type = "pass")
  # filtered
  #filtered <- apply(detrended, 2, function(x) {
  #  filtfilt(bf, x)
  #})
  z_scores <- apply(detrended, 2, scale)
  return(t(z_scores)) # shape : brain region * time points
} 

# 2. calculate  
# brain region * time points

# get FCD 
getFCD <- function(ts, window_length, step_size) {
  ts <- t(ts)
  # get dimensions
  FC_t <- list()
  FC_ave <- cor(ts)
  inter_size <- (nrow(ts) - window_length)%/%step_size
  
  inter <- 1
  while (inter <= inter_size) {
    start <- (inter - 1) * step_size + 1 
    end <- start + window_length
    window_data <- ts[start:end, ]
    FC <- cor(window_data)
    FC_t[[inter]] <- FC # - FC_ave
    inter <- inter + 1 
  }
  
  FCD_p <- matrix(nrow = length(FC[upper.tri(FC_t[[1]])]), ncol = inter_size)
  for(i in 1:inter_size){
    FCD_p[,i] <- FC_t[[i]][upper.tri(FC_t[[i]])]
  }
  FCD <- cor(FCD_p)
  FCD_data <- list(FCD, FCD_p)
  return(FCD_data)
}

# calculate dwell-time
calculate_dwell_time <- function(state_seq) {
  state_seq <- state_seq[state_seq != 0]
  dwell_times <- list()
  current_state <- state_seq[1]
  count <- 1
  
  for (i in 2:length(state_seq)) {
    if (state_seq[i] == current_state) {
      count <- count + 1
    } else {
      dwell_times[[current_state]] <- c(dwell_times[[current_state]], count)
      current_state <- state_seq[i]
      count <- 1
    }
  }
  dwell_times[[current_state]] <- c(dwell_times[[current_state]], count)  # 添加最后一段
  return(dwell_times)
}

#########################
# MAIN ANNALYSIS ########
#########################
# read brain regions
conn_path <- "./"
diected_conn <- read.csv(paste0(conn_path, "/DirectionalConnectome_D99_HCPex_LR_Log10weighted.csv"), header = F)
conn_coor <- read.csv(paste0(conn_path, "/HCPmmp_HCPex_node_coordinate_ranked.csv"), header = T)
region_name <- conn_coor$HCP_MMP_name
# define color list for brain networks
network_colors <- c("#CC3333","#FFFF99", "#339933",
                    "#FFFF00" ,"#66CC33","#E2593D",
                    "#6633CC", "#BC249F", "#DDA036",
                    "gray", "#993399","#6EB4E1","#3366CC")

# get Thr D DS matrix
DTI_conn <- read.csv(paste0(conn_path, "/connectivity_HCP-YA100_streamline_Log10weighted.csv"), header = T)
DTI_conn <- DTI_conn[,-1]
dimnames(DTI_conn) <- list(region_name,region_name)
dimnames(diected_conn) <- list(region_name,region_name)

mask <- matrix(0, nrow = 426, ncol = 426)
mask[diected_conn > 0] <- 1
mask <- mask + t(mask)
mask[mask > 0] <- 1

DS_conn <- matrix(0, nrow = 426, ncol = 426)
DS_conn[mask == 1] <- DTI_conn[mask == 1]
dimnames(DS_conn) <- list(region_name,region_name)

Thr <- sort(as.matrix(DTI_conn))[round(426*426*0.8,0)]
Thr_conn <- matrix(0, nrow = 426, ncol = 426)
Thr_conn[DTI_conn > Thr] <- DTI_conn[DTI_conn > Thr]
dimnames(Thr_conn) <- list(region_name,region_name)


#load trophic levels 
data_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Data"
trophic <- read.csv(paste0(data_path,"/Nodes_informations.csv"))
colnames(trophic)[c(3,4)] <- c("region_names", "Network")
trophic$Network <- factor(trophic$Network, levels = c("Visual1", "Visual2", "Orbito-Affective",
                                                      "Posterior-Multimodal", "Ventral-Multimodal", "Auditory","Language",
                                                      "Dorsal-Attention", "Cingulo-Opercular",
                                                      "subcortex", "Default", "Frontoparietal", "Somatomotor"))




################################################
# calculate group averanged FC for 100 subject #
################################################
# BOLD_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/"
# subjects <- list.files(BOLD_path)
# RS_file <- list.files(paste0(BOLD_path,"/",subjects[1]))

# process_subject <- function(BOLD_path, subject){
#  print(paste0("subject: ",subject))
#  RS_file <- list.files(paste0(BOLD_path,"/",subject))
#  RS_file <- RS_file[grep(".csv",RS_file)]
#  file_count <- 0
#  ave_fc <- matrix(0, nrow = 426, ncol = 426)
#  for(file in RS_file){
#    ts <- read_BOLD(file_path = paste0(BOLD_path,"/",subject,"/",file),
#                    region_name = region_name)
#    data_matrix <- t(ts)
#    gsr_data <- perform_gsr(data_matrix)
#    ts_gsr <- preprocess_fmri(t(gsr_data))
#    ts <- t(ts_gsr)
#    FC_f <- cor(ts)
#    if(mean(FC_f) < 0.4){
#      ave_fc <- ave_fc + FC_f
#      file_count <- file_count + 1 
#      print(paste0(file_count))
#    }
#  }
#  ave_fc <- ave_fc/file_count
#  return(ave_fc)
# }

# subject <- subjects[11]
# RS_file <- list.files(paste0(BOLD_path,"/",subject))
# RS_file <- RS_file[grep(".csv",RS_file)]
# file <- RS_file[1]
# ts <- read_BOLD(file_path = paste0(BOLD_path,"/",subject,"/",file),
#                region_name = region_name)
# data_matrix <- t(ts)
# gsr_data <- perform_gsr(data_matrix)
# ts_gsr <- preprocess_fmri(t(gsr_data))
# ts <- t(ts_gsr)
# FC_f <- cor(ts)
# ts <- t(ts)
# ts <- ts_gsr
# FCD_data <- getFCD(ts,window_length = 42, step_size = 7)

# Fig .5b -----------------------------------------------------
# library(viridis)
# diag(FC_f) <- 0
# bk <- seq(0,1,0.01)
# P <- pheatmap(FCD,
#              scale = "none",
#              cluster_rows = F,
#              cluster_cols = F,
#              breaks = bk,
#              color = turbo(length(bk)),,
#              cellheight = 1,
#              cellwidth = 1,
#              show_rownames = F,
#              show_colnames = F,
#              border_color = F,
#              clustering_method = "complete")
# save_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/annalysis/example"
# file <- paste0(save_path,"/subject118225_FC.pdf")
# pdf(file,width=8,height=8)
# print(P)
# dev.off()

# df_b <- data.frame(time = 101:933, region63=ts[63,101:933], 
#                   region70=ts[73,101:933], region53=ts[53,101:933])
# P <- ggplot() + 
#  geom_line(data=df_b[,], aes(x=time, y=scale01(region63)),color="#6EB4E1") +
#  theme_classic()
# P

# P <- ggplot() + 
#   geom_line(data=df_b[,], aes(x=time, y=scale01(region70)),color="#6EB4E1") +
#   theme_classic()
# P

# P <- ggplot() + 
#   geom_line(data=df_b[,], aes(x=time, y=scale01(region53)),color="#6EB4E1") +
#   theme_classic()
# P


# FC_emprical <- matrix(0, nrow = 426, ncol = 426)
# subject_count <- 0
# for (subject in subjects) {
#   FC_s <- process_subject(BOLD_path, subject)
#   if(!is.na(sum(FC_s)) & mean(FC_s < 0.4)){
#     print(subject)
#    FC_emprical_flit <- FC_emprical_flit + FC_s
#     subject_count <- subject_count + 1
#   }
# }
# FC_emprical_flit <- FC_emprical_flit/subject_count
# dimnames(FC_emprical_flit) <- list(region_name, region_name)
# write.csv(FC_emprical_flit,"D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/annalysis/data/averaged_FC_emprical_gsr.csv")

###################################################
# calculate FCD distribution ######################
###################################################
# FCD_D_matrix <- matrix(0, nrow = length(subjects), ncol = 100) 
# dimnames(FCD_D_matrix) <- list(subjects, round(seq(-1,0.98,0.02),3))

# BOLD_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/"
# subjects <- list.files(BOLD_path)
# RS_file <- list.files(paste0(BOLD_path,"/",subjects[1]))

# process_subject <- function(BOLD_path, subject){
#  print(paste0("subject: ",subject))
#   RS_file <- list.files(paste0(BOLD_path,"/",subject))
#  RS_file <- RS_file[grep(".csv",RS_file)]
#  file_count <- length(RS_file)
#  FCD_D <- numeric(100)
#  for(file in RS_file){
#    ts <- read_BOLD(file_path = paste0(BOLD_path,"/",subject,"/",file),
#                    region_name = region_name)
#    data_matrix <- t(ts)
#    gsr_data <- perform_gsr(data_matrix)
#    ts_gsr <- preprocess_fmri(t(gsr_data))
#    ts <- ts_gsr
#    FCD_data <- getFCD(ts,window_length = 42, step_size = 7)
#    P <-hist(extract_upper_triangle(FCD_data[[1]],k=6),breaks = round(seq(-1,1,0.02),3))
#    FCD_D <- P$counts + FCD_D
#  }
#  return(FCD_D)
#}

# for (subject in subjects) {
#  FC_D_t <- process_subject(BOLD_path, subject)
#  FCD_D_matrix[subject,] <- FC_D_t
#}
# write.csv(FCD_D_matrix, "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/annalysis/data/FCD_Distributions_gsr.csv")

FCD_D_matrix <- read.csv("D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/annalysis/data/FCD_Distributions_gsr.csv", header = T, row.names = 1)
colnames(FCD_D_matrix) <- seq(-1,0.98,0.02)
FCD_D_Empirical <- apply(FCD_D_matrix, 2, sum)
FCD_D_Empirical <- FCD_D_Empirical/sum(FCD_D_Empirical) # FCD values distribution

########################################
# Process simulation data ##############
########################################
# read result of parameter sweeps
# D_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/synthetic/20251202_D"
# DS_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/synthetic/20251215_DS"
# Thr_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/synthetic/20251229_Thr"


process_subject <- function(BOLD){
  subject_list <- list()
  ts <- preprocess_fmri(BOLD)
  FC <- cor(t(ts))
  FCD <- getFCD(ts,window_length = 42, step_size = 7)
  subject_list[["FC"]] <- FC
  subject_list[["FCD"]] <- FCD
  return(subject_list)
}


# Fig .5c data --------------------------------------------------------------------
# plot showed in paper used averaged SC weighted in simulatino 
# ---------------------------------------------------------------------------------
sweep_list_conn <- list()
for(conn in c("D","DS", "Thr")){
  if(conn == "D"){data_path <- D_path}
  if(conn == "DS"){data_path <- DS_path}
  if(conn == "Thr"){data_path <- Thr_path}
  print(data_path)
  files <- list.files(paste0(data_path,"/results"))
  sigma_list <- paste0("sigma", c(0.04,0.043,0.045,0.048,0.05,0.052,0.055,0.057,0.06))
  if(conn == "D"){G_list<- paste0("G",seq(0.04,0.05,0.0005))}
  if(conn == "DS"){G_list <- paste0("G",seq(0.03,0.04,0.0005))}
  if(conn == "Thr"){G_list <- paste0("G",seq(0.03,0.04,0.0005))}
  cor_SC_FC_mat <- matrix(0, nrow = length(sigma_list), ncol = length(G_list))
  dimnames(cor_SC_FC_mat) <- list(sigma_list ,G_list)
  sweep_list <- list(KS_FCD=cor_SC_FC_mat,
                     var_FCD=cor_SC_FC_mat,
                     KS_FC=cor_SC_FC_mat,
                     cor_FC=cor_SC_FC_mat,
                     cor_FSC=cor_SC_FC_mat)
  
  for(file in files){
    G <- strsplit(file,"_")[[1]][2]
    sigma <- strsplit(file,"_")[[1]][3]
    sigma <- gsub(".csv","",sigma)
    if(G %in% dimnames(cor_SC_FC_mat)[[2]] & sigma %in% dimnames(cor_SC_FC_mat)[[1]]){
    BOLD <- read.csv(paste0(data_path,"/results/",file),header = F)
    colnames(BOLD) <- conn_coor$HCP_MMP_name
    BOLD <- BOLD[49:848,]
    subject_list <- process_subject(BOLD)
    FC <- subject_list[["FC"]]
    diag(FC) <- 0
    
    FCD <- subject_list[["FCD"]][[1]]
    cor_FC <- cor(extract_upper_triangle(FC[c(1:180,214:393),c(1:180,214:393)],k=1),extract_upper_triangle(FC_emprical[c(1:180,214:393),c(1:180,214:393)],k=1))
    ks_FC <- ks.test(extract_upper_triangle(FC,k=1),extract_upper_triangle(FC_emprical,k=1))
    ks_FCD <- ks.test(d2,extract_upper_triangle(FCD,k=6))
    var_FCD <- var(extract_upper_triangle(FCD,k=6))
    if(conn == "D"){conn_mat <- diected_conn}
    if(conn == "DS"){conn_mat <- DS_conn}
    if(conn == "Thr"){conn_mat <- Thr_conn}
    cor_FSC <- cor(flatten_mat(conn_mat), flatten_mat(FC))
      sweep_list$KS_FCD[sigma, G] <- ks_FCD$statistic
      sweep_list$var_FCD[sigma, G] <- var_FCD
      sweep_list$cor_FC[sigma, G] <- cor_FC
      sweep_list$KS_FC[sigma, G] <- ks_FC$statistic
      sweep_list$cor_FSC[sigma, G] <- cor_FSC
      print(paste0(file," finished"))
    }
  }
  sweep_list_conn[[conn]] <- sweep_list
}


library(viridis)
library(ggplot2)
sweep_list_conn <- readRDS("./Sweep_results.RDS")

conn <- "Thr"
bk <- seq(0.0,0.05,0.001)
P <- pheatmap(sweep_list_conn[[conn]]$cor_FC[9:1,],
              scale = "none",
              cluster_rows = F,
              cluster_cols = F,
              breaks = bk,
              color = colorRampPalette(cividis(10))(length(bk)),
              cellheight = 10,
              cellwidth = 10,
              show_rownames = T,
              show_colnames = T,
              border_color = F,
              clustering_method = "complete")
P

bk <- seq(0.0,1.0,0.01)
P <- pheatmap(sweep_list_conn[[conn]]$KS_FCD[9:1,],
              scale = "none",
              cluster_rows = F,
              cluster_cols = F,
              breaks = bk,
              color = colorRampPalette(cividis(10))(length(bk)),
              cellheight = 10,
              cellwidth = 10,
              show_rownames = T,
              show_colnames = T,
              border_color = F,
              clustering_method = "complete")
P


bk <- seq(0,0.001,0.00001)
P <- pheatmap(sweep_list_conn[[conn]]$var_FCD[9:1,],
              scale = "none",
              cluster_rows = F,
              cluster_cols = F,
              breaks = bk,
              color = colorRampPalette(cividis(10))(length(bk)),
              cellheight = 10,
              cellwidth = 10,
              show_rownames = T,
              show_colnames = T,
              border_color = F,
              clustering_method = "complete")


bk <- seq(0,0.35,0.01)
P <- pheatmap(sweep_list_conn[[conn]]$cor_FSC[9:1,],
              scale = "none",
              cluster_rows = F,
              cluster_cols = F,
              breaks = bk,
              color = colorRampPalette(cividis(10))(length(bk)),
              cellheight = 10,
              cellwidth = 10,
              show_rownames = T,
              show_colnames = T,
              border_color = F,
              clustering_method = "complete")


#################################################################################################
# process individual levels information #########################################################
#################################################################################################

# data calculate process 
# emprical_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/rfmri_ts_1"
# synthetic_path <- "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/synthetic/individual/selected_BOLD"


process_empirical <- function(BOLD_path, subject){
  print(paste0("subject: ",subject))
  RS_file <- list.files(paste0(BOLD_path,"/",subject))
  RS_file <- RS_file[grep(".csv",RS_file)]
  subject_inf <- list()
  file_count <- 0
  FCD_D <- 0
  ave_fc <- matrix(0, nrow = 426, ncol = 426)
  for(file in RS_file){
    ts <- read_BOLD(file_path = paste0(BOLD_path,"/",subject,"/",file),
                    region_name = region_name)
    data_matrix <- t(ts)
    gsr_data <- perform_gsr(data_matrix)
    ts_gsr <- preprocess_fmri(t(gsr_data))
    ts <- t(ts_gsr)
    FC_f <- cor(ts)
    if(mean(FC_f) < 0.4){
      ave_fc <- ave_fc + FC_f
      file_count <- file_count + 1 
      print(paste0(file_count))
    }
    ts_gsr <- preprocess_fmri(t(gsr_data))
    ts <- ts_gsr
    FCD_data <- getFCD(ts,window_length = 42, step_size = 7)
    var_FCD_tmp <- var(extract_upper_triangle(FCD_data[[1]],k=6))
    FCD_D <- FCD_D + var_FCD_tmp
  }
  subject_inf[["FC"]] <- ave_fc/file_count
  subject_inf[["var_FCD"]]  <- FCD_D/file_count
  subject_inf[["FCD"]] <- FCD_data[[1]]
  return(subject_inf)
}


process_synthetic <- function(BOLD){
  subject_list <- list()
  ts <- preprocess_fmri(BOLD)
  FC <- cor(t(ts))
  FCD_data <- getFCD(ts,window_length = 42, step_size = 7)
  subject_list[["FC"]] <- FC
  subject_list[["var_FCD"]] <- var(extract_upper_triangle(FCD_data[[1]],k=6))
  subject_list[["FCD"]] <- FCD_data[[1]]
  return(subject_list)
}

# read subject list
syn_file_list <- list.files(synthetic_path)
Inf_t <- t(as.data.frame(strsplit(syn_file_list,"_")))
subjects <- names(table(Inf_t[,2]))[table(Inf_t[,2]) == 3]

subjects_Info <- list()
for (subject in subjects) {
  print(paste0(subject, " start"))
  emp_Inf <- process_empirical(emprical_path, subject = subject)
  subjects_Info[[subject]][["empirical"]] <- emp_Inf
  for (conn in c("D", "DS", "T")) {
    file_name <- paste0(conn,"_",subject,"_BOLD.csv")
    BOLD <- read.csv(paste0(synthetic_path ,"/",file_name),header = F)
    colnames(BOLD) <- conn_coor$HCP_MMP_name
    BOLD <- BOLD[49:848,]
    syn_Inf <- process_synthetic(BOLD)
    subjects_Info[[subject]][[conn]] <- syn_Inf
  }
}

# saveRDS(subjects_Info, "D:/7_Research/3_projects/20250104_MacaqueHuman/Results/rfmri_ts/annalysis/data/FC_subjects_emp_syn.RDS")
subjects_Info <- readRDS("D:/7_Research/3_projects/20250104_MacaqueHuman/Results/Materials_for_publication/Scripts/fig.5/FC_subjects_emp_syn.RDS")
subjects <- names(subjects_Info)
subjects
df <- data.frame(subject = "A", type = "A", var_FCD = 0, mean_FC = 0)
for(type in c("empirical", "D", "DS", "T")) {
  for(subject in subjects){
    FC <- subjects_Info[[subject]][[type]][["FC"]]
    var_FCD <- subjects_Info[[subject]][[type]][["var_FCD"]]
    mean_FC <- mean(extract_upper_triangle(FC, k=1))
    df_tmp <- data.frame(subject = subject, type = type, var_FCD = var_FCD, mean_FC = mean_FC)
    df <- rbind(df,df_tmp)
  }
}
df <- df[-1,]

df$type <- factor(df$type, levels = c("empirical","T","DS","D"))
library(ggplot2)
P <- ggplot() + 
  # geom_violin(data = df[! df$subject %in% c("155938", "156334", "157336", "159239"),], aes(x=type, y=var_FCD, fill = type),width = 1.0, alpha = 0.5) +
  geom_boxplot(data = df, aes(x=type, y=var_FCD, fill = type),width = 0.5, alpha = 0.5) +
  geom_point(data = df, aes(x=type, y=var_FCD, color = type),size =2) +
  #geom_line(data = df, aes(x=type, y=var_FCD, group = subject)) +
  scale_fill_manual(values = c("#f39800","#b1afa2","#82ad51","#4c542c")) +
  scale_color_manual(values = c("#f39800","#b1afa2","#82ad51","#4c542c")) +
  scale_y_continuous(limits = c(0,0.01),expand = c(0,0)) +
  theme_classic()


subjects_Info





