# Load packages
library(scISR)
library(ALRA)
#Intsall Rmagic
library(reticulate)
use_python("/usr/bin/python3") 
library(Rmagic)
library(VIPER)
library(scImpute)
library(ScRNAIMM)
library(mclust)
library(flexclust)
library(ggplot2)
library(ggpubr)
#helper function for clustering 
evaluate_ARI <- function(dataset, label, method){
  Imputed <- irlba::prcomp_irlba(t(dataset), n = 50)$x
  cluster_Imputed <- kmeans(Imputed, length(unique(label)),
                                 nstart = 2000, iter.max = 2000)$cluster
  
  ARI <- round(adjustedRandIndex(cluster_Imputed, label),3)
  Jaccard <- unname(comPart(cluster_Imputed, label, type="J"))
  evaluation_frame <- data.frame(method = method, ARI = ARI, Jaccard = Jaccard)
  return(evaluation_frame)
}
evaluation_visualization <- function(evaluation_frame){
  ARI <- ggplot(data = evaluation_frame, aes(x = reorder(method, ARI), y = ARI, fill= method)) +
    geom_bar(stat="identity", width = 0.5) + theme_minimal() + xlab("") + 
    theme(legend.position = "none",
           axis.text.x = element_text(angle = 90))
  
  
  Jaccard <- ggplot(data = evaluation_frame, aes(x = reorder(method, Jaccard), y = Jaccard, fill= method)) +
    geom_bar(stat="identity", width = 0.5) + theme_minimal() + xlab("") + 
    ylab("Jaccard Index") + theme(axis.text.x = element_text(angle = 90))
  
  evaluation_plot <- ggarrange(ARI, Jaccard)
  return(evaluation_plot)
}
#Get the input data and output path 
args <- commandArgs(trailingOnly = TRUE)
#Absolute data set of the data set 
Input_path <- args[1]
#Get dataset name as the file name 
File_name <- tools::file_path_sans_ext(basename(Input_path))
#Absolute path for the output directory 
Output_path <- args[2]
Log10Normalize <- as.logical(args[3])
cores <- as.numeric(args[4])
method <- as.numeric(args[5])
k_clusers <- as.numeric(args[6]) #Optional 
is_10x <- as.logical(args[7])
cluster_lbls <- args[8] #Optional 
  
Input_data <- read.csv(Input_path, row.names = 1)
if (method == 1){ #scISR
  scISR_imputed <- scISR(data = Input_data, ncores = cores)
  #Write the output to csv 
  write.csv(x = scISR_imputed, file = paste0(Output_path, File_name, "_scISR_Imputed.csv"))
}else if (method == 2){ #MAGIC
    MAGIC_data <- magic(t(Input_data), seed = 1, n.jobs = -1)
    #Write the output to csv 
    write.csv(x = t(MAGIC_data$result), file = paste0(Output_path, File_name, "_MAGIC_Imputed.csv"))
}else if (method == 3){ #ALRA
  if (!is.null(k_clusers)){
    ALRA_data <- alra(t(Input_data), k = k_clusers)
  }else{
    ALRA_data <- alra(t(Input_data))
  }
  ALRA_Imputed <- data.frame(t(ALRA_data[[3]]))
  write.csv(x = ALRA_Imputed, file = paste0(Output_path, File_name, "_ALRA_Imputed.csv"))
  
}else if (method == 4){ #VIPER
  print("VIPER Imputation Started...")
  if (is_10x){
    VIPER_res <- VIPER(t(Input_data), num = 1000, percentage.cutoff = 0.5, minbool = FALSE, alpha = 0.5, 
                       report = FALSE, outdir = NULL, prefix = NULL)
  }else{
    VIPER_res <- VIPER(Input_data, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                       report = FALSE, outdir = NULL, prefix = NULL)
  }
  Imputed_VIPER <- VIPER_res$imputed_log
  write.csv(x = Imputed_VIPER, file = paste0(Output_path, File_name, "_VIPER_Imputed.csv"))
  print("VIPER Imputation Ended...")
  
}else if (method == 5){ #scImpute
  print("scImpute Impuation Started...")
  ScImputed <- scimpute(count_path = Input_path, infile = 'csv', outfile = 'csv',
                        out_dir = paste0(Output_path, File_name), 
                        ncores = cores, Kcluster = k_clusers)
  print("scImpute Impuation Ended...")
}else if (method == 6){ #ScRNAIMM-cells
  torch::install_torch(reinstall = F)
  print("ScRNAIMM Impuation Started...")
  if (!is.null(k_clusers)){
    ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,
                                     genes = F, outdir = Output_path,
                                     dataset = File_name,
                                     k = k_clusers)
    
  }else{
    ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,genes = T, outdir = Output_path)
  }
  print("ScRNAIMM Impuation Ended...")
}else if (method == 7){ #ScRNAIMM-genes
  torch::install_torch(reinstall = F)
  print("ScRNAIMM Impuation Started...")
  if(!is.null(k_clusers)){
    ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,genes = T, 
                                     k = k_clusers,
                                     outdir = Output_path,
                                     dataset = paste0(File_name, "_genes"))
  }else{
    ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,genes = T, outdir = Output_path,
                                     dataset = paste0(File_name, "_genes"))
  }

  print("ScRNAIMM Impuation Ended...")
}else{#apply all methods 
  if (!is.null(cluster_lbls))
  {
    print("ARI Evaluation Started...")
    cluster_lbls <- as.numeric(as.factor(read.csv(file = cluster_lbls, header = T)[,1]))
    dataset_frame <- evaluate_ARI(dataset = Input_data, label = cluster_lbls, method = "dataset")
    print("ARI Evaluation Ended...")
    #scISR
    print("scISR Impuation Started...")
    scISR_imputed <- scISR(data = Input_data, ncores = cores)
    print("scISR Impuation Ended...")
    ScISR_frame <- evaluate_ARI(dataset = scISR_imputed, label = cluster_lbls, method = "ScISR")
    ScISR_log <- rbind(dataset_frame, ScISR_frame)
    #MAGIC
    print("MAGIC Impuation Started...")
    MAGIC_data <- magic(t(Input_data), seed = 1, n.jobs = -1)
    print("MAGIC Impuation Ended...")
    MAGIC_frame <- evaluate_ARI(dataset = t(MAGIC_data$result), label = cluster_lbls, method = "MAGIC")
    MAGIC_frame <- rbind(ScISR_log, MAGIC_frame)
    #ALRA
    print("ALRA Impuation Started...")
    if (!is.null(k_clusers)){
      ALRA_data <- alra(t(Input_data), k = k_clusers)
    }else{
      ALRA_data <- alra(t(Input_data))
    }
    ALRA_Imputed <- data.frame(t(ALRA_data[[3]]))
    print("ALRA Impuation Ended...")
    ALRA_frame <- evaluate_ARI(dataset = ALRA_Imputed, label = cluster_lbls, method = "ALRA")
    ALRA_frame <- rbind(ALRA_frame, MAGIC_frame)
    #VIPER
    print("VIPER Impuation Started...")
    if (is_10x)
    {
      VIPER_res <- VIPER(t(Input_data), num = 1000, percentage.cutoff = 0.5, minbool = FALSE, alpha = 0.5, 
                         report = FALSE, outdir = NULL, prefix = NULL)
      Imputed_VIPER <- VIPER_res$imputed_log
      print("VIPER Impuation Ended...")
      VIPER_frame <- evaluate_ARI(dataset = t(Imputed_VIPER), label = cluster_lbls, method = "VIPER")
    }else{
      VIPER_res <- VIPER(Input_data, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                         report = FALSE, outdir = NULL, prefix = NULL)
      Imputed_VIPER <- VIPER_res$imputed_log
      print("VIPER Impuation Ended...")
      VIPER_frame <- evaluate_ARI(dataset = Imputed_VIPER, label = cluster_lbls, method = "VIPER")
    }
    VIPER_frame <- rbind(VIPER_frame, ALRA_frame)
    #scImpute
    print("scImpute Impuation Started...")
    ScImputed <- scimpute(count_path = Input_path, infile = 'csv', outfile = 'csv',
                        out_dir = paste0(Output_path, File_name), 
                        ncores = cores, Kcluster = k_clusers)
    print("scImpute Impuation Ended...")
    ScI_outpath <- list.files(path = Output_path, full.names = T, pattern = "csv")
    ScImputed <- read.csv(ScI_outpath, row.names = 1)
    ScImputed_frame <- evaluate_ARI(dataset = ScImputed, label = cluster_lbls, method = "ScImpute")
    ScImputed_frame <- rbind(ScImputed_frame, ALRA_frame)
    
    print("ScRNAIMM-cells Impuation Started...")
    torch::install_torch(reinstall = F)
    if (!is.null(k_clusers)){
      ScRNAIMM_cells_Imputed <- run_pipeline(Input_data, cells = T,genes = F, k = k_clusers,
                                             outdir = Output_path,
                                             dataset = paste0(File_name, "_cells"))
    }else{
      ScRNAIMM_cells_Imputed <- run_pipeline(Input_data, cells = T,genes = F, 
                                             outdir = Output_path,
                                             dataset = paste0(File_name, "_cells")) 
    }
    MMI_cells_frame <- evaluate_ARI(dataset = ScRNAIMM_cells_Imputed, label = cluster_lbls, method = "ScRNAIMM-cells")
    MMI_cells_frame <- rbind(MMI_cells_frame, ScImputed_frame)
    
    print("ScRNAIMM-genes Impuation Started...")
    torch::install_torch(reinstall = F)
    if(!is.null(k_clusers)){
      ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,genes = T, k = k_clusers,
                                       outdir = Output_path,
                                       dataset = paste0(File_name, "_genes"))
    }else{
      ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,genes = T,
                                       outdir = Output_path,
                                       dataset = paste0(File_name, "_genes")) 
    }
    evaluation_frame <- evaluate_ARI(dataset = ScRNAIMM_Imputed, label = cluster_lbls, method = "ScRNAIMM-genes")
    evaluation_frame <- rbind(evaluation_frame, MMI_cells_frame)
    evaluation_frame <- evaluation_frame[order(evaluation_frame$ARI, decreasing = F),]
    
    write.csv(x = evaluation_frame, file = paste0(Output_path, File_name, "_evaluation.csv"))
    evaluation_plot <- evaluation_visualization(evaluation_frame)
    ggsave(plot = evaluation_plot, filename = paste0(Output_path, File_name, "_plot.jpeg"),
           width = 9, height = 6, dpi = 300)
    
    print("ScRNAIMM Impuation Ended...")
    print("Writing Impuation Files...")
    write.csv(x = scISR_imputed, file = paste0(Output_path, File_name, "_scISR_Imputed.csv"))
    write.csv(x = t(MAGIC_data$result), file = paste0(Output_path, File_name, "_MAGIC_Imputed.csv"))
    write.csv(x = ALRA_Imputed, file = paste0(Output_path, File_name, "_ALRA_Imputed.csv"))
    write.csv(x = Imputed_VIPER, file = paste0(Output_path, File_name, "_VIPER_Imputed.csv"))
    print("Writing Impuation Files Done!...")
  }else{
    #scISR
    print("scISR Impuation Started...")
    scISR_imputed <- scISR(data = Input_data, ncores = cores)
    print("scISR Impuation Ended...")
    #MAGIC
    print("MAGIC Impuation Started...")
    MAGIC_data <- magic(Input_data, seed = 1, n.jobs = -1)
    print("MAGIC Impuation Ended...")
    #ALRA
    print("ALRA Impuation Started...")
    if(!is.null(k_clusers)){
      ALRA_data <- alra(t(Input_data), k = k_clusers)
    }else{
      ALRA_data <- alra(t(Input_data))
    }
    ALRA_Imputed <- data.frame(t(ALRA_data[[3]]))
    print("ALRA Impuation Ended...")
    #VIPER
    print("VIPER Impuation Started...")
    if (is_10x){
      VIPER_res <- VIPER(t(Input_data), num = 1000, percentage.cutoff = 0.5, minbool = FALSE, alpha = 0.5, 
                         report = FALSE, outdir = NULL, prefix = NULL)
    }else{
      VIPER_res <- VIPER(Input_data, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
                         report = FALSE, outdir = NULL, prefix = NULL)
    }
    Imputed_VIPER <- VIPER_res$imputed_log
    print("VIPER Impuation Ended...")
    #scImpute
    print("scImpute Impuation Started...")
    ScImputed <- scimpute(count_path = Input_path, infile = 'csv', outfile = 'csv',
                          out_dir = paste0(Output_path, File_name), 
                          ncores = cores, Kcluster = k_clusers)
    print("scImpute Impuation Ended...")
    print("ScRNAIMM-genes Impuation Started...")
    torch::install_torch(reinstall = F)
    if(!is.null(k_clusers)){
      ScRNAIMM__genes_Imputed <- run_pipeline(Input_data, cells = T,genes = T, k = k_clusers,
                                       outdir = Output_path,
                                       dataset = paste0(File_name, "_genes"))
      
    }else{
      ScRNAIMM__genes_Imputed <- run_pipeline(Input_data, cells = T,genes = T,
                                              outdir = Output_path,
                                              dataset = paste0(File_name, "_genes"))
    }
    print("ScRNAIMM-genes Impuation Ended...")
    print("ScRNAIMM-cells Impuation Started...")
    if(!is.null(k_clusers)){
      ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,genes = F, k = k_clusers,
                                       outdir = Output_path,
                                       dataset = paste0(File_name,"_cells"))
      
    }else{
      ScRNAIMM_Imputed <- run_pipeline(Input_data, cells = T,genes = F, 
                                       outdir = Output_path,
                                       dataset = paste0(File_name, "_cells"))
    }
    print("ScRNAIMM-cells Impuation Ended...")
    print("Writing Impuation Files...")
    write.csv(x = scISR_imputed, file = paste0(Output_path, File_name, "_scISR_Imputed.csv"))
    write.csv(x = t(MAGIC_data$result), file = paste0(Output_path, File_name, "_MAGIC_Imputed.csv"))
    write.csv(x = ALRA_Imputed, file = paste0(Output_path, File_name, "_ALRA_Imputed.csv"))
    write.csv(x = Imputed_VIPER, file = paste0(Output_path, File_name, "_VIPER_Imputed.csv"))
    print("Writing Impuation Files Done!...")
  }

}
