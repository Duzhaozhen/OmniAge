## -----------------------------------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)
data("LungInv")
my_comparisons <- list(c("N\nN=21", "LCIS\nN=13"), c("LCIS\nN=13", "LCIS->LC\nN=22"))
table(df$Group)
## Check available epigenetic clocks
listEpiAge()

## -----------------------------------------------------------------------------
EpiAge.o<-EpiAge(data.m = bmiq.m,clock_names = "mitotic",ages.v = df$Age)
rm(bmiq.m)

## ----make-large-plot-1, fig.width=12, fig.height=10, out.width="100%"---------
df$epiTOC1 <-EpiAge.o$epiTOC1
df$epiTOC2 <-EpiAge.o$epiTOC2$irS
df$epiTOC3 <-EpiAge.o$epiTOC3$irS
df$stemTOCvitro <-EpiAge.o$stemTOCvitro
df$stemTOC<-EpiAge.o$stemTOC
df$HypoClock<-EpiAge.o$HypoClock
df$RepliTali<-EpiAge.o$RepliTali
df$EpiCMIT_Hyper<-EpiAge.o$EpiCMIT_Hyper
df$EpiCMIT_Hypo<-EpiAge.o$EpiCMIT_Hypo
g<-list()
for (i in 1:9){
  temp.df<-data.frame('Group'=df$Group,'Age'=df$Age,'num'=df$num,'score'=df[,i+3])
  y<-summary(lm(temp.df$score~temp.df$num+temp.df$Age))$coefficients[2,4]
  g[[i]]<-ggplot(temp.df,aes(x=Group,y=score,fill=Group))+guides(fill='none')+
  geom_boxplot()+theme_bw()+
  scale_x_discrete(limits=c('N\nN=21','LCIS\nN=13','LCIS->LC\nN=22'))+
  stat_compare_means(method='wilcox.test',comparisons = my_comparisons, size = 4)+
  xlab('')+ylab(colnames(df)[i+3])+
  annotate('text',x=2,y=max(temp.df$score) * 0.9,label=paste0('P(Age-adjusted)=',signif(y,2)), size = 4)+
  theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 11)
    )

}
    
ggpubr::ggarrange(plotlist = g, nrow = 3, ncol = 3)

## -----------------------------------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)
data("Hannum_example", overwrite = TRUE)
age <- PhenoTypesHannum_lv$Age
sex <- ifelse(PhenoTypesHannum_lv$Sex=="F","Female","Male")
EpiAge.out <- EpiAge(hannum_bmiq_m,clock_names = c("Horvath2013","ZhangClock"))

## ----make-large-plot-2, fig.width=8, fig.height=4, out.width="100%"-----------
gp<-list()

for (i in seq_along(EpiAge.out)){
  plot_df <- data.frame(ActualAge = PhenoTypesHannum_lv$Age, PredictedAge = EpiAge.out[[i]])
  cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
  correlation <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
  p_value_formatted <- ifelse(p_value < 0.001,
                            formatC(p_value, format = "e", digits = 2),
                            round(p_value, 3))
  annotation_text <- paste0(
  "R = ", round(correlation, 3), "\n",
  "P = ", p_value_formatted, "\n",
  "MAE = ", round(mae, 2)
  )
  gp[[i]]<-ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm",color = "black", se = FALSE) +
  annotate("text",
           x = min(plot_df$ActualAge, na.rm = TRUE),
           y = max(plot_df$PredictedAge, na.rm = TRUE),
           label = annotation_text,
           hjust = 0,  
           vjust = 1,             size = 5) +

  labs(
    title = NULL,
    x = "Chronological Age",
    y = paste0("Predicted Age(", c("Horvath2013","ZhangClock")[i],")")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

}

ggpubr::ggarrange(plotlist = gp, nrow = 1, ncol = 2)


## -----------------------------------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)
## Check available epigenetic clocks
listEpiGA()
data("GA_example", overwrite = TRUE)
EpiGA.o <- EpiGA(GA_m,clock_names = c("Knight_GA","Mayne_GA"))

## ----make-large-plot-3, fig.width=8, fig.height=4, out.width="100%"-----------
gp<-list()

for (i in seq_along(EpiGA.o)){
  plot_df <- data.frame(Gestational_Age = phenotypeGA$age, PredictedAge = EpiGA.o[[i]])
  cor_test_result <- cor.test(plot_df$Gestational_Age, plot_df$PredictedAge)
  correlation <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  mae <- mean(abs(plot_df$PredictedAge - plot_df$Gestational_Age))
  p_value_formatted <- ifelse(p_value < 0.001,
                            formatC(p_value, format = "e", digits = 2),
                            round(p_value, 3))
  annotation_text <- paste0(
  "R = ", round(correlation, 3), "\n",
  "P = ", p_value_formatted, "\n",
  "MAE = ", round(mae, 2)
  )
  gp[[i]]<-ggplot(data = plot_df, aes(x = Gestational_Age, y = PredictedAge)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm",color = "black", se = FALSE) +
  annotate("text",
           x = min(plot_df$Gestational_Age, na.rm = TRUE),
           y = max(plot_df$PredictedAge, na.rm = TRUE),
           label = annotation_text,
           hjust = 0,  
           vjust = 1,             size = 5) +

  labs(
    title = NULL,
    x = "Gestational Age(weeks)",
    y = paste0("Predicted Age(", c("Knight_GA","Mayne_GA")[i],", weeks)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

}

ggpubr::ggarrange(plotlist = gp, nrow = 1, ncol = 2)


## -----------------------------------------------------------------------------
library(OmniAgeR)
library(Seurat)
library(glmnet)
library(ggplot2)
library(patchwork)
library(ggpubr)
data("Yazar_CD4T_CD8T_example")
scImmuAging.out <- scImmuAging(seurat_obj,c("CD4T","CD8T"))

## ----make-large-plot-4, fig.width=8, fig.height=4, out.width="100%"-----------
sc_gp<-list()

for (i in seq_along(scImmuAging.out)){
  plot_df <- data.frame(ActualAge = scImmuAging.out[[i]]$Donor$age, PredictedAge = scImmuAging.out[[i]]$Donor$predicted)
  cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
  correlation <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
  p_value_formatted <- ifelse(p_value < 0.001,
                            formatC(p_value, format = "e", digits = 2),
                            round(p_value, 3))
  annotation_text <- paste0(
  "R = ", round(correlation, 3), "\n",
  "P = ", p_value_formatted, "\n",
  "MAE = ", round(mae, 2)
  )
  sc_gp[[i]]<-ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm",color = "black", se = FALSE) +
  annotate("text",
           x = min(plot_df$ActualAge, na.rm = TRUE),
           y = max(plot_df$PredictedAge, na.rm = TRUE),
           label = annotation_text,
           hjust = 0,  
           vjust = 1, 
           size = 5) +

  labs(
    title = c("CD4T","CD8T")[i],
    x = "Chronological Age",
    y = "Predicted Age(sc-ImmuAging)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

}

ggpubr::ggarrange(plotlist = sc_gp, nrow = 1, ncol = 2)

## -----------------------------------------------------------------------------
library(OmniAgeR)
library(Seurat)
library(glmnet)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
data("brain_frohlich_control_example_15donors")

brain_clock_results <- Brain_CT_clock(
   seurat_object = brain_seurat,
   cell_types = c("Oligodendrocytes"),
   model_name = "all" ##  c('SC', 'Pseudobulk', 'Bootstrap')
 )

# Use lapply to iterate over each data frame in the list
average_predictions_by_donor <- lapply(brain_clock_results, function(df) {
  df %>%
    group_by(donors,ages,sample_type,celltype) %>%
    summarise(mean_prediction = mean(predictions, na.rm = TRUE))
})


## ----make-large-plot-5, fig.width=8, fig.height=3, out.width="100%"-----------
brain_gp<-list()

for (i in seq_along(average_predictions_by_donor)){
  
  plot_df <- data.frame(ActualAge = average_predictions_by_donor[[i]]$ages, PredictedAge = average_predictions_by_donor[[i]]$mean_prediction)
  cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
  correlation <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
  p_value_formatted <- ifelse(p_value < 0.001,
                            formatC(p_value, format = "e", digits = 2),
                            round(p_value, 3))
  annotation_text <- paste0(
  "R = ", round(correlation, 3), "\n",
  "P = ", p_value_formatted, "\n",
  "MAE = ", round(mae, 2)
  )
  brain_gp[[i]]<-ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm",color = "black", se = FALSE) +
  annotate("text",
           x = min(plot_df$ActualAge, na.rm = TRUE),
           y = max(plot_df$PredictedAge, na.rm = TRUE),
           label = annotation_text,
           hjust = 0,  
           vjust = 1, 
           size = 4) +

  labs(
    title = "Oligodendrocytes",
    x = "Chronological Age",
    y = paste0("Predicted Age(",c('SC', 'Pseudobulk', 'Bootstrap')[i],")")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

}

ggpubr::ggarrange(plotlist = brain_gp, nrow = 1, ncol = 3)


## -----------------------------------------------------------------------------
# Load required libraries
library(OmniAgeR)
library(Seurat)
library(glmnet)
library(magrittr) # For %>% pipe
library(ggplot2)
library(patchwork)
library(ggpubr)    # For simplified plot annotations
# Load the pre-filtered Seurat object
data("seu_gabitto_2024_filtered")

# Extract and clean chronological age from metadata
seu$age <- seu$development_stage %>%
  gsub("-year.*", "", .) %>%             # Remove suffixes
  gsub("-", " ", .) %>%                 # Remove hyphens
  gsub("80 year old and over stage", "85", .) # Standardize 80+ category

set.seed(42)

# Create pseudobulk samples.
# Cells are grouped by 'cell_type' and 'age', then aggregated into
# chunks of 512 cells to simulate bulk RNA-seq samples.
seu_bulk <- making_pseudobulks_from_seurat(seu,
                                           pool_by = c("cell_type", "age"),
                                           chunk_size = 512,
                                           verbose = FALSE)

# Extract the log-normalized expression matrix for prediction
lognorm_matrix <- GetAssayData(seu_bulk, assay = "RNA", layer = "data")
lognorm_matrix <- as.matrix(lognorm_matrix)

# Extract the corresponding metadata for the new pseudobulk samples
seu_bulk_meta <- seu_bulk[[c('chunk_size', 'cell_type', 'age')]]
seu_bulk_meta$age <- as.numeric(seu_bulk_meta$age)

# Apply the PASTA clock
# filter_genes = TRUE: Selects only the genes required by the PASTA model.
# rank_norm = TRUE: Applies the rank-normalization required by PASTA.
PASTA_res <- PASTA_Scores(lognorm_matrix, filter_genes = TRUE, rank_norm = TRUE)

## ----make-large-plot-6, fig.width=4, fig.height=4, out.width="100%"-----------
# Prepare the data frame for plotting
plot_df <- data.frame(ActualAge = seu_bulk_meta$age, 
                        Prediction = PASTA_res$PASTA)

# Create the scatter plot
ggplot(data = plot_df, aes(x = ActualAge, y = Prediction)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  ggpubr::stat_cor(method = "pearson", 
                   label.x.npc = "left",
                   label.y.npc = "top", 
                   hjust = 0,            
                   size = 5) +
  labs(
    title = NULL,
    x = "Chronological Age",
    y = "PASTA Score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )


## -----------------------------------------------------------------------------
library(OmniAgeR)
data("Hannum_example")

CompCRP.o <- CompCRP(hannum_bmiq_m)
print(lapply(CompCRP.o, head, 5))

CompCHIP.o <- CompCHIP(hannum_bmiq_m)
print(lapply(CompCHIP.o, head, 5))



## -----------------------------------------------------------------------------
library(OmniAgeR)
# 1. Get probes for a specific chronological aging predictor
chronological_probes <- getClockProbes(clock_names = "Horvath2013")
print(head(chronological_probes$Horvath2013))
# 2. Get probes for a specific biological aging clocks
bio_probes <- getClockProbes(clock_names = "GrimAge1")
print(head(bio_probes$GrimAge1))
# 3. Get probes for all biomaker
PCClocks_RData <- load_OmniAgeR_data(object_name = "PCClocks_data")
SystemsAge_RData <- load_OmniAgeR_data(object_name = "SystemsAge_data")
all_probes <- getClockProbes(clock_names = "all",PCClocks_RData,SystemsAge_RData)
names(all_probes)[1:10]


## -----------------------------------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)

data("TZH_example_CTF")
DNAm_CTF_Clock_o<-DNAm_CTF_Clock(CTF_m = TZH_Frac_m)



## ----make-large-plot-7, fig.width=4, fig.height=4, out.width="100%"-----------

plot_df <- data.frame(ActualAge = PhenoTypes_TZH_df$Age, PredictedAge = DNAm_CTF_Clock_o)
cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
correlation <- cor_test_result$estimate
p_value <- cor_test_result$p.value
mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
p_value_formatted <- ifelse(p_value < 0.001,
                          formatC(p_value, format = "e", digits = 2),
                          round(p_value, 3))
annotation_text <- paste0(
"R = ", round(correlation, 3), "\n",
"P = ", p_value_formatted, "\n",
"MAE = ", round(mae, 2)
)
ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm",color = "black", se = FALSE) +
  annotate("text",
         x = min(plot_df$ActualAge, na.rm = TRUE),
         y = max(plot_df$PredictedAge, na.rm = TRUE),
         label = annotation_text,
         hjust = 0,  
         vjust = 1,             size = 5) +

  labs(
    title = NULL,
    x = "Chronological Age",
    y = paste0("Predicted Age(DNAm_CTF)")
    ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
    )



#ggpubr::ggarrange(plotlist = gp, nrow = 1, ncol = 1)


## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

