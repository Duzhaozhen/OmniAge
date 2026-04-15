## ----mototic-clock-example----------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)
lungInv <- loadOmniAgeRdata(
    "omniager_lung_inv",
    verbose = FALSE
)
lungInvM <- lungInv$bmiq_m
phenoDf <- lungInv$PhenoTypes

my_comparisons <- list(c("N\nN=21", "LCIS\nN=13"), c("LCIS\nN=13", "LCIS->LC\nN=22"))
table(phenoDf$Group)
## Check available epigenetic clocks
listEpiMarker()

## ----run-mototic-clock--------------------------------------------------------
epiMarkerRes <- epiMarker(
    betaM = lungInvM,
    clockNames = "mitotic",
    chronAge = phenoDf$Age,
    minCoverage = 0,
    verbose = FALSE
)
rm(lungInvM)

## ----mototic-clock-plot, fig.width=12, fig.height=10, out.width="100%"--------
phenoDf$epiTOC1 <- epiMarkerRes$epiTOC1
phenoDf$epiTOC2 <- epiMarkerRes$epiTOC2$irS
phenoDf$epiTOC3 <- epiMarkerRes$epiTOC3$irS
phenoDf$stemTOCvitro <- epiMarkerRes$stemTOCvitro
phenoDf$stemTOC <- epiMarkerRes$stemTOC
phenoDf$HypoClock <- epiMarkerRes$HypoClock
phenoDf$RepliTali <- epiMarkerRes$RepliTali
phenoDf$EpiCMIT_Hyper <- epiMarkerRes$EpiCMIT_Hyper
phenoDf$EpiCMIT_Hypo <- epiMarkerRes$EpiCMIT_Hypo
g <- list()
for (i in 1:9) {
    tempDf <- data.frame("Group" = phenoDf$Group, "Age" = phenoDf$Age, "num" = phenoDf$num, "score" = phenoDf[, i + 3])
    y <- summary(lm(tempDf$score ~ tempDf$num + tempDf$Age))$coefficients[2, 4]
    g[[i]] <- ggplot(tempDf, aes(x = Group, y = score, fill = Group)) +
        guides(fill = "none") +
        geom_boxplot() +
        theme_bw() +
        scale_x_discrete(limits = c("N\nN=21", "LCIS\nN=13", "LCIS->LC\nN=22")) +
        stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, size = 4) +
        xlab("") +
        ylab(colnames(phenoDf)[i + 3]) +
        annotate("text", x = 2, y = max(tempDf$score) * 0.9, label = paste0("P(Age-adjusted)=", signif(y, 2)), size = 4) +
        theme(
            axis.title = element_text(size = 11),
            axis.text = element_text(size = 11)
        )
}

ggpubr::ggarrange(plotlist = g, nrow = 3, ncol = 3)

## ----chronological-age-clock--------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)
hannumExample <- loadOmniAgeRdata(
    "omniager_hannum_example",
    verbose = FALSE
)
hannumBmiqM <- hannumExample[[1]]
phenoTypesHannum <- hannumExample[[2]]
age <- phenoTypesHannum$Age
sex <- ifelse(phenoTypesHannum$Sex == "F", "Female", "Male")
epiMarkerOut <- epiMarker(hannumBmiqM,
    clockNames = c("Horvath2013", "ZhangClock"),
    minCoverage = 0, verbose = FALSE
)

## ----chronological-age-clock-plot, fig.width=8, fig.height=4, out.width="100%"----
gp <- list()

for (i in seq_along(epiMarkerOut)) {
    plot_df <- data.frame(ActualAge = phenoTypesHannum$Age, PredictedAge = epiMarkerOut[[i]])
    cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
    correlation <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
    mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
    p_value_formatted <- ifelse(p_value < 0.001,
        formatC(p_value, format = "e", digits = 2),
        round(p_value, 3)
    )
    annotation_text <- paste0(
        "R = ", round(correlation, 3), "\n",
        "P = ", p_value_formatted, "\n",
        "MAE = ", round(mae, 2)
    )
    gp[[i]] <- ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
        geom_point(color = "steelblue", size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        annotate("text",
            x = min(plot_df$ActualAge, na.rm = TRUE),
            y = max(plot_df$PredictedAge, na.rm = TRUE),
            label = annotation_text,
            hjust = 0,
            vjust = 1, size = 5
        ) +
        labs(
            title = NULL,
            x = "Chronological Age",
            y = paste0("Predicted Age(", c("Horvath2013", "ZhangClock")[i], ")")
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

## ----gestational-age----------------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)
gaExample <- loadOmniAgeRdata(
    "omniager_ga_example",
    verbose = FALSE
)
gaM <- gaExample[[1]]
phenoTypesGa <- gaExample[[2]]
epiGA <- epiMarker(gaM,
    clockNames = c("KnightGA", "MayneGA"),
    minCoverage = 0, verbose = FALSE
)

## ----gestational-age-plot, fig.width=8, fig.height=4, out.width="100%"--------
gp <- list()

for (i in seq_along(epiGA)) {
    plot_df <- data.frame(Gestational_Age = phenoTypesGa$age, PredictedAge = epiGA[[i]])
    cor_test_result <- cor.test(plot_df$Gestational_Age, plot_df$PredictedAge)
    correlation <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
    mae <- mean(abs(plot_df$PredictedAge - plot_df$Gestational_Age))
    p_value_formatted <- ifelse(p_value < 0.001,
        formatC(p_value, format = "e", digits = 2),
        round(p_value, 3)
    )
    annotation_text <- paste0(
        "R = ", round(correlation, 3), "\n",
        "P = ", p_value_formatted, "\n",
        "MAE = ", round(mae, 2)
    )
    gp[[i]] <- ggplot(data = plot_df, aes(x = Gestational_Age, y = PredictedAge)) +
        geom_point(color = "steelblue", size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        annotate("text",
            x = min(plot_df$Gestational_Age, na.rm = TRUE),
            y = max(plot_df$PredictedAge, na.rm = TRUE),
            label = annotation_text,
            hjust = 0,
            vjust = 1, size = 5
        ) +
        labs(
            title = NULL,
            x = "Gestational Age(weeks)",
            y = paste0("Predicted Age(", c("Knight GA", "Mayne GA")[i], ", weeks)")
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

## ----scImmuAging--------------------------------------------------------------
library(OmniAgeR)
library(Seurat)
library(glmnet)
library(ggplot2)
library(patchwork)
library(ggpubr)
#' library(Seurat)
seuratObj <- loadOmniAgeRdata(
    "omniager_yazar_cd4t_cd8t_example",
    verbose = FALSE
)
scImmuAgingOut <- scImmuAging(seuratObj, c("CD4T", "CD8T"))

## ----scImmuAging-plot, fig.width=8, fig.height=4, out.width="100%"------------
sc_gp <- list()

for (i in seq_along(scImmuAgingOut)) {
    plot_df <- data.frame(ActualAge = scImmuAgingOut[[i]]$donor$age, PredictedAge = scImmuAgingOut[[i]]$donor$predicted)
    cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
    correlation <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
    mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
    p_value_formatted <- ifelse(p_value < 0.001,
        formatC(p_value, format = "e", digits = 2),
        round(p_value, 3)
    )
    annotation_text <- paste0(
        "R = ", round(correlation, 3), "\n",
        "P = ", p_value_formatted, "\n",
        "MAE = ", round(mae, 2)
    )
    sc_gp[[i]] <- ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
        geom_point(color = "steelblue", size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        annotate("text",
            x = min(plot_df$ActualAge, na.rm = TRUE),
            y = max(plot_df$PredictedAge, na.rm = TRUE),
            label = annotation_text,
            hjust = 0,
            vjust = 1,
            size = 5
        ) +
        labs(
            title = c("CD4T", "CD8T")[i],
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

## ----brain-ct-clock-----------------------------------------------------------
library(OmniAgeR)
library(Seurat)
library(glmnet)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)

brainSeurat <- loadOmniAgeRdata(
    "omniager_brain_frohlich_control_example_15donors",
    verbose = FALSE
)

# Define cell types of interest
cellTypes <- c("Oligodendrocytes")

# Run all three models for the specified cell types
clockResults <- brainCtClock(
    seuratObj = brainSeurat,
    cellTypes = cellTypes
)

# Use lapply to iterate over each data frame in the list
averagePredictionsDonor <- lapply(clockResults, function(df) {
    df %>%
        group_by(donorId, age, sampleType, celltype) %>%
        summarise(mean_prediction = mean(prediction, na.rm = TRUE))
})

## ----brain-ct-clock-plot, fig.width=8, fig.height=3, out.width="100%"---------
brainPlot <- list()

for (i in seq_along(averagePredictionsDonor)) {
    plot_df <- data.frame(
        ActualAge = averagePredictionsDonor[[i]]$age,
        PredictedAge = averagePredictionsDonor[[i]]$mean_prediction
    )
    cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
    correlation <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
    mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
    p_value_formatted <- ifelse(p_value < 0.001,
        formatC(p_value, format = "e", digits = 2),
        round(p_value, 3)
    )
    annotation_text <- paste0(
        "R = ", round(correlation, 3), "\n",
        "P = ", p_value_formatted, "\n",
        "MAE = ", round(mae, 2)
    )
    brainPlot[[i]] <- ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
        geom_point(color = "steelblue", size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        annotate("text",
            x = min(plot_df$ActualAge, na.rm = TRUE),
            y = max(plot_df$PredictedAge, na.rm = TRUE),
            label = annotation_text,
            hjust = 0,
            vjust = 1,
            size = 4
        ) +
        labs(
            title = "Oligodendrocytes",
            x = "Chronological Age",
            y = paste0("Predicted Age(", c("SC", "Pseudobulk", "Bootstrap")[i], ")")
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12)
        )
}

ggpubr::ggarrange(plotlist = brainPlot, nrow = 1, ncol = 3)

## ----pasta--------------------------------------------------------------------
# Load required libraries
library(OmniAgeR)
library(Seurat)
library(glmnet)
library(magrittr) # For %>% pipe
library(ggplot2)
library(patchwork)
library(ggpubr) # For simplified plot annotations

seu <- loadOmniAgeRdata(
    "omniager_seu_gabitto_2024_filtered",
    verbose = FALSE
)
# Extract and clean chronological age from metadata
seu$age <- seu$development_stage %>%
    gsub("-year.*", "", .) %>%
    gsub("-", " ", .) %>%
    gsub("80 year old and over stage", "85", .)

# Create pseudobulk samples, Pseudobulk construction may follow Salignon et al.
# or be performed using alternative approaches as appropriate
set.seed(42)
seuBulk <- makePseudobulksPasta(
    seu,
    poolBy = c("cell_type", "age"),
    chunkSize = 512,
    verbose = FALSE
)

# Extract the log-normalized expression matrix for prediction
lognormMatrix <- GetAssayData(seuBulk, assay = "RNA", layer = "data")
lognormMatrix <- as.matrix(lognormMatrix)

# Extract the corresponding metadata for the new pseudobulk samples
seuBulkMeta <- seuBulk[[c("chunkSize", "cell_type", "age")]]
seuBulkMeta$age <- as.numeric(seuBulkMeta$age)

# Apply the PASTA clock
# filter_genes = TRUE: Selects only the genes required by the PASTA model.
# rank_norm = TRUE: Applies the rank-normalization required by PASTA.
pastaRes <- pastaScores(lognormMatrix, filterGenes = TRUE, rankNorm = TRUE)

## ----pasta-plot, fig.width=4, fig.height=4, out.width="100%"------------------
# Prepare the data frame for plotting
plot_df <- data.frame(
    ActualAge = as.numeric(seuBulk$age),
    Prediction = as.numeric(pastaRes$PASTA)
)

# Create the scatter plot
ggplot(data = plot_df, aes(x = ActualAge, y = Prediction)) +
    geom_point(color = "steelblue", size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    ggpubr::stat_cor(
        method = "pearson",
        label.x.npc = "left",
        label.y.npc = "top",
        hjust = 0,
        size = 5
    ) +
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

## ----crp-chip-----------------------------------------------------------------
library(OmniAgeR)
hannumBmiqM <- loadOmniAgeRdata(
    "omniager_hannum_example",
    verbose = FALSE
)[[1]]
crpRes <- compCRP(hannumBmiqM)
print(lapply(crpRes, head, 5))

chipRes <- compCHIP(hannumBmiqM)
print(lapply(chipRes, head, 5))

## ----calculateProteomicOrganAges----------------------------------------------
library(OmniAgeR)

# ---------------------------------------------------------
# 1. WyssCorayOrganAge for SomaScan protein data
# ---------------------------------------------------------
protExampleSomaScan <- loadOmniAgeRdata(
    "omniager_proteomic_somascan_v41_example",
    verbose = FALSE
)

# Extract metadata and expression matrix
metadataSoma <- protExampleSomaScan[[1]]
protExpSoma <- protExampleSomaScan[[2]]

# Estimate ages
resultsSomaScan <- wyssCorayOrganAge(protExp = protExpSoma, 
                                     metadata = metadataSoma, 
                                     assayVersion = "v4.1",
                                     reference="cohort")

head(resultsSomaScan$Adipose)


# ---------------------------------------------------------
# 2. GladyshevOrganAge for Olink protein data
# ---------------------------------------------------------
protExampleOlink <- loadOmniAgeRdata(
    "omniager_proteomic_olink1536_example",
    verbose = FALSE
)

# Extract metadata and expression matrix
protExpOlink <- protExampleOlink[[1]]

# Estimate ages
resultsOlink <- gladyshevOrganAge(protExp = protExpOlink,
                                  modelType = "gen2",
                                  platform = "reduced")

head(resultsOlink$Adipose)

## ----getMarkerWeights---------------------------------------------------------
library(OmniAgeR)
# 1. Get probes for a specific chronological aging predictor
chronologicalProbes <- getMarkerWeights(clockNames = "Horvath2013")
print(head(chronologicalProbes$Horvath2013))
# 3. Get probes for all biomaker
allFeatures<- getMarkerWeights(clockNames = "all")
names(allFeatures)[1:10]

## ----dnam-ctf-clock-----------------------------------------------------------
library(OmniAgeR)
library(ggplot2)
library(patchwork)
library(ggpubr)

tzhExample <- loadOmniAgeRdata(
    "omniager_tzh_example_ctf",
    verbose = FALSE
)
phenoTypesTzh <- tzhExample[[1]]
tzhFracM <- tzhExample[[2]]
dnamCTFClockOut <- dnamCTFClock(ctfM = tzhFracM)

## ----dnam-ctf-clock-plot, fig.width=4, fig.height=4, out.width="100%"---------
plot_df <- data.frame(ActualAge = phenoTypesTzh$Age, PredictedAge = dnamCTFClockOut)
cor_test_result <- cor.test(plot_df$ActualAge, plot_df$PredictedAge)
correlation <- cor_test_result$estimate
p_value <- cor_test_result$p.value
mae <- mean(abs(plot_df$PredictedAge - plot_df$ActualAge))
p_value_formatted <- ifelse(p_value < 0.001,
    formatC(p_value, format = "e", digits = 2),
    round(p_value, 3)
)
annotation_text <- paste0(
    "R = ", round(correlation, 3), "\n",
    "P = ", p_value_formatted, "\n",
    "MAE = ", round(mae, 2)
)
ggplot(data = plot_df, aes(x = ActualAge, y = PredictedAge)) +
    geom_point(color = "steelblue", size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    annotate("text",
        x = min(plot_df$ActualAge, na.rm = TRUE),
        y = max(plot_df$PredictedAge, na.rm = TRUE),
        label = annotation_text,
        hjust = 0,
        vjust = 1, size = 5
    ) +
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


# ggpubr::ggarrange(plotlist = gp, nrow = 1, ncol = 1)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

