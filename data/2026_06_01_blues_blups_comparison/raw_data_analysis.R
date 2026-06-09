library(dplyr)
library(lme4)
library(emmeans)

df <- read.csv("ciat_gas_for_shiny_app_and_blues.csv")

str(df)

# table(df$requisitioner)
# length(unique(df$id_lab))
# table(table(df$id_lab))
# table(df$genus, df$species)

df$id_lab  <- factor(df$id_lab)
df$batch   <- factor(df$batch)

traits <- colnames(df)[17:26]

blues <- list()
blups <- list()

write("Analysis Report:", file = "analysis_report.txt")

for (trait in traits) {
  model.1 <- lm(formula(paste0("`", trait, "` ~ id_lab + batch")), data = df)
  
  write("\n========================\n", file = "analysis_report.txt", append = TRUE)
  capture.output(anova(model.1), file = "analysis_report.txt", append = TRUE)
  
  model.2 <- lmer(formula(paste0("`", trait, "` ~ (1 | id_lab) + (1 | batch)")), data = df, REML = TRUE)
  
  vc <- as.data.frame(VarCorr(model.2))
  vc$variation_source <- vc$grp
  vc$variance_ratio <- sprintf("%.2f%%", 100 * vc$vcov / sum(vc$vcov))
  
  write("\nVariance components for random effects:", file = "analysis_report.txt", append = TRUE)
  capture.output(vc[,c("variation_source", "variance_ratio")], file = "analysis_report.txt", append = TRUE)
  
  ranef_info <- ranef(model.2, condVar = TRUE)
  fixed_eff  <- fixef(model.2)["(Intercept)"]
  blups_eff  <- ranef_info$id_lab
  blups_se   <- sqrt(attr(ranef_info$id_lab, "postVar")[1, , ])
  blups_df   <- data.frame(id_lab = rownames(blups_eff), BLUP = blups_eff$`(Intercept)` + fixed_eff, SE = blups_se)
  
  colnames(blups_df) <- c("id_lab", trait, paste0(trait, "_SE"))
  blups[[trait]] <- blups_df
  
  model.3 <- lmer(formula(paste0("`", trait, "` ~ id_lab + (1 | batch)")), data = df, REML = TRUE)
  
  blues_emm <- emmeans(model.3, ~ id_lab)
  blues_df  <- summary(blues_emm)

  pairwise_lsd <- as.data.frame(summary(pairs(blues_emm, adjust = "none")))
  pairwise_lsd$LSD <- qt(0.975, df = pairwise_lsd$df) * pairwise_lsd$SE
  
  avg_LSD <- signif(mean(pairwise_lsd$LSD, na.rm = TRUE), 4)
  write(paste("\nAverage LSD:", avg_LSD), file = "analysis_report.txt", append = TRUE)
  
  blues_df <- blues_df[, c("id_lab", "emmean", "SE")]
  colnames(blues_df) <- c("id_lab", trait, paste0(trait, "_SE"))
  
  blues[[trait]] <- blues_df
}

blue <- Reduce(function(x, y) merge(x, y, by = names(x)[1], all = TRUE), blues)
blup <- Reduce(function(x, y) merge(x, y, by = names(x)[1], all = TRUE), blups)

write.csv(blue, file = "blues.csv", row.names = FALSE)
write.csv(blup, file = "blups.csv", row.names = FALSE)
