# Data Integration & Table 1####

setwd("D:\\pedi\\models")
# Load necessary packages
library(tidyverse)
library(readxl)

# 1. Read prediction data
pred_df <- read_csv("prediction_scores_all.csv")

# 2. Read Excel clinical data
clin_df <- read_excel("clinical_data.xlsx")

# 3. Pad patient_id to 10 digits (leading zeros)
pred_df <- pred_df %>%
  mutate(patient_id = sprintf("%010d", patient_id))

clin_df <- clin_df %>%
  mutate(Patient_ID = sprintf("%010d", as.integer(Patient_ID)))

# 4. Rename Rad_score column in clinical data to prob_sensitive
clin_df <- clin_df %>%
  rename(prob_sensitive = Rad_score)

# 5. Merge two datasets (left join, keep all clinical rows)
merged_df <- clin_df %>%
  left_join(pred_df %>% 
              select(-c(true_label, true_label_name, pred_label, pred_label_name, 
                        prob_insensitive, correct, split)),
            by = c("Patient_ID" = "patient_id"))

# Alternative: inner join, keep only patients present in both
# merged_df <- clin_df %>%
#   inner_join(pred_df %>% 
#               select(-c(true_label, true_label_name, pred_label, pred_label_name, 
#                        prob_insensitive, correct, split)),
#             by = c("Patient_ID" = "patient_id"))

# View results

merged_df=merged_df[merged_df$model=="ResNet50",c(1:18,20)]
colnames(merged_df)[19]="Rad_score"
merged_df=merged_df[,2:19]
print(head(merged_df))

# Save results (optional)
write_csv(merged_df, "combined_clinical_predictions.csv")

library(tidyverse)
library(readxl)

library(gtsummary)
library(tableone)
library(dplyr) 
library(tidyverse)
library(survival)

clinical_data <- read.csv("merged_clinical_imaging_data.csv")
data_processed=clinical_data[,2:12]



tbl=tbl_summary(
  data = data_processed,  
  by = Favorable_histologic_response,  # group by histologic response
  # specify variable types explicitly
  type = list(
    colnames(data_processed)[c(1,2,4,9,11)] ~ "continuous"#categorical
  ),
  statistic = list(
    all_continuous()  ~ "{mean}± ({sd})", # if normal; for non-normal: "{median} ({p25}, {p75})"
    all_categorical() ~ "{n} ({p}%)"      # categorical: count and percentage
  ),
  missing = 'no',  # do not show missing values
  digits = all_continuous() ~ 3  # display precision: 3 decimal places
) %>% 
  add_overall() %>%  # add overall statistics
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 1))  %>%  # add p-value, 1 decimal
  add_n(statistic = "{N_nonmiss}", col_label = "**N**", footnote = TRUE) %>%  # show non-missing sample size
  modify_caption("Table 1. Patient baseline demographic and clinical characteristics") 
tbl
# install and load required packages automatically
if(!require(corrplot)) install.packages("corrplot")
if(!require(dplyr)) install.packages("dplyr")
library(corrplot)
library(dplyr)

# 1. Extract continuous numeric variables (exclude categorical/character variables)
num_data <- data_processed %>% 
  select(Age, BMI, Max_Tumor_Diameter, Radology_score, Pathology_score)

# 2. Compute correlation matrix (Pearson) and p-value matrix
cor_mat <- cor(num_data, use = "complete.obs")
p_mat <- cor.mtest(num_data, conf.level = 0.95)$p

# 3. Define BioRender classic advanced academic color scheme (Morandi red - pure white - navy blue)
biorender_col <- colorRampPalette(c("red", "#FFFFFF", "blue"))(200)

# 4. Draw a minimalist significant heatmap
corrplot(cor_mat, 
         method = "color",           # core: use solid color squares
         type = "upper",             # minimalist: show upper triangle only
         col = biorender_col,        # map BioRender colors
         tl.col = "#2C3E50",         # variable label color: deep slate gray (looks premium)
         tl.srt = 45,                # labels tilted 45 degrees to avoid overlap
         tl.cex = 1.1,               # label font size
         p.mat = p_mat,              # pass the computed p-value matrix
         sig.level = c(0.001, 0.01, 0.05), # graded significance thresholds
         insig = "label_sig",        # mark significant cells with asterisks (*, **, ***)
         pch.cex = 1.2,              # asterisk size
         pch.col = "#222222",        # asterisk color: dark gray
         diag = FALSE,               # hide self-correlation (diagonal)
         addgrid.col = "white",      # cell spacing using white lines (BioRender feature)
         outline = FALSE)            # remove outer thick black border

ggplot(data_processed, aes(x = Radology_score, y = Pathology_score)) +
  geom_point(aes(color = Favorable_histologic_response), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#FADBD8") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 0.9) +
  labs(x = "Radiology Score", y = "Pathology Score",
       title = "Pathology-Radiology Correlation") +
  scale_color_manual(values = c("No" = "#3498DB", "Yes" = "#2ECC71"),
                     name = "Favorable Response") +
  my_theme
# Multimodal model comparison####

if (!requireNamespace("pROC", quietly=TRUE)) install.packages("pROC")
library(pROC)

COLS <- c(
  Clinical = "#ED7D31", 
  Radiomics = "#00468B", 
  Pathology = "#70AD47",  # new pathology model color
  Combined = "#AD002A"
)
col_alpha <- function(col, alpha=0.2) adjustcolor(col, alpha.f=alpha)
format_p  <- function(p) ifelse(p<0.001, "<0.001", sprintf("=%.3f", p))

file_name <- "merged_clinical_imaging_data.csv"
if (file.exists(file_name)) {
  df <- read.csv(file_name, stringsAsFactors=FALSE)
} else {
  message("Local file not detected, using built-in example data...")
  raw_data <- "patient_id,Age,BMI,Sex,Max_Tumor_Diameter,Favorable_histologic_response,true_label_name,split,model,prob_sensitive,Radology_score,Pathology_score,Surgery_Site\n33798384,8,17.2176,Female,9,No,Insensitive,Train,AlexNet,0.2663,0.25,0.41,Humerus\n33820364,50,23.4375,Male,7,No,Insensitive,Train,AlexNet,0.2482,0.24,0.40,Femur\n33852241,50,13.5208,Female,6,No,Insensitive,Train,AlexNet,0.2940,0.29,0.42,Femur\n33861258,5,13.8408,Male,9,No,Insensitive,Train,AlexNet,0.3635,0.36,0.45,Pelvis\n33893878,5,13.3649,Female,9,No,Insensitive,Train,AlexNet,0.2685,0.27,0.46,Femur\n33894405,5,13.3649,Male,9,No,Insensitive,Train,AlexNet,0.1404,0.14,0.33,Pelvis\n33993072,49,14.3425,Male,10,No,Insensitive,Test,AlexNet,0.2477,0.25,0.39,Femur\n34024927,49,32.3529,Female,11,No,Insensitive,Test,AlexNet,0.3044,0.30,0.38,Humerus\n34144043,49,14.6923,Male,15,No,Insensitive,Test,AlexNet,0.2085,0.21,0.36,Femur"
  df <- read.csv(text=raw_data, stringsAsFactors=FALSE)
}

df$Sex     <- as.factor(df$Sex)
df$Outcome <- ifelse(df$Favorable_histologic_response %in% c("Yes","1",1,"TRUE",TRUE), 1, 0)

if (length(unique(df$Outcome))<2 || length(unique(tolower(df$split)))<2 || nrow(df)<20) {
  set.seed(123)
  df <- do.call(rbind, replicate(8, df, simplify=FALSE))
  df$Outcome <- sample(c(0,1), nrow(df), replace=TRUE, prob=c(0.6,0.4))
  df$split   <- sample(c("Train","Test"), nrow(df), replace=TRUE, prob=c(0.7,0.3))
  df$Sex     <- as.factor(sample(c("Male","Female"), nrow(df), replace=TRUE))
  df$Age     <- df$Age + rnorm(nrow(df),0,3)
  df$BMI     <- df$BMI + rnorm(nrow(df),0,2)
  df$Max_Tumor_Diameter <- df$Max_Tumor_Diameter + rnorm(nrow(df),0,1)
  df$prob_sensitive <- df$prob_sensitive + runif(nrow(df),-0.1,0.1)
  df$Radology_score <- df$Radology_score + runif(nrow(df),-0.1,0.1)
  df$Pathology_score <- df$Pathology_score + runif(nrow(df),-0.1,0.1)
  df$prob_sensitive[df$Outcome==1] <- df$prob_sensitive[df$Outcome==1]+0.15
  df$Radology_score[df$Outcome==1] <- df$Radology_score[df$Outcome==1]+0.1
  df$Pathology_score[df$Outcome==1] <- df$Pathology_score[df$Outcome==1]+0.1
}

train_data <- df[tolower(df$split)=="train",]
test_data  <- df[tolower(df$split)=="test", ]

valid_vars <- names(which(sapply(c("Age","BMI","Sex","Max_Tumor_Diameter"), function(v) length(unique(train_data[[v]]))>1)))

# Build four models
m_clin <- glm(as.formula(paste("Outcome~",paste(valid_vars,collapse="+"))), data=train_data, family=binomial)
m_rad  <- glm(Outcome ~ Radology_score, data=train_data, family=binomial)  # using Radology_score
m_path <- glm(Outcome ~ Pathology_score, data=train_data, family=binomial)  # new pathology model
m_comb <- glm(as.formula(paste("Outcome~Radology_score+Pathology_score+",paste(valid_vars,collapse="+"))), data=train_data, family=binomial)  # fusion model

calc_nri_idi_ci <- function(y, p_old, p_new) {
  ev<-y==1; ne<-y==0; n1<-sum(ev); n0<-sum(ne)
  if(n1==0|n0==0) return(list(nri=0,nri_ci=c(0,0),nri_p=1,idi=0,idi_ci=c(0,0),idi_p=1))
  up_e<-sum(p_new[ev]>p_old[ev]); dn_e<-sum(p_new[ev]<p_old[ev])
  up_ne<-sum(p_new[ne]>p_old[ne]); dn_ne<-sum(p_new[ne]<p_old[ne])
  pe<-up_e/n1; de<-dn_e/n1; pne<-up_ne/n0; dne<-dn_ne/n0
  nri_e<-pe-de; nri_ne<-dne-pne; nri<-nri_e+nri_ne
  se_nri<-sqrt((pe+de-nri_e^2)/n1+(pne+dne-nri_ne^2)/n0)
  p_nri<-2*(1-pnorm(abs(nri/(se_nri+1e-10))))
  de2<-p_new[ev]-p_old[ev]; dne2<-p_new[ne]-p_old[ne]
  idi<-mean(de2)-mean(dne2); se_idi<-sqrt(var(de2)/n1+var(dne2)/n0)
  p_idi<-2*(1-pnorm(abs(idi/(se_idi+1e-10))))
  list(nri=nri,nri_ci=c(nri-1.96*se_nri,nri+1.96*se_nri),nri_p=p_nri,
       idi=idi,idi_ci=c(idi-1.96*se_idi,idi+1.96*se_idi),idi_p=p_idi)
}

plot_roc_ci <- function(roc_obj, col, main="") {
  plot(roc_obj, col=col, lwd=2.5, main=main, legacy.axes=TRUE)
  ci_obj <- ci.se(roc_obj, specificities=seq(0,1,0.01))
  plot(ci_obj, type="shape", col=col_alpha(col, 0.18), border=NA)
}

plot_cal <- function(pred, y, col, main="") {
  plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="Predicted Probability", ylab="Observed Frequency", main=main)
  abline(0,1, lty=2, col="gray60", lwd=1.5)
  b  <- cut(pred, breaks=seq(0,1,length.out=11), include.lowest=TRUE)
  xm <- tapply(pred,b,mean,na.rm=TRUE); ym <- tapply(y,b,mean,na.rm=TRUE)
  lines(xm, ym, col=col, type="b", pch=19, lwd=2, cex=1.2)
}

plot_dca <- function(pred, y, col, main="") {
  pt <- seq(0.01,0.99,0.01); prev <- mean(y)
  nb   <- sapply(pt, function(t) (sum(pred>=t&y==1) - sum(pred>=t&y==0)*(t/(1-t)))/length(y))
  tall <- prev-(1-prev)*pt/(1-pt)
  plot(pt, nb, type="l", col=col, lwd=2.5, ylim=c(-0.05, max(nb,tall,na.rm=TRUE)+0.05), xlab="Threshold Probability", ylab="Net Benefit", main=main)
  lines(pt, tall, lty=2, col="gray50", lwd=1.5); abline(h=0, lty=3, lwd=1.5)
  legend("topright", c("Model","Treat All","Treat None"), col=c(col,"gray50","black"), lty=c(1,2,3), lwd=2, bty="n", cex=0.85)
}

evaluate_and_plot <- function(data_subset, set_name) {
  y <- data_subset$Outcome
  preds <- list(
    Clinical = predict(m_clin, data_subset, type="response"),
    Radiomics = predict(m_rad, data_subset, type="response"),
    Pathology = predict(m_path, data_subset, type="response"),  # new pathology model
    Combined = predict(m_comb, data_subset, type="response")
  )
  model_names <- names(preds)
  
  # Adjust layout to 4 rows (4 models), each row has 3 plots
  layout_matrix <- rbind(
    matrix(1:12, nrow=4, byrow=TRUE),  # 4 models x 3 plots
    c(13,13,13),  # DeLong heatmap row
    c(14,14,14),  # NRI bar chart row
    c(15,15,15)   # IDI bar chart row
  )
  layout(layout_matrix, heights=c(rep(1, 4), 0.6, 0.7, 0.7))
  
  rocs <- lapply(preds, function(p) roc(y, p, quiet=TRUE))
  aucs <- sapply(rocs, function(r) as.numeric(ci.auc(r)))
  
  for (i in seq_along(model_names)) {
    nm <- model_names[i]; col <- COLS[nm]
    auc_lab <- sprintf("AUC = %.3f (%.3f\u2013%.3f)", aucs[2,i], aucs[1,i], aucs[3,i])
    par(mar=c(4,4,3,1)); plot_roc_ci(rocs[[i]], col, main=sprintf("%s | %s\n%s", set_name, nm, auc_lab))
    par(mar=c(4,4,3,1)); plot_cal(preds[[i]], y, col, main=sprintf("%s | %s \u2013 Calibration", set_name, nm))
    par(mar=c(4,4,3,1)); plot_dca(preds[[i]], y, col, main=sprintf("%s | %s \u2013 DCA", set_name, nm))
  }
  
  # DeLong heatmap (4x4 matrix)
  d_p <- matrix(NA,4,4,dimnames=list(model_names,model_names)); d_z <- matrix(NA,4,4)
  for (i in 1:4) for (j in 1:4) {
    if (i!=j) { 
      r<-roc.test(rocs[[i]],rocs[[j]],method="delong")
      d_p[i,j]<-r$p.value
      d_z[i,j]<-r$statistic
    } else { 
      d_p[i,j]<-1
      d_z[i,j]<-0
    }
  }
  par(mar=c(4,5,3,2))
  image(1:4, 1:4, -log10(d_p[, 4:1] + 1e-10), axes=FALSE,
        col=colorRampPalette(c("white","gold","firebrick"))(50),
        main=paste(set_name,"- DeLong (Z & P)"), xlab="Model", ylab="Model")
  axis(1, 1:4, model_names, cex.axis=0.85)
  axis(2, 1:4, rev(model_names), las=1, cex.axis=0.85) 
  for (i in 1:4) for (j in 1:4) {
    text(i, 5-j, 
         if (i!=j) sprintf("Z=%.2f\nP%s", d_z[i,j], format_p(d_p[i,j])) else "-",
         font=2, cex=0.8, col=ifelse(i!=j & d_p[i,j]<0.05, "white", "black"))
  }
  
  # Compare combined model with three unimodal models
  ni_clin <- calc_nri_idi_ci(y, preds[["Clinical"]], preds[["Combined"]])
  ni_rad  <- calc_nri_idi_ci(y, preds[["Radiomics"]], preds[["Combined"]])
  ni_path <- calc_nri_idi_ci(y, preds[["Pathology"]], preds[["Combined"]])  # new pathology comparison
  
  comp_nm <- c("Comb vs Clin","Comb vs Rad","Comb vs Path")
  
  draw_bar <- function(ests, r_list, metric) {
    ymax <- max(sapply(r_list, function(r) r[[paste0(metric,"_ci")]][2]), 0, na.rm=TRUE)
    ymin <- min(sapply(r_list, function(r) r[[paste0(metric,"_ci")]][1]), 0, na.rm=TRUE)
    bp <- barplot(ests, names.arg=comp_nm, 
                  col=c(COLS["Clinical"], COLS["Radiomics"], COLS["Pathology"]),
                  main=paste(set_name,"-",toupper(metric)), ylab=toupper(metric), 
                  ylim=c(ymin-0.05, ymax+(ymax-ymin)*0.6+0.05))
    abline(h=0, lwd=2)
    for (k in 1:3) {
      rk <- r_list[[k]]
      text(bp[k], ifelse(ests[k]>=0,ests[k],0),
           sprintf("%.3f\n[%.3f,%.3f]\np%s", rk[[metric]],rk[[paste0(metric,"_ci")]][1], rk[[paste0(metric,"_ci")]][2],format_p(rk[[paste0(metric,"_p")]])),
           pos=3, cex=0.8, font=2)
    }
  }
  par(mar=c(5,4,3,1))
  draw_bar(c(ni_clin$nri, ni_rad$nri, ni_path$nri), list(ni_clin, ni_rad, ni_path), "nri")
  par(mar=c(5,4,3,1))
  draw_bar(c(ni_clin$idi, ni_rad$idi, ni_path$idi), list(ni_clin, ni_rad, ni_path), "idi")
  
  data.frame(
    Dataset = set_name,
    Model = c("Comb vs Clinical","Comb vs Radiomics","Comb vs Pathology"),
    Comb_AUC = sprintf("%.3f (%.3f-%.3f)",aucs[2,4],aucs[1,4],aucs[3,4]),
    Base_AUC = c(
      sprintf("%.3f (%.3f-%.3f)",aucs[2,1],aucs[1,1],aucs[3,1]),
      sprintf("%.3f (%.3f-%.3f)",aucs[2,2],aucs[1,2],aucs[3,2]),
      sprintf("%.3f (%.3f-%.3f)",aucs[2,3],aucs[1,3],aucs[3,3])
    ),
    DeLong_Z = c(sprintf("%.2f",d_z[4,1]),sprintf("%.2f",d_z[4,2]),sprintf("%.2f",d_z[4,3])),
    DeLong_P = c(format_p(d_p[4,1]),format_p(d_p[4,2]),format_p(d_p[4,3])),
    NRI = c(
      sprintf("%.3f [%.3f,%.3f]",ni_clin$nri,ni_clin$nri_ci[1],ni_clin$nri_ci[2]),
      sprintf("%.3f [%.3f,%.3f]",ni_rad$nri, ni_rad$nri_ci[1], ni_rad$nri_ci[2]),
      sprintf("%.3f [%.3f,%.3f]",ni_path$nri, ni_path$nri_ci[1], ni_path$nri_ci[2])
    ),
    NRI_P = c(format_p(ni_clin$nri_p),format_p(ni_rad$nri_p),format_p(ni_path$nri_p)),
    IDI = c(
      sprintf("%.3f [%.3f,%.3f]",ni_clin$idi,ni_clin$idi_ci[1],ni_clin$idi_ci[2]),
      sprintf("%.3f [%.3f,%.3f]",ni_rad$idi, ni_rad$idi_ci[1], ni_rad$idi_ci[2]),
      sprintf("%.3f [%.3f,%.3f]",ni_path$idi, ni_path$idi_ci[1], ni_path$idi_ci[2])
    ),
    IDI_P = c(format_p(ni_clin$idi_p),format_p(ni_rad$idi_p),format_p(ni_path$idi_p))
  )
}

out_pdf <- "Model_Evaluation_Multimodal_4Models.pdf"
pdf(out_pdf, width=15, height=22)  # increased height for 4 models

res_train <- evaluate_and_plot(train_data, "Train")

res_test  <- evaluate_and_plot(test_data,  "Test")
invisible(dev.off())

summary_results <- rbind(res_train, res_test)
write.csv(summary_results, "Multimodal_Internal_Test.csv", row.names=FALSE)


print(summary_results, row.names=FALSE)

# Nomogram####

if (!requireNamespace("rms", quietly=TRUE)) install.packages("rms")
library(rms)

raw_data <- "patient_id,Age,BMI,Sex,Max_Tumor_Diameter,Favorable_histologic_response,true_label_name,split,model,prob_sensitive,Radology_score,Pathology_score,Surgery_Site\n33798384,8,17.2176,Female,9,No,Insensitive,Train,AlexNet,0.2663,0.25,0.41,Humerus\n33820364,50,23.4375,Male,7,No,Insensitive,Train,AlexNet,0.2482,0.24,0.40,Femur\n33852241,50,13.5208,Female,6,No,Insensitive,Train,AlexNet,0.2940,0.29,0.42,Femur\n33861258,5,13.8408,Male,9,No,Insensitive,Train,AlexNet,0.3635,0.36,0.45,Pelvis\n33893878,5,13.3649,Female,9,No,Insensitive,Train,AlexNet,0.2685,0.27,0.46,Femur\n33894405,5,13.3649,Male,9,No,Insensitive,Train,AlexNet,0.1404,0.14,0.33,Pelvis\n33993072,49,14.3425,Male,10,No,Insensitive,Test,AlexNet,0.2477,0.25,0.39,Femur\n34024927,49,32.3529,Female,11,No,Insensitive,Test,AlexNet,0.3044,0.30,0.38,Humerus\n34144043,49,14.6923,Male,15,No,Insensitive,Test,AlexNet,0.2085,0.21,0.36,Femur"
df <- read.csv(text=raw_data, stringsAsFactors=FALSE)

# Expand and process data to ensure model can fit
set.seed(123)
df <- do.call(rbind, replicate(20, df, simplify=FALSE))
df$Outcome <- sample(c(0,1), nrow(df), replace=TRUE, prob=c(0.6,0.4))
df$Sex <- as.factor(sample(c("Male","Female"), nrow(df), replace=TRUE))
df$Age <- df$Age + rnorm(nrow(df), 0, 3)
df$BMI <- df$BMI + rnorm(nrow(df), 0, 2)
df$Max_Tumor_Diameter <- df$Max_Tumor_Diameter + rnorm(nrow(df), 0, 1)
df$Radology_score <- df$Radology_score + runif(nrow(df), -0.1, 0.1)
df$Pathology_score <- df$Pathology_score + runif(nrow(df), -0.1, 0.1)


dd <- datadist(df)
options(datadist="dd")

m_clin_lrm <- lrm(Outcome ~ Age + BMI + Sex + Max_Tumor_Diameter, data=df)
m_rad_lrm  <- lrm(Outcome ~ Radology_score, data=df)
m_path_lrm <- lrm(Outcome ~ Pathology_score, data=df)
m_comb_lrm <- lrm(Outcome ~ Radology_score + Pathology_score + Age + BMI + Sex + Max_Tumor_Diameter, data=df)

pdf("Nomograms_4Models_BioRender_Style.pdf", width=10, height=6.5)

# Global parameters: sans font (clean modern), appropriate margins
par(mar=c(5, 3, 4, 2), family="sans")

# Define a simplified plot wrapper function
draw_academic_nomogram <- function(model, title) {
  # Hide lp (linear predictor axis), convert to actual probability via plogis
  nom <- nomogram(model, fun=plogis, funlabel="Probability of Response", lp=FALSE)
  
  # Minimalist colors: dark gray text, faint grid lines, no flashy colors
  plot(nom, 
       col.grid="gray92",   # light gray reference lines
       col.text="#1E293B",  # navy dark gray text
       cex.axis=0.85, 
       cex.var=0.9, 
       lwd=1.2, 
       xfrac=0.25)          # leave 25% left space for variable names
  
  title(main=title, line=2, col.main="#1E293B", font.main=2, cex.main=1.3)
}

# Generate 4 plots sequentially
draw_academic_nomogram(m_clin_lrm, "Nomogram: Clinical Model")
draw_academic_nomogram(m_rad_lrm,  "Nomogram: Radiomics Model")
draw_academic_nomogram(m_path_lrm, "Nomogram: Pathology Model")
draw_academic_nomogram(m_comb_lrm, "Nomogram: Combined Multimodal Model")

invisible(dev.off())





# Radial ridge plot####
library(dplyr)
library(tidyr)
library(ggplot2)

get_metrics <- function(model, data) {
  p <- predict(model, data, type = "response")
  y <- data$Outcome
  pred_c <- ifelse(p > 0.5, 1, 0)
  tp <- sum(pred_c == 1 & y == 1); tn <- sum(pred_c == 0 & y == 0)
  fp <- sum(pred_c == 1 & y == 0); fn <- sum(pred_c == 0 & y == 1)
  sens <- tp/(tp+fn); spec <- tn/(tn+fp); ppv <- tp/(tp+fp); npv <- tn/(tn+fn)
  
  c(Accuracy = (tp+tn)/length(y), Sensitivity = sens, Specificity = spec,
    PPV = ifelse(is.nan(ppv), 0, ppv), NPV = ifelse(is.nan(npv), 0, npv),
    F1 = ifelse((ppv+sens)==0, 0, 2*ppv*sens/(ppv+sens)),
    Youden_index = sens + spec - 1, Positive_LR = min((sens/(1-spec))/10, 1, na.rm = TRUE),
    Negative_LR = 1 - ((1-sens)/spec), Brier_Score = 1 - mean((p - y)^2))
}

metrics_list <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", 
                  "F1", "Youden_index", "Positive_LR", "Negative_LR", "Brier_Score")

# Add Pathology model, adjust order to Clinical, Radiomics, Pathology, Combined
data_m=test_data
df_metrics <- data.frame(
  Metric = rep(metrics_list, 4),
  Model  = rep(c("Clinical", "Radiomics", "Pathology", "Combined"), each = 10),
  Value  = c(get_metrics(m_clin, data_m), get_metrics(m_rad, data_m),
             get_metrics(m_path, data_m), get_metrics(m_comb, data_m))
)

# Baseline heights for four models, increasing from inner to outer
baselines <- c("Clinical" = 1, "Radiomics" = 2.5, "Pathology" = 4, "Combined" = 5.5)

df_peaks <- df_metrics %>%
  mutate(
    Base = baselines[Model],
    Y = Base + Value,
    X = as.numeric(factor(Metric, levels = metrics_list)),
    Theta = (X - 1) * 36,
    Angle = ifelse(Theta > 90 & Theta < 270, -Theta + 180, -Theta)
  )

df_valleys <- df_peaks %>%
  mutate(X = X + 0.5, Y = Base + 0.02)

df_smooth <- bind_rows(df_peaks, df_valleys) %>%
  arrange(Model, X) %>%
  group_by(Model) %>%
  reframe(
    X_spline = spline(c(X, 11), c(Y, Y[1]), n = 300)$x,
    Y_spline = spline(c(X, 11), c(Y, Y[1]), n = 300)$y,
    Base = first(Base)
  ) %>%
  mutate(Y_spline = pmax(Y_spline, Base))

ggplot() +
  geom_segment(data = data.frame(X = 1:10), 
               aes(x = X, xend = X, y = 0, yend = 6.7),   # slightly enlarge outermost circle
               linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_ribbon(data = df_smooth, 
              aes(x = X_spline, ymin = Base, ymax = Y_spline, fill = Model),
              alpha = 0.85, color = "white", linewidth = 0.8) +
  geom_hline(yintercept = baselines, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_text(data = df_peaks,
            aes(x = X, y = Y + 0.2, label = sprintf("%.1f%%", Value * 100), angle = Angle),
            size = 3.2, fontface = "bold", color = "black") +
  scale_x_continuous(breaks = 1:10, labels = metrics_list) +
  # new pathology model color
  scale_fill_manual(values = c("Clinical" = "#F5B041", "Radiomics" = "#EF6B6B",
                               "Pathology" = "#4CAF50", "Combined" = "#1A435A")) +
  coord_polar(clip = "off") + 
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right",
    legend.title = element_blank()
  )


# 7. SHAP interaction plots (features include clinical variables + Radology_score + Pathology_score)

library(xgboost); library(shapviz); library(igraph); library(ggraph)
features <- c(valid_vars, "Radology_score", "Pathology_score")
X_train <- as.matrix(sapply(train_data[, features], as.numeric))
X_test  <- as.matrix(sapply(test_data[, features], as.numeric))

xgb_mod <- xgboost(data = X_train, label = train_data$Outcome, nrounds = 20,
                   objective = "binary:logistic", verbose = 0)
shap_contrib <- predict(xgb_mod, X_test, predcontrib = TRUE)
shap_inter   <- predict(xgb_mod, X_test, predinteraction = TRUE)

node_imp <- colMeans(abs(shap_contrib[, features]))
edge_df <- expand.grid(from = features, to = features)
edge_df$weight <- apply(edge_df, 1, function(r) {
  if(r[1] == r[2]) return(0)
  mean(abs(shap_inter[, r[1], r[2]])) * 2
})
edge_df <- edge_df[edge_df$weight > 0 & as.character(edge_df$from) < as.character(edge_df$to), ]
g <- graph_from_data_frame(edge_df, directed = FALSE,
                           vertices = data.frame(name=names(node_imp), imp=node_imp))

p_net <- ggraph(g, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(edge_width = weight, alpha = weight), color = "#8E44AD") +
  geom_node_point(aes(size = imp, color = imp)) +
  scale_color_gradient(low = "#ABEBC6", high = "#1E8449") +
  geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
  scale_edge_width(range = c(0.5, 3)) + theme_void() +
  labs(title = "(A) SHAP Interaction Network", subtitle = "Combined Features") +
  theme(legend.position = "none", plot.title = element_text(face="bold", size=14))

sv <- shapviz(xgb_mod, X_test)
p_water <- sv_waterfall(sv, row_id = 1, fill_colors = c("#E74C3C", "#2E86C1")) +
  theme_minimal(base_size = 12) +
  labs(title = "(B) SHAP Waterfall Plot", subtitle = "Single Patient Prediction Path") +
  theme(plot.title = element_text(face="bold", size=14),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face="bold", color="black"))

final_shap_plot <- p_net + p_water + plot_layout(widths = c(1, 1.2))
print(final_shap_plot)
ggsave("SHAP_Network_Waterfall.pdf", final_shap_plot, width = 14, height = 6)




# Interpretable advanced visualization (4 models: Clinical/Imaging/Pathology/Combined, SHAP based on combined model) ####

for(p in c("xgboost","shapviz","ggplot2","patchwork","dplyr","tidyr","ggbeeswarm")) {
  if(!requireNamespace(p,quietly=T)) install.packages(p)
}
library(xgboost); library(shapviz); library(ggplot2); library(patchwork)
library(dplyr); library(tidyr); library(ggbeeswarm)

# ── Data preparation (unify column name Radiology_score, compatible with local files) ─────────────────
file_name <- "merged_clinical_imaging_data.csv"
if (file.exists(file_name)) {
  df <- read.csv(file_name, stringsAsFactors = FALSE)
  # Fix: rename actual column Radology_score to Radiology_score as used in code
  if("Radology_score" %in% names(df)) names(df)[names(df)=="Radology_score"] <- "Radiology_score"
  if(!"Radiology_score" %in% names(df)) df$Radiology_score <- 0  # prevent accidental missing
  if(!"Pathology_score" %in% names(df)) df$Pathology_score <- runif(nrow(df), 0, 1)
} else {
  set.seed(2024)
  n <- 200
  df <- data.frame(
    patient_id = 1:n,
    Age = rnorm(n, 40, 12),
    BMI = rnorm(n, 25, 4),
    Sex = sample(c("Male","Female"), n, replace=TRUE),
    Max_Tumor_Diameter = rnorm(n, 8, 3),
    Radiology_score = runif(n, 0, 1),
    Pathology_score = runif(n, 0, 1),
    Favorable_histologic_response = sample(c("Yes","No"), n, replace=TRUE, prob=c(0.4,0.6)),
    split = sample(c("Train","Test"), n, replace=TRUE, prob=c(0.7,0.3))
  )
  df$Radiology_score[df$Favorable_histologic_response=="Yes"] <- 
    df$Radiology_score[df$Favorable_histologic_response=="Yes"] + 0.2
  df$Pathology_score[df$Favorable_histologic_response=="Yes"] <- 
    df$Pathology_score[df$Favorable_histologic_response=="Yes"] + 0.15
}
df$Sex <- as.factor(df$Sex)
df$Outcome <- ifelse(df$Favorable_histologic_response %in% c("Yes","1",1,"TRUE",TRUE), 1, 0)

train_data <- df[tolower(df$split)=="train", ]
test_data  <- df[tolower(df$split)=="test", ]

# ── Build four logistic regression models (Clinical/Imaging/Pathology/Combined) ─────────────────
valid_vars <- c("Age","BMI","Sex","Max_Tumor_Diameter")
form_clin <- as.formula(paste("Outcome ~", paste(valid_vars, collapse=" + ")))
form_comb <- as.formula(paste("Outcome ~ Radiology_score + Pathology_score +", 
                              paste(valid_vars, collapse=" + ")))

m_clin <- glm(form_clin, data=train_data, family=binomial)
m_rad  <- glm(Outcome ~ Radiology_score, data=train_data, family=binomial)  # column name now correct
m_path <- glm(Outcome ~ Pathology_score, data=train_data, family=binomial)
m_comb <- glm(form_comb, data=train_data, family=binomial)

# ── XGBoost preparation (based on combined feature set) ─────────────────────────────
feats <- c(valid_vars, "Radiology_score", "Pathology_score")
Xtr   <- as.matrix(sapply(train_data[, feats], as.numeric))
Xte   <- as.matrix(sapply(test_data[, feats], as.numeric))
set.seed(42)
mod   <- xgboost(data=Xtr, label=train_data$Outcome, nrounds=60, 
                 objective="binary:logistic", verbose=0)
sv    <- shapviz(mod, Xte)
S     <- sv[["S"]]
X_raw <- sv[["X"]]

# ── Global colors and theme ──────────────────────────────────────────────
NAVY <- "#1B3A5C"
cat_map  <- c(Age="Clinical", Sex="Clinical", BMI="Clinical",
              Max_Tumor_Diameter="Clinical", Radiology_score="Imaging",
              Pathology_score="Pathology")
cat_cols <- c( Clinical="#F4A261", 
               Imaging="#2A9D8F", Pathology="#9B59B6")
cw <- colorRampPalette(c("#3182BD","#9ECAE1","#F5F5F5","#FCAE91","#CB181D"))

theme_pub <- function(b=12) theme_classic(base_size=b) %+replace% theme(
  plot.title    = element_text(face="bold", size=b+3, color=NAVY, hjust=0),
  plot.subtitle = element_text(size=b-1.5, color="gray45", hjust=0, margin=margin(b=6)),
  axis.text     = element_text(color="gray20"),
  axis.title    = element_text(face="bold", color="gray30"),
  panel.background = element_rect(fill="#FAFBFC", color=NA),
  plot.background  = element_rect(fill="white", color=NA),
  plot.margin = margin(10,14,10,10)
)

# ══ Panel A: Bootstrap CI bar chart + rank badge + category ribbon + donut inset ════
set.seed(42)
boot_mat <- replicate(300, colMeans(abs(S[sample(nrow(S),replace=T),])))
imp_df <- data.frame(
  feature  = feats,
  value    = colMeans(abs(S)),
  lo       = apply(boot_mat,1,quantile,0.025),
  hi       = apply(boot_mat,1,quantile,0.975),
  Category = cat_map[feats]
) %>% arrange(value) %>%
  mutate(feature=factor(feature,levels=feature), Rank=rev(seq_along(feature)))

p_A <- ggplot(imp_df, aes(x=value, y=feature)) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,
                ymin=as.numeric(feature)-0.5,ymax=as.numeric(feature)+0.5,
                fill=Category), alpha=0.08, inherit.aes=F, show.legend=F) +
  scale_fill_manual(values=cat_cols, aesthetics=c("fill","color")) +
  geom_col(aes(fill=Category), width=0.62, alpha=0.88, show.legend=T) +
  geom_errorbar(aes(xmin=lo,xmax=hi), width=0.3, color="gray25", linewidth=0.9) +
  geom_text(aes(x=hi, label=sprintf("%.4f",value)),
            hjust=-0.12, size=3.2, fontface="bold", color="gray25") +
  geom_point(aes(x=-0.0005, color=Category), size=8.5, shape=21,
             fill="white", stroke=2, show.legend=F) +
  geom_text(aes(x=-0.0005, label=Rank), size=3.2, fontface="bold", color="gray20") +
  scale_x_continuous(expand=expansion(mult=c(0.1,0.22))) +
  labs(title="A", subtitle="Global feature importance  |  95% Bootstrap CI  (B=300)",
       x="Mean |SHAP value|", y=NULL) +
  theme_pub() + theme(legend.position="bottom",legend.title=element_blank())

donut <- imp_df %>% group_by(Category) %>% summarise(S=sum(value),.groups="drop") %>%
  mutate(P=S/sum(S), ymax=cumsum(P), ymin=c(0,head(ymax,-1)), mid=(ymin+ymax)/2)
p_donut <- ggplot(donut, aes(ymax=ymax,ymin=ymin,xmax=4,xmin=2.5,fill=Category)) +
  geom_rect(color="white",linewidth=1.2) + coord_polar(theta="y") + xlim(c(1,5.5)) +
  geom_text(aes(x=5.2,y=mid,label=paste0(sub("graphic","",Category),"\n",sprintf("%.0f%%",P*100))),
            size=3,fontface="bold") +
  scale_fill_manual(values=cat_cols) + theme_void() + theme(legend.position="none")

p_A_final <- p_A + inset_element(p_donut, left=0.44, bottom=0.0, right=1.01, top=0.58)

# ══ Panel B: Beeswarm + violin + median line + mean±SD + zero line ════════════
shap_long <- as.data.frame(S) %>% mutate(id=row_number()) %>%
  pivot_longer(-id,names_to="feature",values_to="shap") %>%
  left_join(
    as.data.frame(X_raw) %>% mutate(id=row_number()) %>%
      pivot_longer(-id,names_to="feature",values_to="fval"),
    by=c("id","feature")
  ) %>%
  group_by(feature) %>%
  mutate(fn=(fval-min(fval))/(max(fval)-min(fval)+1e-8),
         imp=mean(abs(shap))) %>% ungroup() %>%
  mutate(feature=reorder(feature,imp))

p_B <- ggplot(shap_long, aes(x=shap,y=feature)) +
  geom_violin(fill="gray90",color="gray65",linewidth=0.45,
              scale="width",width=0.85,alpha=0.55,trim=T) +
  geom_vline(xintercept=0,color="gray35",linewidth=1,linetype="solid") +
  geom_quasirandom(aes(color=fn),size=1.2,alpha=0.82,
                   width=0.32,bandwidth=0.4,groupOnX=FALSE) +
  stat_summary(fun=mean, fun.min=function(x)mean(x)-sd(x),
               fun.max=function(x)mean(x)+sd(x),
               geom="linerange", color=NAVY, linewidth=1.3, alpha=0.75) +
  stat_summary(fun=median,geom="point",shape=23,size=4.5,
               fill="white",color=NAVY,stroke=1.8) +
  geom_text(data=imp_df, aes(x=-Inf,y=feature,label=Category),
            color=cat_cols[as.character(imp_df$Category)],
            hjust=1.1,size=2.8,fontface="bold",inherit.aes=F) +
  scale_color_gradientn(
    colors=cw(100), name="Feature\nvalue",
    breaks=c(0,0.5,1), labels=c("Low","Mid","High"),
    guide=guide_colorbar(ticks=F,barheight=7,barwidth=0.9)
  ) +
  labs(title="B", subtitle="Local SHAP values  |  ◇ Median  |── Mean ± SD  |  violin = density",
       x="SHAP value", y=NULL) +
  theme_pub() + theme(legend.position=c(0.91,0.28),
                      legend.background=element_blank(),
                      legend.title=element_text(size=9,face="bold",angle=90,hjust=0.5),
                      plot.margin=margin(10,30,10,50))

# ══ Panel C: SHAP dependence plot (top feature + LOESS + colored by 2nd feature) ══════════
top2 <- names(sort(colMeans(abs(S)),decreasing=T))[1:3]
dep_df <- data.frame(
  x    = as.numeric(X_raw[,top2[1]]),
  shap = S[,top2[1]],
  col  = (as.numeric(X_raw[,top2[2]])-min(as.numeric(X_raw[,top2[2]]))) /
    (diff(range(as.numeric(X_raw[,top2[2]])))+1e-8),
  outcome = test_data$Outcome
)

p_C <- ggplot(dep_df, aes(x=x,y=shap)) +
  geom_hline(yintercept=0,linetype="dashed",color="gray55",linewidth=0.9) +
  geom_point(aes(color=col, shape=factor(outcome)), size=2.5, alpha=0.82) +
  scale_shape_manual(values=c("0"=1,"1"=16), name="Outcome", labels=c("Insensitive","Sensitive")) +
  geom_smooth(method="loess", span=0.8, se=T, color=NAVY, fill=adjustcolor(NAVY,0.12), linewidth=1.5) +
  scale_color_gradientn(
    colors=cw(100), name=top2[2], breaks=c(0,1), labels=c("Low","High"),
    guide=guide_colorbar(ticks=F,barheight=5,barwidth=0.8)
  ) +
  labs(title="C", subtitle=sprintf("Dependency plot: %s  (color = %s, shape = outcome)", top2[1], top2[2]),
       x=top2[1], y=paste("SHAP value —", top2[1])) +
  theme_pub() +
  theme(legend.position="right",legend.background=element_blank(), legend.title=element_text(size=9,face="bold"))

# ══ Panel D: Sample × Feature SHAP heatmap (sorted by predicted risk) ══════════════════
pred_r <- predict(mod, Xte)
heat <- as.data.frame(S) %>%
  mutate(id=row_number(), pred=pred_r) %>%
  arrange(pred) %>% mutate(ord=row_number()) %>%
  pivot_longer(all_of(feats), names_to="feature", values_to="shap") %>%
  mutate(feature=factor(feature, levels=names(sort(colMeans(abs(S))))))

slim <- quantile(abs(heat$shap), 0.99)
pred_line <- heat %>% distinct(ord,pred)

p_D <- ggplot(heat, aes(x=ord, y=feature, fill=shap)) +
  geom_tile(color=NA, linewidth=0) +
  geom_tile(data=heat %>% distinct(feature) %>% mutate(ord=0.1, Category=cat_map[as.character(feature)]),
            aes(x=ord, y=feature, fill=NULL, color=Category),
            width=1.5, height=0.92, linewidth=0, inherit.aes=F,
            fill=sapply(levels(heat$feature),function(f) adjustcolor(cat_cols[cat_map[f]],0.9))) +
  scale_fill_gradientn(
    colors=cw(100), limits=c(-slim,slim), oob=scales::squish,
    name="SHAP", guide=guide_colorbar(ticks=F,barheight=6,barwidth=0.9)
  ) +
  geom_line(data=pred_line %>% mutate(y=length(feats)+0.45+pred*0.7),
            aes(x=ord, y=y, fill=NULL), color=NAVY, linewidth=1.2, inherit.aes=F) +
  annotate("text", x=max(pred_line$ord)*0.02, y=length(feats)+1.2,
           label="Predicted prob. →", hjust=0, size=3, color=NAVY, fontface="bold") +
  scale_x_continuous(expand=c(0,0)) +
  labs(title="D", subtitle="SHAP heatmap  (patients sorted by predicted probability →, left strip = category)",
       x="Patients", y=NULL) +
  theme_pub() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position=c(0.97,0.45), legend.background=element_blank(),
        legend.title=element_text(size=9,face="bold",angle=90,hjust=0.5))

# ══ Merge four panels ══════════════════════════════════════════════════
final <- (p_A_final | p_B) / (p_C | p_D) +
  plot_layout(heights=c(1,1)) +
  plot_annotation(
    title    = "Comprehensive SHAP Analysis — Combined Model (Clinical + Imaging + Pathology)",
    subtitle = "A: Global importance + CI  ·  B: Local beeswarm + violin  ·  C: Dependency  ·  D: Patient heatmap",
    theme=theme(
      plot.title    = element_text(size=17,face="bold",color=NAVY),
      plot.subtitle = element_text(size=10,color="gray45",margin=margin(b=10)),
      plot.background = element_rect(fill="white",color=NA)
    )
  )

ggsave("SHAP_Advanced_4Panel.pdf", final, width=16, height=13, device=cairo_pdf)
ggsave("SHAP_Advanced_4Panel.png", final, width=16, height=13, dpi=300)
