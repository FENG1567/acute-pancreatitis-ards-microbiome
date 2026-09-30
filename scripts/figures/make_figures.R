#!/usr/bin/env Rscript
# Submission figures: R is the exclusive drawing, preview, and export backend.
# This script only re-expresses validated results already present in
# analysis_results; it does not reanalyse, subset, or alter any scientific result.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/figures/make_figures.R REPOSITORY_ROOT")
root <- normalizePath(args[1], mustWork = TRUE)
dat <- file.path(root, "data", "published_results")
out <- file.path(root, "results", "figures")
src <- file.path(root, "results", "figure_source_data_generated")
dir.create(out, recursive = TRUE, showWarnings = FALSE); dir.create(src, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(name) read.delim(file.path(dat, name), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, name) write.table(x, file.path(src, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
# Shared visual grammar: ARDS signal = muted red; reference = blue; secondary = teal.
ink <- "#1F2933"; blue <- "#2C6E9E"; red <- "#B64E4A"; teal <- "#2A8C82"; orange <- "#D0843A"; gray <- "#67737D"; lightgray <- "#E9EEF2"; pale <- "#EAF2F7"; mint <- "#E7F3F0"; rose <- "#F8E7E5"; sand <- "#F8F2E8"; white <- "#FFFFFF"
font_family <- "Arial"; journal_width_mm <- 183; body_text_pt <- 6.5
open_dev <- function(path, fmt, width, height) {
  if (fmt == "pdf") grDevices::cairo_pdf(path, width = width, height = height, family = font_family)
  if (fmt == "svg") svglite::svglite(path, width = width, height = height, system_fonts = list(sans = font_family))
  if (fmt == "png") ragg::agg_png(path, width = width, height = height, units = "in", res = 400)
  if (fmt == "tiff") ragg::agg_tiff(path, width = width, height = height, units = "in", res = 600, compression = "lzw")
}
export_fig <- function(stem, draw, width = 7.205, height = 5.25) {
  formats <- "pdf"
  if (requireNamespace("svglite", quietly = TRUE)) formats <- c(formats, "svg")
  if (requireNamespace("ragg", quietly = TRUE)) formats <- c(formats, "png", "tiff")
  for (fmt in formats) { ext <- if (fmt == "tiff") "tif" else fmt; open_dev(file.path(out, paste0(stem, ".", ext)), fmt, width, height); par(family = font_family, fg = ink, col.axis = ink, col.lab = ink, col.main = ink, cex.axis = 0.74, cex.lab = 0.82, xaxs = "i", yaxs = "i"); draw(); dev.off() }
}
panel <- function(letter) { u <- par("usr"); text(u[1] + 0.012 * diff(u[1:2]), u[4] - 0.025 * diff(u[3:4]), letter, adj = c(0, 1), font = 2, cex = 0.78, xpd = NA) }
title_left <- function(txt, cex = 0.86) title(main = txt, adj = 0, cex.main = cex, font.main = 2, line = 0.55)
draw_ci <- function(est, low, high, y, col = ink, pchcol = teal, lwd = 1.7) { segments(low, y, high, y, col = col, lwd = lwd); points(est, y, pch = 21, bg = pchcol, col = white, cex = 1.25) }

# FIGURE 1 | Evidence hierarchy, primary-cohort flow, and replication eligibility.
eligibility <- read_tsv("public_dataset_ards_eligibility.tsv"); write_tsv(eligibility, "figure1_public_ards_eligibility.tsv")
draw1 <- function() {
  layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), heights = c(0.95, 1.05), widths = c(1.02, 0.98))
  par(mar = c(0.7, 0.5, 1.7, 0.5)); plot.new(); plot.window(c(0, 1), c(0, 1)); panel("a"); title_left("Evidence architecture and analytic separation")
  card <- function(x0, y0, x1, y1, fill, accent, tag, heading, detail, detail_top = y1 - 0.24, heading_cex = 0.70, detail_cex = 0.57) {
    rect(x0, y0, x1, y1, col = fill, border = lightgray, lwd = 1.0)
    rect(x0, y1 - 0.045, x1, y1, col = accent, border = NA)
    text(x0 + 0.018, y1 - 0.085, tag, adj = c(0, 1), cex = 0.45, font = 2, col = accent)
    text(x0 + 0.018, y1 - 0.175, heading, adj = c(0, 1), cex = heading_cex, font = 2, col = ink)
    text(x0 + 0.018, detail_top, detail, adj = c(0, 1), cex = detail_cex, col = gray)
  }
  card(0.035, 0.30, 0.235, 0.76, pale, blue, "DERIVATION", "Fixed score", "Two public cohorts\nFive genera; weights locked", detail_top = 0.47)
  card(0.315, 0.23, 0.610, 0.83, rose, red, "PRIMARY OUTCOME TEST", "Recorded ARDS", "65 unique AP patients\n26 ARDS  |  39 non-ARDS", detail_top = 0.46, heading_cex = 0.76, detail_cex = 0.62)
  card(0.705, 0.55, 0.965, 0.83, mint, teal, "SUPPORTIVE CONTEXT", "Paired compartments", "62 blood-neutrophil pairs\nRelative-score evidence only", detail_top = 0.62, heading_cex = 0.64, detail_cex = 0.52)
  card(0.705, 0.18, 0.965, 0.46, "#EEF0F8", blue, "SUPPORTIVE CONTEXT", "Host response", "87 AP transcriptomes\nSeverity background only", detail_top = 0.25, heading_cex = 0.64, detail_cex = 0.52)
  arrows(0.245, 0.53, 0.298, 0.53, length = 0.07, lwd = 1.4, col = gray)
  segments(0.625, 0.53, 0.662, 0.53, col = gray, lwd = 1.0, lty = 2)
  segments(0.662, 0.325, 0.662, 0.690, col = gray, lwd = 1.0, lty = 2)
  segments(0.662, 0.690, 0.695, 0.690, col = gray, lwd = 1.0, lty = 2)
  segments(0.662, 0.325, 0.695, 0.325, col = gray, lwd = 1.0, lty = 2)
  text(0.655, 0.12, "Separate datasets; no participant-level linkage", cex = 0.50, col = gray, adj = c(0.5, 0.5))
  par(mar = c(2.9, 0.8, 2.3, 0.8)); plot.new(); plot.window(c(0, 1), c(0, 1)); panel("b")
  rect(0.05, 0.53, 0.29, 0.78, col = pale, border = gray, lwd = 1); text(0.17, 0.68, "85 public runs", font = 2, cex = 0.75); text(0.17, 0.59, "PRJNA893348", cex = 0.60); arrows(0.30, 0.655, 0.40, 0.655, length = 0.08, col = gray, lwd = 1.2)
  rect(0.41, 0.53, 0.62, 0.78, col = mint, border = gray, lwd = 1); text(0.515, 0.68, "65 unique AP", font = 2, cex = 0.75); text(0.515, 0.59, "patient-level score", cex = 0.60); arrows(0.63, 0.655, 0.72, 0.78, length = 0.08, col = gray, lwd = 1.2); arrows(0.63, 0.655, 0.73, 0.31, length = 0.08, col = gray, lwd = 1.2)
  rect(0.73, 0.67, 0.95, 0.89, col = rose, border = red, lwd = 1); text(0.84, 0.80, "Recorded ARDS", font = 2, cex = 0.72); text(0.84, 0.72, "n = 26", cex = 0.69); rect(0.73, 0.31, 0.95, 0.53, col = pale, border = blue, lwd = 1); text(0.84, 0.44, "Recorded non-ARDS", font = 2, cex = 0.70); text(0.84, 0.36, "n = 39", cex = 0.69); text(0.50, 0.14, "Sampling within 24 h of AP onset; individual sampling-to-ARDS intervals unavailable.", cex = 0.48, col = gray)
  par(mar = c(3.5, 6.1, 2.3, 1.8))
  cohorts <- eligibility$dataset; cols <- c("Human AP", "Patient map", "ARDS label", "Independent\nreplication"); m <- matrix(0, nrow = length(cohorts), ncol = 4); m[, 1] <- 1; m[, 2] <- ifelse(grepl("patients|aliases|mapping", eligibility$patient_mapping, ignore.case = TRUE), 1, 0); m[, 3] <- ifelse(grepl("Direct AP-ARDS", eligibility$public_outcome_label), 1, 0); m[, 4] <- ifelse(eligibility$independent_AP_ARDS_replication_eligible == "Yes", 1, 0)
  image(1:4, seq_along(cohorts), t(m[nrow(m):1, ]), col = c("#F3F5F6", teal), axes = FALSE, xlab = "", ylab = "", zlim = c(0, 1)); axis(1, 1:4, cols, tick = FALSE, cex.axis = 0.57); axis(2, seq_along(cohorts), rev(cohorts), las = 1, tick = FALSE, cex.axis = 0.63)
  for (i in seq_len(nrow(m))) for (j in 1:4) text(j, length(cohorts) - i + 1, if (m[i, j]) "Y" else "N", cex = 0.64, col = if (m[i,j]) white else gray); abline(v = seq(0.5, 4.5, 1), h = seq(0.5, length(cohorts) + 0.5, 1), col = white, lwd = 1.1); mtext("c", side = 3, adj = -0.12, line = 0.20, font = 2, cex = 0.78)
}
export_fig("Figure_1_study_design", draw1, 7.205, 5.35)

# FIGURE 2 | Score construction and score-component architecture.
comp <- read_tsv("primary_component_details.tsv"); write_tsv(comp, "figure2_component_details.tsv")
draw2 <- function() {
  layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.03, 0.97), heights = c(1, 1))
  par(mar = c(2.5, 0.8, 2.0, 0.8)); plot.new(); plot.window(c(0, 1), c(0, 1)); panel("a"); title_left("Fixed five-genus score construction")
  y <- 0.52; segments(0.10, y, 0.91, y, col = gray, lwd = 2); x <- c(0.15, 0.39, 0.63, 0.87); fill <- c(pale, mint, "#EEF0F8", rose); labs <- c("Genus\ncounts", "Sample-wise\nCLR", "Equal-weight\nmean", "Standardised\nscore"); subs <- c("pseudocount 0.5", "five components", "weight 0.2 each", "per 1 SD")
  for (i in 1:4) { points(x[i], y, pch = 21, bg = fill[i], col = if(i == 4) red else gray, cex = 3.7); text(x[i], 0.74, labs[i], font = 2, cex = 0.66); text(x[i], 0.29, subs[i], cex = 0.57, col = gray) }
  par(mar = c(4.0, 7.0, 2.0, 0.6)); ord <- rev(seq_len(nrow(comp))); yy <- seq_along(ord); plot(NA, xlim = c(-6, 105), ylim = c(0.5, nrow(comp) + 0.5), yaxt = "n", xaxt = "n", xlab = "Detection frequency (%)", ylab = "", bty = "n"); axis(1, at = seq(0,100,25)); axis(2, at = yy, labels = comp$genus[ord], las = 1, tick = FALSE, cex.axis = 0.65)
  for (i in seq_along(ord)) { k <- ord[i]; segments(comp$detected_nonards_pct[k], yy[i]-0.16, comp$detected_ards_pct[k], yy[i]+0.16, col = lightgray, lwd = 4); points(comp$detected_nonards_pct[k], yy[i]-0.16, pch=21,bg=blue,col=white,cex=1.05); points(comp$detected_ards_pct[k], yy[i]+0.16, pch=21,bg=red,col=white,cex=1.05) }; legend("bottomright", legend = c("non-ARDS", "ARDS"), pch = 21, pt.bg = c(blue,red), pt.cex=0.9, bty="n", cex=0.59); panel("b"); title_left("Genus detection differs most for Klebsiella", 0.76)
  par(mar = c(4.0, 7.0, 2.0, 0.6)); yy <- seq_len(nrow(comp)); ord2 <- order(comp$clr_effect_SD); plot(NA, xlim = c(-1.05, 1.45), ylim = c(0.5,nrow(comp)+0.5), yaxt="n", xlab="ARDS - non-ARDS CLR effect (SD)", ylab="", bty="n"); axis(2, at=yy, labels=comp$genus[ord2], las=1,tick=FALSE,cex.axis=0.65); abline(v=0,lty=2,col=gray)
  for(i in seq_along(ord2)) { k <- ord2[i]; cc <- if(comp$genus[k] == "Klebsiella") red else if(comp$genus[k] == "Borreliella") white else teal; bc <- if(comp$genus[k] == "Borreliella") orange else white; draw_ci(comp$clr_effect_SD[k], comp$clr_effect_SD_CI_low[k], comp$clr_effect_SD_CI_high[k], yy[i], pchcol=cc); points(comp$clr_effect_SD[k],yy[i],pch=21,bg=cc,col=bc,cex=1.25) }; panel("c"); title_left("Component-level CLR contrasts",0.82)
  par(mar = c(4.0, 7.0, 2.0, 0.6)); ord3 <- order(comp$fraction_of_total_raw_score_difference); yy <- seq_along(ord3); xx <- comp$fraction_of_total_raw_score_difference[ord3]; plot(NA,xlim=c(-0.10,0.68),ylim=c(0.5,nrow(comp)+0.5),yaxt="n",xlab="Fraction of raw-score group difference",ylab="",bty="n"); axis(2,at=yy,labels=comp$genus[ord3],las=1,tick=FALSE,cex.axis=0.65); abline(v=0,lty=2,col=gray)
  for(i in seq_along(xx)) { k <- ord3[i]; barcol <- if(comp$genus[k] == "Klebsiella") red else if(comp$genus[k] == "Borreliella") gray else teal; rect(min(0,xx[i]),yy[i]-0.24,max(0,xx[i]),yy[i]+0.24,col=barcol,border=NA); text(xx[i],yy[i],sprintf("%.2f",xx[i]),pos=4,cex=0.57,col=ink) }; panel("d"); title_left("Weighted contribution to score separation",0.78)
}
export_fig("Figure_2_score_construction", draw2, 7.205, 5.35)

# FIGURE 3 | Primary ARDS association and supporting adjusted estimates.
pt <- read_tsv("patient_signature_components.tsv"); firth <- read_tsv("firth_logistic_models.tsv"); primary <- read_tsv("primary_effects.tsv"); write_tsv(pt, "figure3_patient_scores.tsv"); write_tsv(firth, "figure3_firth_models.tsv"); write_tsv(primary,"figure3_primary_effect.tsv")
ecdf_xy <- function(x) { z <- sort(x); data.frame(x = c(min(z)-0.02, z), y = c(0, seq_along(z)/length(z))) }
draw3 <- function() {
  layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1,1),heights=c(1,1))
  par(mar=c(4.1,4.3,2.0,0.7)); xj <- ifelse(pt$ARDS==1,2,1) + 0.10*sin(seq_len(nrow(pt))*2.399963); yr3a <- range(pt$signature_score_sd, na.rm = TRUE); pad3a <- max(0.22, 0.08 * diff(yr3a)); plot(xj,pt$signature_score_sd,pch=21,bg=ifelse(pt$ARDS==1,red,blue),col=white,cex=1.02,xaxt="n",xlab="",ylab="Fixed five-component score (SD)",bty="n",xlim=c(0.62,2.38),ylim=yr3a+c(-pad3a,pad3a)); axis(1,at=1:2,labels=c("non-ARDS\nn=39","ARDS\nn=26"),tick=FALSE)
  for(g in 0:1){ q <- quantile(pt$signature_score_sd[pt$ARDS==g],c(.25,.5,.75)); rect(g+0.82,q[1],g+1.18,q[3],col=adjustcolor(if(g==1)red else blue,.17),border=if(g==1)red else blue); segments(g+0.82,q[2],g+1.18,q[2],lwd=1.5,col=ink) }; panel("a"); title_left("Patient-level score distribution",0.84)
  par(mar=c(4.1,4.1,2.0,0.7)); e0 <- ecdf_xy(pt$signature_score_sd[pt$ARDS==0]); e1 <- ecdf_xy(pt$signature_score_sd[pt$ARDS==1]); plot(e0,type="s",lwd=2,col=blue,xlim=range(pt$signature_score_sd)+c(-.12,.12),ylim=c(0,1),xlab="Fixed score (SD)",ylab="Cumulative proportion",bty="n"); lines(e1,type="s",lwd=2,col=red); abline(v=0,lty=2,col=gray); legend("topleft",legend=c("non-ARDS (n=39)","ARDS (n=26)"),col=c(blue,red),lwd=2,bty="n",cex=.61); panel("b"); title_left("Distributional shift without dichotomisation",0.78)
  par(mar=c(4.1,8.4,2.0,0.7)); plot(NA,xlim=c(-.15,1.25),ylim=c(.5,2.5),yaxt="n",xlab="Effect estimate",ylab="",bty="n"); abline(v=0,lty=2,col=gray); draw_ci(primary$estimate,primary$CI_low,primary$CI_high,1.75,pchcol=red,lwd=2.2); points(primary$AUROC_descriptive,0.95,pch=21,bg=blue,col=white,cex=1.25); segments(primary$AUROC_CI_low,.95,primary$AUROC_CI_high,.95,lwd=1.8,col=ink); axis(2,at=c(.95,1.75),labels=c("Descriptive AUROC","Primary mean difference"),las=1,tick=FALSE,cex.axis=.55); text(primary$estimate,2.13,"0.59 SD | 95% CI 0.11-1.07 | P=0.009299",cex=.59,col=red); text(primary$AUROC_descriptive,.62,"0.676 | 95% CI 0.536-0.804",cex=.57,col=blue); panel("c"); title_left("Primary effect and AUROC",0.76)
  par(mar=c(4.1,6.8,2.0,0.7)); yy <- rev(seq_len(nrow(firth))); labs <- c("Unadjusted","Age + sex + BMI","+ log-depth"); plot(firth$OR_firth,yy,log="x",xlim=c(.8,4.2),ylim=c(.5,nrow(firth)+.5),pch=21,bg=teal,col=white,yaxt="n",xlab="Odds ratio per score SD",ylab="",bty="n",axes=FALSE); axis(1,at=c(1,1.5,2,3,4),labels=c("1","1.5","2","3","4")); axis(2,at=yy,labels=labs,las=1,tick=FALSE,cex.axis=.62,line=-.35); abline(v=1,lty=2,col=gray); for(i in seq_len(nrow(firth))) draw_ci(firth$OR_firth[i],firth$profile_CI_low[i],firth$profile_CI_high[i],yy[i],pchcol=teal,lwd=1.8); panel("d"); title_left("Supportive Firth models",0.84)
}
export_fig("Figure_3_primary_ards_association", draw3, 7.205, 5.35)

# FIGURE 4 | Stability checks, kept separate by inferential status.
sens <- read_tsv("sensitivity_with_ci.tsv"); alt <- read_tsv("alternative_compositional_sensitivity.tsv"); lopo <- read_tsv("leave_one_patient_out.tsv"); write_tsv(sens,"figure4_sensitivity_with_ci.tsv"); write_tsv(alt,"figure4_posthoc_compositional.tsv"); write_tsv(lopo,"figure4_leave_one_patient_out.tsv")
forest_effect <- function(d, labs, main, xlim=c(-.35,1.35), colors=NULL, note=NULL) { yy <- rev(seq_len(nrow(d))); plot(NA,xlim=xlim,ylim=c(.5,nrow(d)+.5),yaxt="n",xlab="ARDS - non-ARDS (SD)",ylab="",bty="n"); abline(v=0,lty=2,col=gray); axis(2,at=yy,labels=labs,las=1,tick=FALSE,cex.axis=.58,line=-.35); for(i in seq_len(nrow(d))) draw_ci(d$effect_SD[i],d$bootstrap_CI_low[i],d$bootstrap_CI_high[i],yy[i],pchcol=if(is.null(colors))teal else colors[i]); title_left(main,.76) }
draw4 <- function() {
  layout(matrix(1:4,2,2,byrow=TRUE),widths=c(1.03,.97),heights=c(1,1))
  par(mar=c(4.0,7.7,2.0,.6)); ss <- sens[match(c("primary_five_component_fixed","log_depth_residualized","pseudocount_1.0","leave_out_Borreliella_fixed_weight"),sens$analysis),]; forest_effect(ss,c("Permutation-seed sensitivity","Log-depth residualised","Pseudocount 1.0","Leave Borreliella out"),"Sensitivity analyses on primary scale",note="The primary P=0.009299 is shown in Fig. 3."); panel("a")
  par(mar=c(4.0,6.8,2.0,.8)); keep <- grep("leave_out_(Escherichia|Staphylococcus|Enterococcus|Klebsiella)",sens$analysis); ss2 <- sens[keep,]; labs<-c("Leave Escherichia out","Leave Staphylococcus out","Leave Enterococcus out","Leave Klebsiella out"); forest_effect(ss2,labs,"Fixed-weight component removal",xlim=c(-.25,1.18),colors=ifelse(grepl("Klebsiella",ss2$analysis),red,teal),note="Klebsiella removal attenuates the association: 0.32 SD (95% CI -0.17 to 0.82). "); panel("b")
  par(mar=c(4.0,7.7,2.0,.6)); labs3<-c("4 genera: PC 0.5","4 genera: PC 1.0","4 genera vs background"); forest_effect(alt,labs3,"Post hoc compositional analyses",xlim=c(-.20,1.35),colors=rep(orange,nrow(alt)),note="Exploratory analyses; they do not replace the primary score."); panel("c")
  par(mar=c(4.0,4.4,2.0,.6)); hist(lopo$effect_mean_difference,breaks=12,col=pale,border=white,xlab="Leave-one-patient-out effect (SD)",ylab="Number of omissions",main="",xlim=c(.35,.82)); abline(v=.594743654,lwd=2,col=red); rug(lopo$effect_mean_difference,col=teal,lwd=1.2); legend("topright",legend=c("Full cohort: 0.59 SD","Each omitted patient"),col=c(red,teal),lwd=c(2,1.2),bty="n",cex=.58); panel("d");title_left("Patient-omission influence check",.78)
}
export_fig("Figure_4_robustness", draw4,7.205,5.35)

# FIGURE 5 | Paired blood-neutrophil relative-score/CLR context, not absolute load.
bp <- read_tsv("paired_scores.tsv"); pc <- read_tsv("paired_correlation_ci.tsv"); pg <- read_tsv("paired_genus_effects.tsv"); pr <- read_tsv("paired_results.tsv"); gc <- read_tsv("genus_pair_concordance.tsv"); write_tsv(bp,"figure5_paired_scores.tsv");write_tsv(pc,"figure5_correlation_ci.tsv");write_tsv(pg,"figure5_genus_effects.tsv");write_tsv(gc,"figure5_genus_detection_concordance.tsv")
draw5 <- function() {
  layout(matrix(1:6,2,3,byrow=TRUE),widths=c(1.05,.93,1.05),heights=c(1,1))
  par(mar=c(3.8,3.9,2.0,.5)); yr<-range(c(bp$blood_score,bp$neutrophil_score)); plot(NA,xlim=c(.75,2.25),ylim=yr+c(-.2,.2),xaxt="n",xlab="",ylab="Relative score (SD)",bty="n"); for(i in seq_len(nrow(bp))) segments(1,bp$blood_score[i],2,bp$neutrophil_score[i],col=adjustcolor(gray,.28),lwd=.75); points(rep(1,nrow(bp)),bp$blood_score,pch=21,bg=blue,col=white,cex=.68);points(rep(2,nrow(bp)),bp$neutrophil_score,pch=21,bg=teal,col=white,cex=.68); axis(1,1:2,c("Blood","Neutrophil"),tick=FALSE); panel("a");title_left("Paired participant shifts",.77)
  par(mar=c(3.8,3.9,2.0,.5)); dd<-bp$neutrophil_score-bp$blood_score; hist(dd,breaks=11,col=mint,border=white,xlab="Neutrophil - blood (SD)",ylab="Pairs",main="");abline(v=0,lty=2,col=gray);abline(v=pr$neutrophil_minus_blood_mean[1],col=red,lwd=2);panel("b");title_left("Shift distribution",.82)
  par(mar=c(3.8,4.0,2.0,.5)); y5<-range(bp$neutrophil_score); pad5<-0.10*diff(y5); plot(bp$blood_score,bp$neutrophil_score,pch=21,bg=teal,col=white,cex=.75,xlab="Blood score (SD)",ylab="Neutrophil score (SD)",ylim=y5+c(-pad5,pad5),bty="n");abline(lm(neutrophil_score~blood_score,data=bp),col=ink,lwd=1.3);abline(h=0,v=0,lty=3,col=lightgray);panel("c");title_left("Weak participant-level concordance",.70)
  par(mar=c(3.8,4.8,2.0,.5)); plot(NA,xlim=c(-.5,.6),ylim=c(.5,1.5),yaxt="n",xlab="Spearman rank correlation",ylab="",bty="n");abline(v=0,lty=2,col=gray);draw_ci(pc$spearman_rho,pc$patient_bootstrap_CI_low,pc$patient_bootstrap_CI_high,1,pchcol=blue,lwd=2);text(pc$spearman_rho,1.23,"rho = 0.10",cex=.63,col=blue);panel("d");title_left("Correlation uncertainty",.82)
  par(mar=c(3.8,6.0,2.0,.5)); ord<-order(pg$mean_neutrophil_minus_blood_CLR);yy<-seq_along(ord);plot(NA,xlim=c(-.9,2.3),ylim=c(.5,nrow(pg)+.5),yaxt="n",xlab="Neutrophil - blood CLR",ylab="",bty="n");axis(2,at=yy,labels=pg$genus[ord],las=1,tick=FALSE,cex.axis=.60,line=-.35);abline(v=0,lty=2,col=gray);for(i in seq_along(ord)){k<-ord[i];cc<-if(pg$genus[k]=="Borreliella")white else orange;bc<-if(pg$genus[k]=="Borreliella")red else white;draw_ci(pg$mean_neutrophil_minus_blood_CLR[k],pg$paired_bootstrap_CI_low[k],pg$paired_bootstrap_CI_high[k],yy[i],pchcol=cc);points(pg$mean_neutrophil_minus_blood_CLR[k],yy[i],pch=21,bg=cc,col=bc,cex=1.25)};panel("e");title_left("Genus-level relative shifts",.78)
  par(mar=c(4.5,4.8,2.0,.5)); mm<-rbind(gc$blood_detected/gc$n_pairs,gc$neutrophil_detected/gc$n_pairs,gc$both_detected/gc$n_pairs);image(1:nrow(gc),1:3,t(mm),col=colorRampPalette(c("#F3F5F6",teal))(100),axes=FALSE,xlab="",ylab="",zlim=c(0,1));axis(1,1:nrow(gc),gc$genus,tick=FALSE,cex.axis=.53,las=2);axis(2,1:3,c("Blood","Neutrophil","Both"),tick=FALSE,las=1,cex.axis=.57);for(i in 1:nrow(gc))for(j in 1:3)text(i,j,sprintf("%d%%",round(mm[j,i]*100)),cex=.52,col=if(mm[j,i]>.55)white else ink);abline(v=seq(.5,nrow(gc)+.5,1),h=seq(.5,3.5,1),col=white);panel("f");title_left("Detection pattern across pairs",.76)
}
export_fig("Figure_5_blood_neutrophil",draw5,7.205,5.35)

# FIGURE 6 | Independent host-severity background, explicitly not a microbe-host or ARDS test.
hr <- read_tsv("host_module_results.tsv"); hd <- read_tsv("host_module_dictionary.tsv"); ho <- read_tsv("host_module_overlap.tsv"); hs <- read_tsv("host_module_scores.tsv"); write_tsv(hr,"figure6_host_results.tsv");write_tsv(hd,"figure6_module_dictionary.tsv");write_tsv(ho,"figure6_module_overlap.tsv");write_tsv(hs,"figure6_host_scores.tsv")
draw6 <- function() {
  layout(matrix(c(1,1,2,3,4,4),2,3,byrow=TRUE),widths=c(1.25,.85,.90),heights=c(1.05,.95)); mods<-hr$module; short<-c("Degranulation","Innate NF-kB","NETosis / oxidative","Endothelial barrier")
  par(mar=c(4.1,4.1,2.0,.7)); plot(NA,xlim=c(.55,3.45),ylim=range(as.matrix(hs[,mods]),na.rm=TRUE)+c(-.25,.25),xaxt="n",xlab="AP severity",ylab="Module score (z)",bty="n"); cols<-c(blue,teal,orange,red); lev<-c("Mild AP","Moderately-severe AP","Severe AP"); for(k in seq_along(mods)){ xx<-match(hs$pathology,lev)+(-.18+(k-1)*.12); points(xx,hs[[mods[k]]],pch=21,bg=adjustcolor(cols[k],.62),col=white,cex=.55); med<-tapply(hs[[mods[k]]],match(hs$pathology,lev),median,na.rm=TRUE);lines((1:3)+(-.18+(k-1)*.12),med,col=cols[k],lwd=1.4,type="b",pch=16,cex=.55) }; axis(1,1:3,c("Mild","Moderate","Severe"),tick=FALSE);legend("topleft",legend=short,col=cols,lwd=1.5,pch=16,cex=.54,bty="n",ncol=2);panel("a");title_left("Independent AP-severity module distributions (n=87)",.76)
  par(mar=c(4.1,6.4,2.0,.5)); yy<-rev(seq_len(nrow(hr)));plot(NA,xlim=c(0,.72),ylim=c(.5,nrow(hr)+.5),yaxt="n",xlab="Spearman rho vs severity",ylab="",bty="n");axis(2,at=yy,labels=short,las=1,tick=FALSE,cex.axis=.55,line=-.35);abline(v=0,lty=2,col=gray);for(i in seq_len(nrow(hr))){points(hr$spearman_rho[i],yy[i],pch=21,bg=cols[i],col=white,cex=1.2);text(hr$spearman_rho[i],yy[i],paste0("  q=",formatC(hr$BH_FDR[i],format="g",digits=2)),pos=4,cex=.53)};panel("b");title_left("Severity trends",.80)
  par(mar=c(4.1,6.2,2.0,.5)); total<-table(hd$module);found<-tapply(hd$found_in_filtered_matrix=="true",hd$module,sum);ord<-mods;yy<-rev(seq_along(ord));plot(NA,xlim=c(0,max(total)+3),ylim=c(.5,length(ord)+.5),yaxt="n",xlab="Genes in module dictionary",ylab="",bty="n");axis(2,at=yy,labels=short,las=1,tick=FALSE,cex.axis=.55,line=-.35);for(i in seq_along(ord)){rect(0,yy[i]-.22,total[ord[i]],yy[i]+.22,col=lightgray,border=NA);rect(0,yy[i]-.22,found[ord[i]],yy[i]+.22,col=cols[i],border=NA);text(total[ord[i]]+.4,yy[i],paste0(found[ord[i]]," / ",total[ord[i]]),pos=4,cex=.57)};panel("c");title_left("Predefined module dictionary coverage",.68)
  par(mar=c(4.6,5.6,2.0,1.5)); ov<-diag(1,4);rownames(ov)<-colnames(ov)<-short;for(i in seq_len(nrow(ho))){a<-match(ho$module_a[i],mods);b<-match(ho$module_b[i],mods);ov[a,b]<-ov[b,a]<-ho$jaccard[i]};image(1:4,1:4,t(ov[4:1,]),col=colorRampPalette(c("#F7F8F9",blue))(100),axes=FALSE,zlim=c(0,1),xlab="",ylab="");axis(1,1:4,short,las=2,tick=FALSE,cex.axis=.54);axis(2,1:4,rev(short),las=1,tick=FALSE,cex.axis=.54);for(i in 1:4)for(j in 1:4){val<-ov[5-j,i];text(i,j,if(i==5-j)"1.00" else sprintf("%.2f",val),cex=.56,col=if(val>.45)white else ink)};abline(v=seq(.5,4.5,1),h=seq(.5,4.5,1),col=white,lwd=1.1);panel("d");title_left("Module overlap is explicit",.80)
}
export_fig("Figure_6_host_response",draw6,7.205,5.35)

# Source-data manifest and submission-facing QA boundary note.
files <- list.files(out, recursive=TRUE, full.names=TRUE); manifest <- data.frame(path=basename(files), bytes=file.info(files)$size, md5=unname(tools::md5sum(files)),stringsAsFactors=FALSE); write_tsv(manifest,"figure_manifest.tsv")
writeLines(c("Figure contracts: F1 evidence hierarchy and replication eligibility; F2 score construction and component architecture; F3 primary ARDS association; F4 sensitivity, component-removal, exploratory, and patient-influence checks; F5 paired relative-score/CLR evidence only; F6 independent AP-severity context only.","All figures were drawn, previewed and exported in R. No image panels, cropping, brightness, contrast, or pseudo-colour alterations are used.","Primary inferential anchor: 0.59 SD, 95% CI 0.11 to 1.07, P=0.009299. The additional random-seed sensitivity result is not substituted for the primary result. Klebsiella-excluded estimate is 0.32 SD (95% CI -0.17 to 0.82; P=0.0984).","Borreliella is all-zero in the primary cohort and absent in both paired compartments; CLR effects involving it are denominator-induced. Four-genus and Aitchison results are post hoc. Paired results do not establish absolute load, viability, localisation, transport, or causal direction.","Final visual QA: Figure 1a uses a card hierarchy that distinguishes the primary outcome test from separate supportive datasets; Figure 1c label sits in margin; Figure 2b includes left-axis padding so zero-frequency markers are fully visible; Figure 2d negative label is in-bounds; Figure 3a includes vertical padding for extreme observations; Figure 4b has a wider plotting region; Figure 5c has a 10% vertical data-range margin."), file.path(out,"FIGURE_QA_NOTE.txt"))
cat("REVISION_FIGURES_COMPLETE\n")
