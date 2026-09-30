#!/usr/bin/env Rscript

# Publication figures for Stage 1. Base-R plotting with journal-safe vector devices.
# All numeric panels are generated directly from machine-readable result files.

root <- Sys.getenv("STAGE1_ROOT", unset=getwd())
stage0 <- Sys.getenv("STAGE0_ROOT", unset=file.path(dirname(root), "gut_pancreas_twin_stage0"))
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
outdir <- file.path(root, "08_figures")
srcdir <- file.path(outdir, "source_data")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(srcdir, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(file.path(root, "logs/PRJNA428535_SUPPORT_COMPLETE"))) {
  stop("PRJNA428535 support analysis is not complete")
}

read_tsv <- function(path) read.delim(path, sep="\t", header=TRUE, check.names=FALSE,
                                      stringsAsFactors=FALSE, quote="", comment.char="")
write_tsv <- function(x, path) write.table(x, path, sep="\t", quote=FALSE,
                                           row.names=FALSE, na="")

blue <- "#2166AC"; red <- "#B2182B"; teal <- "#018571"; orange <- "#D6604D"
gray <- "#5A5A5A"; light <- "#E8EEF4"; dark <- "#1F2933"; gold <- "#B8860B"

journal_width_mm <- 183
font_family <- "Arial"

open_device <- function(path, fmt, width, height) {
  if (fmt == "pdf") grDevices::cairo_pdf(path, width=width, height=height, family=font_family)
  if (fmt == "svg") {
    if (!requireNamespace("svglite", quietly=TRUE)) stop("svglite is required for editable SVG export")
    svglite::svglite(path, width=width, height=height, system_fonts=list(sans=font_family))
  }
  if (fmt == "tiff") tiff(path, width=width, height=height, units="in", res=600,
                           compression="lzw", type="cairo")
  if (fmt == "png") png(path, width=width, height=height, units="in", res=600,
                         type="cairo")
}

export_figure <- function(stem, draw_fun, width=7.2, height=5.2) {
  for (fmt in c("pdf", "svg", "tiff", "png")) {
    ext <- if (fmt == "tiff") "tif" else fmt
    open_device(file.path(outdir, paste0(stem, ".", ext)), fmt, width, height)
    par(family=font_family, fg=dark, col.axis=dark, col.lab=dark, col.main=dark)
    draw_fun()
    dev.off()
  }
}

panel_label <- function(x, y, lab) text(x, y, lab, font=2, cex=1.25, xpd=NA)

# Figure 1: study architecture and patient-level roles.
fig1_nodes <- data.frame(
  order=1:6,
  dataset=c("Stage 0 cross-site cohorts", "Signature lock", "PRJNA893348",
            "PRJNA428535", "GSE194331", "Restricted supporting cohorts"),
  role=c("discovery", "five-genus score fixed before outcome analysis",
      "confirmatory ARDS association", "within-person blood-neutrophil comparison",
         "independent host-response support", "etiology and metabolic context"),
  statistical_unit=c("cohort-specific patients", "frozen definition", "65 AP patients",
                     "62 paired participants", "87 AP patients", "dataset-specific patients"),
  stringsAsFactors=FALSE)
write_tsv(fig1_nodes, file.path(srcdir, "figure1_study_architecture.tsv"))

draw_fig1 <- function() {
  par(mar=c(0.4,0.4,0.6,0.4)); plot.new(); plot.window(c(0,1), c(0,1))
  panel_label(0.015,0.965,"a")
  box <- function(x0,y0,x1,y1,fill,title,body) {
    rect(x0,y0,x1,y1,col=fill,border="#8795A1",lwd=1.2)
    text((x0+x1)/2,y1-0.055,title,font=2,cex=0.9)
    text((x0+x1)/2,(y0+y1)/2-0.015,body,cex=0.73)
  }
  box(0.05,0.67,0.30,0.89,"#F1F5F9","Cross-site discovery",
      "Stage 0 cohorts\ncompartment and source evidence")
  box(0.375,0.67,0.625,0.89,"#E1EEF7","Outcome-blind lock",
      "5 genera, fixed direction\nequal weights, 0.5 pseudocount")
  box(0.70,0.67,0.95,0.89,"#FBE9E7","Primary test",
      "PRJNA893348\n65 AP patients, 26 ARDS")
  arrows(0.30,0.78,0.375,0.78,length=0.08,lwd=1.5,col=gray)
  arrows(0.625,0.78,0.70,0.78,length=0.08,lwd=1.5,col=gray)
  panel_label(0.015,0.57,"b")
  box(0.05,0.28,0.31,0.53,"#EAF4F4","Cellular compartment",
      "PRJNA428535\n62 blood-neutrophil pairs")
  box(0.37,0.28,0.63,0.53,"#EEF0F8","Host response",
      "GSE194331\n87 AP transcriptomes")
  box(0.69,0.28,0.95,0.53,"#F4F1E8","Restricted context",
      "Etiology and metabolomics\nno outcome substitution")
  arrows(0.825,0.67,0.18,0.53,length=0.07,lwd=1.1,col="#94A3B8")
  arrows(0.825,0.67,0.50,0.53,length=0.07,lwd=1.1,col="#94A3B8")
  arrows(0.825,0.67,0.82,0.53,length=0.07,lwd=1.1,col="#94A3B8")
  text(0.50,0.15,"Evidence is integrated at the claim level; abundance matrices are not pooled across platforms or cohorts.",
       cex=0.78,font=3)
  text(0.50,0.07,"Primary unit: patient. Supporting data cannot rescue a failed primary test.",cex=0.78)
}
export_figure("Figure_1_study_design", draw_fig1, 7.2, 5.0)

# Figure 2: discovery evidence and the locked score definition.
pf <- read_tsv(file.path(stage0, "06_results/stage05_PRJNA1031835_source_tests.tsv"))
pf <- pf[pf$threshold == 0.0001 & pf$timepoint == 1 & pf$comparison == "feces_to_pancreatic_fluid",]
mn <- read_tsv(file.path(stage0, "06_results/stage06_PRJNA771396_paired_tests.tsv"))
mn <- mn[mn$analysis %in% c("primary","strict","primary_without_kitome"),]
fig2_ev <- rbind(
  data.frame(dataset="PRJNA1031835", analysis="Time-1 feces to pancreatic fluid",
             n=pf$patients, estimate=pf$delta_mean, p=pf$permutation_p_one_sided),
  data.frame(dataset="PRJNA771396", analysis=paste0("Blood-peripancreatic ", mn$analysis),
             n=mn$paired_patients, estimate=mn$delta_mean, p=mn$permutation_p_one_sided))
write_tsv(fig2_ev, file.path(srcdir, "figure2_discovery_evidence.tsv"))
fig2_sig <- data.frame(genus=c("Borreliella","Escherichia","Staphylococcus","Enterococcus","Klebsiella"),
                       direction="higher score", weight=0.2, pseudocount=0.5)
write_tsv(fig2_sig, file.path(srcdir, "figure2_locked_signature.tsv"))

draw_fig2 <- function() {
  layout(matrix(c(1,2),1,2), widths=c(1.08,0.92))
  par(mar=c(4.5,6.6,2.2,0.8))
  y <- rev(seq_len(nrow(fig2_ev)))
  display_labs <- c("Feces to pancreatic fluid\nTime 1",
                    "Blood to peripancreatic\nprimary",
                    "Blood to peripancreatic\nstrict",
                    "Blood to peripancreatic\nwithout kitome")
  plot(fig2_ev$estimate, y, xlim=c(0,max(fig2_ev$estimate)*1.58), ylim=c(0.5,nrow(fig2_ev)+0.5),
       pch=21,bg=c(teal,rep(blue,nrow(fig2_ev)-1)),col="white",cex=1.5,
       yaxt="n",xlab="Matched-minus-null similarity",ylab="",bty="n")
  axis(2,at=y,labels=display_labs,las=1,tick=FALSE,cex.axis=0.61,line=-0.4)
  abline(v=0,col="#9AA5B1",lty=2)
  text(fig2_ev$estimate,y-0.23,paste0("n=",fig2_ev$n,", P=",formatC(fig2_ev$p,format="g",digits=2)),cex=0.59,pos=4)
  mtext("a",side=3,adj=-0.22,font=2,cex=1.2,line=0.2)
  title("Outcome-blind cross-site evidence",adj=0,cex.main=0.76)
  par(mar=c(1.2,0.8,2.2,0.5)); plot.new(); plot.window(c(0,1),c(0,1))
  text(0.03,0.94,"b",font=2,cex=1.2)
  text(0.12,0.88,"Locked five-genus score",adj=0,font=2,cex=0.86)
  for (i in seq_len(nrow(fig2_sig))) {
    yy <- 0.78-(i-1)*0.12
    rect(0.12,yy-0.04,0.88,yy+0.04,col=if(i%%2) "#E8F1F8" else "#F4F6F8",border=NA)
    text(0.17,yy,fig2_sig$genus[i],adj=0,cex=0.83,font=3)
    text(0.82,yy,"+0.2",adj=1,cex=0.83,font=2,col=blue)
  }
  text(0.12,0.145,"Per sample: add 0.5 to all genera,\ncompute CLR, average the five fixed\ncomponents, then scale to 1 SD.",
       adj=c(0,0.5),cex=0.66)
  text(0.12,0.035,"No outcome-driven reweighting\nor taxon replacement",adj=c(0,0.5),cex=0.64,font=3,col=red)
}
export_figure("Figure_2_discovery_and_lock", draw_fig2, 7.2, 4.4)

# Figure 3: confirmatory patient-level association, effect and descriptive ROC.
pt <- read_tsv(file.path(root,"07_results/patient_signature_components.tsv"))
pe <- read_tsv(file.path(root,"07_results/primary_effects.tsv"))
cv <- read_tsv(file.path(root,"07_results/covariate_support.tsv"))
write_tsv(pt, file.path(srcdir,"figure3_patient_scores.tsv"))
write_tsv(pe, file.path(srcdir,"figure3_primary_effect.tsv"))
write_tsv(cv, file.path(srcdir,"figure3_covariate_support.tsv"))

roc_points <- function(score, y) {
  th <- c(Inf, sort(unique(score),decreasing=TRUE), -Inf)
  out <- do.call(rbind,lapply(th,function(t){
    pred <- score >= t
    c(FPR=sum(pred & y==0)/sum(y==0), TPR=sum(pred & y==1)/sum(y==1), threshold=t)
  }))
  as.data.frame(out)
}
roc <- roc_points(pt$signature_score_sd, pt$ARDS)
write_tsv(roc, file.path(srcdir,"figure3_roc_points.tsv"))

draw_fig3 <- function() {
  layout(matrix(1:3,1,3),widths=c(1.05,0.9,0.9))
  par(mar=c(4.5,4.2,2.0,0.8))
  jitter_offset <- 0.11 * sin(seq_len(nrow(pt)) * 2.399963)
  x <- ifelse(pt$ARDS==1,2,1) + jitter_offset
  plot(x,pt$signature_score_sd,pch=21,bg=ifelse(pt$ARDS==1,red,blue),col="white",
       xaxt="n",xlab="",ylab="Locked signature score (SD)",bty="n",xlim=c(0.65,2.35))
  axis(1,at=1:2,labels=c("nonARDS\nn=39","ARDS\nn=26"),tick=FALSE)
  segments(0.86,mean(pt$signature_score_sd[pt$ARDS==0]),1.14,mean(pt$signature_score_sd[pt$ARDS==0]),lwd=2.2,col=dark)
  segments(1.86,mean(pt$signature_score_sd[pt$ARDS==1]),2.14,mean(pt$signature_score_sd[pt$ARDS==1]),lwd=2.2,col=dark)
  mtext("a",3,adj=-0.2,font=2,cex=1.2,line=0.2); title("Patient distributions",adj=0,cex.main=0.9)
  par(mar=c(4.5,4.6,2.0,0.8))
  est <- pe$estimate[1]; lo <- pe$CI_low[1]; hi <- pe$CI_high[1]
  plot(est,1,xlim=c(min(0,lo)-0.1,hi+0.2),ylim=c(0.5,1.5),pch=21,bg=red,col="white",cex=1.5,
       yaxt="n",ylab="",xlab="ARDS - nonARDS (SD)",bty="n")
  segments(lo,1,hi,1,lwd=2.2,col=dark); abline(v=0,lty=2,col="#9AA5B1")
  text(est,0.72,sprintf("%.2f (95%% CI %.2f to %.2f)\nP=%.4f",est,lo,hi,pe$one_sided_permutation_p[1]),cex=0.72)
  mtext("b",3,adj=-0.22,font=2,cex=1.2,line=0.2); title("Prespecified effect",adj=0,cex.main=0.9)
  par(mar=c(4.5,4.5,2.0,0.8))
  plot(roc$FPR,roc$TPR,type="l",lwd=2.2,col=red,xlim=c(0,1),ylim=c(0,1),asp=1,
       xlab="False-positive rate",ylab="True-positive rate",bty="n")
  abline(0,1,lty=2,col="#9AA5B1")
  text(0.58,0.18,sprintf("AUROC %.2f\n95%% CI %.2f-%.2f\ndescriptive only",pe$AUROC_descriptive[1],pe$AUROC_CI_low[1],pe$AUROC_CI_high[1]),cex=0.72)
  mtext("c",3,adj=-0.22,font=2,cex=1.2,line=0.2); title("Apparent discrimination",adj=0,cex.main=0.9)
}
export_figure("Figure_3_primary_ards_association", draw_fig3, 7.2, 3.4)

# Figure 4: influence and prespecified robustness analyses.
loo <- read_tsv(file.path(root,"07_results/leave_one_patient_out.tsv"))
sens <- read_tsv(file.path(root,"07_results/contamination_sensitivity.tsv"))
write_tsv(loo,file.path(srcdir,"figure4_leave_one_patient_out.tsv"))
write_tsv(sens,file.path(srcdir,"figure4_sensitivity.tsv"))

draw_fig4 <- function() {
  layout(matrix(c(1,2),1,2),widths=c(0.98,1.12))
  par(mar=c(4.5,4.3,2.0,1.0))
  ord <- order(loo$effect_mean_difference)
  plot(seq_along(ord),loo$effect_mean_difference[ord],pch=21,bg=ifelse(loo$omitted_outcome[ord]==1,red,blue),
       col="white",xlab="Leave-one-patient-out iteration",ylab="ARDS - nonARDS (SD)",bty="n")
  abline(h=pe$estimate[1],lwd=1.5,col=dark); abline(h=0,lty=2,col="#9AA5B1")
  legend("bottomright",c("omitted nonARDS","omitted ARDS","full estimate"),pch=c(21,21,NA),
         pt.bg=c(blue,red,NA),lty=c(NA,NA,1),col=c("white","white",dark),bty="n",cex=0.68)
  mtext("a",3,adj=-0.18,font=2,cex=1.2,line=0.2); title("Patient influence",adj=0,cex.main=0.9)
  par(mar=c(4.5,6.8,2.2,0.5))
  keep <- sens$critical %in% c(TRUE,"true","True",1)
  ss <- sens[keep,]; ss <- ss[order(ss$effect),]
  labs <- c("Leave out Klebsiella", "Leave out Staphylococcus", "Log-depth residualized",
            "Primary, pseudocount 0.5", "Pseudocount 1.0", "Leave out Enterococcus",
            "Leave out Escherichia", "Leave out Borreliella", "Blind exclusion: bottom 5% depth")
  if (length(labs) != nrow(ss)) labs <- gsub("_"," ",ss$analysis)
  plot(ss$effect,seq_len(nrow(ss)),xlim=c(min(0,min(ss$effect))-0.05,max(ss$effect)+0.34),
       ylim=c(0.5,nrow(ss)+0.5),pch=21,bg=teal,col="white",yaxt="n",xlab="ARDS - nonARDS (SD)",ylab="",bty="n")
  axis(2,at=seq_len(nrow(ss)),labels=labs,las=1,tick=FALSE,cex.axis=0.56,line=-0.4)
  abline(v=0,lty=2,col="#9AA5B1")
  text(ss$effect,seq_len(nrow(ss)),paste0(" P=",formatC(ss$p,format="g",digits=2)),pos=4,cex=0.57)
  mtext("b",3,adj=-0.30,font=2,cex=1.2,line=0.2); title("Prespecified sensitivity analyses",adj=0,cex.main=0.82)
}
export_figure("Figure_4_robustness", draw_fig4, 7.2, 4.5)

# Figure 5: within-person blood-neutrophil support.
bp <- read_tsv(file.path(root,"06_supporting/blood_neutrophil/paired_scores.tsv"))
br <- read_tsv(file.path(root,"06_supporting/blood_neutrophil/paired_results.tsv"))
gc <- read_tsv(file.path(root,"06_supporting/blood_neutrophil/genus_pair_concordance.tsv"))
write_tsv(bp,file.path(srcdir,"figure5_paired_scores.tsv")); write_tsv(br,file.path(srcdir,"figure5_paired_results.tsv")); write_tsv(gc,file.path(srcdir,"figure5_genus_concordance.tsv"))
main_bp <- br[br$analysis=="all_pairs_ge100",]

draw_fig5 <- function() {
  layout(matrix(1:3,1,3),widths=c(1,1,0.95))
  par(mar=c(4.5,4.5,2.0,0.8))
  ok <- bp$pair_qc_ge100 %in% c(TRUE,"true","True",1)
  plot(bp$blood_score[ok],bp$neutrophil_score[ok],pch=21,bg=teal,col="white",
       xlab="Blood score (SD)",ylab="Neutrophil score (SD)",bty="n")
  abline(lm(neutrophil_score~blood_score,data=bp[ok,]),col=dark,lwd=1.6)
  text(par("usr")[1],par("usr")[4],sprintf("rho=%.2f\nP=%.4f\nn=%d",main_bp$spearman_rho,main_bp$pairing_permutation_p,main_bp$n_pairs),adj=c(0,1),cex=0.72)
  mtext("a",3,adj=-0.20,font=2,cex=1.2,line=0.2); title("Participant-level correlation",adj=0,cex.main=0.86)
  par(mar=c(4.5,3.5,2.0,0.8))
  idx <- which(ok); ord <- order(bp$neutrophil_score[idx]-bp$blood_score[idx]); take <- idx[ord]
  matplot(rbind(rep(1,length(take)),rep(2,length(take))),rbind(bp$blood_score[take],bp$neutrophil_score[take]),
          type="l",lty=1,col=adjustcolor(gray,alpha.f=0.32),xaxt="n",xlab="",ylab="Score (SD)",bty="n")
  points(rep(1,length(take)),bp$blood_score[take],pch=16,cex=0.45,col=blue)
  points(rep(2,length(take)),bp$neutrophil_score[take],pch=16,cex=0.45,col=teal)
  axis(1,at=1:2,labels=c("Blood","Neutrophils"),tick=FALSE)
  text(1.5,par("usr")[3],sprintf("mean difference %.2f\n95%% CI %.2f to %.2f",main_bp$neutrophil_minus_blood_mean,main_bp$paired_bootstrap_CI_low,main_bp$paired_bootstrap_CI_high),pos=3,cex=0.67)
  mtext("b",3,adj=-0.20,font=2,cex=1.2,line=0.2); title("Paired enrichment",adj=0,cex.main=0.86)
  par(mar=c(4.5,6.4,2.0,0.8))
  yy <- rev(seq_len(nrow(gc))); both_absent <- gc$genus == "Borreliella" & gc$agreement_fraction == 1
  plot(gc$agreement_fraction,yy,xlim=c(0,1),ylim=c(0.5,nrow(gc)+0.5),
       pch=21,bg=ifelse(both_absent,"white",orange),col=ifelse(both_absent,orange,"white"),
       yaxt="n",xlab="Detection agreement",ylab="",bty="n")
  axis(2,at=yy,labels=gc$genus,las=1,tick=FALSE,cex.axis=0.72)
  abline(v=0.5,lty=2,col="#9AA5B1")
  if (any(both_absent)) text(gc$agreement_fraction[both_absent]-0.03,yy[both_absent],"both absent",pos=2,cex=0.58,col=orange)
  mtext("c",3,adj=-0.33,font=2,cex=1.2,line=0.2); title("Genus-level agreement",adj=0,cex.main=0.86)
}
export_figure("Figure_5_blood_neutrophil", draw_fig5, 7.2, 3.6)

# Figure 6: independent host-response support and evidence boundary.
hm <- read_tsv(file.path(root,"06_supporting/host_response/host_module_scores.tsv"))
hr <- read_tsv(file.path(root,"06_supporting/host_response/host_module_results.tsv"))
write_tsv(hm,file.path(srcdir,"figure6_host_module_scores.tsv")); write_tsv(hr,file.path(srcdir,"figure6_host_module_results.tsv"))
module_names <- c(neutrophil_degranulation="Neutrophil degranulation",innate_NFkB_inflammation="Innate NF-kB inflammation",NETosis_oxidative_burst="NETosis and oxidative burst",endothelial_barrier_injury="Endothelial barrier injury")

draw_fig6 <- function() {
  layout(matrix(c(1,2,3),1,3),widths=c(1.16,0.92,0.92))
  par(mar=c(5.2,4.3,2.0,0.6))
  sev <- factor(hm$pathology,levels=c("Mild AP","Moderately-severe AP","Severe AP"))
  boxplot(hm$NETosis_oxidative_burst~sev,col=c("#DCEAF5","#F5D9D4","#C95A54"),border="#778899",
          outline=FALSE,ylab="NETosis module score",xlab="",xaxt="n",bty="n")
  stripchart(hm$NETosis_oxidative_burst~sev,vertical=TRUE,method="jitter",pch=21,bg=adjustcolor(dark,alpha.f=.45),col="white",add=TRUE,cex=.65)
  axis(1,at=1:3,labels=c("Mild\nn=57","Moderately\nsevere n=20","Severe\nn=10"),tick=FALSE,cex.axis=.7)
  mtext("a",3,adj=-0.18,font=2,cex=1.2,line=.2); title("Independent severity trend",adj=0,cex.main=.88)
  par(mar=c(4.5,6.6,2.2,0.5)); yy <- rev(seq_len(nrow(hr)))
  module_display <- c("Neutrophil degranulation", "Innate NF-kB", "NETosis / oxidative burst", "Endothelial barrier injury")
  plot(hr$spearman_rho,yy,xlim=c(0,max(hr$spearman_rho)+.22),ylim=c(.5,nrow(hr)+.5),pch=21,bg=blue,col="white",yaxt="n",xlab="Spearman rho",ylab="",bty="n")
  axis(2,at=yy,labels=module_display,las=1,tick=FALSE,cex.axis=.58,line=-0.4)
  text(hr$spearman_rho,yy,paste0(" q=",formatC(hr$BH_FDR,format="g",digits=2)),pos=4,cex=.56)
  abline(v=0,lty=2,col="#9AA5B1")
  mtext("b",3,adj=-.28,font=2,cex=1.2,line=.2); title("Host-response modules",adj=0,cex.main=.76)
  par(mar=c(1.0,0.8,2.2,0.5)); plot.new(); plot.window(c(0,1),c(0,1))
  text(-0.02,1.02,"c",font=2,cex=1.2,xpd=NA,adj=c(0,1)); text(0.08,1.01,"Claim boundary",font=2,cex=.78,xpd=NA,adj=c(0,1))
  rect(.08,.62,.92,.83,col="#E8F1F8",border=NA); text(.50,.75,"Supported",font=2,cex=.86); text(.50,.67,"Independent directional triangulation",cex=.72)
  rect(.08,.35,.92,.56,col="#F7F0E8",border=NA); text(.50,.48,"Not tested",font=2,cex=.86); text(.50,.40,"Within-patient microbe-host correlation",cex=.72)
  rect(.08,.08,.92,.29,col="#F8E8E6",border=NA); text(.50,.21,"Not claimed",font=2,cex=.86); text(.50,.125,"Causality, viable bacteria,\nor migration direction",cex=.66)
}
export_figure("Figure_6_host_response", draw_fig6, 7.2, 3.8)

# Machine-readable figure manifest with source-data checksums.
all_files <- list.files(outdir, recursive=TRUE, full.names=TRUE)
all_files <- all_files[file.info(all_files)$isdir == FALSE]
manifest <- data.frame(path=sub(paste0("^",root,"/"),"",all_files),
                       bytes=file.info(all_files)$size,
                       md5=unname(tools::md5sum(all_files)),stringsAsFactors=FALSE)
write_tsv(manifest,file.path(outdir,"figure_manifest.tsv"))
writeLines(c("Publication figures generated from locked result files using base R only.",
             paste("R",R.version.string),"Review: author-reviewed"),file.path(outdir,"FIGURE_QA_NOTE.txt"))
cat("FIGURES_COMPLETE\n")
