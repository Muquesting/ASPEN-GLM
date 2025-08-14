
calc_auc <- function(pvals, actual, sequence){

  tpr <- vector()
  fpr <- vector()
  prec <- vector()
  recall <- vector()
  tp <- vector()
  fn <- vector()
  fp <- vector()
  tn <- vector()

  len <- length(sequence)
  tpr  <- fpr <- prec <- recall <- tp <- fn <- fp <- tn <- numeric(len)

  for (i in seq_along(sequence)){
    prediction <- ifelse(pvals <= sequence[i], 1, 0)
    t <- table(prediction, actual)
    tp[i] <- t[4] #true positives, significant by the test and overlapped PCHiC
    fn[i] <- t[3] #false negative, not significant by the test and overlapped PCHiC
    fp[i] <- t[2] #false positive, significant by the test and did not overlap PCHiC
    tn[i] <- t[1] #true negative, not significant by the test and did not overlap PCHiC
  }

  #tp[length(tp)] <- tp[length(tp) - 1]
  #fn[length(fn)] <- fn[length(fn) - 1]
  #fp[length(fp)] <- fp[length(fp) - 1]
  #tn[length(tn)] <- tn[length(tn) - 1]

  for (i in seq_along(sequence)){
  tpr[i] <- tp[i]/(tp[i]+fn[i])
  fpr[i] <- fp[i]/(fp[i]+tn[i])
  prec[i] <- tp[i]/(tp[i]+fp[i])
  recall[i] <- tp[i]/(tp[i]+fn[i])
  }

  #tpr[is.na(tpr)] <- 1
  #fpr[is.na(fpr)] <- 1

  tpr <- c(0, tpr, 1)
  fpr <- c(0, fpr, 1)

  #prec[is.na(prec)] <- 1
  #recall[is.na(recall)] <- 1

  prec <- c(1, prec, 0.5, 0)
  recall <- c(0, recall, 1, 1)

  height <- (tpr[-1]+tpr[-length(tpr)])/2
  width <- diff(fpr) # = diff(rev(omspec))
  auc <- sum(height*width, na.rm = T)

  height2 <- (recall[-1]+recall[-length(recall)])/2
  width2 <- -diff(prec) # = diff(rev(omspec))
  aupr <- sum(height2*width2, na.rm = T)

  list(tp = tp, fn = fn, fp = fp,
       tn = tn, tpr = tpr, fpr = fpr,
       prec = prec, recall = recall, auc = auc, aupr = aupr)

}


plot_fpr <- function(res, sequence, fdr_orig = FALSE, by_ge = FALSE){

  if ( by_ge == TRUE) {
    res_list <- split(res, f = res$ge_group)
      } else {
    res_list <- split(res, f = res$cell_group)
  }

  if (fdr_orig == TRUE){
    fdr_vals <- lapply(res_list, function(q) q$fdr_orig)
      } else {
    fdr_vals <- lapply(res_list, function(q) q$fdr_shrunk)
  }

  actual <- lapply(res_list, function(q) q$true_AI)

  for (i in seq_along(sequence)){
    prediction <- lapply(fdr_vals, function(q) ifelse(q <= sequence[i], 1, 0))
    t <- mapply(function(p, q) table(p, q), prediction, actual, SIMPLIFY = F)

    if (i == 1){

      tp <- lapply(t, function(q) q[2])
      fn <- lapply(t, function(q) q[1])
      tpr <- mapply(function(p, q) p/(p+q), tp, fn, SIMPLIFY = F)

    } else {

      tp <- Map(c, tp, lapply(t, function(q) q[2]))
      fn <- Map(c, fn, lapply(t, function(q) q[1]))
      tpr <- Map(c, tpr, mapply(function(p, q) p/(p+q), tp, fn, SIMPLIFY = F))
    }

    tp[i] <- lapply(t, function(q) q[2]) #true positives, significant by the test and overlapped PCHiC
    fn[i] <- t[1] #false negative, not significant by the test and overlapped PCHiC
    tpr[i] <- tp[i]/(tp[i]+fn[i])
 }



}

plot_prcurve <- function(auc_res_test1, auc_res_test2){

  require(ggsci)

  data.ppv <- data.frame(rbind(cbind(round(auc_res_test1$prec, 4), round(auc_res_test1$recall, 4), Mode = paste0("Original (aupr - ", round(auc_res_test1$aupr, 2), ")")),
                               cbind(round(auc_res_test2$prec, 4), round(auc_res_test2$recall, 4), Mode = paste0("Emp Bayes (aupr - ", round(auc_res_test2$aupr, 2), ")"))))
  colnames(data.ppv) <- c("Precision", "Recall", "Mode")
  data.ppv$Precision<- as.numeric(data.ppv$Precision)
  data.ppv$Recall <- as.numeric(data.ppv$Recall)

  ggplot(data.ppv, aes(x=Recall, y=Precision, colour=Mode)) +
    geom_line(linewidth = 1) +
    theme_classic(base_size = 20) +
    theme(panel.grid.minor = element_blank(), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 25)) +
    theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14), legend.title = element_blank())+
    labs(title = "Precision Recall Curves") +
    ylim(0,1)+
    geom_segment(aes(x = 0, xend = 1, y = 0.5, yend = 0.5), color="darkgrey", linetype="dashed")+
    scale_color_jama() +
    theme(legend.position=c(0.75, 0.25), legend.box="vertical", legend.text = element_text(size = 14))

}


plot_aucurve <- function(auc_res_test1, auc_res_test2){

  data.auc <- data.frame(rbind(cbind(round(auc_res_test1$tpr, 4), round(auc_res_test1$fpr, 4), Mode = paste0("Original (auc - ", round(auc_res_test1$auc, 2), ")")),
                               cbind(round(auc_res_test2$tpr, 4), round(auc_res_test2$fpr, 4), Mode = paste0("Emp Bayes (auc - ", round(auc_res_test2$auc, 2), ")"))))
  colnames(data.auc) <- c("Sensitivity", "1 - Specificity", "Mode")
  data.auc$Sensitivity <- as.numeric(data.auc$Sensitivity)
  data.auc$`1 - Specificity` <- as.numeric(data.auc$`1 - Specificity`)

  ggplot(data.auc, aes(x=`1 - Specificity`, y=Sensitivity, colour=Mode)) +
    geom_line(linewidth =1) +
    theme_classic(base_size = 20) +
    theme(panel.grid.minor = element_blank(), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 25)) +
    theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14), legend.title = element_blank())+
    labs(title = "AUROC curve") +
    ylim(0,1)+
    geom_abline(slope=1, intercept=0, colour = "darkgrey", linetype = "dashed", linewidth = 0.35) +
    scale_color_jama() +
    theme(legend.position=c(0.75, 0.25), legend.box="vertical", legend.text = element_text(size = 14))
}



#' Plotting distribution of allelic ratio over mean expression
#' @param res output file from beta_binom_test
#' @param fdr_var Binary parameter to specify which test fdr's to plot (default: fdr_var = FALSE, allelic imbalance results are used for plotting )
#' @param fdr_cutoff numeric value of the fdr cut-off to use (default: fdr_cutoff = 0.05)
#' @param min_logFC Log FC value to use
#' @param glob_params output from glob_disp()
#' @keywords
#' @export
#' @examples
#' plot_MA()
plot_MA <- function(res, fdr_var = FALSE, fdr_cutoff = 0.05, min_logFC = 1, glob_params){

  cols_select <- c("AR", "tot_gene_mean", "log2FC", "fdr_var", "fdr_mean")
  MAplot <- data.frame(res[,colnames(res) %in% cols_select])
  MAplot$diffASE <- "NO"
  if (fdr_var == FALSE){
    #MAplot$diffASE[MAplot$fdr_shrunk < 0.05 & abs(MAplot$log2FC) >= min_logFC] <- "YES"
    MAplot$diffASE[MAplot$fdr_mean < fdr_cutoff] <- "YES"
  } else {
    #MAplot$diffASE[MAplot$fdr_orig < 0.05 & abs(MAplot$log2FC) >= min_logFC] <- "YES"
    MAplot$diffASE[MAplot$fdr_var < fdr_cutoff] <- "YES"
  }

  x_pos <- unname(summary(log(res$tot_gene_mean))[5])

  init <- if (fdr_var == FALSE){
    ggplot(data=MAplot, aes(x = log(tot_gene_mean), y = AR, color = -log10(fdr_mean + 1e-10))) +
      geom_point(size = 1)
  } else {
    ggplot(data=MAplot, aes(x = log(tot_gene_mean), y = AR, color = -log10(fdr_var + 1e-10))) +
      geom_point(size = 1)
  }

  p <- init +
    theme(legend.position = c(0.2, 0.8)) +
    theme_classic(base_size = 15) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1)) +
    labs(x = "Log expression", y = "Allelic Ratio") +
    geom_hline(yintercept = glob_params[1], linewidth = 1, colour = "#a3cce9") +
    #geom_hline(yintercept=1, linewidth = 1, colour = "#a3cce9", linetype = "dashed") +
    #geom_hline(yintercept=-1, linewidth = 1, colour = "#a3cce9", linetype = "dashed") +
    #scale_alpha(range = c(1, .1)) + # making the centre dots appear more muted
    scale_colour_gradient_tableau(palette = "Classic Blue", name = "-log10(fdr)") +
    #scale_color_gradient2(low = "#003C67FF", mid = "#7AA6DCFF", high = "#003C67FF", midpoint = 0.5) +
    theme(legend.direction = "horizontal", legend.position = "bottom", #legend.text=element_text(size=15), legend.title=element_text(size=15),
          legend.key.size = unit(0.5, "cm"), legend.box.spacing = unit(0, "pt")) +
    annotate("text", x=x_pos, y=1, label= paste("N signif genes allele1 skew:", nrow(MAplot[MAplot$AR > glob_params[1] & MAplot$diffASE != 'NO',])), size = 4) +
    annotate("text", x=x_pos, y=0, label= paste("N signif genes allele2 skew:", nrow(MAplot[MAplot$AR < glob_params[1] & MAplot$diffASE != 'NO',])), size = 4) +
    geom_point(data = subset(MAplot, diffASE != 'NO'), colour = "#A73030FF", size = 1) +
    #geom_point(data = subset(MAplot, diffASE != 'NO'), colour = "#0d4474") +
    ylim(0,1)
  return(p)

}

#' preparing data frame for plotting
#' @param a1_counts Reference allele counts
#' @param tot_counts Total gene counts counts
#' @param gene gene symbol or other id, the same style that is used for the row names in the count matrices
#' @keywords
#' @export
#' @examples
#' makedf()
makedf <- function(a1_counts, tot_counts, gene, metadata = NULL, order.by = "time", split.var = "group"){

    a1_counts <- as.matrix(a1_counts)
    tot_counts <- as.matrix(tot_counts)

    df <- as.data.frame(cbind(a1 = a1_counts[gene,], tot = tot_counts[gene,]))
    df$a1 <- as.numeric(df$a1)
    df$tot <- as.numeric(df$tot)
    df$AR <- df$a1/df$tot
    #df <- df[order(df$tot),]
    df <- na.omit(df)

      if (is.null(metadata)){
        df$Index <- 1:nrow(df)
      } else {
        df$group <- metadata[match(rownames(df), rownames(metadata)), colnames(metadata) == split.var]
        df$time <- metadata[match(rownames(df), rownames(metadata)), colnames(metadata) == order.by]
        df <- df[order(df$time),]
        df$Index <- 1:nrow(df)
      }

    return(df)
}

#' plotting allelic ratio distribution
#' @param plot_data output of makedf
#' @param gene gene symbol or other id, must match the gene used to generate data frame for plotting
#' @keywords
#' @export
#' @examples
#' plot_distr()
plot_distr <- function(plot_data, gene, add.density = TRUE, min_counts = 0){

  plot_data <- plot_data[plot_data$tot >= min_counts,]

  p <- plot_data %>% ggplot(aes(AR, Index, colour = log(tot))) +
    theme_classic(base_size = 20) +
    geom_pointdensity(size = 0.9) +
    #geom_smooth() +
    theme(plot.subtitle = element_text(hjust = 1)) +
    scale_colour_viridis_c(option = "mako", direction = -1, name =  "log(expression)") +
    #scale_y_log10() +
    #stat_density_2d(aes(color = ..level..), color = "black") +
    theme(legend.position = "bottom", legend.box.spacing = unit(0, "pt"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = "Allelic ratio", caption = gene) +
    xlim(0,1)

  if (add.density == TRUE){
  #geom_vline(xintercept = plot_data$Index[match(unique(plot_data$cell_type), plot_data$cell_type)][-1])
  p2 <-ggExtra::ggMarginal(p, type = "density", margin = "x", groupColour = F, groupFill = F, size = 3)

  p2
  } else{
    p
  }
}


#' plotting mean gene expression vs allelic dispersion
#' @param param_estim output of bb_test_res or correct_theta
#' @param gene gene symbol or other id, the same style that is used for the row names in the count matrices
#' @keywords
#' @export
#' @examples
#' plot_exp_disp()
plot_exp_disp <- function(param_estim, gene = gene){

  p <- ggplot(param_estim, aes(log(tot_gene_mean), log(bb_theta))) +
    #geom_pointdensity(size = 0.7) +
    geom_point(color = "darkgrey", size = 1) +
    #geom_smooth(method = "locfit", color = "#01a2d9", se = T, linewidth = 0.7) +
    #geom_point(aes(log(tot_gene_mean), log(theta_smoothed)), color = "#01a2d9", size = 0.7) +
    geom_line(aes(log(tot_gene_mean), log(theta_common)), color = "#01a2d9", linewidth = 1) +
    #geom_ribbon(
    #  aes(ymin = log(ci_lower + 0.01), ymax = log(ci_upper + 0.01)), fill = "grey70", alpha = 0.4) +
    theme_classic(base_size = 15) +
    geom_point(data = param_estim[gene,], aes(log(tot_gene_mean), log(bb_theta)), color =  "#A73030FF", size = 2) +
    geom_point(data = param_estim[param_estim$gene == gene,], aes(log(tot_gene_mean), log(bb_theta)), color =  "#A73030FF", size = 2) +
    theme(legend.position = "none", legend.title = element_blank()) +
    labs(caption = gene)

    if ("group" %in% colnames(param_estim)){
      p + geom_text_repel(data = param_estim[param_estim$gene == gene,], aes(label = group),
                         size = 3, show.legend = FALSE)
    } else {
      p  + geom_text_repel(data = param_estim[gene,], aes(label = gene),
                           size = 3, show.legend = FALSE)
    }
  #scale_color_gradient2(low = "#003C67FF", mid = "#EFC000FF", high = "#A73030FF", midpoint = midpoint) +
  #annotate("text", x=2, y=1, label= paste("N genes:", nrow(param_reestim)))
}


#' posterior counts simulated under the expected dispersion
#' @param a1_counts Reference allele counts
#' @param tot_counts Total gene counts counts
#' @param gene gene symbol or other id, the same style that is used for the row names in the count matrices
#' @param metadata Metadata object containing cell level information.
#' Metadata rownames must match the column names of the count matrices.
#' @param split.var Name of the variable which will be used to split the cells.
#' @param order.by variable to order allelic ratio by
#' @param estimates_group output of estim_bbparams_bygroup()
#' @param init.seed set seed value
#' @keywords
#' @export
#' @examples
#' make_plotdf_simul()
#instead of calculating posterior counts, the counts are simulated using dispersion
make_plotdf_simul <- function(a1_counts, tot_counts, gene, metadata, split.var = "group", order.by = "time", estimates_group, init.seed = 10001101){

    a1_counts <- as.matrix(a1_counts)
    tot_counts <- as.matrix(tot_counts)

    df <- as.data.frame(cbind(a1 = a1_counts[gene,], tot = tot_counts[gene,]))
    df$a2 <- df$tot - df$a1
    df$AR <- df$a1/df$tot
    df$group <- metadata[match(rownames(df), rownames(metadata)), colnames(metadata) == split.var]
    df$time <- metadata[match(rownames(df), rownames(metadata)), colnames(metadata) == order.by]

    #splitting counts by group
    df_split <- split(df, f = df$group)
    #extracting mean AR estimates per group
    mu_orig <- as.list(estimates_group[[gene]]$bb_mu)
    names(mu_orig) <- as.list(estimates_group[[gene]]$group)
    mu_orig <- mu_orig[levels(df$group)]
    #extracting theta estimates per group
    theta_orig <- as.list(estimates_group[[gene]]$bb_theta)
    names(theta_orig) <- as.list(estimates_group[[gene]]$group)
    theta_orig <- theta_orig[levels(df$group)]

    #extracting common theta values
    theta_common <- as.list(estimates_group[[gene]]$theta_common)
    names(theta_common) <- as.list(estimates_group[[gene]]$group)
    theta_common <- theta_common[levels(df$group)]
    #calculating alpha and beta for background
    alpha_common <- mapply(function(p,q) p/q, mu_orig, theta_common, SIMPLIFY = F)
    beta_common <- mapply(function(p,q) (1-p)/q, mu_orig, theta_common, SIMPLIFY = F)
    #adding posterior counts to df
    set.seed(init.seed)
    df_split <- mapply(function(p,q,r) {p$a1_common <- rbetabinom.ab(length(p$tot), p$tot,
                                                                     shape1 = q,
                                                                     shape2 = r);
                                            p$a2_common <- p$tot - p$a1_common;
                                            p$AR_common <- p$a1_common/(p$a1_common + p$a2_common);
    return(p)}, df_split, alpha_common, beta_common, SIMPLIFY = F)

    df <- do.call("rbind", df_split)
    #df$pc1 <- metadata$PC1_afterfilt_norm[match(rownames(df), rownames(metadata))]
    #df$pc1_resc <- rescale(df$pc1, to = c(1,0))
    #df <- df[order(df$pc1_resc),]
    df <- df[order(df$time),]
    df$Index <- 1:nrow(df)
    df <- na.omit(df)
    return(df)
}


#' Density plot comparing allelic ratio distribution under observed and expected dispersion
#' @param plot_data output of make_plotdf_simul
#' @keywords
#' @export
#' @examples
#' plot_exp_disp()
plot_theta_density <- function(plot_data){

  df <- as.data.frame(rbind(cbind(AR = plot_data$AR,
                                  group = as.character(plot_data$group),
                                  type = "Observed"),
                            #cbind(AR = plot_data$AR_post,
                            cbind(AR = plot_data$AR_common,
                                  group = as.character(plot_data$group),
                                  type = "Expected")))
  df$AR <- as.numeric(df$AR)
  #df$group <- factor(df$group, levels = c("cd8_naive", "cd8_d07_lcmv_arms", "cd8_d40_lcmv_arms"))
  df$type <- factor(df$type, levels = c("Observed", "Expected"))
  ggplot(df, aes(x=AR, color=type)) +
    geom_density() +
    theme_classic(base_size = 15) +
    xlim(0,1) +
    #scale_color_excel_new(theme = "Badge", labels = c("Observed", "Expected")) +
    scale_color_manual(values = ggthemes_data$excel$themes$Badge$accents[c(1,2)]) +
    labs(x = "Allelic Ratio") +
    theme(legend.title = element_blank(), legend.position = "top") +
    #theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
    theme(legend.box.spacing = unit(0, "pt")) +
    #facet_grid(.~ group) #'.~' before variable indicates that the plots will be arranged horizontally
    facet_grid(group  ~. )
}


#' Boxplot of allele-specific expression
#' @param plot_data output of make_plotdf_simul
#' @param allele1 names of the reference strain
#' @param allele1 names of the non-reference strain
#' @keywords
#' @export
#' @examples
#' geBoxplot()
geBoxplot <- function(plot_data, allele1 = "allele1", allele2 = "allele2"){
  #data for the boxplot
  a1_counts <- data.frame(plot_data[,c("a1","group")], allele = allele1)
  colnames(a1_counts)[1] <- "counts"
  a2_counts <- data.frame(plot_data[,c("a2","group")], allele = allele2)
  colnames(a2_counts)[1] <- "counts"
  boxplot_dt <- data.frame(rbind(a1_counts, a2_counts))
  boxplot_dt <- na.omit(boxplot_dt)
  boxplot_dt$log2counts <- log2(boxplot_dt$counts)

  ggboxplot(boxplot_dt, x = "group", y = "log2counts", ylab = "Log2 counts",
            ggtheme=theme_pubr(base_size = 15), fill = "allele", #palette = "Dark2",
            palette = ggthemes_data$excel$themes$Badge$accents, #unname(jcolors('pal8')),
            #orientation = "horizontal") + #, caption = "Mann-Whitney U test") +
            orientation = "vertical") +
    theme(legend.position = "top", legend.title = element_blank(), axis.title.x = element_blank()) +
    stat_compare_means(aes(group = allele), method = "wilcox.test", label = "p.signif", label.y = 8) +
    theme(plot.margin = margin(0.3,0.3,0,0.3, "cm")) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1)) +
    ylim(c(0,9))

}

