library(data.table)

n_list <- c(500, 1000, 2000, 3000, 4000, 5000)
pos_list <- c(0.5, 1, 2, 3)
misp_list <- 1:4
results <-
  rbindlist(unlist(lapply(n_list, function(n){
    unlist(lapply(pos_list, function(pos) {
  lapply(misp_list, function(misp) {
    try({
      dat <-  get_data(100000, pos)
      pi <-dat$pi
      ATE <- dat$ATE
      pos_const <- signif(min(pi, 1-pi),1)
    key <- paste0("DR_iter=", "5000", "_n=", n, "_pos=", pos, "_mode=", misp )
   data <- fread(paste0("simResultsDR/sim_results_", key, ".csv"))
   dat1 <- data[, grep( "audrie", colnames(data)), with = F ]
   colnames(dat1) <- c("estimate", "CI_left", "CI_right")
   dat1$method <- "auDRI"
   dat1$iter <- 1:nrow(dat1)
   dat2 <- data[, grep( "AIPW", colnames(data)), with = F]
   colnames(dat2) <- c("estimate", "CI_left", "CI_right")
   dat2$method <- "AIPW"
   dat2$iter <- 1:nrow(dat2)
    data <- rbind(dat1, dat2)
   data$n <- n
   data$pos <- pos_const
   data$misp <- misp
   data$ATE <- ATE
   return(data)
    })
    return(data.table())
  })
}), recursive = F)
}), recursive = F))


results[, coverage := CI_left <= ATE & CI_right >= ATE]
results <- unique(results[, .( se = sd(estimate), bias = abs(estimate - ATE), coverage = mean(coverage)), by = c("pos", "n", "misp", "method")])
results[, rmse := sqrt(bias^2 + se^2)]
results <- results[!duplicated(results[, c("pos", "n", "misp", "method"), with = F]),]



results$pos <- factor(results$pos, levels = paste0("Overlap: ", rev(sort(unique(results$pos_const)))))
method_names <- c("Rlearner", "Tlearner", "intercept", "AIPW")
names(method_names) <- c("AMLE-partially linear (*)", "AMLE-plugin (*)", "semiparametric", "AIPW")
results$method <- names(method_names)[match(results$method , method_names)]
unique(results[, max(abs(coef)), by =  c("n", "pos", "method" )])

results_summary <- unique(results[,
                                  .(bias = abs(mean(coef)-ATE), se = sd(coef), pos_const = pos_const[1], coverage =  mean(CI_left <= ATE & CI_right >= ATE))
                                  , by = c("n", "pos", "method" )])
results_summary[, mse := sqrt(se^2 + bias^2)]

#results_summary[, se := se/mse]
# results_summary[, bias := bias/mse]

library(ggplot2)
w <- 5
h <- 4
for(pos_const in unique(results_summary$pos_const)) {

  tmp <- as.data.frame(results_summary)[results_summary$pos_const == pos_const,]



  if(local){

    maxval <- min(2*min(tmp$mse[tmp$n ==500]), max(tmp$mse))
    tmp$bias <- pmin(tmp$bias, maxval)
    tmp$se <- pmin(tmp$se, maxval)
    tmp$mse <- pmin(tmp$mse, maxval)
  } else {
    maxval <- 3*min(tmp$mse[tmp$n ==500])
    tmp$bias <- pmin(tmp$bias, maxval)
    tmp$se <- pmin(tmp$se, maxval)
    tmp$mse <- pmin(tmp$mse, maxval)
  }
  limits <- c(0, maxval + .01)
  ggplot(tmp, aes(x= n, y = bias, color = method, linetype = method)) + geom_line(size = 0.8) + facet_wrap(~pos, ncol =1)  + theme_bw() +     theme(  text = element_text(size=18), axis.text.x = element_text(size = 14 , hjust = 1, vjust = 0.5), legend.position = "none",   legend.box = "horizontal"
  )  + labs(x = "Sample Size (n)", y = "Bias", color = "Method", group = "Method", linetype = "Method") +  scale_y_continuous( limits = limits)

  ggsave(file = paste0("simResults/AdaptSimsHard=", complexity,"_",pos_const, "_local_", local, "Bias.pdf") , width = w, height = h)
