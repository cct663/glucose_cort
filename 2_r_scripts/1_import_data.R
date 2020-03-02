######################################################################
## Script to analyze samples of glucose and corticosterone from     ##
## tree swallow adult and nestlings collected in NY, TN, WY, & AK   ##
## from 2016-2019. The main purpose is to ask if increase in cort   ##
## is directly correlated with increase in glucose. Different       ##
## subsets of data are used for different questions becaues not     ##
## all samples were always collected and some birds were part of    ##
## manipulative experiments that could have influenced measures.    ##
##                                                                  ##
## Code by Conor Taff, last updated 25 February 2020                ##
######################################################################

## Some differences in samples collected by age, year, and location
  # AK, TN, WY: no dex or acth glucose; no nestlings
  # NY: dex & acth glucose in some years; acth only in 2019 and only
      # on 3rd caputure of females. Some years second or third capture
      # did not have glucose measures.
  # Nestlings: only from NY in 2019. acth is 3 days after b/s/d series

## Load in packages to be used for various purposes.
  pacman::p_load(plyr, lme4, ggplot2, here, scales)

## Load in main data file that includes all glucose & cort measures
  d <- read.delim(here("/1_raw_data/data_glucose_cort.txt"))
  
## Substituting 0.47 as the minimum detectable amount of corticosterone.
  # In some cases this was done manually when running ELISA, but this is
  # checking to make sure all are set the same way here.
  
  d[which(d$b_cort < 0.47), "b_cort"] <- 0.47
  d[which(d$s_cort < 0.47), "s_cort"] <- 0.47
  d[which(d$d_cort < 0.47), "d_cort"] <- 0.47
  
## Calculate delta values for cort and glucose
  d$s_resp <- d$s_cort - d$b_cort
  d$n_feed <- d$d_cort - d$s_cort
  d$a_inc  <- d$a_cort - d$s_cort
  
  d$gluc_resp <- d$s_gluc - d$b_gluc
  d$gluc_feed <- d$d_gluc - d$s_gluc
  d$gluc_ainc <- d$a_gluc - d$s_gluc
  
## Make a separate data set for adults and nestlings and within adults
  # make one that excludes all post-treatment birds and one that excludes
  # birds from all but NY. These will be used for different analyses
    da <- subset(d, d$class == "adult")   
    da2 <- subset(da, da$post_trt == "no") # excludes all post treatment measures
    da2n <- subset(da2, da2$state == "NY") # New York only
    dn <- subset(d, d$class == "nestling")
    
    da2n <- subset(da, da$state == "NY")
    
## set colors to be used for all plotting
    b_col <- "#56B4E9"  # a light blue
    s_col <- "#E69F00"  # an orange
    d_col <- "#009E73"  # a dark green
    a_col <- "#F0E442"  # yellow
    
## Make a basic plot for NY adults showing change in cort and change in glucose
    
    par(mfrow = c(1, 2))
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
        ylab = expression(paste("Corticosterone (ng/", mu, "l)")), bty = "n",
        xlim = c(1.5, 5.5), ylim = c(0, 103), xaxs = "i", yaxs = "i",
        main = "Adult New York Corticosterone")
    axis(1, c(-10, 2, 3, 4, 5, 10), c("", "Baseline", "Stress", "Negative", "ACTH", ""))  
    mtext("Induced", 1, at = 3, line = 2)
    mtext("Feedback", 1, at = 4, line = 2)
    mtext("Challenge", 1, at = 5, line = 2)
    axis(2, seq(-20, 160, 20), las = 2)
    points(rep(2, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$b_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$s_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$d_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(5, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$a_cort, pch = 16,
           col = alpha("gray40", 0.6), cex = .5)
    boxplot(da2n$b_cort, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(b_col, 0.6))
    boxplot(da2n$s_cort, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(da2n$d_cort, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(da2n$a_cort, add = TRUE, at = 5, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = "Blood Glucose (mg/dl)", bty = "n",
         xlim = c(1.5, 5.5), ylim = c(80, 420), xaxs = "i", yaxs = "i",
         main = "Adult New York Glucose")
    axis(1, c(-10, 2, 3, 4, 5, 10), c("", "Baseline", "Stress", "Negative", "ACTH", ""))  
    mtext("Induced", 1, at = 3, line = 2)
    mtext("Feedback", 1, at = 4, line = 2)
    mtext("Challenge", 1, at = 5, line = 2)
    axis(2, seq(-50, 800, 50), las = 2)
    points(rep(2, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$b_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$s_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$d_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(5, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$a_gluc, pch = 16,
           col = alpha("gray40", 0.6), cex = 0.5)
    boxplot(da2n$b_gluc, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(b_col, 0.6))
    boxplot(da2n$s_gluc, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(da2n$d_gluc, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(da2n$a_gluc, add = TRUE, at = 5, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
    par(mfrow = c(1, 2))
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste("Corticosterone (ng/", mu, "l)")), bty = "n",
         xlim = c(1.5, 5.5), ylim = c(0, 85), xaxs = "i", yaxs = "i",
         main = "Nestling New York Corticosterone")
    axis(1, c(-10, 2, 3, 4, 5, 10), c("", "Baseline", "Stress", "Negative", "ACTH", ""))  
    mtext("Induced", 1, at = 3, line = 2)
    mtext("Feedback", 1, at = 4, line = 2)
    mtext("Challenge", 1, at = 5, line = 2)
    axis(2, seq(-20, 160, 20), las = 2)
    points(rep(2, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$b_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$s_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$d_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(5, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$a_cort, pch = 16,
           col = alpha("gray40", 0.6), cex = .5)
    boxplot(dn$b_cort, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(b_col, 0.6))
    boxplot(dn$s_cort, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(dn$d_cort, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(dn$a_cort, add = TRUE, at = 5, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = "Blood Glucose (mg/dl)", bty = "n",
         xlim = c(1.5, 5.5), ylim = c(80, 420), xaxs = "i", yaxs = "i",
         main = "Nestling New York Glucose")
    axis(1, c(-10, 2, 3, 4, 5, 10), c("", "Baseline", "Stress", "Negative", "ACTH", ""))  
    mtext("Induced", 1, at = 3, line = 2)
    mtext("Feedback", 1, at = 4, line = 2)
    mtext("Challenge", 1, at = 5, line = 2)
    axis(2, seq(-50, 800, 50), las = 2)
    points(rep(2, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$b_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$s_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$d_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(5, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$a_gluc, pch = 16,
           col = alpha("gray40", 0.6), cex = 0.5)
    boxplot(dn$b_gluc, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(b_col, 0.6))
    boxplot(dn$s_gluc, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(dn$d_gluc, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(dn$a_gluc, add = TRUE, at = 5, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
## Make a basic plot for NY adults showing change in cort and change in glucose
    
    par(mfrow = c(1, 2))
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")), bty = "n",
         xlim = c(1.5, 4.5), ylim = c(-100, 100), xaxs = "i", yaxs = "i",
         main = "Adult New York Delta Corticosterone")
    axis(1, c(-10, 2, 3, 4, 10), c("", "Base to", "Stress to", "Stress to", ""))  
    mtext("Stress", 1, at = 2, line = 2)
    mtext("Dex", 1, at = 3, line = 2)
    mtext("ACTH", 1, at = 4, line = 2)
    axis(2, seq(-500, 500, 50), las = 2)
    abline(h = 0, lty = 2)
    points(rep(2, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$s_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$n_feed, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$a_inc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    boxplot(da2n$s_resp, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(da2n$n_feed, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(da2n$a_inc, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), bty = "n",
         xlim = c(1.5, 4.5), ylim = c(-100, 150), xaxs = "i", yaxs = "i",
         main = "Adult New York Delta Glucose")
    axis(1, c(-10, 2, 3, 4, 10), c("", "Base to", "Stress to", "Stress to", ""))  
    mtext("Stress", 1, at = 2, line = 2)
    mtext("Dex", 1, at = 3, line = 2)
    mtext("ACTH", 1, at = 4, line = 2)
    axis(2, seq(-500, 500, 50), las = 2)
    abline(h = 0, lty = 2)
    points(rep(2, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$gluc_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$gluc_feed, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(da2n)) + runif(nrow(da2n), -.1, .1), da2n$gluc_ainc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    boxplot(da2n$gluc_resp, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(da2n$gluc_feed, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(da2n$gluc_ainc, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
    par(mfrow = c(1, 2))
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")), bty = "n",
         xlim = c(1.5, 4.5), ylim = c(-100, 100), xaxs = "i", yaxs = "i",
         main = "Nestling New York Delta Corticosterone")
    axis(1, c(-10, 2, 3, 4, 10), c("", "Base to", "Stress to", "Stress to", ""))  
    mtext("Stress", 1, at = 2, line = 2)
    mtext("Dex", 1, at = 3, line = 2)
    mtext("ACTH", 1, at = 4, line = 2)
    axis(2, seq(-500, 500, 50), las = 2)
    abline(h = 0, lty = 2)
    points(rep(2, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$s_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$n_feed, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$a_inc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    boxplot(dn$s_resp, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(dn$n_feed, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(dn$a_inc, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), bty = "n",
         xlim = c(1.5, 4.5), ylim = c(-160, 130), xaxs = "i", yaxs = "i",
         main = "Nestling New York Delta Glucose")
    axis(1, c(-10, 2, 3, 4, 10), c("", "Base to", "Stress to", "Stress to", ""))  
    mtext("Stress", 1, at = 2, line = 2)
    mtext("Dex", 1, at = 3, line = 2)
    mtext("ACTH", 1, at = 4, line = 2)
    axis(2, seq(-500, 500, 50), las = 2)
    abline(h = 0, lty = 2)
    points(rep(2, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$gluc_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$gluc_feed, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4, nrow(dn)) + runif(nrow(dn), -.1, .1), dn$gluc_ainc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    boxplot(dn$gluc_resp, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(s_col, 0.6))
    boxplot(dn$gluc_feed, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(d_col, 0.6))
    boxplot(dn$gluc_ainc, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(a_col, 0.6))
    
    
    
  ## Scatter plots of cort values and glucose values for adults
    par(mfrow = c(2, 2))
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 38), ylim = c(80, 370), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Adults Baseline")
    axis(1, seq(-20, 200, 5))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(da2n$b_cort, da2n$b_gluc, pch = 21, bg = alpha(b_col, 0.6), cex = 0.8)
    m <- lm(b_gluc ~ b_cort, data = da2n)
    newx = seq(min(na.omit(da2n$b_cort)), max(na.omit(da2n$b_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(b_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 110), ylim = c(80, 370), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Adults Stress Induced")
    axis(1, seq(-20, 200, 20))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(da2n$s_cort, da2n$s_gluc, pch = 21, bg = alpha(s_col, 0.6), cex = 0.8)
    m <- lm(s_gluc ~ s_cort, data = da2n)
    newx = seq(min(na.omit(da2n$s_cort)), max(na.omit(da2n$s_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(s_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 48), ylim = c(80, 370), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Adults Post Dex")
    axis(1, seq(-20, 200, 10))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(da2n$d_cort, da2n$d_gluc, pch = 21, bg = alpha(d_col, 0.6), cex = 0.8)
    m <- lm(d_gluc ~ d_cort, data = da2n)
    newx = seq(min(na.omit(da2n$d_cort)), max(na.omit(da2n$d_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(d_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 110), ylim = c(80, 420), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Adults Post ACTH")
    axis(1, seq(-20, 200, 10))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(da2n$a_cort, da2n$a_gluc, pch = 21, bg = alpha(a_col, 0.6), cex = 0.8)
    m <- lm(a_gluc ~ a_cort, data = da2n)
    newx = seq(min(na.omit(da2n$a_cort)), max(na.omit(da2n$a_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(a_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
## Scatter plots of cort values and glucose values for nestlings
    par(mfrow = c(2, 2))
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 32), ylim = c(80, 380), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Nestling Baseline")
    axis(1, seq(-20, 200, 5))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(dn$b_cort, dn$b_gluc, pch = 21, bg = alpha(b_col, 0.6), cex = 0.8)
    m <- lm(b_gluc ~ b_cort, data = dn)
    newx = seq(min(na.omit(dn$b_cort)), max(na.omit(dn$b_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(b_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 110), ylim = c(80, 380), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Nestling Stress Induced")
    axis(1, seq(-20, 200, 20))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(dn$s_cort, dn$s_gluc, pch = 21, bg = alpha(s_col, 0.6), cex = 0.8)
    m <- lm(s_gluc ~ s_cort, data = dn)
    newx = seq(min(na.omit(dn$s_cort)), max(na.omit(dn$s_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(s_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 40), ylim = c(80, 370), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Nestling Post Dex")
    axis(1, seq(-20, 200, 10))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(dn$d_cort, dn$d_gluc, pch = 21, bg = alpha(d_col, 0.6), cex = 0.8)
    m <- lm(d_gluc ~ d_cort, data = dn)
    newx = seq(min(na.omit(dn$d_cort)), max(na.omit(dn$d_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(d_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste("Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(0, 82), ylim = c(80, 420), xaxs = "i", yaxs = "i", 
         ylab = "Blood Glucose (mg/dl)", main = "New York Nestling Post ACTH")
    axis(1, seq(-20, 200, 10))    
    axis(2, seq(0, 500, 50), las = 2)    
    points(dn$a_cort, dn$a_gluc, pch = 21, bg = alpha(a_col, 0.6), cex = 0.8)
    m <- lm(a_gluc ~ a_cort, data = dn)
    newx = seq(min(na.omit(dn$a_cort)), max(na.omit(dn$a_cort)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(a_cort = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
## Delta cort vs. delta glucose in adults and nestlings
    par(mfrow = c(3, 2))
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(-50, 100), ylim = c(-100, 200), xaxs = "i", yaxs = "i", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), main = "New York Adults Stress Response")
    axis(1, seq(-200, 200, 50))    
    axis(2, seq(-500, 500, 50), las = 2)    
    points(da2n$s_resp, da2n$gluc_resp, pch = 21, bg = alpha("gray40", 0.6), cex = 0.8)
    m <- lm(gluc_resp ~ s_resp, data = da2n)
    newx = seq(min(na.omit(da2n$s_resp)), max(na.omit(da2n$s_resp)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(s_resp = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(-30, 85), ylim = c(-100, 200), xaxs = "i", yaxs = "i", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), main = "New York Nestling Stress Response")
    axis(1, seq(-500, 500, 20))    
    axis(2, seq(-500, 500, 50), las = 2)    
    points(dn$s_resp, dn$gluc_resp, pch = 21, bg = alpha("gray40", 0.6), cex = 0.8)
    m <- lm(gluc_resp ~ s_resp, data = dn)
    newx = seq(min(na.omit(dn$s_resp)), max(na.omit(dn$s_resp)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(s_resp = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(-80, 5), ylim = c(-100, 100), xaxs = "i", yaxs = "i", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), main = "New York Adults Neg Feedback")
    axis(1, seq(-200, 200, 25))    
    axis(2, seq(-500, 500, 50), las = 2) 
    points(da2n$n_feed, da2n$gluc_feed, pch = 21, bg = alpha("gray40", 0.6), cex = 0.8)
    m <- lm(gluc_feed ~ n_feed, data = da2n)
    newx = seq(min(na.omit(da2n$n_feed)), max(na.omit(da2n$n_feed)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(n_feed = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(-100, 25), ylim = c(-100, 104), xaxs = "i", yaxs = "i", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), main = "New York Nestling Neg Feedback")
    axis(1, seq(-200, 200, 50))    
    axis(2, seq(-200, 200, 50), las = 2)    
    points(dn$n_feed, dn$gluc_feed, pch = 21, bg = alpha("gray40", 0.6), cex = 0.8)
    m <- lm(gluc_feed ~ n_feed, data = dn)
    newx = seq(min(na.omit(dn$n_feed)), max(na.omit(dn$n_feed)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(n_feed = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(-25, 100), ylim = c(-90, 90), xaxs = "i", yaxs = "i", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), main = "New York Adults ACTH Response")
    axis(1, seq(-200, 200, 25))    
    axis(2, seq(-500, 500, 25), las = 2)    
    points(da2n$a_inc, da2n$gluc_ainc, pch = 21, bg = alpha("gray40", 0.6), cex = 0.8)
    m <- lm(gluc_ainc ~ a_inc, data = da2n)
    newx = seq(min(na.omit(da2n$a_inc)), max(na.omit(da2n$a_inc)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(a_inc = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")),
         bty = "n", xlim = c(-50, 70), ylim = c(-210, 200), xaxs = "i", yaxs = "i", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), main = "New York Nestling ACTH Response")
    axis(1, seq(-200, 200, 25))    
    axis(2, seq(-500, 500, 50), las = 2)    
    points(dn$a_inc, dn$gluc_ainc, pch = 21, bg = alpha("gray40", 0.6), cex = 0.8)
    m <- lm(gluc_ainc ~ a_inc, data = dn)
    newx = seq(min(na.omit(dn$a_inc)), max(na.omit(dn$a_inc)), 0.05)
    conf_interval <- predict(m, newdata = data.frame(a_inc = newx), interval = "confidence",
                             level = 0.95)
    polygon(c(rev(newx), newx), c(rev(conf_interval[, 3]), 
                                  conf_interval[, 2]), col = alpha("gray40", 0.4), border = NA)
    lines(newx, conf_interval[, 1], lwd = 2)
    
## Between population comparisons, just base and stress (no dex glucose)
    
    par(mfrow = c(1, 2))
    dany <- subset(da2, da2$state == "NY")
    dawy <- subset(da2, da2$state == "WY")
    daak <- subset(da2, da2$state == "AK")
    datn <- subset(da2, da2$state == "TN")
    
    col_ny <- "red"
    col_tn <- "orange"
    col_wy <- "yellow"
    col_ak <- "slateblue"
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste("Corticosterone (ng/", mu, "l)")), bty = "n",
         xlim = c(1.5, 5.5), ylim = c(0, 103), xaxs = "i", yaxs = "i",
         main = "Population Comparison Corticosterone")
    axis(1, c(-10, 2.5, 4.5, 10), c("", "Baseline", "Stress", ""))  
    mtext("Induced", 1, at = 4.5, line = 2)
    axis(2, seq(-20, 160, 20), las = 2)
    points(rep(2, nrow(dany)) + runif(nrow(dany), -.1, .1), dany$b_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.33, nrow(datn)) + runif(nrow(datn), -.1, .1), datn$b_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.67, nrow(dawy)) + runif(nrow(dawy), -.1, .1), dawy$b_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(daak)) + runif(nrow(daak), -.1, .1), daak$b_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    
    points(rep(4, nrow(dany)) + runif(nrow(dany), -.1, .1), dany$s_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4.33, nrow(datn)) + runif(nrow(datn), -.1, .1), datn$s_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4.67, nrow(dawy)) + runif(nrow(dawy), -.1, .1), dawy$s_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(5, nrow(daak)) + runif(nrow(daak), -.1, .1), daak$s_cort, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    
    boxplot(dany$b_cort, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(col_ny, 0.6), boxwex = .5)
    boxplot(datn$b_cort, add = TRUE, at = 2.33, axes = FALSE, outline = FALSE, col = alpha(col_tn, 0.6), boxwex = .5)
    boxplot(dawy$b_cort, add = TRUE, at = 2.67, axes = FALSE, outline = FALSE, col = alpha(col_wy, 0.6), boxwex = .5)
    boxplot(daak$b_cort, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(col_ak, 0.6), boxwex = .5)
    
    
    boxplot(dany$s_cort, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(col_ny, 0.6), boxwex = .5)
    boxplot(datn$s_cort, add = TRUE, at = 4.33, axes = FALSE, outline = FALSE, col = alpha(col_tn, 0.6), boxwex = .5)
    boxplot(dawy$s_cort, add = TRUE, at = 4.67, axes = FALSE, outline = FALSE, col = alpha(col_wy, 0.6), boxwex = .5)
    boxplot(daak$s_cort, add = TRUE, at = 5, axes = FALSE, outline = FALSE, col = alpha(col_ak, 0.6), boxwex = .5)
    
    legend("topleft", c("NY", "TN", "WY", "AK"), bty = "n", pch = 22, pt.bg = c(col_ny,
        col_tn, col_wy, col_ak), pt.cex = 2)
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = "Blood Glucose (mg/dl)", bty = "n",
         xlim = c(1.5, 5.5), ylim = c(80, 390), xaxs = "i", yaxs = "i",
         main = "Population Comparison Glucose")
    axis(1, c(-10, 2.5, 4.5, 10), c("", "Baseline", "Stress", ""))  
    mtext("Induced", 1, at = 4.5, line = 2)
    axis(2, seq(0, 500, 50), las = 2)
    points(rep(2, nrow(dany)) + runif(nrow(dany), -.1, .1), dany$b_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.33, nrow(datn)) + runif(nrow(datn), -.1, .1), datn$b_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.67, nrow(dawy)) + runif(nrow(dawy), -.1, .1), dawy$b_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(daak)) + runif(nrow(daak), -.1, .1), daak$b_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    
    points(rep(4, nrow(dany)) + runif(nrow(dany), -.1, .1), dany$s_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4.33, nrow(datn)) + runif(nrow(datn), -.1, .1), datn$s_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(4.67, nrow(dawy)) + runif(nrow(dawy), -.1, .1), dawy$s_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(5, nrow(daak)) + runif(nrow(daak), -.1, .1), daak$s_gluc, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    
    boxplot(dany$b_gluc, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(col_ny, 0.6), boxwex = .5)
    boxplot(datn$b_gluc, add = TRUE, at = 2.33, axes = FALSE, outline = FALSE, col = alpha(col_tn, 0.6), boxwex = .5)
    boxplot(dawy$b_gluc, add = TRUE, at = 2.67, axes = FALSE, outline = FALSE, col = alpha(col_wy, 0.6), boxwex = .5)
    boxplot(daak$b_gluc, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(col_ak, 0.6), boxwex = .5)
    
    
    boxplot(dany$s_gluc, add = TRUE, at = 4, axes = FALSE, outline = FALSE, col = alpha(col_ny, 0.6), boxwex = .5)
    boxplot(datn$s_gluc, add = TRUE, at = 4.33, axes = FALSE, outline = FALSE, col = alpha(col_tn, 0.6), boxwex = .5)
    boxplot(dawy$s_gluc, add = TRUE, at = 4.67, axes = FALSE, outline = FALSE, col = alpha(col_wy, 0.6), boxwex = .5)
    boxplot(daak$s_gluc, add = TRUE, at = 5, axes = FALSE, outline = FALSE, col = alpha(col_ak, 0.6), boxwex = .5)
        
    
## Population delta glucose and delta cort
    
    
    par(mfrow = c(1, 2))
    dany <- subset(da2, da2$state == "NY")
    dawy <- subset(da2, da2$state == "WY")
    daak <- subset(da2, da2$state == "AK")
    datn <- subset(da2, da2$state == "TN")
    
    col_ny <- "red"
    col_tn <- "orange"
    col_wy <- "yellow"
    col_ak <- "slateblue"
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste(Delta, "Corticosterone (ng/", mu, "l)")), bty = "n",
         xlim = c(1.5, 3.5), ylim = c(-20, 200), xaxs = "i", yaxs = "i",
         main = "Population Comparison Corticosterone")
    axis(1, c(-10, 2.5, 10), c("", "Base to", ""))  
    mtext("Stress", 1, at = 2.5, line = 2)
    axis(2, seq(-50, 260, 50), las = 2)
    points(rep(2, nrow(dany)) + runif(nrow(dany), -.1, .1), dany$s_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.33, nrow(datn)) + runif(nrow(datn), -.1, .1), datn$s_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.67, nrow(dawy)) + runif(nrow(dawy), -.1, .1), dawy$s_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(daak)) + runif(nrow(daak), -.1, .1), daak$s_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    
    boxplot(dany$s_resp, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(col_ny, 0.6), boxwex = .5)
    boxplot(datn$s_resp, add = TRUE, at = 2.33, axes = FALSE, outline = FALSE, col = alpha(col_tn, 0.6), boxwex = .5)
    boxplot(dawy$s_resp, add = TRUE, at = 2.67, axes = FALSE, outline = FALSE, col = alpha(col_wy, 0.6), boxwex = .5)
    boxplot(daak$s_resp, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(col_ak, 0.6), boxwex = .5)
    
    legend("topleft", c("NY", "TN", "WY", "AK"), bty = "n", pch = 22, pt.bg = c(col_ny,
                                                                                col_tn, col_wy, col_ak), pt.cex = 2)
    
    abline(h = 0, lty = 2)
    
    
    plot(1, 1, type = "n", yaxt = "n", xaxt = "n", xlab = "", 
         ylab = expression(paste(Delta, "Blood Glucose (mg/dl)")), bty = "n",
         xlim = c(1.5, 3.5), ylim = c(-150, 225), xaxs = "i", yaxs = "i",
         main = "Population Comparison Glucose")
    axis(1, c(-10, 2.5, 10), c("", "Base to", ""))  
    mtext("Stress", 1, at = 2.5, line = 2)
    axis(2, seq(-500, 5000, 50), las = 2)
    points(rep(2, nrow(dany)) + runif(nrow(dany), -.1, .1), dany$gluc_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.33, nrow(datn)) + runif(nrow(datn), -.1, .1), datn$gluc_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(2.67, nrow(dawy)) + runif(nrow(dawy), -.1, .1), dawy$gluc_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    points(rep(3, nrow(daak)) + runif(nrow(daak), -.1, .1), daak$gluc_resp, pch = 16, 
           col = alpha("gray40", 0.6), cex = .5)
    
    boxplot(dany$gluc_resp, add = TRUE, at = 2, axes = FALSE, outline = FALSE, col = alpha(col_ny, 0.6), boxwex = .5)
    boxplot(datn$gluc_resp, add = TRUE, at = 2.33, axes = FALSE, outline = FALSE, col = alpha(col_tn, 0.6), boxwex = .5)
    boxplot(dawy$gluc_resp, add = TRUE, at = 2.67, axes = FALSE, outline = FALSE, col = alpha(col_wy, 0.6), boxwex = .5)
    boxplot(daak$gluc_resp, add = TRUE, at = 3, axes = FALSE, outline = FALSE, col = alpha(col_ak, 0.6), boxwex = .5)
    
    abline(h = 0, lty = 2)
    
    
## Breeding stage