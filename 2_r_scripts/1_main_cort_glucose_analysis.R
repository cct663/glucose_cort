## Script to analyze samples of glucose and corticosterone from     ##
## tree swallow adult and nestlings collected in NY, TN, WY, & AK   ##
## from 2016-2019. The main purpose is to ask if increase in cort   ##
## is directly correlated with increase in glucose. Different       ##
## subsets of data are used for different questions becaues not     ##
## all samples were always collected and some birds were part of    ##
## manipulative experiments that could have influenced measures.    ##
##                                                                  ##
## Code by Conor Taff, last updated 8 November 2021                 ##

# Notes ----
## Some differences in samples collected by age, year, and location
# AK, TN, WY: no dex or acth glucose; no nestlings
# NY: dex & acth glucose in some years; acth only in 2019 and only
# on 3rd caputure of females. Some years second or third capture
# did not have glucose measures.
# Nestlings: only from NY in 2019. acth is 3 days after b/s/d series

## Load packages ----
    pacman::p_load(plyr, lme4, ggplot2, here, scales, lmerTest, sjPlot, scales, emmeans, emojifont,
                   tidyverse, raincloudplots, viridis, ggExtra, MASS, rethinking, rptR, DHARMa)

## Load & clean data ----
      # ACTH validation experiment data
            d_acth_fem <- read.delim(here::here("1_raw_data", "adult_acth_validation.txt"))
            d_acth_nest <- read.delim(here::here("1_raw_data", "nestling_acth_validation.txt"))


      # Main data
          d <- read.delim(here::here("1_raw_data/data_glucose_cort.txt"))

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

      ## Make a separate dataset for adults and exclude any post treatment measures.
            # Also make a separate dataset just for New York.
              da <- subset(d, d$class == "adult")   
            # Remove post treatment samples from corticosterone dosed birds
              da$post_trt2 <- paste(da$post_trt, da$treatment, sep = "_")
              da2 <- subset(da, da$post_trt2 != "yes_CORT_3d")
              da2 <- subset(da2, da2$post_trt2 != "yes_CORT_6d")
              da2n <- subset(da, da$state == "NY") # New York only
      
      # Make a separate dataframe for nestlings
              dn <- subset(d, d$class == "nestling")
              dn$post_trt2 <- paste(dn$post_trt, dn$treatment, sep = "_")
              
## Glucose repeatability ----
    da2nx <- subset(da2, is.na(da2$b_gluc) == FALSE & is.na(da2$s_gluc) == FALSE & is.na(da2$gluc_resp) == FALSE)          
    for(i in 1:nrow(da2nx)){
      da2nx$count[i] <- nrow(subset(da2nx, da2nx$band == da2nx$band[i]))
    }          
    da2x <- subset(da2nx, da2nx$count > 1)
              
    r_base <- rpt(b_gluc ~ (1|band), grname = "band", data = da2x, npermut = 0, datatype = "Gaussian")
    r_str <- rpt(s_gluc ~ (1|band), grname = "band", data = da2x, npermut = 0, datatype = "Gaussian")
    r_resp <- rpt(gluc_resp ~ (1|band), grname = "band", data = da2x, npermut = 0, datatype = "Gaussian")
    
    r_basec <- rpt(b_cort ~ (1|band), grname = "band", data = da2x, npermut = 0, datatype = "Gaussian")
    r_strc <- rpt(s_cort ~ (1|band), grname = "band", data = da2x, npermut = 0, datatype = "Gaussian")
    r_respc <- rpt(s_resp ~ (1|band), grname = "band", data = da2x, npermut = 0, datatype = "Gaussian")
    
## Model group means NY ----
              recode <- data.frame(type = c("b_cort", "s_cort", "d_cort", "a_cort", 
                                            "b_gluc", "s_gluc", "d_gluc", "a_gluc"),
                                   timepoint = rep(c("base", "induce", "dex", "acth"), 2))
              recode$timepoint <- factor(recode$timepoint, levels = c("base", "induce", "dex", "acth"))
              
            #Adults
              d1 <- da2n %>%
                pivot_longer(cols = c("b_cort", "s_cort", "a_cort"), names_to = "type", values_to = "cort")
              d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "a_cort"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m1 <- lmer(cort ~ timepoint + sex + (1|band), data = d1)
              m1_em <- as.data.frame(emmeans(m1, "timepoint", lmer.df = "satterthwaite"))
              m1_em$type <- factor(c("b_cort", "s_cort", "a_cort"), levels = c("b_cort", "s_cort", "a_cort"))
              m1_eml <- pivot_longer(m1_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              
              
              d1 <- da2n %>%
                pivot_longer(cols = c("b_gluc", "s_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
              d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "a_gluc"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m2 <- lmer(glucose ~ timepoint + sex + (1|band), data = d1)
              m2_em <- as.data.frame(emmeans(m2, "timepoint", lmer.df = "satterthwaite"))
              m2_em$type <- factor(c("b_gluc", "s_gluc", "a_gluc"), levels = c("b_gluc", "s_gluc", "a_gluc"))
              m2_eml <- pivot_longer(m2_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              
            #Nestlings
              d1 <- dn %>%
                pivot_longer(cols = c("b_cort", "s_cort", "a_cort"), names_to = "type", values_to = "cort")
              d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "a_cort"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1$nest <- paste(d1$site, d1$box, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m3 <- lmer(cort ~ timepoint + (1|nest), data = d1)
              m3_em <- as.data.frame(emmeans(m3, "timepoint", lmer.df = "satterthwaite"))
              m3_em$type <- factor(c("b_cort", "s_cort", "a_cort"), levels = c("b_cort", "s_cort", "a_cort"))
              m3_eml <- pivot_longer(m3_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              d1 <- dn %>%
                pivot_longer(cols = c("b_gluc", "s_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
              d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "a_gluc"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1$nest <- paste(d1$site, d1$box, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m4 <- lmer(glucose ~ timepoint + (1|nest), data = d1)
              m4_em <- as.data.frame(emmeans(m4, "timepoint", lmer.df = "satterthwaite"))
              m4_em$type <- factor(c("b_gluc", "s_gluc", "a_gluc"), levels = c("b_gluc", "s_gluc", "a_gluc"))
              m4_eml <- pivot_longer(m4_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              t1 <- tab_model(m1, m2, m3, m4, show.re.var = FALSE, show.p = FALSE,
                        dv.labels = c("Adult Corticosterone", "Adult Glucose", "Nestling Corticosterone", "Nestling Glucose"),
                        pred.labels = c("Intercept (Base / Female)", "Induced", "Post-Cortrosyn", "Sex (Male)"))
              
              saveRDS(t1,
                      here::here("2_r_scripts/table_s3.rds"))
              
              
## Plot group means NY ----
    # NY nestling       
        # Corticosterone
            d1 <- dn %>%
            pivot_longer(cols = c("b_cort", "s_cort", "a_cort"), names_to = "type", values_to = "cort")
          d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "a_cort"))
          d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
          p1 <- ggplot(data = d1, mapping = aes(x = type, y = cort, fill = type, color = type)) + 
            geom_boxplot(alpha = 0.2, outlier.shape = NA, position = position_nudge(x = -0.3), width = 0.25) +
            geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
            theme_classic() +
            scale_fill_viridis(discrete = TRUE) + 
            scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
            xlab("") + ylab("Corticosterone (ng/mL)") +
            scale_x_discrete(labels = c("Base", "Induced", "Cortrosyn")) +
            annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5) +
            theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), axis.text.x = element_text(angle = 30, hjust = 1)) 
          
          # add in confidence intevals from emmeans model
            p1 <- p1 + geom_line(data = m3_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
              geom_point(data = m3_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
          
                  
       # Glucose           
          d1 <- dn %>%
            pivot_longer(cols = c("b_gluc", "s_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
          d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "a_gluc"))
          d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
            p2 <- ggplot(data = d1, mapping = aes(x = type, y = glucose, fill = type, color = type)) + 
              geom_boxplot(alpha = 0.2, outlier.shape = NA, position = position_nudge(x = -0.3), width = 0.25) +
              geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
              theme_classic() +
              scale_fill_viridis(discrete = TRUE) + 
              scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
              xlab("") + ylab("Glucose (mg/dl)") +
              scale_x_discrete(labels = c("Base", "Induced", "Cortrosyn"))+
              annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5) +
              theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), axis.text.x = element_text(angle = 30, hjust = 1))
            
            # add in emmeans intervals
            p2 <- p2 + geom_line(data = m4_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
              geom_point(data = m4_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
            
          ggsave(here::here("2_r_scripts/figure_2b.pdf"), 
                 ggpubr::ggarrange(p1, p2),
              device = "pdf", width = 5, height = 4, units = "in")    
          
      
      # NY adult       
          # Corticosterone
          d1 <- da2n %>%
            pivot_longer(cols = c("b_cort", "s_cort", "a_cort"), names_to = "type", values_to = "cort")
          d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "a_cort"))
          d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
          p1 <- ggplot(data = d1, mapping = aes(x = type, y = cort, fill = type, color = type)) + 
            geom_boxplot(alpha = 0.2, position = position_nudge(x = -0.3), width = 0.25, outlier.shape = NA) +
            geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
            theme_classic() +
            scale_fill_viridis(discrete = TRUE) + 
            scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
            xlab("") + ylab(paste("Corticosterone (ng/mL)")) +
            scale_x_discrete(labels = c("Base", "Induced", "Cortrosyn")) +
            ylim(0, 110) +
            annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5) +
            theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), axis.text.x = element_text(angle = 30, hjust = 1))
          
                # add in confidence intevals from emmeans model
                    p1 <- p1 + geom_line(data = m1_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                      geom_point(data = m1_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
          
          # Glucose           
              d1 <- da2n %>%
                pivot_longer(cols = c("b_gluc", "s_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
              d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "a_gluc"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              p2 <- ggplot(data = d1, mapping = aes(x = type, y = glucose, fill = type, color = type)) + 
                geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.25, position = position_nudge(x = -0.3)) +
                geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
                theme_classic() +
                scale_fill_viridis(discrete = TRUE) + 
                scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
                xlab("") + ylab("Glucose (mg/dl)") +
                scale_x_discrete(labels = c("Base", "Induced", "Cortrosyn"))+
                annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5) +
                theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), axis.text.x = element_text(angle = 30, hjust = 1))
              
              # add in emmeans intervals
                  p2 <- p2 + geom_line(data = m2_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                    geom_point(data = m2_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
                
          
          ggsave(here::here("2_r_scripts/figure_1b.pdf"), 
                 ggpubr::ggarrange(p1, p2),
                 device = "pdf", width = 5, height = 4, units = "in")      
                  
              
          
         
          
          
## Plot individual variation NY ----
    dny <- rbind(da2n, dn)        
   p1 <- ggplot(data = dny, mapping = aes(x = b_cort, y = b_gluc, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     #guides(fill = FALSE, color = FALSE) +
     theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
     annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5) +
     xlab("Baseline cort \n (log ng/mL)") +
     ylab("Baseline glucose (mg/dl)") +
     theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13),
           legend.text = element_text(size = 12)) +
     coord_cartesian(xlim = c(0, 20))
   
   p2 <- ggplot(data = dny, mapping = aes(x = s_resp, y = gluc_resp, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5) +
     xlim(-50, 100) +
     xlab("Induced - baseline \n corticosterone (ng/mL)") +
     ylab("Induced - baseline \n glucose (mg/dl)") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
     theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13))
   
   p3 <- ggplot(data = dny, mapping = aes(x = n_feed, y = gluc_feed, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5) +
     xlim(-100, 25) +
     xlab("Induced - post-dex \n corticosterone (ng/mL)") +
     ylab("Induced - post-dex \n glucose (mg/dl)") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
     theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13))
   
   p4 <- ggplot(data = dny, mapping = aes(x = a_inc, y = gluc_ainc, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5) +
     xlab("Post-cortosyn - induced \n corticosterone (ng/mL)") +
     ylab("Post-cortrosyn - induced \n glucose (mg/dl)") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
     theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13))
   
   #p2m <- ggMarginal(p2, type = "boxplot", margins = "y", groupColour = TRUE, groupFill = TRUE) 
   #p3m <- ggMarginal(p3, type = "boxplot", margins = "y", groupColour = TRUE, groupFill = TRUE) 
   #p4m <- ggMarginal(p4, type = "boxplot", margins = "y", groupColour = TRUE, groupFill = TRUE, xparams = list(varwidth = FALSE))
   
   ggsave(here::here("2_r_scripts/figure_3b.pdf"),
    ggpubr::ggarrange(p1, p2, p4, nrow = 1, ncol = 3),
    device = "pdf", width = 10.5, height = 3.75, units = "in")
   
## Modeling individual variation NY ----
   
   # Adults
        # Baseline
            da2n$predictor <- da2n$b_cort
            mb <- lmer(b_gluc ~ scale(predictor) + scale(mass) + sex + (1|band), data = da2n)
            res_mb <- simulateResiduals(mb)
            plot(res_mb)
            # plotQQunif(mb)
            # plotResiduals(mb)
          
        # Induced
            da2n$predictor <- da2n$s_resp
            ms <- lmer(gluc_resp ~ scale(predictor) * scale(mass) + sex + (1|band), data = da2n)
            # plotQQunif(ms)
            # plotResiduals(ms)
            
        # Post-dex
            da2n$predictor <- da2n$n_feed
            md <- lmer(gluc_feed ~ scale(predictor) + sex + (1|band), data = da2n)
            # plotQQunif(md)
            # plotResiduals(md)
            
        # Post-cortrosyn
            da2n$predictor <- da2n$a_inc
            ma <- lm(gluc_ainc ~ scale(predictor) + scale(mass), data = da2n)
            # plotQQunif(ma)
            # plotResiduals(ma)
            
            ta <- tab_model(mb, ms, ma, show.re.var = FALSE,
                      dv.labels = c("Baseline Glucose", "Induced - Base Glucose", "Post-Cortrosyn - Induced Glucose"),
                      pred.labels = c("Intercept", "Corticosterone", "Mass", "Sex (male)", "Corticosterone * Mass"))
            saveRDS(ta, here::here("2_r_scripts/table_s4.rds"))
            
            #emmeans(ms, "sex", lmer.df = "Satterthwaite")
            
      # interaction for induced by base
           post <- mvrnorm(n = 1e6, mu = fixef(ms), vcov(ms))   
           r <- seq(-3, 3, 0.1)
           mu_neg1 <- sapply(r, function(z)mean(post[, 1] + post[, 2]*z + post[, 3] * -1 + post[, 4] * z * -1))
           ci_neg1 <- sapply(r, function(z)HPDI(post[, 1] + post[, 2]*z + post[, 3] * -1 + post[, 4] * z * -1))
           
           mu_zero <- sapply(r, function(z)mean(post[, 1] + post[, 2]*z + post[, 3] * 0 + post[, 4] * z * 0))
           ci_zero <- sapply(r, function(z)HPDI(post[, 1] + post[, 2]*z + post[, 3] * 0 + post[, 4] * z * 0))
           
           mu_pos1 <- sapply(r, function(z)mean(post[, 1] + post[, 2]*z + post[, 3] * 1 + post[, 4] * z * 1))
           ci_pos1 <- sapply(r, function(z)HPDI(post[, 1] + post[, 2]*z + post[, 3] * 1 + post[, 4] * z * 1))
           
           colss <- viridis(n = 5, option = "C")
           
           pdf(here::here("2_r_scripts/figure_4b.pdf"), width = 6.5, height = 6.5)
             plot(r, mu_neg1, lwd = 2, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlim = c(-2.1, 2.1), 
                  ylim = c(-40, 100), bty = "n", xlab = "Induced - baseline corticosterone (SD units)", ylab = "Induced - baseline glucose (mg/dl)",
                  cex.lab = 1.5)
             axis(1, seq(-5, 5, 1)) 
             axis(2, seq(-500, 500, 20), las = 2)
             abline(h = 0, lty = 2, col = "gray60")
             shade(ci_neg1, r, col = alpha(colss[2], 0.3))
             lines(r, mu_neg1, lwd = 2, col = colss[2])
             
             shade(ci_zero, r, col = alpha(colss[3], 0.3))
             lines(r, mu_zero, lwd = 2, col = colss[3])
             
             shade(ci_pos1, r, col = alpha(colss[4], 0.3))
             lines(r, mu_pos1, lwd = 2, col = colss[4])
             
             legend(-0.7, -12, c("mass -1 SD", "mass  0 SD", "mass  1 SD"), bty = "n", col = colss[2:4], lwd = 2, cex = 1.2)
            dev.off()
           
           
           
            
      # Nestling
            dn$nest <- paste(dn$site, dn$box, sep = "_")
          # Baseline
            dn$predictor <- dn$b_cort
            mb <- lmer(b_gluc ~ scale(predictor) * scale(mass) + (1|nest) + (1|ID), data = dn)
            # plotQQunif(mb)
            # plotResiduals(mb)
            
          # Induced
            dn$predictor <- dn$s_resp
            ms <- lmer(gluc_resp ~ scale(predictor) + scale(mass) + (1|nest), data = dn)
            
          # Post-dex
            dn$predictor <- dn$n_feed
            md <- lmer(gluc_feed ~ scale(predictor) * scale(mass) + (1|nest), data = dn)
            
          # Post-cortrosyn
            dn$predictor <- dn$a_inc
            ma <- lmer(gluc_ainc ~ scale(predictor) + scale(mass) + (1|nest), data = dn)
            
            tn <- tab_model(mb, ms, ma, show.re.var = FALSE, 
                      dv.labels = c("Baseline Glucose", "Induced - Base Glucose", 
                                    "Post-Cortrosyn - Induced Glucose"),
                      pred.labels = c("Intercept", "Corticosterone", "Mass", "Corticosterone * Mass"))
            saveRDS(tn, here::here("2_r_scripts/table_s5.rds"))
   
            
## Modeling population comparison ----
    # models for each state
        dwy <- subset(da2, da2$state == "WY")
        dak <- subset(da2, da2$state == "AK")
        dtn <- subset(da2, da2$state == "TN")
        
        dwy$cort_pred <- dwy$b_cort
        mwyb <- lmer(b_gluc ~ scale(cort_pred) * scale(mass) + (1|band), data = subset(dwy, dwy$sex == "F"))
        dwy$cort_pred <- dwy$s_resp
        mwyr <- lm(gluc_resp ~ scale(cort_pred) * scale(mass), data = subset(dwy, dwy$sex == "F"))
        
        dak$cort_pred <- dak$b_cort
        makb <- lmer(b_gluc ~ scale(cort_pred) * scale(mass) + (1|band), data = subset(dak, dak$sex == "F"))
        dak$cort_pred <- dak$s_resp
        makr <- lmer(gluc_resp ~ scale(cort_pred) * scale(mass) + (1|band), data = subset(dak, dak$sex == "F"))
        
        dtn$cort_pred <- dtn$b_cort
        mtnb <- lmer(b_gluc ~ scale(cort_pred) * scale(mass) + (1|band), data = subset(dtn, dtn$sex == "F"))
        dtn$cort_pred <- dtn$s_resp
        mtnr <- lm(gluc_resp ~ scale(cort_pred) * scale(mass), data = subset(dtn, dtn$sex == "F"))
        
        da2n$cort_pred <- da2n$b_cort
        mnyb <- lmer(b_gluc ~ scale(cort_pred) * scale(mass) + (1|band), data = subset(da2n, da2n$sex == "F"))
        da2n$cort_pred <- da2n$s_resp
        mnyr <- lmer(gluc_resp ~ scale(cort_pred) * scale(mass) + (1|band), data = subset(da2n, da2n$sex == "F"))
        
        tc1 <- tab_model(makb, mnyb, mtnb, mwyb, show.re.var = FALSE,
                         dv.labels = c("AK Base Glucose", "NY Base Glucose", "TN Base Glucose", "WY Base Glucose"),
                         pred.labels = c("Intercept", "Base Corticosterone", "Mass", "Corticosterone * Mass"))
        tc2 <- tab_model(makr, mnyr, mtnr, mwyr, show.re.var = FALSE,
                         dv.labels = c("AK Induced - Base Glucose", "NY Induced - Base Glucose", "TN Induced - Base Glucose", "WY Induced - Base Glucose"),
                         pred.labels = c("Intercept", "Induced - Base Corticosterone", "Mass", "Corticosterone * Mass"))
        
        saveRDS(tc1,
                here::here("2_r_scripts/table_s6.rds"))
        
        saveRDS(tc2,
                here::here("2_r_scripts/table_s7.rds"))
        

      # Do states differ in glucose levels
          md <- lmer(b_gluc ~ state + (1|band), data = da2)
          md2 <- lmer(s_gluc ~ state + (1|band), data = da2)
          md3 <- lmer(s_gluc - b_gluc ~ state + (1|band), data = da2)
          
          em_md <- as.data.frame(emmeans(md, "state"))
          em_md$type <- factor(c("AK", "NY", "TN", "WY"), levels = c("AK", "NY", "TN", "WY"))
          em_mdl <- pivot_longer(em_md, cols = c("lower.CL", "upper.CL"), values_to = "y")
          
          em_md2 <- as.data.frame(emmeans(md2, "state"))
          em_md2$type <- factor(c("AK", "NY", "TN", "WY"), levels = c("AK", "NY", "TN", "WY"))
          em_md2l <- pivot_longer(em_md2, cols = c("lower.CL", "upper.CL"), values_to = "y")
          
          em_md3 <- as.data.frame(emmeans(md3, "state"))
          em_md3$type <- factor(c("AK", "NY", "TN", "WY"), levels = c("AK", "NY", "TN", "WY"))
          em_md3l <- pivot_longer(em_md3, cols = c("lower.CL", "upper.CL"), values_to = "y")
          
          tab_model(md, md2, md3)
            
## Plotting population comparison ----
          
    # base glucose
          
          # Base glucose           
          
          
              pop1 <- ggplot(data = da2, mapping = aes(x = state, y = b_gluc, fill = state, color = state)) + 
                geom_boxplot(alpha = 0.4, outlier.shape = NA, position = position_nudge(x = -0.3), width = 0.25) +
                geom_jitter(width = 0.1, alpha = 0.6, size = 0.2) +
                theme_classic() +
                scale_fill_manual(values = c("slateblue", "coral3", "goldenrod", "purple")) + 
                scale_color_manual(values = c("slateblue", "coral3", "goldenrod", "purple")) +
                xlab("") + ylab("Baseline glucose (mg/dl)") +
                guides(color = FALSE, fill = FALSE) +
                #scale_x_discrete(labels = c("Base", "Induced", "Dex.", "Cortrosyn"))+
                annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5) +
              theme(axis.title = element_text(size = 14), axis.text.x = element_text(size = 12))
              
              # add in emmeans intervals
              pop1 <- pop1 + geom_line(data = em_mdl, mapping = aes(x = state, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                geom_point(data = em_md, mapping = aes(x = state, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
              
              
          # Stress glucose
              pop2 <- ggplot(data = da2, mapping = aes(x = state, y = s_gluc, fill = state, color = state)) + 
                geom_boxplot(alpha = 0.4, outlier.shape = NA, position = position_nudge(x = -0.3), width = 0.25) +
                geom_jitter(width = 0.1, alpha = 0.6, size = 0.2) +
                theme_classic() +
                scale_fill_manual(values = c("slateblue", "coral3", "goldenrod", "purple")) + 
                scale_color_manual(values = c("slateblue", "coral3", "goldenrod", "purple")) +
                xlab("") + ylab("Induced glucose (mg/dl)") +
                guides(color = FALSE, fill = FALSE) +
                #scale_x_discrete(labels = c("Base", "Induced", "Dex.", "Cortrosyn"))+
                annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5) +
                theme(axis.title = element_text(size = 14), axis.text.x = element_text(size = 12))
              
              # add in emmeans intervals
              pop2 <- pop2 + geom_line(data = em_md2l, mapping = aes(x = state, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                geom_point(data = em_md2, mapping = aes(x = state, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
              
          # delta glucose
              pop3 <- ggplot(data = da2, mapping = aes(x = state, y = s_gluc - b_gluc, fill = state, color = state)) + 
                geom_boxplot(alpha = 0.4, outlier.shape = NA, position = position_nudge(x = -0.3), width = 0.25) +
                geom_jitter(width = 0.1, alpha = 0.6, size = 0.2) +
                theme_classic() +
                scale_fill_manual(values = c("slateblue", "coral3", "goldenrod", "purple")) + 
                scale_color_manual(values = c("slateblue", "coral3", "goldenrod", "purple")) +
                xlab("") + ylab("Induced - base glucose (mg/dl)") +
                guides(color = FALSE, fill = FALSE) +
                #scale_x_discrete(labels = c("Base", "Induced", "Dex.", "Cortrosyn"))+
                annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5) +
                theme(axis.title = element_text(size = 14), axis.text.x = element_text(size = 12))
              
              # add in emmeans intervals
              pop3 <- pop3 + geom_line(data = em_md3l, mapping = aes(x = state, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                geom_point(data = em_md3, mapping = aes(x = state, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2)) +
                geom_hline(yintercept = 0, linetype = "dashed", col = "gray60")
              
          # save figure    
              
              ggsave(here::here("2_r_scripts/figure_5b.pdf"), 
                     ggpubr::ggarrange(pop1, pop2, pop3, nrow = 1),
                     device = "pdf", width = 7.5, height = 4, units = "in")    
          
                    
## Within-individual covariance ----
              library(standardize)
              
              
   da2n$byd <- paste(da2n$band, da2n$year, da2n$date, sep = "_")
   da2nw <- subset(da2n, is.na(da2n$b_cort) == FALSE & is.na(da2n$s_cort) == FALSE &
                     is.na(da2n$b_gluc) == FALSE & is.na(da2n$s_gluc) == FALSE)
   for(i in 1:nrow(da2nw)){
     da2nw$num_samps[i] <- nrow(subset(da2nw, da2nw$band == da2nw$band[i]))
   }
   da2nw2 <- subset(da2nw, da2nw$num_samps > 3)
   
   #Standardize within individuals
      da2nw2$b_gluc_s <- scale_by(b_gluc ~ as.factor(band), data = da2nw2)
      da2nw2$b_cort_s <- scale_by(b_cort ~ as.factor(band), data = da2nw2)
      da2nw2$s_gluc_s <- scale_by(s_gluc ~ as.factor(band), data = da2nw2)
      da2nw2$s_cort_s <- scale_by(s_cort ~ as.factor(band), data = da2nw2)
      da2nw2$s_resp_s <- scale_by(s_resp ~ as.factor(band), data = da2nw2)
      da2nw2$gluc_resp_s <- scale_by(gluc_resp ~ as.factor(band), data = da2nw2)
   
   p1 <- ggplot(data = da2nw2, mapping = aes(x = b_cort_s, y = b_gluc_s, by = as.factor(band))) +
     geom_point(alpha = 0.5) +
     geom_line(stat = "smooth", method = "lm", alpha = 0.9, color = "lightblue") +
     theme_classic() +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5) +
     geom_smooth(data = da2nw2, mapping = aes(x = b_cort_s, y = b_gluc_s, by = state), method = "lm", color = "coral3",
                 fill = "coral3") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
     ylim(-2.5, 2.5) +
     xlab("Within-individual \n base corticosterone") +
     ylab("Within-individual \n base glucose") +
     theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
   
   
   p2 <- ggplot(data = da2nw2, mapping = aes(x = s_cort_s, y = s_gluc_s, by = as.factor(band))) +
     geom_point(alpha = 0.5) +
     geom_line(stat = "smooth", method = "lm", alpha = 0.9, color = "lightblue") +
     theme_classic() +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5) +
     geom_smooth(data = da2nw2, mapping = aes(x = s_cort_s, y = s_gluc_s, by = state), method = "lm", color = "coral3",
                 fill = "coral3") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
     ylim(-2.5, 2.5) +
     xlab("Within-individual \n induced corticosterone") +
     ylab("Within-individual \n induced glucose") +
     theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
   
   
   p3 <- ggplot(data = da2nw2, mapping = aes(x = s_resp_s, y = gluc_resp_s, by = as.factor(band))) +
     geom_point(alpha = 0.5) +
     geom_line(stat = "smooth", method = "lm", alpha = 0.9, color = "lightblue") +
     theme_classic() +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5) +
     geom_smooth(data = da2nw2, mapping = aes(x = s_resp_s, y = gluc_resp_s, by = state), method = "lm", color = "coral3",
                 fill = "coral3") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
     ylim(-2.5, 2.1) +
     xlab("Within-individual \n \u0394 corticosterone") +
     ylab("Within-individual \n \u0394 glucose") +
     theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
   
   ggpubr::ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
   
   
   wi_b <- lm(b_gluc_s ~ 0 + b_cort_s, data = da2nw2)
   wi_s <- lm(s_gluc_s ~ 0 + s_cort_s, data = da2nw2)
   wi_sr <- lm(s_resp_s ~ 0 + gluc_resp_s, data = da2nw2)
   
   tab_model(wi_b, wi_s, wi_sr)
   
   ggsave(here::here("2_r_scripts/figure_6b.pdf"),
          ggpubr::ggarrange(p1, p2, p3, nrow = 1, ncol = 3),
          device = "pdf", width = 10.5, height = 3.75, units = "in")
          
## ACTH Validation Models ----
        
        # Nestlings
            nest_l <- d_acth_nest %>%
              pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                           values_to = "cort", values_drop_na = TRUE)
            nest_l <- as.data.frame(nest_l)
            nest_l$treatment <- as.factor(nest_l$treatment)
            nest_l$treatment <- relevel(nest_l$treatment, ref = "Saline")
            m_n_acth <- lmer(cort ~ timepoint*treatment + (1|band) + (1|unit_box), data = nest_l)
        
        # Adults
            fem_l <- d_acth_fem %>%
              pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                           values_to = "cort", values_drop_na = TRUE)
            fem_l <- as.data.frame(fem_l)
            fem_l$treatment <- gsub("BCC", "Saline", fem_l$treatment)
            fem_l$treatment <- gsub("BCA", "ACTH", fem_l$treatment)
            fem_l$treatment <- as.factor(fem_l$treatment)
            fem_l$treatment <- relevel(fem_l$treatment, ref = "Saline")
            m_a_acth <- lmer(cort ~ timepoint*treatment + (1|band), data = fem_l)
        
        # Make a table
            t1 <- tab_model(m_n_acth, m_a_acth, 
                      pred.labels = c("Intercept (Baseline)", "Timepoint 2", "Timepoint 3",
                                      "Cortrosyn at Baseline", "Cortrosyn at Timepoint 2", "Cortrosyn at Timepoint 3"),
                      dv.labels = c("Nestling Corticosterone", "Adult Corticosterone"))
            
            saveRDS(t1, here::here("2_r_scripts/table_s1.rds"))
        
        # Note that I want to put this table in the pdf output supplementary materials, but there is no direct way to 
        # use sjplot to put html tables into a markdown pdf. I manually saved the html as a pdf to put it in. That
        # means if this code is modified the saved pdf file needs to be overwritten with a new version.
        
## ACTH Validation Plots ----
        # These are plotted and saved to file here then pasted into the rmarkdown file for the supplemental materials
        
        # Nestlings
            nest_a <- d_acth_nest %>%
              pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                           values_to = "cort", values_drop_na = TRUE) %>%
              ggplot(aes(x = timepoint, y = cort, fill = treatment)) + 
              geom_line(mapping = aes(x = timepoint, y = cort, group = band, color = treatment), alpha = 0.45) +
              geom_boxplot(width = 0.25) + theme_classic() + xlab("Minutes After Capture") + ylab("Corticosterone ng/mL") +
              scale_x_discrete(labels = c("<3", "15", "30")) + geom_vline(xintercept = 1.15, lty = 2, col = "gray40") + 
              annotate("text", x = 1.1, y = 40, label = "Cortrosyn or Saline Injection", angle = 90) + labs(fill = "Treatment") + 
              ggtitle("15 Day Old Nestlings") +
              scale_fill_discrete(name = "Treatment", labels = c("Cortrosyn", "Saline")) + guides(color = FALSE) + 
              theme(legend.position = c(0.12, 0.9))
            ggsave(here::here("2_r_scripts/figure_s2b.pdf"), plot = nest_a, width = 6, height = 5, units = "in", device = "pdf")
        
        # Adults
            fem_a <- d_acth_fem %>%
              pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                           values_to = "cort", values_drop_na = TRUE) %>%
              ggplot(aes(x = timepoint, y = cort, fill = treatment)) + 
              geom_line(mapping = aes(x = timepoint, y = cort, group = band, color = treatment), alpha = 0.45) +
              geom_boxplot(width = 0.25, outlier.size = 0, outlier.stroke = 0) + theme_classic() + xlab("Minutes After Capture") + 
              ylab("Corticosterone ng/mL") +
              scale_x_discrete(labels = c("<3", "30", "60")) + geom_vline(xintercept = 1.15, lty = 2, col = "gray40") + 
              annotate("text", x = 1.1, y = 45, label = "Saline Injection", angle = 90) + 
              geom_vline(xintercept = 2.15, lty = 2, col = "gray40") +
              annotate("text", x = 2.1, y = 50, label = "Cortrosyn or Saline Injection", angle = 90) +
              labs(fill = "Treatment") + ggtitle("Adult Females") +
              scale_fill_discrete(name = "Treatment", labels = c("Cortrosyn", "Saline")) + guides(color = FALSE) + 
              theme(legend.position = c(0.12, 0.9))
            ggsave(here::here("2_r_scripts/figure_s1b.pdf"), plot = fem_a, width = 6, height = 5, units = "in", device = "pdf")
        
        
        
        
        
        
              
              