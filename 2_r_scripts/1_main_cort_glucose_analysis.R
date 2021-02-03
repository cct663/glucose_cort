######################################################################
## Script to analyze samples of glucose and corticosterone from     ##
## tree swallow adult and nestlings collected in NY, TN, WY, & AK   ##
## from 2016-2019. The main purpose is to ask if increase in cort   ##
## is directly correlated with increase in glucose. Different       ##
## subsets of data are used for different questions becaues not     ##
## all samples were always collected and some birds were part of    ##
## manipulative experiments that could have influenced measures.    ##
##                                                                  ##
## Code by Conor Taff, last updated 1 February 2021                 ##
######################################################################

## Some differences in samples collected by age, year, and location
# AK, TN, WY: no dex or acth glucose; no nestlings
# NY: dex & acth glucose in some years; acth only in 2019 and only
# on 3rd caputure of females. Some years second or third capture
# did not have glucose measures.
# Nestlings: only from NY in 2019. acth is 3 days after b/s/d series

## Load packages ----
    pacman::p_load(plyr, lme4, ggplot2, here, scales, lmerTest, sjPlot, scales,
                   tidyverse, raincloudplots, viridis, ggExtra, MASS, rethinking)

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
              da2 <- subset(da, da$post_trt == "no") # excludes all post treatment measures
              da2n <- subset(da2, da2$state == "NY") # New York only
      
      # Make a separate dataframe for nestlings
              dn <- subset(d, d$class == "nestling")
              

## Model group means NY ----
              recode <- data.frame(type = c("b_cort", "s_cort", "d_cort", "a_cort", 
                                            "b_gluc", "s_gluc", "d_gluc", "a_gluc"),
                                   timepoint = rep(c("base", "induce", "dex", "acth"), 2))
              recode$timepoint <- factor(recode$timepoint, levels = c("base", "induce", "dex", "acth"))
              
            #Adults
              d1 <- da2n %>%
                pivot_longer(cols = c("b_cort", "s_cort", "d_cort", "a_cort"), names_to = "type", values_to = "cort")
              d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "d_cort", "a_cort"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m1 <- lmer(cort ~ timepoint + sex + as.factor(year) + (1|band), data = d1)
              m1_em <- as.data.frame(emmeans(m1, "timepoint", lmer.df = "satterthwaite"))
              m1_em$type <- factor(c("b_cort", "s_cort", "d_cort", "a_cort"), levels = c("b_cort", "s_cort", "d_cort", "a_cort"))
              m1_eml <- pivot_longer(m1_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              
              
              d1 <- da2n %>%
                pivot_longer(cols = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
              d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m2 <- lmer(glucose ~ timepoint + sex + (1|band), data = d1)
              m2_em <- as.data.frame(emmeans(m2, "timepoint", lmer.df = "satterthwaite"))
              m2_em$type <- factor(c("b_gluc", "s_gluc", "d_gluc", "a_gluc"), levels = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"))
              m2_eml <- pivot_longer(m2_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              
            #Nestlings
              d1 <- dn %>%
                pivot_longer(cols = c("b_cort", "s_cort", "d_cort", "a_cort"), names_to = "type", values_to = "cort")
              d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "d_cort", "a_cort"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1$nest <- paste(d1$site, d1$box, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m3 <- lmer(cort ~ timepoint + (1|nest) + (1|band), data = d1)
              m3_em <- as.data.frame(emmeans(m3, "timepoint", lmer.df = "satterthwaite"))
              m3_em$type <- factor(c("b_cort", "s_cort", "d_cort", "a_cort"), levels = c("b_cort", "s_cort", "d_cort", "a_cort"))
              m3_eml <- pivot_longer(m3_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              d1 <- dn %>%
                pivot_longer(cols = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
              d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              d1$nest <- paste(d1$site, d1$box, sep = "_")
              d1 <- plyr::join(d1, recode, "type")
              
              m4 <- lmer(glucose ~ timepoint + (1|nest) + (1|band), data = d1)
              m4_em <- as.data.frame(emmeans(m4, "timepoint", lmer.df = "satterthwaite"))
              m4_em$type <- factor(c("b_gluc", "s_gluc", "d_gluc", "a_gluc"), levels = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"))
              m4_eml <- pivot_longer(m4_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
              
              t1 <- tab_model(m1, m2, m3, m4, show.re.var = FALSE, show.p = FALSE,
                        dv.labels = c("Adult Corticosterone", "Adult Glucose", "Nestling Corticosterone", "Nestling Glucose"),
                        pred.labels = c("Intercept (Base / Female)", "Induced", "Post-Dexamethasone", "Post-Cortrosyn", "Sex (Male)"))
              
              saveRDS(t1,
                      here::here("2_r_scripts/ny_ad_nestling_basic_model.rds"))
              
              
## Plot group means NY ----
    # NY nestling       
        # Corticosterone
            d1 <- dn %>%
            pivot_longer(cols = c("b_cort", "s_cort", "d_cort", "a_cort"), names_to = "type", values_to = "cort")
          d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "d_cort", "a_cort"))
          d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
          p1 <- ggplot(data = d1, mapping = aes(x = type, y = cort, fill = type, color = type)) + 
            geom_boxplot(alpha = 0.2, outlier.shape = NA, position = position_nudge(x = -0.3), width = 0.25) +
            geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
            theme_classic() +
            scale_fill_viridis(discrete = TRUE) + 
            scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
            xlab("") + ylab(paste("Corticosterone (ng/\u03BCl)")) +
            scale_x_discrete(labels = c("Base", "Induced", "Dex.", "Cortrosyn")) +
            annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5)
          
          # add in confidence intevals from emmeans model
            p1 <- p1 + geom_line(data = m3_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
              geom_point(data = m3_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
          
                  
       # Glucose           
          d1 <- dn %>%
            pivot_longer(cols = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
          d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"))
          d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
            p2 <- ggplot(data = d1, mapping = aes(x = type, y = glucose, fill = type, color = type)) + 
              geom_boxplot(alpha = 0.2, outlier.shape = NA, position = position_nudge(x = -0.3), width = 0.25) +
              geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
              theme_classic() +
              scale_fill_viridis(discrete = TRUE) + 
              scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
              xlab("") + ylab("Glucose (mg/dl)") +
              scale_x_discrete(labels = c("Base", "Induced", "Dex.", "Cortrosyn"))+
              annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5)
            
            # add in emmeans intervals
            p2 <- p2 + geom_line(data = m4_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
              geom_point(data = m4_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
            
          ggsave(here::here("2_r_scripts/NY_basic_comparison_nestling.png"), 
              ggpubr::ggarrange(p1, p2),
              device = "png", width = 5, height = 4, units = "in")    
      
      # NY adult       
          # Corticosterone
          d1 <- da2n %>%
            pivot_longer(cols = c("b_cort", "s_cort", "d_cort", "a_cort"), names_to = "type", values_to = "cort")
          d1$type <- factor(d1$type, levels = c("b_cort", "s_cort", "d_cort", "a_cort"))
          d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
          p1 <- ggplot(data = d1, mapping = aes(x = type, y = cort, fill = type, color = type)) + 
            geom_boxplot(alpha = 0.2, position = position_nudge(x = -0.3), width = 0.25, outlier.shape = NA) +
            geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
            theme_classic() +
            scale_fill_viridis(discrete = TRUE) + 
            scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
            xlab("") + ylab(paste("Corticosterone (ng/\u03BCl)")) +
            scale_x_discrete(labels = c("Base", "Induced", "Dex.", "Cortrosyn")) +
            ylim(0, 110) +
            annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5)
          
                # add in confidence intevals from emmeans model
                    p1 <- p1 + geom_line(data = m1_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                      geom_point(data = m1_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
          
          # Glucose           
              d1 <- da2n %>%
                pivot_longer(cols = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"), names_to = "type", values_to = "glucose")
              d1$type <- factor(d1$type, levels = c("b_gluc", "s_gluc", "d_gluc", "a_gluc"))
              d1$unique <- paste(d1$band, d1$year, d1$date, sep = "_")
              p2 <- ggplot(data = d1, mapping = aes(x = type, y = glucose, fill = type, color = type)) + 
                geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.25, position = position_nudge(x = -0.3)) +
                geom_jitter(width = 0.1, alpha = 0.4, size = 0.2) +
                theme_classic() +
                scale_fill_viridis(discrete = TRUE) + 
                scale_color_viridis(discrete = TRUE) + guides(fill = FALSE, color = FALSE) +
                xlab("") + ylab("Glucose (mg/dl)") +
                scale_x_discrete(labels = c("Base", "Induced", "Dex.", "Cortrosyn"))+
                annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5)
              
              # add in emmeans intervals
                  p2 <- p2 + geom_line(data = m2_eml, mapping = aes(x = type, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                    geom_point(data = m2_em, mapping = aes(x = type, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2))
                
          
          ggsave(here::here("2_r_scripts/NY_basic_comparison.png"), 
                 ggpubr::ggarrange(p1, p2),
                 device = "png", width = 5, height = 4, units = "in")      
                  
              
          
         
          
          
## Plot individual variation NY ----
    dny <- rbind(da2n, dn)        
   p1 <- ggplot(data = dny, mapping = aes(x = log(b_cort), y = b_gluc, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     #guides(fill = FALSE, color = FALSE) +
     theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
     annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5) +
     xlab("Baseline corticosterone (log ng/\u03BCl)") +
     ylab("Baseline glucose (mg/dl)")
   
   p2 <- ggplot(data = dny, mapping = aes(x = s_resp, y = gluc_resp, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5) +
     xlim(-50, 100) +
     xlab("Induced - base corticosterone (ng/\u03BCl)") +
     ylab("Induced - base glucose (mg/dl)") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray60")
   
   p3 <- ggplot(data = dny, mapping = aes(x = n_feed, y = gluc_feed, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5) +
     xlim(-100, 25) +
     xlab("Induced - post-dex corticsterone (ng/\u03BCl)") +
     ylab("Induced - post-dex glucose (mg/dl)") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray60")
   
   p4 <- ggplot(data = dny, mapping = aes(x = a_inc, y = gluc_ainc, color = class, fill = class)) +
     geom_point(alpha = 0.5, size = 0.7) +
     geom_smooth(method = "lm") +
     theme_classic() +
     scale_fill_manual(values = c("slateblue", "orange")) +
     scale_color_manual(values = c("slateblue", "orange")) +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1.5) +
     xlab("Post-cotrosyn - induced corticosterone (ng/\u03BCl)") +
     ylab("Post-cortrosyn - induced glucose (mg/dl)") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray60")
   
   #p2m <- ggMarginal(p2, type = "boxplot", margins = "y", groupColour = TRUE, groupFill = TRUE) 
   #p3m <- ggMarginal(p3, type = "boxplot", margins = "y", groupColour = TRUE, groupFill = TRUE) 
   #p4m <- ggMarginal(p4, type = "boxplot", margins = "y", groupColour = TRUE, groupFill = TRUE, xparams = list(varwidth = FALSE))
   
   ggsave(here::here("2_r_scripts/delta_plots.png"),
    ggpubr::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2),
    device = "png", width = 8, height = 7, units = "in")
   
## Modeling individual variation NY ----
   
   # Adults
        # Baseline
            da2n$predictor <- da2n$b_cort
            mb <- lmer(b_gluc ~ scale(predictor) + sex + (1|band), data = da2n)
          
        # Induced
            da2n$predictor <- da2n$s_resp
            ms <- lmer(gluc_resp ~ scale(predictor) * scale(mass) + sex + (1|band), data = da2n)
            
        # Post-dex
            da2n$predictor <- da2n$n_feed
            md <- lmer(gluc_feed ~ scale(predictor) + sex + (1|band), data = da2n)
            
        # Post-cortrosyn
            da2n$predictor <- da2n$a_inc
            ma <- lm(gluc_ainc ~ scale(predictor), data = da2n)
            
            ta <- tab_model(mb, ms, md, ma, show.re.var = FALSE,
                      dv.labels = c("Baseline Glucose", "Induced - Base Glucose", "Induced - Post-Dex Glucose", "Post-Cortrosyn - Induced Glucose"),
                      pred.labels = c("Intercept", "Corticosterone", "Mass", "Sex (male)", "Corticosterone * Mass"))
            saveRDS(ta, here::here("2_r_scripts/adult_covariation.rds"))
            
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
           
           png(here::here("2_r_scripts/mass_interaction.png"), width = 6.5, height = 6.5, units = "in", res = 300)
             plot(r, mu_neg1, lwd = 2, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlim = c(-2.1, 2.1), 
                  ylim = c(-40, 100), bty = "n", xlab = "Induced - baseline corticosterone (SD units)", ylab = "Induced - baseline glucose (mg/dl)")
             axis(1, seq(-5, 5, 1)) 
             axis(2, seq(-500, 500, 20), las = 2)
             abline(h = 0, lty = 2, col = "gray60")
             shade(ci_neg1, r, col = alpha(colss[2], 0.3))
             lines(r, mu_neg1, lwd = 2, col = colss[2])
             
             shade(ci_zero, r, col = alpha(colss[3], 0.3))
             lines(r, mu_zero, lwd = 2, col = colss[3])
             
             shade(ci_pos1, r, col = alpha(colss[4], 0.3))
             lines(r, mu_pos1, lwd = 2, col = colss[4])
             
             legend(-0.55, -15, c("mass -1 SD", "mass  0 SD", "mass  1 SD"), bty = "n", col = colss[2:4], lwd = 2)
            dev.off()
           
           
           
            
      # Nestling
            dn$nest <- paste(dn$site, dn$box, sep = "_")
          # Baseline
            dn$predictor <- dn$b_cort
            mb <- lmer(b_gluc ~ scale(predictor) + scale(mass) + (1|nest), data = dn)
            
          # Induced
            dn$predictor <- dn$s_resp
            ms <- lmer(gluc_resp ~ scale(predictor) + scale(mass) + (1|nest), data = dn)
            
          # Post-dex
            dn$predictor <- dn$n_feed
            md <- lmer(gluc_feed ~ scale(predictor) + scale(mass) + (1|nest), data = dn)
            
          # Post-cortrosyn
            dn$predictor <- dn$a_inc
            ma <- lmer(gluc_ainc ~ scale(predictor) + scale(mass) + (1|nest), data = dn)
            
            tn <- tab_model(mb, ms, md, ma, show.re.var = FALSE, 
                      dv.labels = c("Baseline Glucose", "Induced - Base Glucose", "Induced - Post-Dex Glucose", 
                                    "Post-Cortrosyn - Induced Glucose"),
                      pred.labels = c("Intercept", "Corticosterone", "Mass"))
            saveRDS(tn, here::here("2_r_scripts/nestling_covariation.rds"))
   
            
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
                here::here("2_r_scripts/pop_base_glucose.rds"))
        
        saveRDS(tc2,
                here::here("2_r_scripts/pop_change_glucose.rds"))
        

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
                annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5)
              
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
                annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5)
              
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
                annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5)
              
              # add in emmeans intervals
              pop3 <- pop3 + geom_line(data = em_md3l, mapping = aes(x = state, y = y), col = "black", size = 1, position = position_nudge(x = 0.2)) +
                geom_point(data = em_md3, mapping = aes(x = state, y = emmean), color = "black", shape = 23, position = position_nudge(x = 0.2)) +
                geom_hline(yintercept = 0, linetype = "dashed", col = "gray60")
              
          # save figure    
              
              ggsave(here::here("2_r_scripts/pop_comparison.png"), 
                     ggpubr::ggarrange(pop1, pop2, pop3, nrow = 1),
                     device = "png", width = 7.5, height = 4, units = "in")    
          
                    
## Within-individual covariance DELETE ----
   ## DELETE THIS SECTION
              library(standardize)
              
              
   da2n$byd <- paste(da2n$band, da2n$year, da2n$date, sep = "_")
   ggplot(data = da2n, mapping = aes(x = s_resp, y = gluc_resp, color = as.factor(band))) +
     geom_point(alpha = 0.5) +
     geom_smooth(method = "lm", se = FALSE) +
     theme_classic() +
     guides(fill = FALSE, color = FALSE) +
     annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5) +
     xlim(-5, 100) + ylim(-75, 150)
   
   
   da2nx <- subset(da2n, is.na(da2n$s_resp) == FALSE & is.na(da2n$gluc_resp) == FALSE)
   slopes <- data.frame(band = unique(da2nx$band))
   for(i in 1:nrow(slopes)){
     sub <- subset(da2nx, da2nx$band == slopes$band[i])
     if(nrow(sub) > 1){
       mx <- lm(b_gluc ~ b_cort, data = sub)
       slopes$intercept[i] <- coef(mx)[1]
       slopes$slope[i] <- coef(mx)[2]
       slopes$count[i] <- nrow(sub)
     }
   }
   
   ggplot(data = slopes, mapping = aes(y = slope, x = 1)) + geom_boxplot() + geom_jitter(width = 0.1) + ylim(-5, 5)
   
   da2n$b_gluc_s <- scale_by(b_gluc ~ as.factor(band), data = da2n, scale = 0)
   da2n$b_cort_s <- scale_by(b_cort ~ as.factor(band), data = da2n, scale = 0)
   da2nx <- subset(da2n, da2n$b_cort_s < 3 & da2n$b_cort_s > -5)
   ggplot(data = da2nx, mapping = aes(x = b_cort_s, y = b_gluc_s)) +
     geom_smooth(method = "lm", se = TRUE) +
     guides(fill = FALSE, color = FALSE) +
     geom_point()
   
   da2n$s_resp_s <- scale_by(s_resp ~ as.factor(band), data = da2n, scale = 0)
   da2n$g_resp_s <- scale_by(gluc_resp ~ as.factor(band), data = da2n, scale = 0)
   da2nx <- subset(da2n, da2n$s_resp_s > -40 & da2n$s_resp_s < 75)
   ggplot(data = da2nx, mapping = aes(x = s_resp_s, y = g_resp_s)) +
     geom_smooth(method = "lm", se = TRUE) +
     guides(fill = FALSE, color = FALSE) +
     geom_point()
   
          
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
            
            saveRDS(t1, here::here("2_r_scripts/acth_table.rds"))
        
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
              geom_boxplot(width = 0.25) + theme_classic() + xlab("Minutes After Capture") + ylab(expression(paste("Corticosterone ng/", mu, "l"))) +
              scale_x_discrete(labels = c("<3", "15", "30")) + geom_vline(xintercept = 1.15, lty = 2, col = "gray40") + 
              annotate("text", x = 1.1, y = 40, label = "Cortrosyn or Saline Injection", angle = 90) + labs(fill = "Treatment") + 
              ggtitle("15 Day Old Nestlings") +
              scale_fill_discrete(name = "Treatment", labels = c("Cortrosyn", "Saline")) + guides(color = FALSE) + 
              theme(legend.position = c(0.12, 0.9))
            ggsave(here::here("2_r_scripts/cortrosyn_nestlings.png"), plot = nest_a, width = 6, height = 5, units = "in", device = "png")
        
        # Adults
            fem_a <- d_acth_fem %>%
              pivot_longer(cols = c("cort1", "cort2", "cort3"), names_to = "timepoint",
                           values_to = "cort", values_drop_na = TRUE) %>%
              ggplot(aes(x = timepoint, y = cort, fill = treatment)) + 
              geom_line(mapping = aes(x = timepoint, y = cort, group = band, color = treatment), alpha = 0.45) +
              geom_boxplot(width = 0.25, outlier.size = 0, outlier.stroke = 0) + theme_classic() + xlab("Minutes After Capture") + 
              ylab(expression(paste("Corticosterone ng/", mu, "l"))) +
              scale_x_discrete(labels = c("<3", "30", "60")) + geom_vline(xintercept = 1.15, lty = 2, col = "gray40") + 
              annotate("text", x = 1.1, y = 45, label = "Saline Injection", angle = 90) + 
              geom_vline(xintercept = 2.15, lty = 2, col = "gray40") +
              annotate("text", x = 2.1, y = 50, label = "Cortrosyn or Saline Injection", angle = 90) +
              labs(fill = "Treatment") + ggtitle("Adult Females") +
              scale_fill_discrete(name = "Treatment", labels = c("Cortrosyn", "Saline")) + guides(color = FALSE) + 
              theme(legend.position = c(0.12, 0.9))
            ggsave(here::here("2_r_scripts/cortrosyn_adults.png"), plot = fem_a, width = 6, height = 5, units = "in", device = "png")
        
        
        
        
        
        
              
              