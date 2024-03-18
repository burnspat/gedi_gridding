
# Title: gedi_gridding_validation_als_plots.R
# Authors: Patrick Burns [pb463@nau.edu] and Chris Hakkenberg
# Purpose: calculate summary stats and produce plots for GEDI-ALS validation
###############################################################################


#####
# Libraries
#####
library(Metrics)
library(ggplot2)
library(dplyr)
library(Polychrome)
library(data.table)
library(ggh4x)



stat_labs <- c('mean'='Mean', 'med'='Median', "sd"="SD", "iqr"="IQR", "p95"="95th Perc.", "shan"="Shannon's H")



# 1. NEON validation ------------------------------------------------------------
# validation data extracted per site by Chris Hakkenberg
validation_mat <- read.csv("C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\neon_validation_samples.csv")
agg_methods <- c("mean","median","sd","IQR","p95","_H")
strct_metrics <- c("RH98","PAI","RH50","FHD") 

index<-1
results_list<-list()
for (strct_metric in strct_metrics){
  for (agg_method in agg_methods){
    stat_labs <- c('mean'='Mean', 'median'='Median', "sd"="SD", "IQR"="IQR", "p95"="95th Perc.", "_H"="Shannon's H")
    
    
    tmp <- select(validation_mat,contains(c(agg_method,"site")))
    tmp <- select(tmp,contains(c(strct_metric,"site")))
    tmp <- tmp[complete.cases(tmp),]
    names(tmp)[1:2] <- c("ALS","GEDI")
    
    dt_as <- as.data.table(tmp)
    dt_as$rmse <- rmse(dt_as$ALS,dt_as$GEDI)
    dt_as$prmse <- 100*dt_as$rmse/mean(dt_as$ALS)
    dt_as$mae <- mae(dt_as$ALS,dt_as$GEDI)
    
    tmp_lm <- lm(dt_as$ALS ~ dt_as$GEDI)
    dt_as$lm_rmse <- sqrt(mean(tmp_lm$residuals^2))
    dt_as$lm_mae <- mean(abs(tmp_lm$residuals))
    dt_as$lm_ar2 <- summary(tmp_lm)$adj.r.squared
    dt_as$lm_slp <- tmp_lm$coef[[2]]
    dt_as$lm_int <- tmp_lm$coef[[1]]
    
    dt_as$aggstat <- agg_method
    dt_as$aggstat_lab <- as.character(stat_labs[agg_method])
    dt_as$gedimet <- strct_metric
    dt_as$n_1km2_samp <- nrow(dt_as)
    
    results_list[[index]] <- dt_as
    index <- index +1
  }
}
dt_neon_c <- rbindlist(results_list)

# plots
# separate plots per metrics, faceted by stat

# color palette
pal <- palette36.colors(31)

for (met in unique(dt_neon_c$gedimet)){
  
  if (met == "PAI"){
    unit <- "(m2/m2)"
  } else if (met == "FHD"){
    unit <- ""
  } else {
    unit <- "(m)"
  }
  
  dt_neon_met <- dt_neon_c[gedimet == met]
  minv <- dt_neon_met[,.(min_v = floor(min(min(ALS), min(GEDI)))), by=aggstat_lab]
  maxv <- dt_neon_met[,.(max_v = ceiling(max(max(ALS), max(GEDI)))), by=aggstat_lab]
  
  ggplot(dt_neon_met, aes(x = GEDI, y = ALS, group = site)) +
    geom_abline(aes(slope = 1, intercept = 0), col = alpha("black", 0.6), size=1.5) +
    geom_point(alpha = 0.5, size = 2, stroke = 0, aes(col=site)) +
    scale_color_manual(values = as.vector(pal)) +
    geom_smooth(method = "lm", se=FALSE, color="purple", aes(group=NULL)) +
    geom_text(mapping = aes(x = -Inf, y = Inf, 
                            label = paste0("R2 = ", signif(lm_ar2,2), "\n", 
                                           "RMSE = ", signif(rmse,2)),
                            group = NULL),
              hjust = -0.1, vjust = 1.1, color = 'purple') +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(plot.title = element_text(size=16, hjust = 0.5),
          legend.title = element_text(size=12),
          legend.text=element_text(size=11),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
    labs(title =  paste0("USA, NEON Sites, 2020+2021 \n",
                         met, " gridded at 1 km"),
         x = paste0("GEDI ", met, " ", unit),
         y= paste0("ALS ", met, " ", unit)) +
    #coord_cartesian(xlim = c(0, maxv), ylim = c(0, maxv), expand = FALSE) +
    facet_wrap(factor(dt_neon_met$aggstat_lab, 
                      levels = c("Mean", "Median", "SD", "IQR", "95th Perc.", "Shannon's H")), 
               nrow = 3, scales = c("free")) +
    ggh4x::facetted_pos_scales(
      x = list(
        scale_x_continuous(limits = c(as.numeric(minv[1,2]), as.numeric(maxv[1,2]))),
        scale_x_continuous(limits = c(as.numeric(minv[2,2]), as.numeric(maxv[2,2]))),
        scale_x_continuous(limits = c(as.numeric(minv[3,2]), as.numeric(maxv[3,2]))),
        scale_x_continuous(limits = c(as.numeric(minv[4,2]), as.numeric(maxv[4,2]))),
        scale_x_continuous(limits = c(as.numeric(minv[5,2]), as.numeric(maxv[5,2]))),
        scale_x_continuous(limits = c(as.numeric(minv[6,2]), as.numeric(maxv[6,2])))
      ),
      y = list(
        scale_y_continuous(limits = c(as.numeric(minv[1,2]), as.numeric(maxv[1,2]))),
        scale_y_continuous(limits = c(as.numeric(minv[2,2]), as.numeric(maxv[2,2]))),
        scale_y_continuous(limits = c(as.numeric(minv[3,2]), as.numeric(maxv[3,2]))),
        scale_y_continuous(limits = c(as.numeric(minv[4,2]), as.numeric(maxv[4,2]))),
        scale_y_continuous(limits = c(as.numeric(minv[5,2]), as.numeric(maxv[5,2]))),
        scale_y_continuous(limits = c(as.numeric(minv[6,2]), as.numeric(maxv[6,2])))
      )
    )
  
  
  ggsave(filename = paste0("NEON_val_",met, ".png"), plot = last_plot(), device = "png", width = 8, height = 10, units = "in", path = getwd())
}



# 2. Other ALS validation ------------------------------------------------------------
# data exported from GEE

gee_val_files <- c("C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_borneocms.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_borneosafe.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_suma.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_6000m_val_als_98p_sonoma.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_6000m_val_als_98p_coco.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_sonoma.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_coco.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH50_1000m_val_als_98p_suma.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_PAI_1000m_val_als_98p_suma.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_PAI_1000m_val_als_98p_borneosafe.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_FHD_1000m_val_als_98p_suma.csv",
                   "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_FHD_1000m_val_als_98p_borneosafe.csv")

region_lab_dict <- c('sonoma'='USA, CA, Sonoma Cty., 2013',
                     'coco'='USA, AZ, Coconino NF, 2019',
                     'borneocms'='Indonesia, Kalimantan, 2014',
                     'borneosafe'='Malaysia, Sabah, 2014',
                     'suma'='Indonesia, Jambi, 2020+2022')

stat_labs <- c('mean'='Mean', 'med'='Median', "sd"="SD", "iqr"="IQR", "p95"="95th Perc.", "shan"="Shannon's H")

# per-pixel shot thresholds to try
shot_ns <- c(2,5,10,15,20,25,30,35,40,45,50)
dt_s_list <- list()
j<-1

# run validation for each site separately
for(f in gee_val_files){
    base <- tools::file_path_sans_ext(basename(f))
    met <- unlist(strsplit(x = base, split = "_", fixed = TRUE))[2]
    met <- gsub(pattern = "rh", replacement = "RH", x = met)
    resm <- unlist(strsplit(x = base, split = "_", fixed = TRUE))[3]
    res_n <- as.numeric(gsub(pattern = "m", replacement = "", x = resm))
    res_n_km <- res_n/1000
    reg <- unlist(strsplit(x = base, split = "_", fixed = TRUE))[7]
    region_lab <- as.character(region_lab_dict[reg])
    
    dt_met <- fread(f)
    cat("Validation file", f, "\n")
    
    for(s in 1:length(shot_ns)){
    # filter by shot count
    sc <- paste0('gedi_',met,'_shotcount')
    dt_met <- dt_met[get(sc)>=shot_ns[s]]
    if(nrow(dt_met)>=20){
    
    cat("working on GEDI metric", met, "\n")
    cat("shot threshold", shot_ns[s], "\n")
    
    # reformat original DT
    dt_list<-list()
    i <- 1
    for (stat in c('mean', 'med', 'sd', 'iqr', 'p95', 'shan')){
      keep <- c(colnames(dt_met)[grepl(pattern = glob2rx(pattern = paste0("als*_", stat)), x = colnames(dt_met))],
                colnames(dt_met)[grepl(pattern = glob2rx(pattern = paste0("gedi*_", stat)), x = colnames(dt_met))])
      dt_stat <- dt_met[,..keep]
      colnames(dt_stat) <- c('ALS', 'GEDI')
      dt_stat <- dt_stat[ALS >= -10 & GEDI >= -10]
      dt_stat$strct_metric <- met
      dt_stat$agg_method <- stat
      dt_stat$aggstat_lab <- as.character(stat_labs[stat])
      
      # stats
      dt_stat$rmse <- rmse(dt_stat$ALS,dt_stat$GEDI)
      dt_stat$prmse <- 100*dt_stat$rmse/mean(dt_stat$ALS)
      dt_stat$mae <- mae(dt_stat$ALS,dt_stat$GEDI)
      
      # linear model stats
      lm <- lm(formula = "ALS ~ GEDI", data = dt_stat)
      dt_stat$lm_slp <- lm$coef[[2]]
      dt_stat$lm_int <- lm$coef[[1]]
      dt_stat$lm_ar2 <- summary(lm)$adj.r.squared
      dt_stat$lm_rmse <- sqrt(mean(lm$residuals^2))
      dt_stat$lm_rrmse <- 100*dt_stat$rmse/mean(lm$model[,1])
      mae <- mean(abs(lm$residuals))
      dt_list[[i]] <- dt_stat
      i <- i+1
    }
    dt_r <- rbindlist(dt_list)
    
    if(res_n == 1000){
      psize <- 1.5
      ptrans <- 0.5
    } else if (res_n > 1000){
      psize <- 3
      ptrans <- 0.75
    }
    
    if (met == "PAI"){
      unit <- "(m2/m2)"
    } else if (met == "FHD"){
      unit <- ""
    } else {
      unit <- "(m)"
    }
    
    minv <- dt_r[,.(min_v = floor(min(min(ALS), min(GEDI)))), by=aggstat_lab]
    maxv <- dt_r[,.(max_v = ceiling(max(max(ALS), max(GEDI)))), by=aggstat_lab]
    
    ggplot(dt_r, aes(x = GEDI, y = ALS)) +
      geom_abline(aes(slope = 1, intercept = 0), col = alpha("black", 0.6), size=1.5) +
      geom_point(alpha = ptrans, size = psize, color = '#d95f02', stroke = 0) +
      geom_smooth(method = "lm", se=FALSE, color="purple", linewidth = 1.5) +
      geom_text(mapping = aes(x = -Inf, y = Inf,
                              label = paste0("R2 = ", signif(lm_ar2,2), "\n",
                                             "RMSE = ", signif(rmse,2)),
                              group = NULL),
                hjust = -0.1, vjust = 1.1, color = 'purple') +
      theme_bw() +
      theme(aspect.ratio = 1) +
      theme(plot.title = element_text(size=16, hjust = 0.5),
            legend.title = element_text(size=12),
            legend.text=element_text(size=11),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
      labs(title =  paste0(region_lab, "\n",
                           met, " gridded at ", res_n_km, " km"),
           x = paste0("GEDI ", met, " ", unit),
           y= paste0("ALS ", met, " ", unit)) +
      #coord_cartesian(xlim = c(0, maxv), ylim = c(0, maxv), expand = FALSE) +
      facet_wrap(factor(dt_r$aggstat_lab, 
                        levels = c("Mean", "Median", "SD", "IQR", "95th Perc.", "Shannon's H")), 
                 nrow = 3, scales = c("free")) +
      ggh4x::facetted_pos_scales(
        x = list(
          scale_x_continuous(limits = c(as.numeric(minv[1,2]), as.numeric(maxv[1,2]))),
          scale_x_continuous(limits = c(as.numeric(minv[2,2]), as.numeric(maxv[2,2]))),
          scale_x_continuous(limits = c(as.numeric(minv[3,2]), as.numeric(maxv[3,2]))),
          scale_x_continuous(limits = c(as.numeric(minv[4,2]), as.numeric(maxv[4,2]))),
          scale_x_continuous(limits = c(as.numeric(minv[5,2]), as.numeric(maxv[5,2]))),
          scale_x_continuous(limits = c(as.numeric(minv[6,2]), as.numeric(maxv[6,2])))
        ),
        y = list(
          scale_y_continuous(limits = c(as.numeric(minv[1,2]), as.numeric(maxv[1,2]))),
          scale_y_continuous(limits = c(as.numeric(minv[2,2]), as.numeric(maxv[2,2]))),
          scale_y_continuous(limits = c(as.numeric(minv[3,2]), as.numeric(maxv[3,2]))),
          scale_y_continuous(limits = c(as.numeric(minv[4,2]), as.numeric(maxv[4,2]))),
          scale_y_continuous(limits = c(as.numeric(minv[5,2]), as.numeric(maxv[5,2]))),
          scale_y_continuous(limits = c(as.numeric(minv[6,2]), as.numeric(maxv[6,2])))
        )
      ) 
    
    # save individual scatter plots
    #ggsave(filename = paste0(reg, "_val_",met,"_", resm, ".png"), plot = last_plot(), device = "png", width = 8, height = 10, units = "in", path = getwd())
    
    dt_save <- dt_r[,.(rmse = mean(rmse),
                       prmse = mean(prmse),
                       mae =mean(mae),
                       ar2 = mean(lm_ar2),
                       n_samp = .N), 
                    by = c('strct_metric', 'aggstat_lab')]
    dt_save$Region <- region_lab
    dt_save$res <- res_n
    dt_save$shotn_t <- shot_ns[s]
    
    # save individual tables
    #fwrite(x = dt_save, file = paste0(reg, "_val_",met,"_", resm, "_summ_stats.csv"))
    dt_s_list[[j]] <- dt_save
    j <- j+1
    cat("\n")
    } else {
      cat("Not enough shots pixels after shot thresh...\n")
    }
  }

}

# examine the role of shot count threshold on accuracy metrics
dt_c <- rbindlist(dt_s_list)

for (m in c('RH50', 'RH98', 'PAI', 'FHD')){
  if (met == "PAI"){
    unit <- "(m2/m2)"
  } else if (met == "FHD"){
    unit <- ""
  } else {
    unit <- "(m)"
  }
  
# RMSE
dt_c_m <- dt_c[strct_metric == m & res == 1000]
ggplot(dt_c_m, aes(x= shotn_t, y = rmse, color = Region, shape = Region)) +
  geom_point(size=2) + 
  geom_smooth(se = FALSE) +
  facet_wrap(factor(dt_c_m$aggstat_lab, 
                    levels = c("Mean", "Median", "SD", "IQR", "95th Perc.", "Shannon's H")), 
             nrow = 3, scales = c("free")) +
  xlab('GEDI Shot Count Threshold')+
  ylab(paste0('RMSE ', unit)) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text=element_text(size=11),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  labs(title =  paste0(m, " gridded at 1 km")) +
  scale_shape_manual(values=c(0,1,2,3,4)) 

ggsave(filename = paste0("rmse_byshotn_",m,"_1km.png"), plot = last_plot(), device = "png", width = 8, height = 10, units = "in", path = getwd())


# R2
dt_c_m <- dt_c[strct_metric == m & res == 1000]
ggplot(dt_c_m, aes(x= shotn_t, y = ar2, color = Region, shape = Region)) +
  geom_point(size=2) + 
  geom_smooth(se = FALSE) +
  facet_wrap(factor(dt_c_m$aggstat_lab, 
                    levels = c("Mean", "Median", "SD", "IQR", "95th Perc.", "Shannon's H")), 
             nrow = 3) +
  xlab('GEDI Shot Count Threshold')+
  ylim(0,1) + 
  ylab(bquote('Adj.'~R^2)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text=element_text(size=11),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  labs(title =  paste0(m, " gridded at 1 km")) + 
  scale_shape_manual(values=c(0,1,2,3,4))

ggsave(filename = paste0("r2_byshotn_",m,"_1km.png"), plot = last_plot(), device = "png", width = 8, height = 10, units = "in", path = getwd())
}


# combine SE Asia RH98 validation
gee_val_rh98_files <- c("C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_borneocms.csv",
                        "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_borneosafe.csv",
                        "C:\\Users\\pb463\\Desktop\\gedi_gridding\\ancillary\\validation\\gedi_RH98_1000m_val_als_98p_suma.csv")
dt_all_list <- list()
i <- 1
for(f in gee_val_rh98_files){
  base <- tools::file_path_sans_ext(basename(f))
  met <- unlist(strsplit(x = base, split = "_", fixed = TRUE))[2]
  met <- gsub(pattern = "rh", replacement = "RH", x = met)
  resm <- unlist(strsplit(x = base, split = "_", fixed = TRUE))[3]
  res_n <- as.numeric(gsub(pattern = "m", replacement = "", x = resm))
  res_n_km <- res_n/1000
  reg <- unlist(strsplit(x = base, split = "_", fixed = TRUE))[7]
  region_lab <- as.character(region_lab_dict[reg])
  dt_f <- fread(f)
  dt_f$region_lab <- region_lab
  dt_all_list[[i]] <- dt_f
  i <- i+1
}
dt_se_rh98 <- rbindlist(dt_all_list)

# reformat original DT
dt_list<-list()
i <- 1
sc_min <- 2
dt_se_rh98_sc <- dt_se_rh98[gedi_RH98_shotcount >= sc_min]
for (stat in c('mean', 'med', 'sd', 'iqr', 'p95', 'shan')){
  keep <- c(colnames(dt_se_rh98_sc)[grepl(pattern = glob2rx(pattern = paste0("als*_", stat)), x = colnames(dt_se_rh98_sc))],
            colnames(dt_se_rh98_sc)[grepl(pattern = glob2rx(pattern = paste0("gedi*_", stat)), x = colnames(dt_se_rh98_sc))])
  dt_stat <- dt_se_rh98_sc[,..keep]
  colnames(dt_stat) <- c('ALS', 'GEDI')
  dt_stat <- dt_stat[ALS >= -10 & GEDI >= -10]
  dt_stat$strct_metric <- 'RH98'
  dt_stat$agg_method <- stat
  dt_stat$aggstat_lab <- as.character(stat_labs[stat])
  dt_stat$region_lab <- dt_se_rh98_sc$region_lab
  
  # stats
  dt_stat$rmse <- rmse(dt_stat$ALS,dt_stat$GEDI)
  dt_stat$prmse <- 100*dt_stat$rmse/mean(dt_stat$ALS)
  dt_stat$mae <- mae(dt_stat$ALS,dt_stat$GEDI)
  
  # linear model stats
  lm <- lm(formula = "ALS ~ GEDI", data = dt_stat)
  dt_stat$lm_slp <- lm$coef[[2]]
  dt_stat$lm_int <- lm$coef[[1]]
  dt_stat$lm_ar2 <- summary(lm)$adj.r.squared
  dt_stat$lm_rmse <- sqrt(mean(lm$residuals^2))
  dt_stat$lm_rrmse <- 100*dt_stat$rmse/mean(lm$model[,1])
  mae <- mean(abs(lm$residuals))
  dt_list[[i]] <- dt_stat
  i <- i+1
}
dt_r <- rbindlist(dt_list)
minv <- dt_r[,.(min_v = floor(min(min(ALS), min(GEDI)))), by=aggstat_lab]
maxv <- dt_r[,.(max_v = ceiling(max(max(ALS), max(GEDI)))), by=aggstat_lab]

ggplot(dt_r, aes(x = GEDI, y = ALS, color = region_lab)) +
  geom_abline(aes(slope = 1, intercept = 0), col = alpha("black", 0.6), size=1.5) +
  geom_point(alpha = 0.5, size = 2.5, stroke = 0) +
  geom_smooth(method = "lm", se=FALSE, color="purple", linewidth = 1.5) +
  geom_text(mapping = aes(x = -Inf, y = Inf,
                          label = paste0("R2 = ", signif(lm_ar2,2), "\n",
                                         "RMSE = ", signif(rmse,2)),
                          group = NULL),
            hjust = -0.1, vjust = 1.1, color = 'purple') +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        legend.title = element_text(size=12),
        legend.text=element_text(size=11),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  labs(title =  paste0("Southeast Asia Validation \n",
                       met, " gridded at ", res_n_km, " km \n",
                       ">= ", sc_min, "shots per pixel"),
       x = paste0("GEDI ", met, " ", unit),
       y= paste0("ALS ", met, " ", unit)) +
  #coord_cartesian(xlim = c(0, maxv), ylim = c(0, maxv), expand = FALSE) +
  facet_wrap(factor(dt_r$aggstat_lab, 
                    levels = c("Mean", "Median", "SD", "IQR", "95th Perc.", "Shannon's H")), 
             nrow = 3, scales = c("free")) +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(as.numeric(minv[1,2]), as.numeric(maxv[1,2]))),
      scale_x_continuous(limits = c(as.numeric(minv[2,2]), as.numeric(maxv[2,2]))),
      scale_x_continuous(limits = c(as.numeric(minv[3,2]), as.numeric(maxv[3,2]))),
      scale_x_continuous(limits = c(as.numeric(minv[4,2]), as.numeric(maxv[4,2]))),
      scale_x_continuous(limits = c(as.numeric(minv[5,2]), as.numeric(maxv[5,2]))),
      scale_x_continuous(limits = c(as.numeric(minv[6,2]), as.numeric(maxv[6,2])))
    ),
    y = list(
      scale_y_continuous(limits = c(as.numeric(minv[1,2]), as.numeric(maxv[1,2]))),
      scale_y_continuous(limits = c(as.numeric(minv[2,2]), as.numeric(maxv[2,2]))),
      scale_y_continuous(limits = c(as.numeric(minv[3,2]), as.numeric(maxv[3,2]))),
      scale_y_continuous(limits = c(as.numeric(minv[4,2]), as.numeric(maxv[4,2]))),
      scale_y_continuous(limits = c(as.numeric(minv[5,2]), as.numeric(maxv[5,2]))),
      scale_y_continuous(limits = c(as.numeric(minv[6,2]), as.numeric(maxv[6,2])))
    )
  ) 
ggsave(filename = paste0("rh98_val_seasia_1km_", sc_min, "shots.png"), plot = last_plot(), device = "png", width = 8, height = 10, units = "in", path = getwd())
