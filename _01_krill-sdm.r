# load libraries ###############################################################
# stats stuff
library(sdmTMB)
# library(sdmTMBextra)
library(tidyverse)
library(mgcv)
# plotting 
# library(INLA)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lubridate)

# helper functions for sdmTMB ###################################################

source('~/GitHub/Utility-code/NT-sdmTMB-helper-functions.r')

# file locations ###############################################################
home_dir = getwd()
data_dir = paste0(home_dir,'/data/')
fig_dir = paste0(home_dir,'/figures/')
# results_dir = paste0(home_dir,'/results-krill/')

# run after running models to accumulate aic output
# AIC_TABLE = get_aic_table(results_dir); AIC_TABLE

# data import, cleanup, prep ###################################################


# df = data.frame(read.csv( paste0(data_dir,"PB_NASC_v5_with_depth.csv"), header=TRUE))
df = data.frame(read.csv( paste0(data_dir,"PB_NASC_sizeclass2_with_depth.csv"), header=TRUE))

head(df)

df = df %>% rename(lat = LAT, lon = LON, temp=TEMP, sal=SAL, depth = matched_depth)

## combine groups for some analyses ############################################

# Cerataulina,Dactyliosolen,Detonula,Guinardia
# Lauderia

# df = df %>% mutate(big1 = Cera_Dact_Deto_Guin + Lauderia,
#                    all_big = big1+ Chaetoceros+Eucampia+Thalassiosira)


# add year and julian day ######################################################
df$date = lubridate::as_date(df$DT)
df$year = lubridate::year(df$date)
df$fyear = as.factor(df$year)
df$day = lubridate::yday(df$date)

# add easting / northing ####
df = sdmTMB::add_utm_columns(df, c('lon', 'lat'))
head(df)

# data prep ####################################################################
# diat_diatDino=log10((PB.diatom+1./(PB.diatom+PB.dino+1)))

# df = df %>% mutate(diat_diatDino = log10((diatom+1/(diatom+dino+1))),
#                    Dino_diatDino = log10((dino+1/(diatom+dino+1))),
#                    diat_dino = diatom/dino
#                    )
# ##############################################################################
# FIT SDM MODELS ###############################################################
################################################################################

# select spp of interest #######################################################
colnames(df)

# main species = response variable
spp = c("avNASC", "diat_diatDino", "diatom", "dino", "Pseudonitzschia")[1]
spp

df_sdm = df[!is.na(df[,spp]),]
df_sdm = df_sdm[df_sdm[,spp] >= 0,]
dim(df_sdm)

# add covarate if appropriate
phyto = c("diat_diatDino", "diatom", "dino", "Pseudonitzschia",
         'Cera_Dact_Deto_Guin', 'big1', 'all_big','NBSS_slope',NA)[8] 
phyto

# scale phyto for fitting or not
# df_sdm$phyto = scale(df_sdm[,phyto])
df_sdm$phyto = df_sdm[,phyto]

ggplot(df_sdm, aes(phyto)) + 
  geom_histogram()

# remove negative and zero values for diat_diatDino ration
# total of about 14 points
nrow(df_sdm)
if(phyto == "diat_diatDino"){df_sdm = df_sdm[df_sdm$phyto>0,]}
nrow(df_sdm)

ggplot(df_sdm, aes(phyto)) + 
  geom_histogram()
ggsave( paste0(results_dir,"phyto-hist.png"))

# make results dir #############################################################
# file for testing distributions etc
# results_dir = paste0(home_dir,'/results-',spp,'-base/') 

# file for complete analysis
if(is.na(phyto)==TRUE){results_dir = paste0(home_dir,'/results-',spp,'/')}else{
  results_dir = paste0(home_dir,'/results-',spp,'-',phyto,'/')}
results_dir
dir.create(results_dir)

################################################################################

mesh <- make_mesh(df_sdm, xy_cols = c("X", "Y"), cutoff = 10)

# scale some predictors ####
df_sdm$scale_depth = scale(-1*df_sdm$depth)
df_sdm$scale_day = scale(df_sdm$day)
df_sdm$scale_sal = scale(df_sdm$sal)
df_sdm$scale_temp = scale(df_sdm$temp)
df_sdm$scale_fl = scale(df_sdm$FL)
df_sdm$f_year = as.factor(df_sdm$year)

# how many zeros in the DV? ####
(rng = round(range(df_sdm[,spp]),3))
df0 = df_sdm[ df_sdm[,spp]== 0,]
nrow(df0)
nrow(df_sdm)
zeros = paste0("There were ", nrow(df0)," zeros in the data file out of ", 
               nrow(df_sdm), " observations total = ", 100*nrow(df0)/nrow(df_sdm), "%. Data range was: ",rng[1],'-',rng[2])
zeros
capture.output(zeros, file = paste0(results_dir, 'zeros.txt'))


# formulae for iteration below

# set phytoplankton covariate ##################################################


################################################################################

forms = list(
  paste0(spp," ~ 0 + f_year"), #1
  
  paste0(spp," ~ 0 + f_year + scale_day"), #2
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3)"), #3
  
  paste0(spp," ~ 0 + f_year + scale_depth"), #4
  paste0(spp," ~ 0 + f_year + s(scale_depth, k = 3)"), #5
  
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + scale_depth"), #6
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k = 3)"),  #7
  
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + phyto"), #8
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + phyto + I(phyto^2)"), #9
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + s(phyto,k=3)" ), #10

  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + scale_depth + phyto"),#11
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + scale_depth + + s(phyto,k=3)"),#12
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k = 3) + phyto"), #13
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k = 3) + phyto + I(phyto^2)"), #14
  paste0(spp," ~ 0 + f_year + s(scale_day, k = 3) + s(scale_depth, k = 3) + s(phyto,k=3)") #15
  )

(n = length(forms))

for(i in 1:n){ # short to look at distributions and QQ plots
  xform = as.formula(forms[[i]])
  # keep track of which run is going
  write.csv( i, paste0(results_dir,"iteration.csv"))
  print(xform)
  # prepare model name ####
  # rm(dprefix,zprefix,rprefix)
  # set prefix for day
  dprefix = paste0(spp,'')
  if(str_detect(forms[[i]],"day")){dprefix=paste0(spp,'-day')}
  if(str_detect(forms[[i]],"I\\(scale_day")){dprefix=paste0(spp,'-qrd_day')}
  if(str_detect(forms[[i]],"s\\(scale_day")){dprefix=paste0(spp,'-sm_day')}
  # set prefix for depth 
  zprefix= paste0('')
  if(str_detect(forms[[i]],"depth")){zprefix='-depth'}
  if(str_detect(forms[[i]],"I\\(scale_depth")){zprefix='-qrd_depth'}
  if(str_detect(forms[[i]],"s\\(scale_depth")){zprefix='-sm_depth'}
  # set prefix for ratio
  phyto_exists = grep('phyto', as.character(xform))
    if(length(phyto_exists)==0){rprefix=""}else{rprefix=paste0("-",phyto)}
  #rprefix = phyto
  if(str_detect(forms[[i]],"phyto")){rprefix='-linear_phyto'; sv = "~ 0 + phyto "}
  if(str_detect(forms[[i]],"I\\(phyto")){rprefix='-quadratic_phyto'; sv = "~ 0 + phyto + I(phyto^2)"}
  if(str_detect(forms[[i]],"s\\(phyto")){rprefix='-gam_phyto'; sv = "~ 0 + s(phyto , k = 3)"}
  xprefix  = paste0(dprefix,zprefix,rprefix)
  print(xprefix)
  
  if(exists('sv')){svform = as.formula(sv); svform}
  # fit model
  rm(fit)
  fit <- sdmTMB(
    formula = xform,
    # spatial_varying = svform,
    data = df_sdm,
    mesh = mesh,
    time = "year",
    spatial = 'on',
    spatiotemporal = 'iid', 
    anisotropy=TRUE,
    family = tweedie(),
    silent=FALSE
  ) # end fit
  fit
  sanity(fit)
  get_model_name(fit, prefix = xprefix)
  sdm_save_output(fit = fit, results_dir = results_dir, prefix=xprefix)
}

# end run models ###############################################################


# compile sanity results for the tweedie models
x = dir(results_dir)
x = x[grep("sdm-", x)]
x = x[grep("tweedie", x)]
x
for(i in 1:length(x)){
  print(i)
  m1 = readRDS(paste0(results_dir,x[i]))
  aic = AIC(m1)
  san = data.frame(unlist(sanity(m1)))
  convg = san['nlminb_ok',]
  mod_perf = data.frame(t(san))
  mod_perf$model = x[i]
  mod_perf$aic = aic
  mod_perf$converged = convg
  mod_perf = mod_perf[,c('model','converged','all_ok','aic',
                         "hessian_ok", "eigen_values_ok", "nlminb_ok", 
                         "range_ok", "gradients_ok","se_magnitude_ok",
                         "se_na_ok","sigmas_ok" )]
  if(i==1){dfperf = mod_perf}else{dfperf = data.frame(rbind(dfperf,mod_perf))}
}

## NOTE the AIC table isn't really appropriate for comparing 
## different distributions.  Use QQ plots to check best fit for residuals
## then look at GCV output.

dfperf$model = stringr::str_remove(dfperf$model,'sdm-')
dfperf$model = stringr::str_remove(dfperf$model,'.rds')
dfperf$delta = dfperf$aic - min(dfperf$aic, na.rm = TRUE)
dfperf = dfperf[ order(dfperf$delta),]
write.csv(dfperf,paste0(results_dir,"AIC-sanity-table-full.csv"), row.names = FALSE)
# view(dfperf)
dfperfshort = dfperf[dfperf$all_ok == TRUE,c('model','all_ok','aic','delta')]

write.csv(dfperfshort, paste0(results_dir,'AIC-sanity-table-short.csv'))
view(dfperf)
view(dfperfshort)

################################################################################
# cross-validation to select best model #########################################
################################################################################

cross_forms = forms
cross_forms

(ncf = length(cross_forms))

# run twice; once with spatial & spatio temporal on; once with spatio temporal off;
# select species etc above

cross_dir = paste0(results_dir,'cross-validation/')
dir.create(cross_dir)

# library(future)
# plan(multisession)

for(i in 1:7){ # SET AS NEEDED
  xform = as.formula(cross_forms[[i]])
  print(xform)
  if(exists('sv')){rm(sv)}
  # set prefix 
  dprefix = paste0(spp,'')
  if(str_detect(forms[[i]],"day")){dprefix=paste0(spp,'-day')}
  if(str_detect(forms[[i]],"I\\(scale_day")){dprefix=paste0(spp,'-qrd_day')}
  if(str_detect(forms[[i]],"s\\(scale_day")){dprefix=paste0(spp,'-sm_day')}
  # set prefix for depth 
  zprefix= paste0('')
  if(str_detect(forms[[i]],"depth")){zprefix='-depth'}
  if(str_detect(forms[[i]],"I\\(scale_depth")){zprefix='-qrd_depth'}
  if(str_detect(forms[[i]],"s\\(scale_depth")){zprefix='-sm_depth'}
  # set prefix for ratio
  phyto_exists = grep('phyto', as.character(xform))
  if(length(phyto_exists)==0){rprefix=""}else{rprefix=paste0("-",phyto)}
  #rprefix = phyto
  if(str_detect(forms[[i]],"phyto")){rprefix='-linear_phyto'; sv = "~ 0 + phyto "}
  if(str_detect(forms[[i]],"I\\(phyto")){rprefix='-quadratic_phyto'; sv = "~ 0 + phyto + I(phyto^2)"}
  if(str_detect(forms[[i]],"s\\(phyto")){rprefix='-gam_phyto'; sv = "~ 0 + s(phyto , k = 3)"}
  xprefix  = paste0(dprefix,zprefix,rprefix)
  if(exists('sv')){svform = as.formula(sv); svform}
  # fit model
  rm(cv_fit)
  cv_fit <- sdmTMB_cv(
    formula = xform,
    # turn on of manually
    # spatial_varying = svform,
    data = df_sdm,
    mesh = mesh,
    time = "year",
    spatial = 'on',
    spatiotemporal = 'iid', 
    anisotropy=TRUE,
    family = tweedie(),
    silent= FALSE
  ) # end fit
  mname = get_model_name(cv_fit$models[[1]], prefix = xprefix)
  mname = stringr::str_remove(mname, "sdm-")
  mname = paste0("GCV_", mname)
  mname
  saveRDS(cv_fit, paste0(cross_dir, mname))
}

### get log likelihoods ########################################################

x = dir(cross_dir)
x = x[grep("GCV_", x)]
x

for(i in 1:length(x)){
  print(i)
  xfit = readRDS( paste0(cross_dir,x[i]))
  LL = xfit$sum_loglik
  modname = x[i]
  modname = stringr::str_remove(modname,"GCV_")
  modname = stringr::str_remove(modname,".rds")
  print(modname)
  print(xfit$models[[1]])
  print(LL)
  print(xfit$fold_loglik)
  dfLL = data.frame(mname = modname,LL = LL)
  if(i == 1){df_loglik = dfLL}else(df_loglik = rbind(df_loglik, dfLL))
}
# more positive; more better!
df_loglik = df_loglik[order(df_loglik$LL, decreasing = TRUE),]
df_loglik$mname = stringr::str_remove(df_loglik$mname,"GCV_")
df_loglik$mname = stringr::str_remove(df_loglik$mname,".rds")
df_loglik = df_loglik %>% rename(model = mname)
# remove sm_day because never fits
san_tbl = data.frame(read.csv( paste0(results_dir,"AIC-sanity-table-full.csv")))
san_tbl = san_tbl[,c("model", "all_ok","aic",'delta')]
df_loglik = full_join(df_loglik, san_tbl)

write.csv(df_loglik, paste0(cross_dir,"GCV-log-likelihoods.csv"), 
          row.names = FALSE)
view(df_loglik)
df_good = df_loglik %>% filter(all_ok==TRUE)
view(df_good)
write.csv(df_good, paste0(cross_dir,"GCV-log-likelihoods-good-sanity.csv"), 
          row.names = FALSE)

################################################################################
# RELOAD BEST FIT MODEL ########################################################
################################################################################

spp = 'avNASC'
phyto = 'diat_diatDino'
results_dir = paste0(home_dir,'/results-',spp,"-", phyto,"/")
results_dir
cross_dir = paste0(results_dir,'cross-validation/')

df_good = read.csv(paste0(cross_dir,"GCV-log-likelihoods-good-sanity.csv"))

# double check row here
best_fit_name = df_good[1,'model']
best_fit_name
bfit = readRDS( paste0(results_dir,'sdm-',best_fit_name,'.rds'))
# save out best fit for future reference
capture.output(best_fit_name, file = paste0(results_dir,"best-fit-model-name.txt"))
capture.output(bfit, file = paste0(results_dir,"best-fit-model-results.txt"))
saveRDS(bfit, paste0(results_dir, "Best-fit-model.rds"))
bfit
sanity(bfit)

# # alternate reload models 
# q = dir(results_dir)
# q = q[grep("sdm",q)]
# q
# 
# best_fit_name = q[6]
# best_fit_name
# bfit = readRDS( paste0(results_dir,best_fit_name))
################################################################################
# predictions & index ##########################################################
################################################################################

x = as.character(unlist(bfit$formula)[1])
z = stringr::str_locate(x, "\\~")
spp = substring(x,1,z[1]-2)
spp
# load grid ####
wc_grid = data.frame(read.delim("Tolimieri_et_al_2015_Ecosphere_2km_grid.txt", header=TRUE))

wc_grid = sdmTMB::add_utm_columns(wc_grid, c('LON','LAT'),units = 'km' )
head(wc_grid)
# Note, depth is a negative number in both files before scaling
# need to transfer over the scaled depth from the df_sdm file NOT
# scale the depths here. Math would be wrong.

wc_grid = wc_grid %>% 
  mutate(depth_m = SRTM_M) %>%
  rename(lat = LAT, lon=LON) %>%
  filter(depth_m != -9999)

# transfer over known values BUT lots of NAs
xdepth = data.frame(depth_wc = wc_grid$depth_m)
xdepth$scale_depth = df_sdm$scale_depth[ match(xdepth$depth_wc, df_sdm$depth) ]
# quick look, don't save
# plot(scale_depth ~ depth_wc, data=xdepth)
# interpolate missing depths
lm1 = lm(scale_depth ~ depth_wc, data = xdepth)
new_depth = data.frame(depth_wc = wc_grid$depth_m)
p1 = predict(lm1, newdata = new_depth)
p2 = data.frame(depth_m = new_depth$depth_wc, scale_depth = p1)

# double check
plot(p2$scale_depth, p2$depth_m)

# add to wc file
wc_grid$scale_depth = p2$scale_depth[ match(wc_grid$depth_m,p2$depth_m)]
head(wc_grid)  
wc <- replicate_df(wc_grid, "year", unique(df_sdm$year))
wc$f_year = as.factor(wc$year)

# smoothed terms in the grid ##########
# for terms that do not have specific values in the grid (like day or phyto)
# set the value to either the mean/median or zero (if the term was scaled
# to mean = 0)

wc$scale_day <- 0
# wc$scale_depth <- 0
wc$phyto = mean(df_sdm$phyto)

# get predictions ####
p_grid = predict(bfit, newdata = wc, return_tmb_object = TRUE)
# index <- get_index(p_grid, bias_correct = TRUE, 
#                    area = rep(4, nrow(west_coast_grid)))
# # scaled for presentation
# index$index = index$est/ max(index$upr)
# index$i_upr = index$upr/ max(index$upr)
# index$i_lwr = index$lwr/ max(index$upr)
# # save out
saveRDS(p_grid, file = paste0(results_dir,'predictions-',spp,'.rds'))
# write.csv(index, paste0(results_dir,'index-',spp,'.csv'))

################################################################################
# plots ########################################################################
################################################################################

## index #######################################################################

# ggplot(index, aes(x=year,y=est)) + 
#   geom_ribbon( aes(ymin=lwr, ymax=upr), fill='grey80') +
#   geom_line() + geom_point() + ylab(spp) +
#   # scale_x_continuous(breaks = seq(2005,max(index$date),5), 
#   #                    minor_breaks = min(index$date):max(index$date))+
#   theme_bw() + theme(axis.text = element_text(size=8))
# 
# ggsave(paste0(results_dir,'index.png'), width = 4, height = 2)

## map #########################################################################
p = p_grid$data
head(p)
range(p$est)
# model_type = grep('est1',colnames(p))

### back transform est
### may need to modify for different distributions

if(bfit$family$link[1] == 'identity'){p$cpue = p$est}
if(bfit$family$link[1] != 'identity'){
  if(bfit$family$family[1]!="tweedie"){ # lognormal or gamma
    p$prob = (exp(p$est1)/(1+exp(p$est1)))
    p$cpue = (exp(p$est1)/(1+exp(p$est1))) * exp(p$est2)}else{
      p$cpue = exp(p$est) 
  }
}

# grid size is 2 x 2 km; covert value to per km2
p$cpue = p$cpue/4
p$index = p$cpue/(max(p$cpue))

### subsetting ####

min(df_sdm$lat)
latrange = c(min(df_sdm$lat),max(df_sdm$lat))
lonrange = c(NA, -118)
p = p %>% filter(lat > latrange[1])
p$Year = p$year
pshort = p %>% filter(year %in% 2019:2024)
head(p)

# map note #####################################################################
note = "-best"
################################################################################

# abundance from detla model
dist_map <- plot_dist_maps(pshort, z_var = "cpue",
                            Title = spp,
                            lat_range = latrange,
                            lon_range = lonrange,
                            color_opt = 'rbgradient',
                            midp = 50,
                            trans = 'log')
print(dist_map)
ggsave(paste0(fig_dir,'dist-map-abundance-',spp,note,'.png'),
       width = 5, height = 5)
ggsave(paste0(results_dir,'dist-map-abundance-',spp,note,'.png'),
       width = 5, height = 5)

# covarate effect if appropriate ###############################################

dist_map2 <- plot_dist_maps(pshort, z_var = "phyto",
                           Title = "zeta",
                           lat_range = latrange,
                           lon_range = lonrange,
                           color_opt = 'rbgradient',
                           midp = 0.1, 
                           trans = 'identity')
print(dist_map2)
ggsave(paste0(fig_dir,'dist-map-zeta-',spp,note,'.png'), 
       width = 5, height = 5)
ggsave(paste0(results_dir,'dist-map-zeta-',spp,note,'.png'), 
       width = 5, height = 5)

# anisotropy
anis = plot_anisotropy(bfit)
anis + theme_bw()
ggsave(paste0(results_dir,'anisotropy-',spp,'.png'), width = 3.5, height=5)

