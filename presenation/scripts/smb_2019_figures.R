library(tidyverse)
library(R.utils)
library(lme4)
library(car)
library(gtable)
library(grid)
library(raster)
library(RColorBrewer)
library(maps)
library(mapdata)
library(gridBase)
library(mapproj)
library(grDevices)

# source the optimality model
source('optimal_vcmax_R/calc_optimal_vcmax_revised.R')
sourceDirectory('optimal_vcmax_R/functions')

# VPD from T function
vpd_from_t <- function(t, rh){ # t is temperature in °C and rh is relative humidity (%), from http://cronklab.wikidot.com/calculation-of-vapour-pressure-deficit
  
  svp = 610.7 * 10^(7.5*t/(237.3+t))
  vpd = ((100-rh) / 100)*svp
  vpd #Pa
  
}

calc_vcmax_tresp_mult = function(tleaf, tmean, tref){
  
  temp = tleaf + 273.15
  Ha= 71513
  Hd= 200000
  adelS= 668.39
  bdelS= -1.07
  tmeanK=tmean+273.15
  trefK=tref+273.15
  R=8.314
  kbeg=exp(Ha*(temp-trefK)/(trefK*R*temp))
  kend=((1+exp((trefK*(adelS+bdelS*tmean)-Hd)/(trefK*R)))/(1+exp((temp*(adelS+bdelS*tmean)-Hd)/(temp*R))))
  kbeg*kend
  
}

calc_jmax_tresp_mult = function(tleaf, tmean, tref){
  
  temp = tleaf + 273.15
  Ha= 49884
  Hd= 200000
  adelS= 659.7
  bdelS= -0.75
  tmeanK=tmean+273.15
  trefK=tref+273.15
  R=8.314
  kbeg=exp(Ha*(temp-trefK)/(trefK*R*temp))
  kend=((1+exp((trefK*(adelS+bdelS*tmean)-Hd)/(trefK*R)))/(1+exp((temp*(adelS+bdelS*tmean)-Hd)/(temp*R))))
  kbeg*kend
  
}

# build some response matrices
tg_c_respone = calc_optimal_vcmax(tg_c = seq(15, 35, 1), 
                                  z = 0, 
                                  vpdo = 1, 
                                  cao = 400,
                                  paro = 800)
vpdo_respone = calc_optimal_vcmax(tg_c = 25, 
                                  z = 0, 
                                  vpdo = seq(0.5, 9.5, 0.5), 
                                  cao = 400,
                                  paro = 800)
cao_respone = calc_optimal_vcmax(tg_c = 25, 
                                  z = 0, 
                                  vpdo = 1, 
                                  cao = seq(300, 1000, 50),
                                  paro = 800)
paro_respone = calc_optimal_vcmax(tg_c = 25, 
                                 z = 0, 
                                 vpdo = 1, 
                                 cao = 400,
                                 paro = seq(100, 2000, 50))
z_respone = calc_optimal_vcmax(tg_c = 25, 
                                  z = seq(0, 5000, 500), 
                                  vpdo = 1, 
                                  cao = 400,
                                  paro = 800)

chi_tg_c_plot = ggplot(data = tg_c_respone, aes(y = chi, x = tg_c)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab('Temperature (°C)') +
  ylab('χ') +
  ylim(0.5, 1)

chi_vpdo_plot = ggplot(data = vpdo_respone, aes(y = chi, x = vpdo)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab('VPD (kPa)') +
  ylab('χ') +
  ylim(0.5, 1) +
  xlim(0, 10)

chi_z_plot = ggplot(data = z_respone, aes(y = chi, x = z)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab('Elevation (m)') +
  ylab('χ') +
  ylim(0.5, 1) 
  # xlim(0, 10)

chi_cao_plot = ggplot(data = cao_respone, aes(y = chi, x = cao)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab(expression('CO'[2] * '(Pa)')) +
  ylab('χ') +
  ylim(0.5, 1) 
# xlim(0, 10)

vcmax_tg_c_plot = ggplot(data = tg_c_respone, aes(y = vcmax_prime, x = tg_c)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab('Temperature (°C)') +
  ylab(expression(italic('V')[cmax] * ' (µmol m' ^ '-2' * ' s' ^ '-1' * ')')) +
  ylim(0, 200)

vcmax_vpdo_plot = ggplot(data = vpdo_respone, aes(y = vcmax_prime, x = vpdo)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab('VPD (kPa)') +
  ylab(expression(italic('V')[cmax] * ' (µmol m' ^ '-2' * ' s' ^ '-1' * ')')) +
  ylim(0, 200) +
  xlim(0, 10)

vcmax_z_plot = ggplot(data = z_respone, aes(y = vcmax_prime, x = z)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab('Elevation (m)') +
  ylab(expression(italic('V')[cmax] * ' (µmol m' ^ '-2' * ' s' ^ '-1' * ')')) +
  ylim(0, 200)  
# xlim(0, 10)

vcmax_cao_plot = ggplot(data = cao_respone, aes(y = vcmax_prime, x = cao)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab(expression('CO'[2] * '(Pa)')) +
  ylab(expression(italic('V')[cmax] * ' (µmol m' ^ '-2' * ' s' ^ '-1' * ')')) +
  ylim(0, 200)  
# xlim(0, 10)

vcmax_paro_plot = ggplot(data = paro_respone, aes(y = vcmax_prime, x = paro)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(3), colour = 'black'),
        axis.title.x=element_text(size=rel(3), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(colour = 'darkblue', size = 6, linetype = 1) +
  xlab(expression('PAR' * ' (µmol m' ^ '-2' * ' s' ^ '-1' * ')')) +
  ylab(expression(italic('V')[cmax] * ' (µmol m' ^ '-2' * ' s' ^ '-1' * ')')) +
  ylim(0, 200) +
  xlim(0, 2000)

# future map figure


# get global input data

tmp_globe_4model = read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/cru_tmp_climExtract_growingseason_globe.csv')
par_globe_4model = read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/cru_par_climExtract_growingseason_globe.csv')
vpd_globe_4model = read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/cru_vpd_climExtract_growingseason_globe.csv')
z_globe_4model =  read.csv('/Users/nicksmith/Documents/Research/Colimitation/Spatial_Maps/z_globe.csv')

modis_2001 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2001.tif')
modis_2002 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2002.tif')
modis_2003 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2003.tif')
modis_2004 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2004.tif')
modis_2005 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2005.tif')
modis_2006 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2006.tif')
modis_2007 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2007.tif')
modis_2008 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2008.tif')
modis_2009 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2009.tif')
modis_2010 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2010.tif')
modis_2011 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2011.tif')
modis_2012 = raster('/Users/nicksmith/Documents/Research/Colimitation/landCoverMODIS/LC_hd_global_2012.tif')

modis = overlay(modis_2001, modis_2002, modis_2003, modis_2004, modis_2005, modis_2006, modis_2007, modis_2008, modis_2009, modis_2010, modis_2011, modis_2012, fun = mean)
modis[modis == 16] <- 0 #barren
modis[modis > 0] <- 1 # vegetated

# predict Vcmax at each global site
pred_globe=calc_optimal_vcmax(cao= 400, 
                              tg_c = tmp_globe_4model$tmp, 
                              paro = par_globe_4model$par,
                              vpdo = vpd_globe_4model$vpd, 
                              z= z_globe_4model$z)
pred_globe$lon = tmp_globe_4model$lon
pred_globe$lat = tmp_globe_4model$lat
pred_globe$ag = pred_globe$vcmax_prime * pred_globe$mc
pred_globe$an = pred_globe$ag - (pred_globe$vcmax_prime * 0.015)
pred_globe$vcmax25 = pred_globe$vcmax_prime / calc_tresp_mult(pred_globe$tg_c, pred_globe$tg_c, 25)
pred_globe$narea = pred_globe$vcmax25 / (47.3 * 6.26) # from CLM5 tech note

# predict the future
pred_globe_fut=calc_optimal_vcmax(cao= 1000, 
                              tg_c = tmp_globe_4model$tmp + 4, 
                              paro = par_globe_4model$par,
                              vpdo = vpd_globe_4model$vpd, 
                              z= z_globe_4model$z)
pred_globe_fut$lon = tmp_globe_4model$lon
pred_globe_fut$lat = tmp_globe_4model$lat
pred_globe_fut$ag = pred_globe_fut$vcmax_prime * pred_globe_fut$mc
pred_globe_fut$an = pred_globe_fut$ag - (pred_globe_fut$vcmax_prime * 0.015)
pred_globe_fut$vcmax25 = pred_globe_fut$vcmax_prime / calc_tresp_mult(pred_globe_fut$tg_c, pred_globe_fut$tg_c, 25)
pred_globe_fut$narea = pred_globe_fut$vcmax25 / (47.3 * 6.26) # from CLM5 tech note

# get differences
pred_globe_diff_ag = ((pred_globe_fut$ag - pred_globe$ag) / pred_globe$ag) * 100
pred_globe_diff_ag[pred_globe_diff_ag > 100] <- 100
pred_globe_diff_ag[pred_globe_diff_ag < -100] <- -100
pred_globe_diff_narea = ((pred_globe_fut$narea - pred_globe$narea) / pred_globe$narea) * 100
pred_globe_diff_narea[pred_globe_diff_narea > 100] <- 100
pred_globe_diff_narea[pred_globe_diff_narea < -100] <- -100

# create raster
ag_globe_diff = cbind(pred_globe$lon, pred_globe$lat, pred_globe_diff_ag)
ag_globe_diff_ras = rasterFromXYZ(ag_globe_diff)
ag_globe_diff_ras_mod = modis * ag_globe_diff_ras
narea_globe_diff = cbind(pred_globe$lon, pred_globe$lat, pred_globe_diff_narea)
narea_globe_diff_ras = rasterFromXYZ(narea_globe_diff)
narea_globe_diff_ras_mod = modis * narea_globe_diff_ras

# plot
palette = colorRampPalette(c(rev(brewer.pal(10,'RdBu'))[1:5], 'white', 
                             'white', rev(brewer.pal(10,'RdBu'))[6:10]))
cols = palette(10)
arg = list(at = seq(-100, 100, 20), labels = seq(-100, 100, 20))

par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
plot(ag_globe_diff_ras_mod, col=cols, breaks=seq(-100, 100, 20), cex.axis=1.5, 
     yaxt = 'n', xaxt = 'n', lab.breaks = seq(-100, 100, 20), 
     ylim = c(-90, 90), 
     legend.args=list(text=expression('∆Photosyntesis'*' (%)'), line = 4, side = 4, las = 3, cex = 1.5), 
     axis.args = arg)
map('world',col=c('black'),fill=F, add = T)
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)
title(expression('+600 ppm CO'[2] * '; +4°C'), line = -2, cex.main = 2)

par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
plot(narea_globe_diff_ras_mod, col=cols, breaks=seq(-100, 100, 20), cex.axis=1.5, 
     yaxt = 'n', xaxt = 'n', lab.breaks = seq(-100, 100, 20), 
     ylim = c(-90, 90), 
     legend.args=list(text=expression('∆Leaf N'*' (%)'), line = 4, side = 4, las = 3, cex = 1.5), 
     axis.args = arg)
map('world',col=c('black'),fill=F, add = T)
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)
title(expression('+600 ppm CO'[2] * '; +4°C'), line = -2, cex.main = 2)





pale = colorRampPalette(c('white', rev(brewer.pal(10,'Spectral'))))
cols = pale(9)
arg = list(at = seq(0, 45, 5), labels = seq(0, 45, 5))

par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
plot(ag_globe_ras_mod, col=cols, breaks = seq(0, 45, 5), 
     cex.axis=1.5, yaxt = 'n', xaxt = 'n', lab.breaks = seq(0, 45, 5), 
     ylim = c(-90, 90), 
     legend.args=list(text=expression('Photosynthesis'*' (µmol m'^'-2'*' s'^'-1'*')'), 
                      line = 4, side = 4, las = 3, cex = 1.5), 
     legend = T, xlim = c(-180, 180), axis.args = arg)
map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)

# effect of future conditions
predVcmaxJF_globe_present= calc_vcmax_jf_dAjdJmax(ca= 400, tmean = tmp_globe_4model$tmp, par = par_globe_4model$par,vpdo = vpd_globe_4model$vpd, z= z_globe_4model$z, tleaf=tmp_globe_4model$tmp, q0= q0o * abs * ext, theta = theta)
predVcmaxJF_globe_present $An = (predVcmaxJF_globe_present $Vcmax_tmean * predVcmaxJF_globe_present $mc) - (0.015 * predVcmaxJF_globe_present $Vcmax_tmean)
predVcmaxJF_globe_co2= calc_vcmax_jf_dAjdJmax(ca= 1000, tmean = tmp_globe_4model$tmp, par = par_globe_4model$par,vpdo = vpd_globe_4model$vpd, z= z_globe_4model$z, tleaf=tmp_globe_4model$tmp, q0= q0o * abs * ext, theta = theta)
predVcmaxJF_globe_co2 $An = (predVcmaxJF_globe_co2 $Vcmax_tmean * predVcmaxJF_globe_co2 $mc) - (0.015 * predVcmaxJF_globe_co2 $Vcmax_tmean)
predVcmaxJF_globe_temp= calc_vcmax_jf_dAjdJmax(ca= 400, tmean = tmp_globe_4model$tmp+4, par = par_globe_4model$par,vpdo = vpd_globe_4model$vpd, z= z_globe_4model$z, tleaf=tmp_globe_4model$tmp+4, q0= q0o * abs * ext, theta = theta)
predVcmaxJF_globe_temp $An = (predVcmaxJF_globe_temp $Vcmax_tmean * predVcmaxJF_globe_temp $mc) - (0.015 * predVcmaxJF_globe_temp $Vcmax_tmean)
predVcmaxJF_globe_co2_temp= calc_vcmax_jf_dAjdJmax(ca= 1000, tmean = tmp_globe_4model$tmp+4, par = par_globe_4model$par,vpdo = vpd_globe_4model$vpd, z= z_globe_4model$z, tleaf=tmp_globe_4model$tmp+4, q0= q0o * abs * ext, theta = theta)
predVcmaxJF_globe_co2_temp$An = (predVcmaxJF_globe_co2_temp$Vcmax_tmean * predVcmaxJF_globe_co2_temp$mc) - (0.015 * predVcmaxJF_globe_co2_temp$Vcmax_tmean)

# difference
dif_co2 = ((predVcmaxJF_globe_co2-predVcmaxJF_globe_present)/predVcmaxJF_globe_present)*100
dif_co2_4ras = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, dif_co2$Vcmax_tmean)
dif_co2_ras = rasterFromXYZ(dif_co2_4ras)
dif_co2_veg = modis * dif_co2_ras
dif_temp = ((predVcmaxJF_globe_temp-predVcmaxJF_globe_present)/predVcmaxJF_globe_present)*100
dif_temp_4ras = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, dif_temp$Vcmax_tmean)
dif_temp_ras = rasterFromXYZ(dif_temp_4ras)
dif_temp_veg = modis * dif_temp_ras

dif_co2_temp = ((predVcmaxJF_globe_co2_temp-predVcmaxJF_globe_present)/predVcmaxJF_globe_present)*100
dif_co2_temp_4ras = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, dif_co2_temp$Vcmax_tmean)
dif_co2_temp_ras = rasterFromXYZ(dif_co2_temp_4ras)
dif_co2_temp_veg = modis * dif_co2_temp_ras

dif_co2_temp_4ras_an  = cbind(tmp_globe_4model$lon, tmp_globe_4model$lat, dif_co2_temp$An)
dif_co2_temp_ras_an = rasterFromXYZ(dif_co2_temp_4ras_an)
dif_co2_temp_veg_an = modis * dif_co2_temp_ras_an

pale = colorRampPalette(c(rev(brewer.pal(10,'RdBu'))[1:5], 'white', 'white', rev(brewer.pal(10,'RdBu'))[6:10]))
cols = pale(20)
arg = list(at = c(-50, -25, 0, 25, 50), labels = seq(-50,50, 25))

par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
plot(dif_co2_veg, col=cols, breaks=seq(-50,50, 5), cex.axis=1.5, yaxt = 'n', xaxt = 'n', lab.breaks = seq(-50, 50, 5), ylim = c(-90, 90), legend.args=list(text=expression(italic('∆V'[cmax])*' (%)'), line = 4, side = 4, las = 3, cex = 1.5), axis.args = arg)
map('world',col=c('black'),fill=F, add = T)
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)
title('+600 ppm [CO2]', line = -2, cex.main = 2)

par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
plot(dif_temp_veg, col=cols, breaks=seq(-50,50, 5), cex.axis=1.5, yaxt = 'n', xaxt = 'n', lab.breaks = seq(-50, 50, 5), ylim = c(-90, 90), legend.args=list(text=expression(italic('∆V'[cmax])*' (%)'), line = 4, side = 4, las = 3, cex = 1.5), axis.args = arg)
map('world',col=c('black'),fill=F, add = T)
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)
title('+4°C', line = -2, cex.main = 2)

palette = colorRampPalette(c(rev(brewer.pal(10,'RdBu'))[1:5], 'white', 'white', rev(brewer.pal(10,'RdBu'))[6:10]))
cols = palette(20)
arg = list(at = seq(-100, 100, 50), labels = seq(-100, 100, 50))

par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
plot(dif_co2_temp_veg, col=cols, breaks=seq(-100, 100, 10), cex.axis=1.5, yaxt = 'n', xaxt = 'n', lab.breaks = seq(-100, 100, 10), ylim = c(-90, 90), legend.args=list(text=expression(italic('∆V'[cmax])*' (%)'), line = 4, side = 4, las = 3, cex = 1.5), axis.args = arg)
map('world',col=c('black'),fill=F, add = T)
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)
title('elevated [CO2] & Temperature', line = -2, cex.main = 2)

par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
plot(dif_co2_temp_veg_an, col=cols, breaks=seq(-100, 100, 10), cex.axis=1.5, yaxt = 'n', xaxt = 'n', lab.breaks = seq(-100, 100, 10), ylim = c(-90, 90), legend.args=list(text=expression(italic('∆A'[net])*' (%)'), line = 4, side = 4, las = 3, cex = 1.5), axis.args = arg)
map('world',col=c('black'),fill=F, add = T)
axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)
title('elevated [CO2] & Temperature', line = -2, cex.main = 2)








