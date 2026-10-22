# Clip the 1 km SRTM DEM of China to the study box (109.4-111.6 E, 31.2-33.4 N) and write
# data/dem_ShiYan_1km.tif, which scripts/plot_study_area_map.py draws as paper Fig. 1's relief.
#
#   Rscript scripts/clip_dem.R      # from the repository root
#
# The input `f` is a local path on the machine that first ran this; point it at your copy of
# dem_china_1km_SRTM.tif. `st_rect` comes from the sf2 package. That nationwide DEM is not in this
# project; to restore the clipped output instead, run scripts/download_dem_shiyan.jl.
# Tracked copy of data/clip_dem.R, which sits in the gitignored data/ folder.

pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  terra, sf2
)

f = "X:/rpkgs/VICTools.R/inst/database/dem_china_1km_SRTM.tif"
ra = rast(f)

poly = st_rect(c(109.4, 111.6, 31.2, 33.4))
ra2 = crop(ra, poly)

writeRaster(ra2, "./data/dem_ShiYan_1km.tif")
