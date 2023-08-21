#install.packages("devtools") 
library(devtools)
#install_github("sdellicour/seraphim/unix_OS") # for Unix systems
#install_github("sdellicour/seraphim/windows") # for Windows systems
library(seraphim)
#install.packages("diagram")
library(diagram)
library(rgdal)
library(cowplot)
library(treeio)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)


# 1. Extracting spatio-temporal information embedded in posterior trees

localTreesDirectory = "Run_09052023_v1/Run_09052023_v1/Tree_extractions"
allTrees = scan(file="Run_09052023_v1/Run_09052023_v1/Selected_sequences_v1_newcoord.trees", what="", sep="\n", quiet=T)

localTreesDirectory = "Run_09052023_m1/Run_09052023_m1/Tree_extractions"
allTrees = scan(file="Run_09052023_m1/Run_09052023_m1/Selected_sequences_m1_newcoord.trees", what="", sep="\n", quiet=T)

localTreesDirectory = "Run_09052023_ns1/Run_09052023_ns1/Tree_extractions"
allTrees = scan(file="Run_09052023_ns1/Run_09052023_ns1/Selected_sequences_ns1_newcoord.trees", what="", sep="\n", quiet=T)

burnIn = 1000
randomSampling = FALSE
nberOfTreesToSample = 1000
mostRecentSamplingDatum = decimal_date(ymd("2022-02-07"))
coordinateAttributeName = "location"

treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)



# 2. Extracting spatio-temporal information embedded in the MCC tree

source("Script/mccExtractions.r")
mcc_tre = readAnnotatedNexus("Run_09052023_v1/Run_09052023_v1/Selected_sequences_v1_newcoord_MCC.trees")
mcc_data = as.treedata(mcc_tre)
mcc_tab = mccExtractions(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, "Run_09052023_v1/Run_09052023_v1/Selected_sequences_v1_newcoord_MCC.csv", row.names=F, quote=F)

source("Script/mccExtractions.r")
mcc_tre = readAnnotatedNexus("Run_09052023_m1/Run_09052023_m1/Selected_sequences_m1_newcoord_MCC.trees")
mcc_data = as.treedata(mcc_tre)
mcc_tab = mccExtractions(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, "Run_09052023_m1/Run_09052023_m1/Selected_sequences_m1_newcoord_MCC.csv", row.names=F, quote=F)

source("Script/mccExtractions.r")
mcc_tre = readAnnotatedNexus("Run_09052023_ns1/Run_09052023_ns1/Selected_sequences_ns1_newcoord_MCC.trees")
mcc_data = as.treedata(mcc_tre)
mcc_tab = mccExtractions(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, "Run_09052023_ns1/Run_09052023_ns1/Selected_sequences_ns1_newcoord_MCC.csv", row.names=F, quote=F)



# 3. Estimating the HPD region for each time slice

nberOfExtractionFiles = nberOfTreesToSample
prob = 0.95; precision = 0.025
startDatum = min(mcc_tab[,"startYear"])

polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))

### New version
xMin = 9999; xMax = 0; yMin = 9999; yMax = 0; minYears = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
{
  tab = read.csv(paste0("Run_09052023_v1/Run_09052023_v1/Tree_extractions/TreeExtractions_",i,".csv"), head=T)
  if (xMin > min(tab[,"endLon"])) xMin = min(tab[,"endLon"])
  if (xMax < max(tab[,"endLon"])) xMax = max(tab[,"endLon"])
  if (yMin > min(tab[,"endLat"])) yMin = min(tab[,"endLat"])
  if (yMax < max(tab[,"endLat"])) yMax = max(tab[,"endLat"])
  minYears[i] = min(tab[,"startYear"])
}
localTreesDirectory = "Run_09052023_v1/Run_09052023_v1/Tree_extractions"
mcc = read.csv("Run_09052023_v1/Run_09052023_v1/Selected_sequences_v1_newcoord_MCC.csv", head=T)  
mostRecentSamplingDatum = max(mcc[,"endYear"])
minYear = startDatum = quantile(minYears, c(0.025,0.975))[1] # min(mcc[,"startYear"])
maxYear = mostRecentSamplingDatum; mcc = mcc[order(mcc[,"endYear"],mcc[,"startYear"]),]
prob = 0.80; precision = 2; startDatum = quantile(minYears, c(0.025,0.975))[1] # min(mcc[,"startYear"])
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


xMin = 9999; xMax = 0; yMin = 9999; yMax = 0; minYears = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
{
  tab = read.csv(paste0("Run_09052023_m1/Run_09052023_m1/Tree_extractions/TreeExtractions_",i,".csv"), head=T)
  if (xMin > min(tab[,"endLon"])) xMin = min(tab[,"endLon"])
  if (xMax < max(tab[,"endLon"])) xMax = max(tab[,"endLon"])
  if (yMin > min(tab[,"endLat"])) yMin = min(tab[,"endLat"])
  if (yMax < max(tab[,"endLat"])) yMax = max(tab[,"endLat"])
  minYears[i] = min(tab[,"startYear"])
}
localTreesDirectory = "Run_09052023_m1/Run_09052023_m1/Tree_extractions"
mcc = read.csv("Run_09052023_m1/Run_09052023_m1/Selected_sequences_m1_newcoord_MCC.csv", head=T)  
mostRecentSamplingDatum = max(mcc[,"endYear"])
minYear = startDatum = quantile(minYears, c(0.025,0.975))[1] # min(mcc[,"startYear"])
maxYear = mostRecentSamplingDatum; mcc = mcc[order(mcc[,"endYear"],mcc[,"startYear"]),]
prob = 0.80; precision = 2; startDatum = quantile(minYears, c(0.025,0.975))[1] # min(mcc[,"startYear"])
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


xMin = 9999; xMax = 0; yMin = 9999; yMax = 0; minYears = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
{
  tab = read.csv(paste0("Run_09052023_ns1/Run_09052023_ns1/Tree_extractions/TreeExtractions_",i,".csv"), head=T)
  if (xMin > min(tab[,"endLon"])) xMin = min(tab[,"endLon"])
  if (xMax < max(tab[,"endLon"])) xMax = max(tab[,"endLon"])
  if (yMin > min(tab[,"endLat"])) yMin = min(tab[,"endLat"])
  if (yMax < max(tab[,"endLat"])) yMax = max(tab[,"endLat"])
  minYears[i] = min(tab[,"startYear"])
}
localTreesDirectory = "Run_09052023_ns1/Run_09052023_ns1/Tree_extractions"
mcc = read.csv("Run_09052023_ns1/Run_09052023_ns1/Selected_sequences_ns1_newcoord_MCC.csv", head=T)  
mostRecentSamplingDatum = max(mcc[,"endYear"])
minYear = startDatum = quantile(minYears, c(0.025,0.975))[1] # min(mcc[,"startYear"])
maxYear = mostRecentSamplingDatum; mcc = mcc[order(mcc[,"endYear"],mcc[,"startYear"]),]
prob = 0.80; precision = 2; startDatum = quantile(minYears, c(0.025,0.975))[1] # min(mcc[,"startYear"])
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


# 4. Defining the different colour scales to use

colour_scale = colorRampPalette(brewer.pal(11,"RdYlGn"))(141)[21:121]
minYear = min(mcc_tab[,"startYear"]); maxYear = max(mcc_tab[,"endYear"])
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
		polygons_colours[i] = paste0(colour_scale[polygon_index],"40")
	}

### New version 
cols = gsub("FF","",rev(viridis::magma(241))[1:201])
startYear_index = (((min(mcc[,"startYear"])-minYear)/(maxYear-minYear))*100)+1
startYear_colour = cols[startYear_index]
endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_indices[which(endYears_indices<1)] = NA
endYears_colours = cols[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
{
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  if (polygon_index >= 1)
  {
    polygons_colours[i] = paste0(cols[polygon_index],"20")
  }
}


# 5. Co-plotting the HPD regions and MCC tree

template_shp <-readOGR("Data/Vector/countries.shp") 
template_shp <- subset(template_shp, template_shp@data$UNREG2 == "Africa")
#template_raster = raster("YFV_studyArea.asc") #open shapefile instead of raster
#borders = crop(getData("GADM", country="BRA", level=1), extent(template_raster))
tiff("Results/Figure 4_phylomap_08.tiff", width = 20, height = 20, units = "cm",
     compression = "lzw", res = 300)
#dev.new(width=6, height=6.3)
par(mar=c(0,0,0,0)) #, oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
#plot(template_raster, col="white", box=F, axes=F, colNA="grey90", legend=F)
plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "black", lwd=0.5)
#plot(density_poultry, add=T)
for (i in 1:length(polygons))
	{
    plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
	}
#plot(borders, add=T, lwd=0.1, border="gray10")
for (i in 1:dim(mcc_tab)[1])
	{
		curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
	}
for (i in dim(mcc_tab)[1]:1)
	{
		if (i == 1) {
				points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=16, col=colour_scale[1], cex=0.8)
				points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=1, col="gray10", lwd=0.2, cex=0.8)
		}
  if (mcc_tab[i, "node2"]%in%mcc_tab[, "node1"]) 
    { # for internal nodes (dots)
    points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
    points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=1, col="gray10", lwd=0.4, cex=0.8)
  } 
  else 
    { # for the tip nodes (squares)
    points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
    points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=0, col="gray10", lwd=0.4, cex=0.8)
  }
}
#rect(xmin(template_raster), ymin(template_raster), xmax(template_raster), ymax(template_raster), xpd=T, lwd=0.2)
#axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
#axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tab[,"startYear"]); rast[2] = max(mcc_tab[,"endYear"])
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.45,0.1,0.11),
     legend.args=list(text="", cex=1, line=0, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.8, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0)) #, mgp=c(0,-0.02,0), at=seq(2016.4,2017.2,0.2)))
dev.off()


### New version 
template_shp <-readOGR("Data/Vector/countries.shp") 
template_shp <- subset(template_shp, template_shp@data$UNREG2 == "Africa")
#template_raster = raster("YFV_studyArea.asc") #open shapefile instead of raster
#borders = crop(getData("GADM", country="BRA", level=1), extent(template_raster))

tiff("Results/Figure 4_phylomap_v1.tiff", width = 20, height = 20, units = "cm",
     compression = "lzw", res = 300)
#dev.new(width=6, height=6.3)
par(mar=c(0,0,0,0)) #, oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
#plot(template_raster, col="white", box=F, axes=F, colNA="grey90", legend=F)
plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "gray10", lwd=0.5)
#plot(density_poultry, add=T)

for (i in 1:length(polygons))
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
#plot(borders, add=T, lwd=0.1, border="gray10")
for (i in 1:dim(mcc)[1])
{
  curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.4, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
}
for (i in dim(mcc)[1]:1)
{
  if (i == 1) {
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYear_colour, cex=0.8)
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  }
  if (mcc[i, "node2"]%in%mcc[, "node1"]) 
  { # for internal nodes (dots)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  } 
  else 
  { # for the tip nodes (squares)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="black", lwd=0.4, cex=0.8)
  }
}
#rect(xmin(template_raster), ymin(template_raster), xmax(template_raster), ymax(template_raster), xpd=T, lwd=0.2)
#axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
#axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
#rast_breaks <- seq(floor(startDatum/5)*5, ceiling(max(mcc[,"endYear"])/5)*5, by=5)
plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.45,0.1,0.11),
     legend.args=list(text="", cex=0.5, line=0, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-1, col.axis="gray30", line=0))
                    #at=rast_breaks, labels=rast_breaks)) #, mgp=c(0,-0.02,0), at=seq(2016.4,2017.2,0.2)))
dev.off()  ### How can I add more dates on the legend scale? 


tiff("Results/Figure 4_phylomap_m1.tiff", width = 20, height = 20, units = "cm",
     compression = "lzw", res = 300)
#dev.new(width=6, height=6.3)
par(mar=c(0,0,0,0)) #, oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
#plot(template_raster, col="white", box=F, axes=F, colNA="grey90", legend=F)
plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "gray10", lwd=0.5)
#plot(density_poultry, add=T)
title("Reduced", adj = 0.05, line = -1, font.main=1)

for (i in 1:length(polygons))
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
#plot(borders, add=T, lwd=0.1, border="gray10")
for (i in 1:dim(mcc)[1])
{
  curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.4, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
}
for (i in dim(mcc)[1]:1)
{
  if (i == 1) {
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYear_colour, cex=0.8)
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  }
  if (mcc[i, "node2"]%in%mcc[, "node1"]) 
  { # for internal nodes (dots)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  } 
  else 
  { # for the tip nodes (squares)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="black", lwd=0.4, cex=0.8)
  }
}
#rect(xmin(template_raster), ymin(template_raster), xmax(template_raster), ymax(template_raster), xpd=T, lwd=0.2)
#axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
#axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
#rast_breaks <- seq(floor(startDatum/5)*5, ceiling(max(mcc[,"endYear"])/5)*5, by=5)
plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.45,0.1,0.11),
     legend.args=list(text="", cex=0.5, line=0, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-1, col.axis="gray30", line=0))
#at=rast_breaks, labels=rast_breaks)) #, mgp=c(0,-0.02,0), at=seq(2016.4,2017.2,0.2)))
dev.off()  ### How can I add more dates on the legend scale? 


tiff("Results/Figure 4_phylomap_ns1.tiff", width = 20, height = 20, units = "cm",
     compression = "lzw", res = 300)
#dev.new(width=6, height=6.3)
par(mar=c(0,0,0,0)) #, oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
#plot(template_raster, col="white", box=F, axes=F, colNA="grey90", legend=F)
plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "gray10", lwd=0.5)
#plot(density_poultry, add=T)
title("Equivalent", adj = 0.05, line = -1, font.main=1)

for (i in 1:length(polygons))
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
#plot(borders, add=T, lwd=0.1, border="gray10")
for (i in 1:dim(mcc)[1])
{
  curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.4, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
}
for (i in dim(mcc)[1]:1)
{
  if (i == 1) {
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYear_colour, cex=0.8)
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  }
  if (mcc[i, "node2"]%in%mcc[, "node1"]) 
  { # for internal nodes (dots)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  } 
  else 
  { # for the tip nodes (squares)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="black", lwd=0.4, cex=0.8)
  }
}
#rect(xmin(template_raster), ymin(template_raster), xmax(template_raster), ymax(template_raster), xpd=T, lwd=0.2)
#axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
#axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
#rast_breaks <- seq(floor(startDatum/5)*5, ceiling(max(mcc[,"endYear"])/5)*5, by=5)
plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.45,0.1,0.11),
     legend.args=list(text="", cex=0.5, line=0, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-1, col.axis="gray30", line=0))
#at=rast_breaks, labels=rast_breaks)) #, mgp=c(0,-0.02,0), at=seq(2016.4,2017.2,0.2)))
dev.off()  ### How can I add more dates on the legend scale? 



# 6. Analysing the long-distance dispersal event

latitude_cutoff = 21.3 # roughly corresponding to the southern border of Morocco
trees_to_check = c() # trees with more than one long-distance dispersal events #### Why looking at those with more than 1 long-dispersal event? 
for (i in 1:nberOfExtractionFiles)
{
  print(i)
  tab = read.csv(paste0("Run_09052023_v1/Run_09052023_v1/Tree_extractions/TreeExtractions_",i,".csv"), head=T)
  indices = which((tab[,"startLat"]>latitude_cutoff)&(tab[,"endLat"]<latitude_cutoff))
  if (length(indices) != 1)
  {
    print(c(i,indices)) # 349 (57,58), 528 (13,14), 979 (15,16)
    trees_to_check = c(trees_to_check, i)
  }
} #0 trees out of 1000

latitude_cutoff = 21.3 # roughly corresponding to the southern border of Morocco
trees_to_check = c() # trees with more than one long-distance dispersal events #### Why looking at those with more than 1 long-dispersal event? 
for (i in 1:nberOfExtractionFiles)
{
  print(i)
  tab = read.csv(paste0("Run_09052023_m1/Run_09052023_m1/Tree_extractions/TreeExtractions_",i,".csv"), head=T)
  indices = which((tab[,"startLat"]>latitude_cutoff)&(tab[,"endLat"]<latitude_cutoff))
  if (length(indices) != 1)
  {
    print(c(i,indices)) # 349 (57,58), 528 (13,14), 979 (15,16)
    trees_to_check = c(trees_to_check, i)
  }
} #0 trees out of 1000

latitude_cutoff = 21.3 # roughly corresponding to the southern border of Morocco
trees_to_check = c() # trees with more than one long-distance dispersal events #### Why looking at those with more than 1 long-dispersal event? 
for (i in 1:nberOfExtractionFiles)
{
  print(i)
  tab = read.csv(paste0("Run_09052023_ns1/Run_09052023_ns1/Tree_extractions/TreeExtractions_",i,".csv"), head=T)
  indices = which((tab[,"startLat"]>latitude_cutoff)&(tab[,"endLat"]<latitude_cutoff))
  if (length(indices) != 1)
  {
    print(c(i,indices)) # 349 (57,58), 528 (13,14), 979 (15,16)
    trees_to_check = c(trees_to_check, i)
  }
} #15 trees out of 1000

# Plot trees with more than one long-distance dispersal events
for (h in 1:length(trees_to_check))
{
  tre = read.csv(paste0("Run_09052023_v1/Run_09052023_v1/Tree_extractions/TreeExtractions_",trees_to_check[h],".csv"), head=T)
  startYear_index = (((min(tre[,"startYear"])-minYear)/(maxYear-minYear))*100)+1
  startYear_colour = cols[startYear_index]
  endYears_indices = (((tre[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
  endYears_indices[which(endYears_indices<1)] = NA
  endYears_colours = cols[endYears_indices]
  tiff(paste0("Results/Figure 4_phylomap_v1_", trees_to_check[h],".tiff"), width = 20, height = 20, units = "cm",
       compression = "lzw", res = 300)
  par(mar=c(0,0,0,0)) 
  plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "gray10", lwd=0.5)
  for (i in 1:length(polygons))
  {
    plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
  }
  for (i in 1:dim(tre)[1])
  {
    curvedarrow(cbind(tre[i,"startLon"],tre[i,"startLat"]), cbind(tre[i,"endLon"],tre[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=0.4, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
  }
  for (i in dim(tre)[1]:1)
  {
    if (i == 1) {
      points(tre[i,"startLon"], tre[i,"startLat"], pch=16, col=startYear_colour, cex=0.8)
      points(tre[i,"startLon"], tre[i,"startLat"], pch=1, col="black", lwd=0.4, cex=0.8)
    }
    if (tre[i, "node2"]%in%tre[, "node1"]) 
    { # for internal nodes (dots)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=1, col="black", lwd=0.4, cex=0.8)
    } 
    else 
    { # for the tip nodes (squares)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=0, col="black", lwd=0.4, cex=0.8)
    }
  }
  rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(tre[,"endYear"])
  plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.45,0.1,0.11),
       legend.args=list(text="", cex=0.5, line=0, col="gray30"), horizontal=T,
       axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-1, col.axis="gray30", line=0))
  dev.off()
}

for (h in 1:length(trees_to_check))
{
  tre = read.csv(paste0("Run_09052023_m1/Run_09052023_m1/Tree_extractions/TreeExtractions_",trees_to_check[h],".csv"), head=T)
  startYear_index = (((min(tre[,"startYear"])-minYear)/(maxYear-minYear))*100)+1
  startYear_colour = cols[startYear_index]
  endYears_indices = (((tre[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
  endYears_indices[which(endYears_indices<1)] = NA
  endYears_colours = cols[endYears_indices]
  tiff(paste0("Results/Figure 4_phylomap_m1_", trees_to_check[h],".tiff"), width = 20, height = 20, units = "cm",
       compression = "lzw", res = 300)
  par(mar=c(0,0,0,0)) 
  plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "gray10", lwd=0.5)
  for (i in 1:length(polygons))
  {
    plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
  }
  for (i in 1:dim(tre)[1])
  {
    curvedarrow(cbind(tre[i,"startLon"],tre[i,"startLat"]), cbind(tre[i,"endLon"],tre[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=0.4, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
  }
  for (i in dim(tre)[1]:1)
  {
    if (i == 1) {
      points(tre[i,"startLon"], tre[i,"startLat"], pch=16, col=startYear_colour, cex=0.8)
      points(tre[i,"startLon"], tre[i,"startLat"], pch=1, col="black", lwd=0.4, cex=0.8)
    }
    if (tre[i, "node2"]%in%tre[, "node1"]) 
    { # for internal nodes (dots)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=1, col="black", lwd=0.4, cex=0.8)
    } 
    else 
    { # for the tip nodes (squares)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=0, col="black", lwd=0.4, cex=0.8)
    }
  }
  rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(tre[,"endYear"])
  plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.45,0.1,0.11),
       legend.args=list(text="", cex=0.5, line=0, col="gray30"), horizontal=T,
       axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-1, col.axis="gray30", line=0))
  dev.off()
}

for (h in 1:length(trees_to_check)) #15 trees out of 1000
{
  tre = read.csv(paste0("Run_09052023_ns1/Run_09052023_ns1/Tree_extractions/TreeExtractions_",trees_to_check[h],".csv"), head=T)
  startYear_index = (((min(tre[,"startYear"])-minYear)/(maxYear-minYear))*100)+1
  startYear_colour = cols[startYear_index]
  endYears_indices = (((tre[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
  endYears_indices[which(endYears_indices<1)] = NA
  endYears_colours = cols[endYears_indices]
  tiff(paste0("Results/Figure 4_phylomap_ns1_", trees_to_check[h],".tiff"), width = 20, height = 20, units = "cm",
       compression = "lzw", res = 300)
  par(mar=c(0,0,0,0)) 
  plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "gray10", lwd=0.5)
  for (i in 1:length(polygons))
  {
    plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
  }
  for (i in 1:dim(tre)[1])
  {
    curvedarrow(cbind(tre[i,"startLon"],tre[i,"startLat"]), cbind(tre[i,"endLon"],tre[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=0.4, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
  }
  for (i in dim(tre)[1]:1)
  {
    if (i == 1) {
      points(tre[i,"startLon"], tre[i,"startLat"], pch=16, col=startYear_colour, cex=0.8)
      points(tre[i,"startLon"], tre[i,"startLat"], pch=1, col="black", lwd=0.4, cex=0.8)
    }
    if (tre[i, "node2"]%in%tre[, "node1"]) 
    { # for internal nodes (dots)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=1, col="black", lwd=0.4, cex=0.8)
    } 
    else 
    { # for the tip nodes (squares)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
      points(tre[i,"endLon"], tre[i,"endLat"], pch=0, col="black", lwd=0.4, cex=0.8)
    }
  }
  rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(tre[,"endYear"])
  plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.45,0.1,0.11),
       legend.args=list(text="", cex=0.5, line=0, col="gray30"), horizontal=T,
       axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-1, col.axis="gray30", line=0))
  dev.off()
}


# Look at introduction times for trees which have only one long distance dispersal event
introduction_times_in_western_Africa = rep(NA, length(nberOfExtractionFiles))
trees_one_dispersal_event = c() # trees with one long-distance dispersal event
for (i in 1:nberOfExtractionFiles)
{
  if (!i%in%trees_to_check)  ### Does that mean you don't take the trees that have more than one long dipersal event?
  {
    tab = read.csv(paste0("Run_09052023_v1/Run_09052023_v1/Tree_extractions/TreeExtractions_",i,".csv"), head=T)
    indices = which((tab[,"startLat"]>latitude_cutoff)&(tab[,"endLat"]<latitude_cutoff))
    if (length(indices) == 1)
    {
      introduction_times_in_western_Africa[i] = (tab[indices,"startYear"]+tab[indices,"endYear"])/2  ### So you take the middle date?
      print(c(i,indices))
      trees_one_dispersal_event = c(trees_one_dispersal_event, i)
    }
  }
}
introduction_times_in_western_Africa = introduction_times_in_western_Africa[!is.na(introduction_times_in_western_Africa)]
summary(introduction_times_in_western_Africa)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2012    2015    2015    2015    2016    2017


# Plot introduction times on the map 
template_shp <-readOGR("Data/Vector/countries.shp") 
template_shp <- subset(template_shp, template_shp@data$UNREG2 == "Africa")
#template_raster = raster("YFV_studyArea.asc") #open shapefile instead of raster
#borders = crop(getData("GADM", country="BRA", level=1), extent(template_raster))

tiff("Results/Figure 4_phylomap_v1_time.tiff", width = 20, height = 20, units = "cm",
     compression = "lzw", res = 300)
par(mar=c(0,0,0,0)) 
plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "gray10", lwd=0.5)
for (i in 1:length(polygons))
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
for (i in 1:dim(mcc)[1])
{
  curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.4, lty=1, lcol="black", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
}
for (i in dim(mcc)[1]:1)
{
  if (i == 1) {
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYear_colour, cex=0.8)
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  }
  if (mcc[i, "node2"]%in%mcc[, "node1"]) 
  { # for internal nodes (dots)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="black", lwd=0.4, cex=0.8)
  } 
  else 
  { # for the tip nodes (squares)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
    points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="black", lwd=0.4, cex=0.8)
  }
}
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.60,0.05,0.06),
     legend.args=list(text="", cex=0.5, line=0, col="gray30"), horizontal=T,
     axis.args=list(at=seq(1985, 2020, by=5), cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-0.6, col.axis="gray30", line=0, mgp=c(3, 0.5, 0)))
     #axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.5, tck=-1, col.axis="gray30", line=0))

par(fig = c(0.05, 0.28, 0.78, 1), mar=c(0,0,0,0), new=TRUE) #xlim=c(-20, 10), ylim=c(0, 37.5)
xLim = c(2011,2017.2); yLim=c(0,0.9) #; X = 280000; Y = 0.0001
cols1 = list(); cols1[[1]] = rgb(204,0,0,255,maxColorValue=255); cols1[[2]] = rgb(120,120,120,255,maxColorValue=255)
cols2 = list(); cols2[[1]] = rgb(204,0,0,100,maxColorValue=255); cols2[[2]] = rgb(120,120,120,100,maxColorValue=255)
plot(density(introduction_times_in_western_Africa), lwd=0.7, axes=F, ann=F, col=cols1[[1]], xlim=xLim, ylim=yLim) 
polygon(density(introduction_times_in_western_Africa), col=cols2[[1]], border=NA)
axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0.3, mgp=c(0,0.25,0), tck=-0.025, col.tick="gray10", col.axis="gray10", col="gray10", at=c(2010:2018)) #, 
axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0.3, mgp=c(0,-0.02,0), tck=-0.025, col.tick="gray10", col.axis="gray10", col="gray10") # 
mtext("Date of dispersal event to West Africa", side=1, line=1, cex=0.7, col="gray10")
mtext("Density", side=2, line=1, cex=0.7, col="gray10")
# title(xlab="Date", cex.lab=0.7, mgp=c(0,0.26,0), col.lab="gray10") #
# title(ylab="Density", cex.lab=0.7, mgp=c(0,-0.03,0), col.lab="gray10") # 
dev.off() 



# 7. Extracting statistics

nberOfExtractionFiles = 100
timeSlices = 100
onlyTipBranches = FALSE
showingPlots = TRUE
outputName = "stat_v1"
nberOfCores = 1
slidingWindow = 1
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices,
                   onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)


wavefront <- read.csv("Run_09052023_v1/Run_09052023_v1/Statistics/stat_v1_spatial_wavefront_distances.csv", sep=";") 
wavefront_summary <- wavefront %>%
  pivot_longer(cols=starts_with("distance"),
               names_to="trees", values_to="value") %>%
  group_by(time) %>%
  summarise(l95_value = HDInterval::hdi(value)[[1]],
            h95_value = HDInterval::hdi(value)[[2]],
            median_value = median(value))
wavefront_summary$time <- as.Date(date_decimal(wavefront_summary$time))

wavefront_plot <- ggplot(wavefront_summary) + 
  geom_ribbon(aes(time, ymin = l95_value, ymax = h95_value), fill="grey90") +
  geom_step(aes(time, median_value), color="grey10") +
  labs(x="",y="Spatial wavefront distance (km)") +
  #ylim(0,5) +
  theme_bw() +
  scale_x_date(limits = c(ymd("1977-01-01"), ymd("2022-02-07")), breaks= "2 year", date_labels = "%Y") +
  theme(panel.grid = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", linewidth = 0.2, linetype = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size=11, angle = 45, hjust = 1)) 

velocity_summary <- read.csv("Run_09052023_v1/Run_09052023_v1/Statistics/stat_v1_weighted_branch_dispersal_velocity.csv", sep=";") 
velocity_summary$time <- as.Date(date_decimal(velocity_summary$time))
  
vel_summary <- velocity_summary %>% 
  mutate(year = format(time, "%Y")) %>%
  group_by(year) %>%
  filter(year %in% c("2003":"2013")) %>%
  ungroup() %>%
  summarise(median_velocity = median(velocity, na.rm = T),
            l95_value = HDInterval::hdi(velocity)[[1]],
            h95_value = HDInterval::hdi(velocity)[[2]])

vel_summary <- velocity_summary %>% 
  mutate(year = format(time, "%Y")) %>% 
  group_by(year) %>% 
  filter(year %in% c("2017":"2020")) %>% 
  ungroup() %>% 
  summarise(median_velocity = median(velocity, na.rm = T),
            l95_value = HDInterval::hdi(velocity)[[1]],
            h95_value = HDInterval::hdi(velocity)[[2]])

velocity_plot <- ggplot(velocity_summary) + 
  geom_ribbon(aes(time, ymin = X95.HPD_lower_value, ymax = X95.HPD_higher_value), fill="cornflowerblue", alpha = 0.5) +
  geom_step(aes(time, velocity), color="dimgrey") +
  labs(x="",y="Weighted branch velocity (km/year)") +
  #ylim(0,5) +
  theme_bw() +
  scale_x_date(limits = c(ymd("1977-01-01"), ymd("2022-02-07")), breaks= "2 year", date_labels = "%Y") +
  theme(panel.grid = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", linewidth = 0.2, linetype = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size=11, angle = 45, hjust = 1)) 

tiff("Results/Statistics.tiff", width = 20, height = 30, units = "cm",
     compression = "lzw", res = 300)
plot_grid(wavefront_plot, velocity_plot, ncol=1, align="v", rel_heights = c(1,1)) #labels = c("A", "B", "C", "D", "E"), 
dev.off()

diffusion_summary <- read.csv("Run_09052023_v1/Run_09052023_v1/Statistics/stat_v1_estimated_dispersal_statistics.csv", sep=";") 
diffusion_sum <- diffusion_summary %>% 
  summarise(median_velocity = median(weighted_branch_dispersal_velocity),
            l95_value_v = HDInterval::hdi(weighted_branch_dispersal_velocity)[[1]],
            h95_value_v = HDInterval::hdi(weighted_branch_dispersal_velocity)[[2]],
            median_diffusion = median(weighted_diffusion_coefficient),
            l95_value_d = HDInterval::hdi(weighted_diffusion_coefficient)[[1]],
            h95_value_d = HDInterval::hdi(weighted_diffusion_coefficient)[[2]])

# Plot introduction times and statistics 

tab_introduction <- as.data.frame(introduction_times_in_western_Africa)
introduction <- ggplot(tab_introduction, aes(x=introduction_times_in_western_Africa)) +
  geom_density(color=cols1[[1]],
               fill=cols2[[1]]) +
  labs(x="Date of dispersal event to West Africa",y="Density") +
  theme_classic() +
  #scale_x_date(breaks= "1 year", date_labels = "%Y") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size=11)) 
  
wavefront <- ggplot(wavefront_summary) + 
  geom_ribbon(aes(time, ymin = l95_value, ymax = h95_value), fill=cols2[[1]]) +
  geom_step(aes(time, median_value), color=cols1[[1]]) +
  labs(x="Time",y="Spatial wavefront distance (km)") +
  theme_classic() +
  scale_x_date(limits = c(ymd("1990-01-01"), ymd("2022-02-07")), breaks= "10 year", date_labels = "%Y") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size=11)) 

velocity <- ggplot(velocity_summary) + 
  geom_ribbon(aes(time, ymin = X95.HPD_lower_value, ymax = X95.HPD_higher_value), fill=cols2[[1]], alpha = 0.5) +
  geom_step(aes(time, velocity), color=cols1[[1]]) +
  labs(x="Time",y="Weighted branch velocity (km/year)") +
  theme_classic() +
  scale_x_date(limits = c(ymd("1990-01-01"), ymd("2022-02-07")), breaks= "10 year", date_labels = "%Y") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size=11)) 


tiff("Results/Figure 5_phylomap_v1_stat.tiff", width = 30, height = 10, units = "cm",
     compression = "lzw", res = 300)
plot_grid(introduction, wavefront, velocity, ncol=3) #align="h", rel_heights = c(1,1,1)) #labels = c("A", "B", "C", "D", "E"), 
dev.off() 


# 8. Extracting parameters from log file 

log <- read.table("Run_09052023_v1/Run_09052023_v1/Selected_sequences_v1_newcoord.log", header = TRUE)
log <- log %>% 
  mutate(date = date(date_decimal(age.root.)))
summary(log$date)
#Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
#"1389-03-19" "1999-02-25" "2003-11-11" "2002-06-29" "2007-10-08" "2012-04-01" 
#Look at tracer and manually enter 
date(date_decimal(2003.9729))
date(date_decimal(1989.3779))
date(date_decimal(2012.2513))



# 9. Compare parameters from sensitivity analysis

log_v1 <- read.table("Run_09052023_v1/Run_09052023_v1/Selected_sequences_v1_newcoord.log", header = TRUE)
log_m1 <- read.table("Run_09052023_m1/Run_09052023_m1/Selected_sequences_m1_newcoord.log", header = TRUE)
log_ns1 <- read.table("Run_09052023_ns1/Run_09052023_ns1/Selected_sequences_ns1_newcoord.log", header = TRUE)

v1 <- log_v1 %>% 
  select(c(meanRate, age.root.)) %>% 
  mutate(date = date(date_decimal(age.root.))) %>% 
  mutate(deme = "Random")
m1 <- log_m1 %>% 
  select(c(meanRate, age.root.)) %>% 
  mutate(date = date(date_decimal(age.root.))) %>% 
  mutate(deme = "Equivalent")
ns1 <- log_ns1 %>% 
  select(c(meanRate, age.root.)) %>% 
  mutate(date = date(date_decimal(age.root.))) %>% 
  mutate(deme = "Reduced")
sensitivity <- rbind(v1, m1, ns1)
levels(as.factor(sensitivity$deme))
sensitivity$deme <- factor(sensitivity$deme, levels=c("Equivalent", "Reduced", "Random"))

sensitivity_summary <- sensitivity %>%
  group_by(deme) %>%
  summarise(clock_l95_value = HDInterval::hdi(meanRate)[[1]],
            clock_h95_value = HDInterval::hdi(meanRate)[[2]],
            clock_median_value = median(meanRate),
            tmrca_l95_value = HDInterval::hdi(age.root.)[[1]],
            tmrca_h95_value = HDInterval::hdi(age.root.)[[2]],
            tmrca_median_value = median(age.root.))

clock <- ggplot(sensitivity, aes(x=meanRate, y=deme)) + 
  geom_violin(draw_quantiles=0.5, fill="gray") +
  labs(y="Sampling Scheme", x="Clock rate") +
  theme_classic() +
  theme(text = element_text(size = 15),
    legend.position = "none")

root <- ggplot(sensitivity, aes(x=age.root., y=deme)) + 
  geom_violin(draw_quantiles=0.5, fill="gray") +
  scale_x_continuous(limits=c(1980, 2015)) +
  labs(y="Sampling Scheme", x="tMRCA") +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "none")

tiff("Results/Figure 6_phylomap_v1m1ns1_sens.tiff", width = 20, height = 10, units = "cm",
     compression = "lzw", res = 300)
plot_grid(clock, root, ncol=2) #align="h", rel_heights = c(1,1,1)) #labels = c("A", "B", "C", "D", "E"), 
dev.off() 



# Combine tree with summary statistics
tree_figure2 <- ggtree(
  tr = beast_phylo_08_new, 
  mrsd="2022-02-07", 
  color="gray30",
  as.Date = T,
  aes()) +
  theme_tree2() +
  scale_x_date(limits = c(ymd("1991-01-01"), ymd("2022-02-07")), breaks= "2 year", date_labels = "%Y") +
  geom_tippoint(aes(fill = Country), shape=21, color="black", size=2) +
  scale_fill_manual(values = col) +
  #coord_cartesian(clip = 'off') + 
  #theme_tree2(plot.margin=margin(6, 30, 6, 6)) +
  theme(plot.margin = unit(c(0.5,1,0,1), "cm"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
        legend.position=c(.05,.85),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 12),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x =  element_line(colour = "grey85", linewidth = 0.1, linetype = 2))

phylo_08_data_new <- phylo_08_data_new %>% 
  mutate(node_prob_above_threshold = case_when(
    posterior > 0.75 ~ T,
    T ~ F))

tree_figure2 <- tree_figure2 %<+% 
  phylo_08_data_new + 
  geom_nodelab(aes(subset=node_prob_above_threshold==T, x=branch, label=round(as.numeric(posterior), digits=1)), vjust=-0.5, size = 1.9, colour = "gray50") +
  geom_range(range="height_0.95_HPD_dates", color="grey80", alpha=.6, size=2) 

show(tree_figure2)

wavefront_plot2 <- ggplot(wavefront_summary) + 
  geom_ribbon(aes(time, ymin = l95_value, ymax = h95_value), fill="#CC000064") +
  geom_step(aes(time, median_value), color="#CC0000FF") +
  labs(x="",y="Spatial wavefront distance (km)") +
  #ylim(0,5) +
  theme_bw() +
  scale_x_date(limits = c(ymd("1991-01-01"), ymd("2022-02-07")), breaks= "2 year", date_labels = "%Y") +
  theme(panel.grid = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey85", linewidth = 0.1, linetype = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size=11, angle = 45, hjust = 1)) 


tiff("Results/Figure 5_MCC_Stat_08.tiff", width = 20, height = 30, units = "cm",
     compression = "lzw", res = 300)
plot_grid(tree_figure2, wavefront_plot2, ncol=1, align="v", rel_heights = c(3,1)) #labels = c("A", "B", "C", "D", "E"), 
dev.off()


# 7. Plotting poultry density in Africa
density_poultry = raster("Data/dataverse_files/5_Ch_2015_Da.tif")

summary(getValues(density_poultry))
cuts=c(10, 25, 50, 150, 100, 500, 1000, 5000, 10000) #set breaks
pal <- brewer.pal(9, "YlOrBr")


tiff("Results/Figure 5.tiff", width = 20, height = 20, units = "cm",
     compression = "lzw", res = 300)
par(mar=c(0,0,0,0)) #, oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(template_shp, xlim=c(-20, 10), ylim=c(0, 37.5), col="ivory2", border =  "black", lwd=0.5)
plot(density_poultry, xlim=c(-20, 10), ylim=c(0, 37.5), breaks=cuts, col = pal, add=T) # 
dev.off()

