## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE-----------------------------------------------------------
#  install.packages("devtools")

## ---- eval=FALSE-----------------------------------------------------------
#  library(devtools)
#  install_github("jbergenstrahle/SpatialCPie")

## ---- eval=TRUE, message=FALSE, warning=FALSE------------------------------
library(SpatialCPie)

## ---- eval=FALSE-----------------------------------------------------------
#  counts = read.table("countsFromSTpipeline.tsv")
#  counts = t(counts)

## ---- eval=TRUE, echo=FALSE, results='asis'--------------------------------
setwd('..')
load("data/fetalHeartCountsPrepared.RData")
knitr::kable(head(counts[,1:5], 5))

## ---- eval=FALSE, echo=TRUE------------------------------------------------
#  load("SpatialCPieData.RData") #Provided RData file

## ---- eval=FALSE-----------------------------------------------------------
#  #Lets remove the spots that are outside of the tissue
#  D2spotsUnderTissue = paste(D2spots$V1, D2spots$V2, sep="x")
#  D2spotsUnderTissue = D2spotsUnderTissue[2:length(D2spotsUnderTissue)]
#  D2counts = D2counts[, which(colnames(D2counts) %in% D2spotsUnderTissue)]
#  

## ---- eval=FALSE-----------------------------------------------------------
#  library(tidyverse)
#  clusters = 2:5 %>% map(~kmeans(t(D2counts), centers=.)$cluster)
#  

## ---- eval=FALSE-----------------------------------------------------------
#  viz = runCPie(D2counts, clusters)

## ---- eval=FALSE-----------------------------------------------------------
#  viz = runCPie(D2counts, clusters, view="browser")
#  

## ---- eval=FALSE-----------------------------------------------------------
#  coords <- parseSpotFile("spot_data-sel-BC24044_D2.tsv")
#  #In the provided RData file, coords is already loaded

## ---- eval=FALSE-----------------------------------------------------------
#  library(jpeg)
#  HE_img = jpeg::readJPEG( [path to image] )

## ---- eval=FALSE-----------------------------------------------------------
#  viz = runCPie(D2counts, clusters, img = HE_img, pixel.coords = coords)

## ---- eval=FALSE-----------------------------------------------------------
#  viz$piePlots$`4`

## ---- eval=FALSE-----------------------------------------------------------
#  subCounts = D2counts[, spots]

## ---- eval=FALSE-----------------------------------------------------------
#  subClusters =  2:5 %>% map(~kmeans(t(subCounts), centers=.)$cluster)
#  viz = runCPie(subCounts, subClusters)

## ---- eval=FALSE-----------------------------------------------------------
#  library(jpeg)
#  library(grid)
#  
#  #counts = our count matrix
#  xcoord = as.numeric(sapply(strsplit(colnames(counts), "x"), "[[", 1))
#  ycoord = as.numeric(sapply(strsplit(colnames(counts), "x"), "[[", 2))
#  coord_df = as.data.frame(cbind(x=xcoord, y=ycoord))
#  df <- cbind(coord_df, color=as.factor(as.character(viz$clusters$`4`)))
#  
#  #read in the image and make a rasterGrob
#  he.image = readJPEG("HE_BT24044_D2.jpg")
#  grobHE <- rasterGrob(he.image, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
#  
#  #The scales are manually adjusted
#  ggPIE = ggplot(df, aes(x=x, y=36-y, colour=color))+
#            labs(colour="Cluster")+
#            scale_x_continuous(limits = c(1, 33), expand = c(0,0))+
#            scale_y_continuous(limits = c(1, 35), expand = c(0,0))
#  
#  #Specify the output name of the picture
#  out.file = "HE_BT24044_D2_Supp_4clusters.pdf"
#  pdf(out.file, onefile = TRUE, useDingbats = FALSE, height = 35 * 0.5, width = 33 * 0.5)
#  
#  #We add our image as background to our ggplot object
#  ggPIE = ggPIE + annotation_custom(grobHE, -Inf, Inf, -Inf, Inf)
#  
#  #We use viz$arrayPiePlotsPieInfo to extract info regarding the pie chart ratios for cluster dimension 4 that we extracted earlier
#  ggPIE + geom_scatterpie(data=viz$arrayPiePlotsPieInfo[[4]], aes(x=x, y=36-y, r=0.5),
#                                       cols=colnames(viz$arrayPiePlotsPieInfo[[4]][1:4]), alpha=0.5)
#  
#  dev.off()
#  

