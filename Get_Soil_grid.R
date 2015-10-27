########################################################################################
###   This script is provided by the Soil & Landscape Grid of Australia to           ###
###   demonstrate accessing the soil and landscape data sets using OGC web services  ###
###   OGC web services.                                                              ###
###                                                                                  ###
###   URL for Soil and Landscape Grid is www.csiro.au/soil-and-landscape-grid        ###
###                                                                                  ###
###   Author : Ross Searle - email: ross.searle@csiro.au                             ###
###                                                                                  ###
########################################################################################


library(RCurl)
library(rgdal)
library(raster)
library(sp)
library(XML)
library(maptools)
library(hash)
library(plyr)
library(GSIF)
library(aqp)
library(epiR)

###################  Required functions     ############################################

urlRoot = "http://www.asris.csiro.au/ArcGis/services/"

iniAttributeInfo <- function(){
  #cat("Generating attribute info dictionary")
  dict <- hash()
  assign("Depth_to_Rock", c("None", "Meters", "DER"), dict)
  assign("Rooting_Depth", c("None", "Meters", "DPE"), dict)
  assign("Organic_Carbon", c("Log", "%", "SOC"), dict)
  assign("pH_Soil_Water", c("None", "", "PHW"), dict)  
  assign("pH_Soil_CaCl2", c("None", "", "PHC"), dict)
  assign("Clay", c("None", "%", "CLY"), dict)
  assign("Silt", c("None", "%", "SLT"), dict)
  assign("Sand", c("None", "%", "SND"), dict)
  assign("ECEC", c("Log", "meq/100g", "ECE"), dict)
  assign("Bulk_Density", c("None", "g/cm", "BDW"), dict)
  assign("Available_Water_Capacity", c("None", "%", "AWC"), dict)
  assign("EC", c("Log", "dS/m", "ECD"), dict) 
  assign("Total_P", c("Log", "%", "PTO"), dict)
  assign("Total_Nitrogen", c("Log", "%", "NTO"), dict)  
  return (dict)
}

iniDepthInfo <- function(){
  #cat("Generating Depth info dictionary")
  dict <- hash()
  assign("0-5", c(1), dict)
  assign("5-15", c(4), dict)
  assign("15-30", c(7), dict)
  assign("30-60", c(10), dict)
  assign("60-100", c(13), dict)
  assign("100-200", c(16), dict) 
  return (dict)
}

iniComponentInfo <- function(){
  #cat("Generating Depth info dictionary")
  dict <- hash()
  assign("value", c(0), dict)
  assign("lower_ci", c(2), dict)
  assign("upper_ci", c(1), dict)
  return (dict)
}

iniProductInfo <- function(){
  
  dict <- hash()
  assign("National", c("TERN", "ACLEP_AU_NAT_C"), dict)
  assign("WA", c("TERN", "ACLEP_AU_WAT_D"), dict)
  assign("SA", c("TERN", "ACLEP_AU_SAT_D"), dict)
  assign("TAS", c("TERN", "ACLEP_AU_TAS_N"), dict)
  assign("Australian-Wide", c("TERN", "ACLEP_AU_TRN_N"), dict)
  return (dict)
}

getServiceInfo <- function(region, attribute, component, depth){
  #cat("Getting Layer Info")
  
  attributeInfo <- iniAttributeInfo()
  depthInfo <- iniDepthInfo()
  componentInfo <- iniComponentInfo()
  productInfo <- iniProductInfo()
  
  ainfo<- get(attribute, attributeInfo)
  attCode = ainfo[3]
  
  cinfo = get(component, componentInfo)
  compCode = cinfo[1]
  
  pinfo = get(region, productInfo)
  
  
  dinfo = get(depth, depthInfo)
  depthCode = dinfo[1]
  layerId <- depthCode + compCode
  
  URLp1 = paste(urlRoot, pinfo[1],"/", attCode, "_", pinfo[2], "/MapServer/WCSServer?", sep="") 
  return (list(URLp1, layerId))
}

getWCSURL <- function(product, attribute, component, depth, extents){
  
  cols = as.integer((maxX - minX) / 0.000833333);
  rows = as.integer((maxY - minY) / 0.000833333);
  
  id = getServiceInfo(product, attribute, component, depth )
  urlBase = id[[1]][1]
  layer = id[[2]][1]
  
  wcs = paste(urlBase, "REQUEST=GetCoverage&SERVICE=WCS&VERSION=1.0.0&COVERAGE=", layer, "&CRS=EPSG:4283&BBOX=", extents[[1]][1], ",", extents[[3]][1], ",", extents[[2]][1],",", extents[[4]][1], "&WIDTH=", cols, "&HEIGHT=", rows, "&FORMAT=GeoTIFF",sep="")
  
  return (wcs)
}


getANZSoilWFSasDT <- function(url, region, extents, attribute, localData){
  
  xml.request = paste('<?xml version="1.0"?><wfs:GetFeature xmlns:om="http://www.opengis.net/om/2.0" xmlns:gml="http://www.opengis.net/gml/3.2" xmlns:anzsmlss="http://anzsoil.org/ns/soilsample/2.0.0" xmlns:sam="http://www.opengis.net/sampling/2.0" xmlns:wfs="http://www.opengis.net/wfs" xmlns:fes="http://www.opengis.net/fes/2.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" service="WFS" version="1.1.0" outputFormat="application/gml+xml; version=3.2" xsi:schemaLocation="http://www.opengis.net/wfs http://schemas.opengis.net/wfs/1.1.0/wfs.xsd http://anzsoil.org/ns/soilsample/2.0.0 http://www.clw.csiro.au/aclep/ANZSoilML/trunk/ANZSoilML_SoilSample/xsd/anzsoilml-soilsample.xsd" resultType="results">
                      <wfs:Query typeName="om:OM_Observation">
                      <ogc:Filter>
                      <ogc:And>
                      <ogc:PropertyIsEqualTo>
                      <!--<ogc:PropertyName>om:procedure/gsmlla:AnalyticalProcess/gml:name</ogc:PropertyName>-->
                      <ogc:PropertyName>om:observedProperty/@xlink:title</ogc:PropertyName>
                      <ogc:Literal>DENSITY</ogc:Literal>
                      </ogc:PropertyIsEqualTo>
                      <ogc:BBOX>
                      <gml:Envelope srsName="EPSG:4283">
                      <gml:lowerCorner>' , extents[[1]][1]  , ' ' , extents[[3]][1] ,'</gml:lowerCorner>
                      <gml:upperCorner>' , extents[[2]][1]  , ' ' , extents[[4]][1] ,'</gml:upperCorner>
                      </gml:Envelope>
                      </ogc:BBOX>
                      </ogc:And>
                      </ogc:Filter>
                      </wfs:Query>
                      </wfs:GetFeature>', sep='')
  
  # cat(xml.request)
  if(localData)
  {
    xmltext  <- xmlTreeParse('c:/temp/sali.xml',useInternalNodes=T)
  }
  else
  {    
    myheader=c(Connection="close", 'Content-Type' = "application/xml", 'Content-length'=nchar(xml.request))
    data =  getURL(url = url, postfields=xml.request, httpheader=myheader, verbose=TRUE)
    xmltext  <- xmlTreeParse(data, asText = TRUE,useInternalNodes=T)
    sink(paste('c:/temp/', region, '.xml', sep=''))
    print(xmltext)
    sink()
  }
  
  print(xmltext)
  
  #writeLines(xmltext, "c:/temp/natsoil.xml")
  
  valsFilt = '[../om:featureOfInterest/sams:SF_SpatialSamplingFeature/sams:shape/gml:Point/gml:pos]'
  Vals <- as.numeric(unlist(xpathApply(xmltext, paste('//wfs:FeatureCollection/wfs:member/om:OM_Observation/om:result', valsFilt , sep=''),xmlValue)) )
  idsFilt = '[../om:featureOfInterest/sams:SF_SpatialSamplingFeature/sams:shape/gml:Point/gml:pos]'
  ids <- unlist(xpathApply(xmltext, paste('//wfs:FeatureCollection/wfs:member/om:OM_Observation/@gml:id', valsFilt , sep=''))) 
  methodsFilt = '[../../../om:featureOfInterest/sams:SF_SpatialSamplingFeature/sams:shape/gml:Point/gml:pos]'
  methods <- unlist(xpathApply(xmltext, paste('//wfs:FeatureCollection/wfs:member/om:OM_Observation/om:procedure/gsmlla:AnalyticalProcess/gml:name', methodsFilt , sep=''),xmlValue)) 
  datesFilt = '[../om:featureOfInterest/sams:SF_SpatialSamplingFeature/sams:shape/gml:Point/gml:pos]'
  dates <- unlist(xpathApply(xmltext,paste('//wfs:FeatureCollection/wfs:member/om:OM_Observation/om:resultTime', datesFilt , sep=''),xmlValue)) 
  locFilt = '[../om:featureOfInterest/sams:SF_SpatialSamplingFeature/sams:shape/gml:Point/gml:pos]'
  locs <- unlist(xpathApply(xmltext,'//wfs:FeatureCollection/wfs:member/om:OM_Observation/om:featureOfInterest/sams:SF_SpatialSamplingFeature/sams:shape/gml:Point/gml:pos',xmlValue))  
  upper_depths <- as.numeric(unlist(xpathApply(xmltext,'//wfs:FeatureCollection/wfs:member/om:OM_Observation/om:featureOfInterest/sams:SF_SpatialSamplingFeature/sam:parameter/om:NamedValue/om:value/anzsml:DepthQuantityRange/anzsml:upperBoundary/swe:Quantity/swe:value',xmlValue)))
  lower_depths <- as.numeric(unlist(xpathApply(xmltext,'//wfs:FeatureCollection/wfs:member/om:OM_Observation/om:featureOfInterest/sams:SF_SpatialSamplingFeature/sam:parameter/om:NamedValue/om:value/anzsml:DepthQuantityRange/anzsml:lowerBoundary/swe:Quantity/swe:value',xmlValue)))
  
  lats <- numeric()
  longs<- numeric()
  IDcomps<-character()
  regionID<-character()
  for (i in 1:length(locs) ) { 
    s =  strsplit(locs[i], " ")
    lats<- append(lats, as.numeric(s[[1]][2]))
    longs<- append(longs, as.numeric(s[[1]][1]))
    
    idc = strsplit(ids[i], "_")
    idc2 = strsplit(idc[[1]][3], "\\.")
    IDcomps<-append(IDcomps, paste(idc2[[1]][1], idc2[[1]][2], idc2[[1]][3], sep="_"))
    regionID<-append(IDcomps, region)
    
  }
  
  
  # just for checking purposes
  length(methods)
  length(IDcomps)
  length(ids)
  length(dates)
  length(lats)
  length(longs)
  length(upper_depths)
  length(lower_depths)
  length(Vals)
  
  
  
  dt <- data.frame(region, IDcomps,lats, longs, methods, dates, upper_depths, lower_depths, Vals)
  
  return (dt)
}

###################  End of Required functions    ######################################

###################  Set these parameters to meet your requirements  ###################
###                                                                                  ###
###  PLEASE NOTE : the current Web Coverage Service (WCS) is a restricted download   ###
###               size limit. It is limited to 3x3 arc seconds. If you try to        ###
###               a larger area R will hang.                                         ###

for(Attribute in c("Organic_Carbon")){

Region = "National"  #  see 'iniProductInfo' function for choices
Attribute = Attribute   #  see 'iniProductInfo' function for choices
Component = "value"   #  see 'iniDepthInfo' function for choices
DepthInterval = "0-5"  #  see 'iniDepthInfo' function for choices
DownloadDirectory = 'c:/temp'

# extent to extract in WGS 84 coordinates - refer note above about max size
minX = 149
maxX = 152
minY = -33
maxY = -28



###################  End of parameters section   ######################################

###################  Run the commands below to extract and download the raster data ###

setwd(DownloadDirectory)
extents = list(minX, maxX, minY, maxY)
downloadName = paste(DownloadDirectory,'/', Region,'_', Attribute,'_', Component,'_', DepthInterval, '.tif', sep='')


url = getWCSURL(Region, Attribute, Component, DepthInterval, extents)
url

bin = getBinaryURL(url)
con <- file(downloadName, open = "wb")
writeBin(bin, con)
close(con)
}
r <- raster(downloadName)
attributeInfo <- iniAttributeInfo()
ainfo<- get(Attribute, attributeInfo)
attUnits = ainfo[2]
plot(r, main=paste(Region, Attribute, Component, DepthInterval, 'cm (', attUnits ,')', sep=' '))

