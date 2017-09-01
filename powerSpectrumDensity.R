feature_set = list()
for(i in names(tremor.tbl.list)[1:4]){
  spec_features = plyr::dlply(.data = tremor.tbl.list[[i]], .variables = "recordId", .parallel=T,
                              .progress = 'text',.fun = generateSpectrograms, rotate = T) %>% Filter(f = function(x) !is.na(x$recordID) && !is.null(x$recordID))
  if(!grepl("Nose", i) && T){
    rotRemoved = removeRotation_angle(spec_features, rel.threshold = 1.2)
    spec_features[names(rotRemoved)] = rotRemoved
  }
  ##Frequency binned features (averaged over time)
  freqFeatures = frequencyFeatures(spec_features, nbins=20, freq_range = c(0.7,10), logScale = T, 
                                   tremor.activity = i, bin = "linear", wl=256) %>% 
    data.table::rbindlist(fill=T)
  
  freqFeaturesAgg = freqFeatures[,lapply(.SD, mean, na.rm=T), by=.(recordId, dataGroups), .SDcols = grep("Hz", colnames(freqFeatures), value=T)]
  feature_set[[i]] = freqFeaturesAgg

}

combinedData = list()
res = list()
for(i in c("handInLap", "handAtShoulder", "handToNose")){
  combinedData[[i]] = plyr::join_all(feature_set[grep(i, names(feature_set), value=T)], type="inner", by = "recordId", match="all")
  res[[i]] = featureTest(combinedData[[i]], grep("Hz",colnames(combinedData[[i]]), value = T))
  
}

for( i in names(combinedData)){
  colAnno = HeatmapAnnotation(df = data.frame(test = c(rep("Right", sum(grepl("right", colnames(combinedData[[i]])))), 
                                                       rep("Left", sum(grepl("left", colnames(combinedData[[i]])))))), 
                              col = list(test = c("Right" = "red","Left" = "blue")),
                              na_col = "gray", width = unit(1, "cm"))
  hmStatus = Heatmap(combinedData[[i]][,grep("Hz", colnames(combinedData[[i]])), with = F], show_row_dend = F, show_column_dend = F,
                     split = combinedData[[i]]$dataGroups, name = paste(i, "log_power_density", sep="_"), col = rev(heat.colors(10)), cluster_rows = T, cluster_columns = F,
                     top_annotation = colAnno, show_heatmap_legend = F)
  # split=5
  dataGroup = HeatmapAnnotation(df = data.frame(type = combinedData[[i]]$dataGroups), #zAxisMovement = anomalousData),
                                col = list(type = c("control" = "blue", "parkinson" = "red")), #zAxisMovement = c("TRUE" = "red", "FALSE" = "blue")),
                                na_col = "black", which = "row", 
                                width = unit(1, "cm"), show_legend = F)
  png(paste0("plots/", i,"_energyDensity_bandPass.png"))
  draw(hmStatus + dataGroup)
  dev.off()
}

