generateSpectrograms = function(row, TSLICE=c(100, 100), wl=256, remTimeSlice = T, rotate = T){
  #print(row['tremorJSONFileLocation'])
  if(is.na(as.character(row['tremorJSONFileLocation']))){
    feat = list()
    feat$acceleration = NA
    feat$timeStamp = NA
    feat$sampleRate = NA
    feat$spectrogram = NA
    feat$healthCode = NA
    feat$recordID = NA
    feat$dataGroups = NA
    feat$momentInDay = NA
    return(feat)
  }
  
  data = jsonlite::fromJSON(as.character(row[1,'tremorJSONFileLocation']))
  if(nrow(data) <= sum(TSLICE)){
    feat = list()
    feat$acceleration = NA
    feat$timeStamp = NA
    feat$sampleRate = NA
    feat$spectrogram = NA
    feat$healthCode = NA
    feat$recordID = NA
    feat$dataGroups = NA
    feat$momentInDay = NA
    return(feat)
  }
  if(remTimeSlice){
    data = data[-c((1:TSLICE[1]), (nrow(data) - TSLICE[2]):nrow(data)),]
  }

  if(rotate)
    correctedAcc = sqrt(rowSums(mpowertools:::get_quaternary_rotated_userAccel(data)^2))
  else
    correctedAcc = sqrt(rowSums((data$userAcceleration)^2))

  feat = list()
  
  ##Acceleration range feature
  if(diff(range(correctedAcc)) > 10){
    feat$recordID = NA
    return(feat)
  }
  
  feat$acceleration = correctedAcc 
  # negative of gravity gives orientation wrt gravity or ground
  #dot product with x=0,y=0,z=1 gives angle with z axis, similarly for x and y axes
  feat$angleWithZ = -1*data$gravity$z/(rowSums(data$gravity ^ 2, na.rm=T)^0.5)
  
  feat$timeStamp = data$timestamp
  feat$sampleRate = length(unique(data$timestamp))/diff(range(data$timestamp))
 
  #feat$spectrogram = signal::specgram(correctedAcc - mean(correctedAcc, na.rm=T), Fs = feat$sampleRate, window = wl)
  
  
  #calculate fundamental frequency
  # fundamental frequency in kHz, need to multiply by 1000 when finding median and variance
  #feat$funFrequency = plyr::adply(.data = feat$spectrogram$t, .margins = 1, .fun = function(x){
  #  fundfreq = seewave::fund(feat$acceleration - mean(feat$acceleration, na.rm=T), wl=wl, at= (x + 128/feat$sampleRate), plot=F, f = feat$sampleRate)
  #  data.table(fundfreq)
  #})
  
  feat$healthCode = row[1,'healthCode']
  feat$recordID = row[1,'recordId']
  feat$dataGroups = row[1,'dataGroups']
  feat$momentInDay = row[1,'momentInDayFormat.json.choiceAnswers']
  feat
}

frequencyFeatures = function(spectrograms, freq_range = NULL, nbins = 5, normalise = T, logScale=T, bin = NULL, tremor.activity = "", wl = 256, ovlp = 0.5){
  featureSet = Filter(function(obj) !is.na(obj$recordID), spectrograms) %>% plyr::llply(.fun = function(obj){
    if(is.null(freq_range))
      freq_range = c(0.1, obj$sampleRate/2)
    else
      freq_range[1] = max(freq_range[1], 0.1, na.rm=T)
    
    
    tmp = tuneR::Wave(left = obj$acceleration, samp.rate = obj$sampleRate) %>% 
      tuneR::periodogram(width=wl, overlap = ovlp*wl, frqRange = freq_range, normalize = normalise)
    spec = data.frame(tmp@spec)
    #spec = abs(obj$spectrogram$S)
    colnames(spec) = paste0("tBin", 1:ncol(spec))
    if(logScale)
      spec = log(spec)
    
    #f = obj$spectrogram$f
    f= tmp@freq
    
    if(is.null(bin)){
      cutPoints=f
      nbins=length(f)
    }
    else{
      if(bin =="log")
        cutPoints = exp(seq(log(freq_range[1]), log(freq_range[2]), length.out = nbins))
      else
        cutPoints = seq(freq_range[1], freq_range[2], length.out = nbins)
    }
    binDF = data.frame(fbin = paste0("bin", 0:nbins), minFreq = c(0,cutPoints), maxFreq = c(cutPoints, obj$sampleRate/2))
    

    
    if(!is.null(freq_range)){
      idx = (f >= freq_range[1] & f<freq_range[2])
      f = f[idx]
      spec = spec[idx,]
    }
    else
      freq_range = range(f) + 1e-5
    
    bins = findInterval(f, cutPoints)
    bins[is.na(bins)] = 0
    
    #colnames(spec) = paste(tremor.activity, prettyNum(binDF$minFreq, digits=3, format="f"), 
    #                       prettyNum(binDF$maxFreq, digits=3, format="f"), "Hz",sep="_")
    ##Attach spectral power density for different frquencies in different time slices
    features = data.table(recordId = obj$recordID, dataGroups = obj$dataGroups, healthCode = obj$healthCode, medication = obj$momentInDay,
                          fbin = paste0("bin", bins), f=f) %>% cbind(spec) %>% 
      melt(id.vars = c('recordId', 'dataGroups', 'healthCode', 'medication', 'fbin'), measure.vars = colnames(spec), variable.name = "tBin", value.name = "powerDensity")
    ##Average/Sum spectral density for each frequency bin in different time slices
    aggFreq = features[,mean(powerDensity, na.rm=T), 
                       by=.(tBin, fbin, recordId, dataGroups, healthCode, medication)]
    aggFreq = merge(aggFreq, binDF)
    colnames(aggFreq) = c("frequencyBin", "time_slice", "recordId", "dataGroups", "healthCode", "medication", "meanPowerDensity", "minFreq", "maxFreq")
    aggFreq[,freqRange := paste(tremor.activity, prettyNum(minFreq, digits=3, format="f"), prettyNum(maxFreq, digits=3, format="f"), "Hz",sep="_")]
    
    freqFeaturesAgg = dcast(aggFreq, recordId + healthCode+ dataGroups + time_slice~ freqRange, value.var='meanPowerDensity', fill = 0)
  
    freqFeaturesAgg
  }, .parallel = T)
  return(featureSet)
}

ShapeTremorData = function(tremor.tbl, tremor.activity = "", freq_range=c(0.3,10), nbins=20, remove_rotation = T){
  records = plyr::dlply(.data = tremor.tbl, .variables = "recordId", .parallel=T,
                              .progress = 'text',.fun = generateSpectrograms, rotate = T) %>% 
    Filter(f = function(x) !is.na(x$recordID) && !is.null(x$recordID))
  if(remove_rotation && !grepl("Nose", tremor.activity)){
    rotRemove = removeRotation_angle(records, rel.threshold = 1.2)
    records[names(rotRemove)] = rotRemove
  }
  shapedFeatures = frequencyFeatures(records, freq_range = freq_range, tremor.activity = tremor.activity, nbins=nbins, bin = "linear")
  return(shapedFeatures)
}

removeRotation_angle = function(spec_features, abs.threshold = NULL, rel.threshold = NULL){
  angleWithZRange = spec_features %>% lapply(function(x) zoo::rollapply(x$angleWithZ, width=256, by = 128, FUN = function(x) diff(range(x))))
  filteredRecs = list()
  rotationRemoved = list()
  print("range calculated")
  if(!is.null(abs.threshold)){
    filteredRecs = Filter(function(x) max(x) > abs.threshold, angleWithZRange)
    rotationRemoved = plyr::alply(names(filteredRecs), .margins = 1, .fun = function(recId){
      remTS = ind2Loc(which(angleWithZRange[[recId]] > abs.threshold))
      obj = spec_features[[recId]]
      obj$acceleration = obj$acceleration[-remTS]
      obj$angleWithZ = obj$angleWithZ[-remTS]
      return(obj)
    })
    names(rotationRemoved) = names(filteredRecs)
  }
  
  if(!is.null(rel.threshold)){
    filteredRecs = Filter(function(x) max(x) > rel.threshold && sum(x > median(x)+2*sd(x)) > 0, angleWithZRange)
    rotationRemoved = plyr::alply(names(filteredRecs), .margins = 1, .fun = function(recId){
      remTS = ind2Loc(which((angleWithZRange[[recId]] > rel.threshold) & 
                              (angleWithZRange[[recId]] > median(angleWithZRange[[recId]]) + 2*sd(angleWithZRange[[recId]]))))
      obj = spec_features[[recId]]
      obj$acceleration = obj$acceleration[-remTS]
      obj$angleWithZ = obj$angleWithZ[-remTS]
      return(obj)
    })
    names(rotationRemoved) = names(filteredRecs)
  }
  return(rotationRemoved)
}

getSpecDensityFeatures = function(tremor.tbl.list, remove_rotation_slice = T){
  feature_list = list()
  for(i in names(tremor.tbl.list)){
    ### Frequency bins and power density in different time slices for each record
    test = gsub(pattern = "deviceMotion_tremor_",replacement = "",x = i)
    shapedFeatures = ShapeTremorData(tremor.tbl.list[[i]], tremor.activity = strsplit(x = test, split="[\r.]")[[1]][1])
    ## Convert the power density to time series features
    feature_list[[i]] = plyr::llply(shapedFeatures, .fun = function(rec){
        features = c("recordId" = rec$recordId[1])
        for(j in grep("Hz",colnames(rec), value=T))
          features = c(features, SingleAxisFeatures(unlist(rec[,j,with=F]), 
                                                    tmp_time = 128*((1:nrow(rec)) - 1),varName = j))
        data.frame(t(features))
      }, .parallel = T) %>% data.table::rbindlist(fill=T)
  }
  
  
  # by = c("recordId", "healthCode", "dataGroups", "time_slice")
  tremor.spectral.features = plyr::join_all(feature_list, by = c("recordId"), type = "full", match="all")
  return(tremor.spectral.features)
}


SingleAxisFeatures <- function(x, tmp_time, varName) {
  meanX <- mean(x, na.rm = TRUE)
  sdX <- sd(x, na.rm = TRUE)
  modeX <- pracma::Mode(x)
  skewX <- e1071::skewness(x)
  kurX <- e1071::kurtosis(x)
  auxX <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  q1X <- auxX[[2]]
  medianX <- auxX[[3]]
  q3X <- auxX[[4]]
  iqrX <- q3X - q1X
  rangeX <- auxX[[5]] - auxX[[1]]
  acfX <- stats::acf(x, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
  zcrX <- mpowertools:::ZCR(x)
  dfaX <- tryCatch({
    fractal::DFA(x, sum.order = 1)[[1]]
  }, error = function(err){ NA })
  cvX <- mpowertools:::Cv(x)
  tkeoX <- mpowertools:::MeanTkeo(x)

  out <- c(meanX, sdX, modeX, skewX, kurX,
           q1X, medianX, q3X, iqrX, rangeX, acfX,
           zcrX, dfaX, cvX, tkeoX)
  
  nms <- c("mean", "sd", "mode", "skew", "kur", "q1",
           "median", "q3", "iqr", "range",
           "acf", "zcr", "dfa", "cv", "tkeo")
  
  names(out) <- paste(nms, varName, sep = "_")
  return(out)
}
