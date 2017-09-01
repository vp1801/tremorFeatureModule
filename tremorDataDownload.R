TREMOR_TABLE_ID = "syn5734657"
OBJECTIVE_PD_ID = "syn8533708"

FEATURES_PARENT_ID = "syn10154509"

SYNAPSE_CACHE = "~/synapseFeatures/synapseCache/tremorFiles"
WORK_DIR = "~/"
LIBRARY_PATH = "~/mylibs"
ACTIVITY_NAME = "Generate power density features for tremor tests"
thisFile = "https://github.com/vp1801/mPowerAnalysis/blob/tremor_module_dev/featureExtraction/tremorModule_V2.R"
.libPaths(LIBRARY_PATH)
setwd(WORK_DIR)

source("~/featureEval.R")
library('synapseClient')
library('signal')
library("plyr")
library("dplyr")
library("ggplot2")
library("doParallel")
library("parallel")
library("tidyr")
library("lubridate")
library('doMC')
library('Rwave')
library('data.table')
# github based module
library("mpowertools")


synapseLogin()

registerDoMC(detectCores()-10)
#get objectivePD health codes
healthCodesPD = read.table(synGet( OBJECTIVE_PD_ID)@filePath, sep=",", as.is=T, header=T)
#get data for above health codes
dataSet="objectivePD"
tremorDataTablePD = synTableQuery(paste0("SELECT * FROM ",TREMOR_TABLE_ID,
                                         " WHERE healthCode IN (", paste0("'",healthCodesPD$healthcode, "'", collapse=","), ")"))
#get all data
dataSet="allControl"
tremorDataTablePD = synTableQuery(paste0("SELECT * FROM ",TREMOR_TABLE_ID)," WHERE dataGroups == 'control'")
dataFields = c("deviceMotion_tremor_handInLap_right.json.items"           
               #,"accel_tremor_handInLap_left.json.items"                   
               ,"deviceMotion_tremor_handInLap_left.json.items"            
               #,"accel_tremor_handAtShoulderLength_right.json.items"       
               ,"deviceMotion_tremor_handAtShoulderLength_right.json.items"
               #,"accel_tremor_handAtShoulderLength_left.json.items"        
               ,"deviceMotion_tremor_handAtShoulderLength_left.json.items" 
               #,"accel_tremor_handToNose_right.json.items"                 
               ,"deviceMotion_tremor_handToNose_right.json.items"          
               #,"accel_tremor_handToNose_left.json.items"                  
               ,"deviceMotion_tremor_handToNose_left.json.items")
FEATURES_ID = c()
tremor.tbl.list = list()
feature_set = list()
#FEATURES_ID[dataFields[1]] = "syn10163170"
#FEATURES_ID[dataFields[1]] = ""
for (i in dataFields[1:4])  
{
  if(file.exists(paste0(i,"_", dataSet,".tsv"))){
    tremor.tbl = read.table(paste0(i,"_", dataSet,".tsv"), sep="\t", header = T)
  }
  else{
    synapseCacheDir(file.path(SYNAPSE_CACHE, i))
    ## Download all the json files and attach file locations as columns
    tbl.file.ids = synapseClient::synDownloadTableColumns(tremorDataTablePD, i) %>%
      unlist %>% data.frame 
    colnames(tbl.file.ids) = 'tremorJSONFileLocation'
    tbl.file.ids[,i] = rownames(tbl.file.ids)
    tbl.file.ids$tremorJSONFileLocation = as.character(tbl.file.ids$tremorJSONFileLocation)
    
    ## Join tremor file location to the tremor.tbl data frame
    tremor.tbl = tremorDataTablePD@values
    tremor.tbl$row.id = rownames(tremor.tbl)
    tremor.tbl = tremor.tbl %>%
      tidyr::separate(row.id, into = c("source.row.id","source.row.version"), sep = "_") %>%
      dplyr::mutate(source.row.id = as.integer(source.row.id),
                    source.row.version = as.integer(source.row.version),
                    row.id = paste(source.row.id, source.row.version, sep = '_')) %>%
      full_join(tbl.file.ids)
    #tremor.tbl.list[[i]] = tremor.tbl
    write.table(tremor.tbl,paste0(i,".tsv"), sep="\t", col.names = T, row.names = F)
  }
  tremor.tbl.list[[i]] = tremor.tbl
}