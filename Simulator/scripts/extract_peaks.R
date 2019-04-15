# needs xcms3, see https://bioconductor.org/packages/release/bioc/html/xcms.html for installation
library(xcms)

### some useful functions ###

getDesign = function(mzml_files, label) {
  design <- data.frame(sample_name = sub(basename(mzml_files), pattern = '.mzML', replacement = "", fixed = TRUE),
                       sample_group = c(rep(label, length(mzml_files))),
                       stringsAsFactors = FALSE)
  return(design)  
}

readData = function(files, design, mode="onDisk", msLevel=1) {
  raw_data <- MSnbase::readMSData(files = files, pdata = new("NAnnotatedDataFrame", design),
                                  mode = mode, msLevel. = msLevel, centroided. = TRUE)
  return(raw_data)
}

pickPeaks = function(raw_data, ppm, peakwidth, snthresh, prefilter, mzdiff) {
  # https://rdrr.io/bioc/xcms/man/findChromPeaks-centWave.html
  cwp = CentWaveParam(ppm=ppm, peakwidth=peakwidth, snthresh=snthresh,
                      prefilter=prefilter, mzdiff=mzdiff)
  xdata <- findChromPeaks(raw_data, param = cwp)
  return(xdata)
}

extract_peaks <- function(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff, msLevel) {
  # TODO: at the moment, findChromPeaks only works with OnDiskMSnExp, that's why we pass it raw_data_on_disk for peak picking
  raw_data <- readData(mzml_files, design, mode='onDisk', msLevel=msLevel)
  xdata <- pickPeaks(raw_data, ppm=ppm, peakwidth=peakwidth, snthresh=snthresh, prefilter=prefilter, mzdiff=mzdiff)
  return(xdata)
}

get_df = function(xdata, msLevel) {
  df <- data.frame(chromPeaks(xdata))
  df$msLevel <- msLevel
  df$filename <- basename(fileNames(xdata)[df$sample])
  return(df)  
}

write_df <- function(df, filename) {
  # write.csv(df, file = gzfile(filename), row.names = FALSE)  
  write.csv(df, file = filename, row.names = FALSE)  
}

process_dir <- function(mzml_dir, extracted_peaks_out, aligned_features_out, aligned_intensities_out, 
                        msLevel, ppm, bandwidth, snthresh, prefilter, mzdiff, binSize, minFraction, bw) {

  # run for all files in mzml_dir
  mzml_files <- list.files(path=mzml_dir, pattern='*.mzML', full.names=TRUE)
  design <- getDesign(mzml_files, "group1")
  
  xdata <- extract_peaks(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff, msLevel)
  peaks_df <- get_df(xdata, msLevel)
  write_df(peaks_df, paste(mzml_dir, extracted_peaks_out, sep='/'))
  
  # Retention time alignment
  xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = binSize))
  
  ## Correspondence: group peaks across samples
  pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                          minFraction = minFraction, bw=bw)
  xdata <- groupChromPeaks(xdata, param = pdp)
  feature_df <- featureDefinitions(xdata)
  write_df(feature_df, paste(mzml_dir, aligned_features_out, sep='/'))
  
  ## Extract the into column for each feature.
  into_df <- featureValues(xdata, value = "into")
  write_df(into_df, paste(mzml_dir, aligned_intensities_out, sep='/'))
    
}

### processing starts here ###

# ppm value is set fairly large.
# Other parameters for peakwidth, snthresh, prefilter for justin's data, taken from 
# https://www.dropbox.com/home/Meta_clustering/ms2lda/large_study/r/beer_method_3_pos?preview=xcmsPeakPicking.R
ppm <- 10
peakwidth <- c(5, 100)
snthresh <- 3
prefilter <- c(3, 1000)
mzdiff = 0.001

# alignment and grouping parameters
binSize <- 0.6 # ObiWarp bin size
minFraction <- 0.8 # minimum fractions for a group to occur in all samples
bw <- 30 # bandwidth

msLevel <- 1
mzml_dir <- '../models/dda_results_test'
process_dir(mzml_dir, 'extracted_peaks_ms1.csv', 'aligned_features_ms1.csv', 'aligned_intensities_ms1.csv', 
            msLevel, ppm, bandwidth, snthresh, prefilter, mzdiff, binSize, minFraction, bw)
