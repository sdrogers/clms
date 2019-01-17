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

filterDf <- function(df, rt_start, rt_end, min_intensity, min_sn) {
  peak_info <- df[which(df$rt >= rt_start),]
  peak_info <- peak_info[which(peak_info$rt <= rt_end),]
  peak_info <- peak_info[which(peak_info$rtmin >= rt_start),]
  peak_info <- peak_info[which(peak_info$into >= min_intensity),]
  peak_info <- peak_info[which(peak_info$sn >= min_sn),]
  return(peak_info)  
}

get_feature = function(raw_data, row, msLevel) {
  sample <- row$sample
  mz_range_df <- row[, c('mzmin', 'mzmax')]
  rt_range_df <- row[, c('rtmin', 'rtmax')]

  # TODO: we need to add some padding to the mz_range to get the complete chromatogram, why?!
  mz_range <- c(mz_range_df$mzmin-0.0003, mz_range_df$mzmax+0.0003)
  
  # alternatively we have a ppm constant for the bounds
  # mz_range = c(mz * (1-mzppm/1e6), mz * (1+mzppm/1e6))
  
  rt_range <- c(rt_range_df$rtmin, rt_range_df$rtmax)

  # below will not give you the mz values of each data point  
  # chr_raw <- chromatogram(raw_data, mz=mz_range, rt=rt_range)
  # chr_raw <- chr_raw[1, sample]

  # below is the correct code
  one_file = filterFile(raw_data, sample)
  filtered_file = filterRt(filterMz(one_file, mz_range), rt_range)
  chrom_data = as(filtered_file, 'data.frame')
  # plotMsData(chrom_data)

  # TODO: normalise the rt values, first value (start) is 0
  # TODO: normalise the mz values, the average is 0
  # TODO: not sure about the intensities
  mz_values = chrom_data$mz
  rt_values = chrom_data$rt
  intensity_values = chrom_data$i
  
  return(list(id=rownames(row), mz=row$mz, mzmin=row$mzmin, mzmax=row$mzmax, 
              rt=row$rt, rtmin=row$rtmin, rtmax=row$rtmax,
              into=row$into, maxo=row$maxo, sn=row$sn, 
              mz_values = mz_values, rt_values = rt_values, intensity_values = intensity_values,
              msLevel=msLevel, fromFile=sample))    
}

sortDf <- function(df, decreasing=TRUE) {
  speaks = df[order(df[,'into'], decreasing=decreasing),]
  return(speaks)
}

get_all_features <- function(data, filtered_list, msLevel, N) {
  start.time <- Sys.time()
  # all_features <- lapply(filtered_list[1:N], function(row) {
  #   feature <- get_feature(data, row, msLevel)
  # })
  all_features <- vector("list", length = N)
  for (i in 1:N) {
    print(paste(i, '/', N))
    row <- filtered_list[[i]]
    all_features[[i]] <- get_feature(data, row, msLevel)
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(all_features)  
}

get_all_features_par <- function(data, filtered_list, msLevel, N) {
  library(foreach)
  library(doParallel)
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  start.time <- Sys.time()

  # use biocparallel instead?
  # all_features = bplapply(1:N, function(i, filtered_list, msLevel, get_feature) {
  #   row <- filtered_list[[i]]
  #   feature <- get_feature(data, row, msLevel)
  # }, filtered_list, msLevel, get_feature)
  
  all_features <- foreach(i=1:N,
                          .combine = list,
                          .multicombine = FALSE,
                          .export = 'get_feature',
                          .packages = 'xcms') %dopar% {
                            row <- filtered_list[[i]]
                            feature <- get_feature(data, row, msLevel)
                          }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(all_features)
}

### processing starts here ###

mzml_dir <- 'C:\\Users\\joewa\\Work\\docs\\clms\\FourBeers_mzML_files\\POS\\'
mzml_files <- dir(mzml_dir, full.names=TRUE)
design <- getDesign(mzml_files, "BEER")
msLevel = 1
raw_data_on_disk <- readData(mzml_files, design, mode='onDisk', msLevel=msLevel)
raw_data_in_memory <- readData(mzml_files, design, mode='inMemory', msLevel=msLevel)

# ppm value is set fairly large.
# Other parameters for peakwidth, snthresh, prefilter for justin's data, taken from 
# https://www.dropbox.com/home/Meta_clustering/ms2lda/large_study/r/beer_method_3_pos?preview=xcmsPeakPicking.R
# TODO: at the moment, findChromPeaks only works with OnDiskMSnExp, that's why we pass it raw_data_on_disk for peak picking
xdata <- pickPeaks(raw_data_on_disk, ppm=10, peakwidth=c(5, 100), snthresh=3, prefilter=c(3, 1000), mzdiff=0.001)
peaks_df <- data.frame(chromPeaks(xdata))

# Filtering parameters, taken from 
# see https://www.dropbox.com/home/Meta_clustering/ms2lda/large_study/r/config?preview=config_beer_pos_3.yml
#     rt_start: 3                         # too early (top of peak)
#     rt_end: 21                          # too late (top of peak)
#     min_MS1_intensity : 1E5 or 2E5 should do a decent job on Polyomics Q-Exactive data
#     min_MS2_intensity: 5000
# filtered_df <- filterDf(peaks_df, rt_start=3*60, rt_end=21*60, min_intensity=1E6, min_sn=100)
# sorted_df <- sortDf(filtered_df, TRUE)

# sort dataframe by integrated intensity descending
sorted_df <- sortDf(peaks_df, FALSE)
filtered_list <-  setNames(split(sorted_df, seq(nrow(sorted_df))), rownames(sorted_df)) # convert df to list

N <- length(filtered_list)
# all_features_par <- get_all_features_par(xdata, filtered_list, msLevel, N)
all_features <- get_all_features(xdata, filtered_list, msLevel, N)
dfs <- lapply(all_features, data.frame, stringsAsFactors = FALSE) # list of dataframes, one for each peak
df <- do.call("rbind", dfs) # a single big dataframe
write.csv(df, file = 'peaks.csv')
