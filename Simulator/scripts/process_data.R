# needs xcms3, see https://bioconductor.org/packages/release/bioc/html/xcms.html for installation
library(xcms)

### some useful functions ###

getDesign = function(mzml_files, label) {
  design <- data.frame(sample_name = sub(basename(mzml_files), pattern = '.mzML', replacement = "", fixed = TRUE),
                       sample_group = c(rep(label, length(mzml_files))),
                       stringsAsFactors = FALSE)
  return(design)  
}

readData = function(files, design, mode="onDisk") {
  raw_data <- MSnbase::readMSData(files = files, pdata = new("NAnnotatedDataFrame", design),
                                  mode = mode, centroided. = TRUE, msLevel. = 1)
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

get_feature = function(data, row, mzppm=NULL, make_plot=FALSE) {
  sample <- row$sample
  mz_range_df <- row[, c('mzmin', 'mzmax')]
  rt_range_df <- row[, c('rtmin', 'rtmax')]

  # TODO: we need to add some padding to the mz_range to get the complete chromatogram, why?!
  # alternatively we can also have a ppm constant for the bounds
  mz = row$mz
  if (is.null(mzppm)) {
    mz_range <- c(mz_range_df$mzmin-0.0003, mz_range_df$mzmax+0.0003)
  } else {
    mz_range <- c(mz * (1-mzppm/1e6), mz * (1+mzppm/1e6))
  }
  rt_range <- c(rt_range_df$rtmin, rt_range_df$rtmax)

  # below will not give you the mz values of each data point  
  # chr_raw <- chromatogram(data, mz=mz_range, rt=rt_range)
  # chr_raw <- chr_raw[1, sample]

  # below is the correct code
  one_file = filterFile(data, sample)
  sample_name = one_file@phenoData@data$sample_name
  filtered_file = filterRt(filterMz(one_file, mz_range), rt_range)
  chrom_data = as(filtered_file, 'data.frame')
  
  if (make_plot) {
    plotMsData(chrom_data)
  }

  mz_values = chrom_data$mz
  rt_values = chrom_data$rt
  intensity_values = chrom_data$i
  
  return(list(id=rownames(row), mz=mz, mzmin=row$mzmin, mzmax=row$mzmax, 
              rt=row$rt, rtmin=row$rtmin, rtmax=row$rtmax,
              into=row$into, maxo=row$maxo, sn=row$sn, 
              mz_values=mz_values, rt_values=rt_values, intensity_values=intensity_values,
              msLevel=1, fromFile=sample, sample_name=sample_name))    
}

sortDf <- function(df, decreasing=TRUE) {
  speaks = df[order(df[,'into'], decreasing=decreasing),]
  return(speaks)
}

get_all_features <- function(data, filtered_list, N, mzppm=NULL, make_plot=FALSE) {
  start.time <- Sys.time()
  # all_features <- lapply(filtered_list[1:N], function(row) {
  #   feature <- get_feature(data, row)
  # })
  all_features <- vector("list", length = N)
  for (i in 1:N) {
    print(paste(i, '/', N))
    row <- filtered_list[[i]]
    all_features[[i]] <- get_feature(data, row, mzppm=mzppm, make_plot=make_plot)
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(all_features)  
}

get_all_features_par <- function(data, filtered_list, N, mzppm=NULL, make_plot=FALSE) {
  library(foreach)
  library(doParallel)
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  start.time <- Sys.time()

  # use biocparallel instead?
  # all_features = bplapply(1:N, function(i, filtered_list, get_feature) {
  #   row <- filtered_list[[i]]
  #   feature <- get_feature(data, row)
  # }, filtered_list, get_feature)
  
  all_features <- foreach(i=1:N,
                          .combine = list,
                          .multicombine = FALSE,
                          .export = 'get_feature',
                          .packages = 'xcms') %dopar% {
                            row <- filtered_list[[i]]
                            feature <- get_feature(data, row, mzppm=mzppm, make_plot=make_plot)
                          }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(all_features)
}

extract_peaks <- function(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff) {
  # TODO: at the moment, findChromPeaks only works with OnDiskMSnExp, that's why we pass it raw_data_on_disk for peak picking
  raw_data <- readData(mzml_files, design, mode='onDisk')
  xdata <- pickPeaks(raw_data, ppm=ppm, peakwidth=peakwidth, snthresh=snthresh, prefilter=prefilter, mzdiff=mzdiff)
  return(xdata)
}

extract_features <- function(data, mzppm, make_plot, N=NULL) {
  peaks_df <- data.frame(chromPeaks(data))
  sorted_df <- sortDf(peaks_df, TRUE)
  filtered_list <-  setNames(split(sorted_df, seq(nrow(sorted_df))), rownames(sorted_df)) # convert df to list
  
  if (is.null(N)) {
    N <- length(filtered_list)
  }
  # all_features_par <- get_all_features_par(data, filtered_list, N, mzppm=mzppm)
  all_features <- get_all_features(data, filtered_list, N, mzppm=mzppm, make_plot=make_plot)
  
  # remove features that have no chromatographic peak shapes
  selected <- all_features[sapply(all_features, function(f) length(f$mz_values)>0)]
  return(selected)
}

get_df = function(all_features) {
  dfs <- lapply(all_features, data.frame, stringsAsFactors = FALSE) # list of dataframes, one for each peak
  df <- do.call("rbind", dfs) # a single big dataframe
  return(df)  
}

write_df <- function(df, filename) {
  write.csv(df, file = gzfile(filename), row.names = FALSE)  
}

### processing starts here ###

# ppm value is set fairly large.
# Other parameters for peakwidth, snthresh, prefilter for justin's data, taken from 
# https://www.dropbox.com/home/Meta_clustering/ms2lda/large_study/r/beer_method_3_pos?preview=xcmsPeakPicking.R
ppm = 10
peakwidth = c(5, 100)
snthresh = 3
prefilter = c(3, 1000)
mzdiff = 0.001
mzppm <- 10 # the ppm value used to draw the bounding box to retrieve chromatographic peak
make_plot = FALSE

# # run for beer1pos only
# mzml_dir <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Data\\multibeers_urine_data\\beers\\fullscan'
# mzml_files <- list.files(path=mzml_dir, pattern='*.mzML', full.names=TRUE)
# mzml_files = mzml_files[1]
# design <- getDesign(mzml_files, "group1")
# 
# outfile <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Trained Models\\chromatogram_beer1pos.csv.gz'
# data <- extract_peaks(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff)
# ms1_features <- extract_features(data, mzppm, make_plot)
# df = get_df(ms1_features)
# write_df(df, outfile)

# # run for beer2pos only
# mzml_dir <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Data\\multibeers_urine_data\\beers\\fullscan'
# mzml_files <- list.files(path=mzml_dir, pattern='*.mzML', full.names=TRUE)
# mzml_files = mzml_files[2]
# design <- getDesign(mzml_files, "group1")
# 
# outfile <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Trained Models\\chromatogram_beer2pos.csv.gz'
# data <- extract_peaks(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff)
# ms1_features <- extract_features(data, mzppm, make_plot)
# df = get_df(ms1_features)
# write_df(df, outfile)

# # run for all 19 beers
# mzml_dir <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Data\\multibeers_urine_data\\beers\\fullscan'
# mzml_files <- list.files(path=mzml_dir, pattern='*.mzML', full.names=TRUE)
# design <- getDesign(mzml_files, "group1")
# 
# outfile <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Trained Models\\chromatogram_19_beers.csv.gz'
# data <- extract_peaks(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff)
# ms1_features <- extract_features(data, mzppm, make_plot)
# df = get_df(ms1_features)
# write_df(df, outfile)

# run for urine02pos only
mzml_dir <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Data\\multibeers_urine_data\\urines\\fullscan'
mzml_files <- list.files(path=mzml_dir, pattern='*.mzML', full.names=TRUE)
mzml_files = mzml_files[1]
design <- getDesign(mzml_files, "group1")

outfile <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Trained Models\\chromatogram_urine02pos.csv.gz'
data <- extract_peaks(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff)
ms1_features <- extract_features(data, mzppm, make_plot)
df = get_df(ms1_features)
write_df(df, outfile)

# run for urine03pos only
mzml_dir <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Data\\multibeers_urine_data\\urines\\fullscan'
mzml_files <- list.files(path=mzml_dir, pattern='*.mzML', full.names=TRUE)
mzml_files = mzml_files[2]
design <- getDesign(mzml_files, "group1")

outfile <- 'C:\\Users\\joewa\\University of Glasgow\\Vinny Davies - CLDS Metabolomics Project\\Trained Models\\chromatogram_urine03pos.csv.gz'
data <- extract_peaks(mzml_files, design, ppm, peakwidth, snthresh, prefilter, mzdiff)
ms1_features <- extract_features(data, mzppm, make_plot)
df = get_df(ms1_features)
write_df(df, outfile)
