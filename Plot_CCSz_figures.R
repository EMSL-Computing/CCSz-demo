library(rhdf5)
library(tidyverse)
library(ggplot2)

# Set the base and output paths
basePath <- "E:/mza/CCSz-demo" # Update this path
outputPath <- file.path(basePath, "R_OutputFigures")
if (!dir.exists(outputPath)) {
  dir.create(outputPath)
}

mzafiles <- c("DT.mza", "DT-linear_CCSz.mza", "SLIM.mza", "SLIM-2nd-order_CCSz.mza", "SLIM-3rd-order_CCSz.mza")
myYtitles <- c("DTIMS AT (ms)", "DTIMS CCS/z", "SLIM AT (ms)", "SLIM 2nd order CCS/z", "SLIM 3rd order CCS/z")

tics <- list()

for (f in seq_along(mzafiles)) {
  fullPath <- file.path(basePath, "RawDataMza", mzafiles[f])
  print(fullPath)
  mza <- H5Fopen(fullPath)
  
  # Reading Metadata
  metadata <- as.data.frame(h5read(mza, "Metadata"))
  metadata <- subset(metadata, IonMobilityBin > 0)
  
  full_mz <- h5read(mza, "Full_mz_array")
  
  points <- tibble(mz=double(), intensity=double(), at=double())
  
  for (k in 1:(nrow(metadata) - 1)) {
    scan <- metadata$Scan[k]
    mzaPath <- as.character(metadata$MzaPath[k])
    mzbins <- h5read(mza, paste0("Arrays_mzbin", mzaPath, "/", scan))
    mz_array <- full_mz[mzbins + 1] # +1 for R's 1-indexing
    intensity_array <- h5read(mza, paste0("Arrays_intensity", mzaPath, "/", scan))
    
    new_points <- tibble(
      mz = mz_array[-length(mz_array)], # Exclude last element to match Python code
      intensity = intensity_array[-length(intensity_array)], # Exclude last element
      at = rep(metadata$IonMobilityTime[k], length(mz_array) - 1)
    )
    
    points <- bind_rows(points, new_points)
  }
  
  H5Fclose(mza)
  
  # Filtering based on intensity
  threshold <- max(points$intensity) * 0.008
  points <- filter(points, intensity > threshold)
  
  # Plotting
  p <- ggplot(points, aes(x=mz, y=at, color=log10(intensity))) +
    geom_point(size=1.5, stroke=0.2) +
    scale_color_viridis_c() +
    labs(x="m/z", y=myYtitles[f], color="log10(Intensity)") +
    theme_bw()
  
  ggsave(file.path(outputPath, paste0(gsub(".mza", "", mzafiles[f]), ".png")), plot=p)
  
  if (grepl("CCS", myYtitles[f])) {
    # Grouping by 'at' and summing 'intensity'
    dftic <- aggregate(intensity ~ at, data = points, FUN = sum)
    # Normalizing 'intensity'
    dftic$intensity <- dftic$intensity / max(dftic$intensity)
    dftic$method <- myYtitles[f]
    # Appending dftic to 'tics'
    tics[[f]] <- dftic
  }
}

# Combine all TICs and plot
tics_df <- bind_rows(tics)
tics_df <- filter(tics_df, at > 180)

p_tic <- ggplot(tics_df, aes(x=at, y=intensity, color=method)) +
  geom_line() +
  labs(x="CCS/z", y="Normalized total ion intensity", color="Method") +
  theme_bw()

ggsave(file.path(outputPath, "TICs.png"), plot=p_tic)
