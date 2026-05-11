# running aldex3 on all tiyani dataset combinations

# Scott Dos Santos
# Last edited: 8th May 2026

#################################### setup ####################################

library(dplyr)
library(ALDEx3)

# get data directory containing all tiyani data, as well as all count table and
# metadata .Rda objects [CHANGE THIS TO REFLECT YOUR FILE PATH]
path.data <- "~/Documents/GitHub/usri/data/tiyani_pairs/"
path.count <- list.files(path = path.data, pattern = "tiyani_.*\\.ft\\.Rda")
path.mdata <- list.files(path = path.data, pattern = "tiyani_.*\\.meta\\.Rda")

#################################### aldex2 ####################################

# make lists to store aldex3 outputs at increasing gamma values
tiyani.aldex3.0 <- list()
tiyani.aldex3.1 <- list()
tiyani.aldex3.2 <- list()
tiyani.aldex3.3 <- list()
tiyani.aldex3.4 <- list()
tiyani.aldex3.5 <- list()

# loop over all pairwise combinations of the data and run ALDEx3 at a gamma of
# virtually zero, 0.1, 0.2, 0.3, 0.4, and 0.5, not using any randomisation of
# groups or binomial thinning (i.e. just the original data, as-is)
if(length(list.files(path.data, pattern = "aldex3")) != 6){
  for(i in 1:length(path.count)){
    
    load(paste0(path.data, path.count[i]))
    load(paste0(path.data, path.mdata[i]))
    
    combo <- gsub("\\.ft", "", ls(pattern = "ft"))
    
    assign(x = "ft", value = get(ls(pattern = "ft")))
    assign(x = "meta", value = get(ls(pattern = "meta")))
    
    # gamma = 0 (virtually)
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0")
    clr.0 <- aldex(Y = ft, X = ~group, data = meta, nsample = 256, scale = clr.sm, gamma = 1e-3)
    tiyani.aldex3.0[[combo]] <- summary(clr.0)
    message("Finished running ALDEx3 on: ", combo, " at gamma = 0\n\n")
    
    # gamma = 0.1
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.1")
    clr.1 <- aldex(Y = ft, X = ~group, data = meta, nsample = 256, scale = clr.sm, gamma = 0.1)
    tiyani.aldex3.1[[combo]] <- summary(clr.1)
    message("Finished running ALDEx3 on: ", combo, " at gamma = 0.1\n\n")
    
    # gamma = 0.2
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.2")
    clr.2 <- aldex(Y = ft, X = ~group, data = meta, nsample = 256, scale = clr.sm, gamma = 0.2)
    tiyani.aldex3.2[[combo]] <- summary(clr.2)
    message("Finished running ALDEx3 on: ", combo, " at gamma = 0.2\n\n")
    
    # gamma = 0.3
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.3")
    clr.3 <- aldex(Y = ft, X = ~group, data = meta, nsample = 256, scale = clr.sm, gamma = 0.3)
    tiyani.aldex3.3[[combo]] <- summary(clr.3)
    message("Finished running ALDEx3 on: ", combo, " at gamma = 0.3\n\n")
    
    # gamma = 0.4
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.4")
    clr.4 <- aldex(Y = ft, X = ~group, data = meta, nsample = 256, scale = clr.sm, gamma = 0.4)
    tiyani.aldex3.4[[combo]] <- summary(clr.4)
    message("Finished running ALDEx3 on: ", combo, " at gamma = 0.4\n\n")
    
    # gamma = 0.5
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.5")
    clr.5 <- aldex(Y = ft, X = ~group, data = meta, nsample = 256, scale = clr.sm, gamma = 0.5)
    tiyani.aldex3.5[[combo]] <- summary(clr.5)
    message("Finished running ALDEx3 on: ", combo, " at gamma = 0.5\n\n")
    
    # remove all temporary objects
    rm(list = c(ls(pattern = "clr"),
                ls(pattern = "ft"),
                ls(pattern = "meta")))
  }
  
  # save results object
  for(i in ls(pattern = "aldex3")){
    
    j <- gsub("tiyani.aldex3.", "", i)
    save(list = i, file = paste0(path.data, "results_aldex3.", j, ".Rda"))
    
  }
}
