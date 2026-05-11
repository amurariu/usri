# running aldex2 on all tiyani dataset combinations

# Scott Dos Santos
# Last edited: 8th May 2026

#################################### setup ####################################

library(dplyr)
library(ALDEx2)

# get data directory containing all tiyani data, as well as all count table and
# metadata .Rda objects [CHANGE THIS TO REFLECT YOUR FILE PATH]
path.data <- "~/Documents/GitHub/usri/data/tiyani_pairs/"
path.count <- list.files(path = path.data, pattern = "tiyani_.*\\.ft\\.Rda")
path.mdata <- list.files(path = path.data, pattern = "tiyani_.*\\.meta\\.Rda")

#################################### aldex2 ####################################

# make lists to store aldex2 outputs at increasing gamma values
tiyani.aldex2.0 <- list()
tiyani.aldex2.1 <- list()
tiyani.aldex2.2 <- list()
tiyani.aldex2.3 <- list()
tiyani.aldex2.4 <- list()
tiyani.aldex2.5 <- list()

# loop over all pairwise combinations of the data and run ALDEx2 at a gamma of
# virtually zero, 0.1, 0.2, 0.3, 0.4, and 0.5, not using any randomisation of
# groups or binomial thinning (i.e. just the original data, as-is)
if(length(list.files(path.data, pattern = "aldex2")) != 6){
  for(i in 1:length(path.count)){
    
    load(paste0(path.data, path.count[i]))
    load(paste0(path.data, path.mdata[i]))
    
    combo <- gsub("\\.ft", "", ls(pattern = "ft"))
    
    assign(x = "ft", value = get(ls(pattern = "ft")))
    assign(x = "meta", value = get(ls(pattern = "meta")))
    
    
    # gamma = 0 (virtually)
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0")
    clr.0 <- aldex.clr(reads = ft, conds = meta$group, mc.samples = 128,
                       denom = "all", gamma = 1e-3, verbose = TRUE)
    
    clr.0.e <- aldex.effect(clr = clr.0, verbose = TRUE, include.sample.summary = FALSE)
    
    clr.0.t <- aldex.ttest(clr = clr.0, verbose = TRUE)
    
    tiyani.aldex2.0[[combo]] <- cbind(clr.0.e, clr.0.t)
    message("Finished running ALDEx2 on: ", combo, " at gamma = 0\n\n")
    
    
    
    # gamma = 0.1
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.1")
    clr.1 <- aldex.clr(reads = ft, conds = meta$group, mc.samples = 128,
                       denom = "all", gamma = 0.1, verbose = TRUE)
    
    clr.1.e <- aldex.effect(clr = clr.1, verbose = TRUE, include.sample.summary = FALSE)
    
    clr.1.t <- aldex.ttest(clr = clr.1, verbose = TRUE)
    
    tiyani.aldex2.1[[combo]] <- cbind(clr.1.e, clr.1.t)
    message("Finished running ALDEx2 on: ", combo, " at gamma = 0.1\n\n")
    
    
    
    # gamma = 0.2
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.2")
    clr.2 <- aldex.clr(reads = ft, conds = meta$group, mc.samples = 128,
                       denom = "all", gamma = 0.2, verbose = TRUE)
    
    clr.2.e <- aldex.effect(clr = clr.2, verbose = TRUE, include.sample.summary = FALSE)
    
    clr.2.t <- aldex.ttest(clr = clr.2, verbose = TRUE)
    
    tiyani.aldex2.2[[combo]] <- cbind(clr.2.e, clr.2.t)
    message("Finished running ALDEx2 on: ", combo, " at gamma = 0.2\n\n")
    
    
    
    # gamma = 0.3
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.3")
    clr.3 <- aldex.clr(reads = ft, conds = meta$group, mc.samples = 128,
                       denom = "all", gamma = 0.3, verbose = TRUE)
    
    clr.3.e <- aldex.effect(clr = clr.3, verbose = TRUE, include.sample.summary = FALSE)
    
    clr.3.t <- aldex.ttest(clr = clr.3, verbose = TRUE)
    
    tiyani.aldex2.3[[combo]] <- cbind(clr.3.e, clr.3.t)
    message("Finished running ALDEx2 on: ", combo, " at gamma = 0.3\n\n")
    
    
    
    # gamma = 0.4
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.4")
    clr.4 <- aldex.clr(reads = ft, conds = meta$group, mc.samples = 128,
                       denom = "all", gamma = 0.4, verbose = TRUE)
    
    clr.4.e <- aldex.effect(clr = clr.4, verbose = TRUE, include.sample.summary = FALSE)
    
    clr.4.t <- aldex.ttest(clr = clr.4, verbose = TRUE)
    
    tiyani.aldex2.4[[combo]] <- cbind(clr.4.e, clr.4.t)
    message("Finished running ALDEx2 on: ", combo, " at gamma = 0.4\n\n")
    
    
    
    # gamma = 0.5
    set.seed(2026 + i)
    message("Running ALDEx2 on: ", combo, " at gamma = 0.5")
    clr.5 <- aldex.clr(reads = ft, conds = meta$group, mc.samples = 128,
                       denom = "all", gamma = 0.5, verbose = TRUE)
    
    clr.5.e <- aldex.effect(clr = clr.5, verbose = TRUE, include.sample.summary = FALSE)
    
    clr.5.t <- aldex.ttest(clr = clr.5, verbose = TRUE)
    
    tiyani.aldex2.5[[combo]] <- cbind(clr.5.e, clr.5.t)
    message("Finished running ALDEx2 on: ", combo, " at gamma = 0.5\n\n")
    
    
    # remove all temporary objects
    rm(list = c(ls(pattern = "clr"),
                ls(pattern = "ft"),
                ls(pattern = "meta")))
  }
  
  # save results object
  for(i in ls(pattern = "aldex2")){
    
    j <- gsub("tiyani.aldex2.", "", i)
    save(list = i, file = paste0(path.data, "results_aldex2.", j, ".Rda"))
    
  }
}
