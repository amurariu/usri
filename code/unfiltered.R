library(ALDEx2,warn.conflicts = F) #do not load aldex2 and aldex3 at the same time
library(seqgendiff, warn.conflicts=F)

####### yeast #######
##setting up data
url <- "https://raw.githubusercontent.com/amurariu/usri/main/data/transcriptome.tsv"
yeast <- read.table(url, header=T, row.names=1)

# remove the one gene with 0 reads
yeast <- yeast[rownames(yeast) != "YOR072W-B",]

# Gierlinski:2015aa
yeast[,c('SNF2.6', 'SNF2.13','SNF2.25','SNF2.35')] <- NULL
yeast[,c('WT.21','WT.22','WT.25','WT.28','WT.34','WT.36')] <- NULL

##setting up conditions vector
conditions_y <- c(rep('S', 44), rep('W', 42))
yeast.conds <- data.frame(conditions_y)

# convert to matrix and set gamma
yeast <- as.matrix(yeast)
gamma <- 1e-3

# make objects to hold data as in code
thin.data.out.aldex <- list()
data.out.aldex.r <- list()
data.out.aldex.t <- list()

# set seed based on number chosen in fig4 code
i <- 68
seed <- 20 + i
set.seed(seed)

# run code as in function
message(paste0("\nLoop iteration: ", i))

thin <- thin_2group(yeast, prop_null=0.95, alpha=0,
                    signal_fun = stats::rnorm, 
                    signal_params = list(mean = 0, sd = 2))
thin.data.out.aldex[[i]] <- thin
conds_th <- as.vector(thin$designmat)   # permuted and thinned conditions and data
data_th <- thin$mat

#randomized without FP addition 
message("Running ALDEx2 (scale = ", gamma, ") on original data with randomised groups...")
aldex.r <- aldex(yeast, conditions=conds_th, gamma = gamma) #uses original dataset but permuted conditions
data.out.aldex.r[[i]] <- aldex.r

#randomized with FP addition 
message("Running ALDEx2 (scale = ", gamma, ") on thinned data with randomised groups...")
aldex.t <- aldex(data_th, conditions=conds_th, gamma = gamma) #uses new dataset with permuted conditions
data.out.aldex.t[[i]] <- aldex.t

yeastUnf.data_0.aldex2 <- aldex.t
save(yeastUnf.data_0.aldex2, file = "~/Documents/GitHub/usri/data/yeastUnf.A2g0.Rda")







######## immuno ########

#pull the data 
raw_counts_immuno <- 'https://raw.githubusercontent.com/amurariu/usri/main/data/imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv'
immuno<-read.table(file=raw_counts_immuno, header = T, skip=35, sep='\t', row.names = 1)

# conv to matrix and set gamma
immuno <- as.matrix(immuno)
gamma <- 1e-3

# make objects to hold data as in code
thin.data.out.aldex <- list()
data.out.aldex.r <- list()
data.out.aldex.t <- list()

# set seed based on number chosen in fig4 code
i <- 68
seed <- 20 + i
set.seed(seed)

# run code as in function
message(paste0("\nLoop iteration: ", i))

thin <- thin_2group(immuno, prop_null=0.95, alpha=0,
                    signal_fun = stats::rnorm, 
                    signal_params = list(mean = 0, sd = 2))
thin.data.out.aldex[[i]] <- thin
conds_th <- as.vector(thin$designmat)   # permuted and thinned conditions and data
data_th <- thin$mat

#randomized without FP addition 
message("Running ALDEx2 (scale = ", gamma, ") on original data with randomised groups...")
aldex.r <- aldex(immuno, conditions=conds_th, gamma = gamma) #uses original dataset but permuted conditions
data.out.aldex.r[[i]] <- aldex.r

#randomized with FP addition 
message("Running ALDEx2 (scale = ", gamma, ") on thinned data with randomised groups...")
aldex.t <- aldex(data_th, conditions=conds_th, gamma = gamma) #uses new dataset with permuted conditions
data.out.aldex.t[[i]] <- aldex.t

immunoUnf.data_0.aldex2 <- aldex.t
save(immunoUnf.data_0.aldex2, file = "~/Documents/GitHub/usri/data/immunoUnf.A2g0.Rda")
