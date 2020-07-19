### Combined code for Seedling Drought Experiment



# Read data from Github

path <- "https://raw.github.com/bobmuscarella/Luquillo_LTER_Seedling_Drought_Experiment/master/Data/"

growth_url <- paste0(path, "LUQ_DroughtExp_Seedling_growth.csv")
survive_url <- paste0(path, "LUQ_DroughtExp_Seedling_survival.csv")
traits_url <- paste0(path, "LUQ_DroughtExp_Seedling_traits.csv")
photo_url <- paste0(path, "LUQ_DroughtExp_Survey_photosynthesis.csv")

growth <- read.csv(growth_url)
survival <- read.csv(survive_url)
traits <- read.csv(traits_url)
photo <- read.csv(photo_url)









