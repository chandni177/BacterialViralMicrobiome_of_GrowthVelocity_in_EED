# calculating_delta_HAZ.R

library("tidyverse")

# read in list of analyzed samples and generate needed data
analyzedSamples <- read.delim("../documents/all_analyzed_sample_data.txt",
                              header = TRUE, sep = "\t")
analyzedSamples <- analyzedSamples %>% 
  mutate("Manary_Barcode" = str_extract(string = sample_name,
                                        pattern = "NA[:digit:]+"))

# pull in the Patient/Subject ID for each sample from "all_samples_with_HAZ.txt"
allSamplesWithHAZ <- read.delim("../documents/all_samples_with_HAZ.txt",
                                header = TRUE, sep = "\t")

subjectByBarcode <- select(allSamplesWithHAZ, Manary_Barcode, Subject_ID)

analyzedSamples <- merge(analyzedSamples, subjectByBarcode,
                         by = "Manary_Barcode", all = FALSE)

# add timepoint column
analyzedSamples <- analyzedSamples %>% 
  mutate("timepoint" = str_extract(string = sample_name,
                                   pattern = "[:alpha:]+(_|-)NA[:digit:]+"),
         "timepoint" = str_remove(string = timepoint,
                                  pattern = "(_|-)NA[:digit:]+"))

# We will later loop over analyzedSamples and get each one's delta HAZ for table 1

#----- Build Subject Class -----#
# S4 class will hold the Subject_ID, the subject's inital, midstudy, and endstudy
#   HAZ

allSamplesWithHAZ$Subject_ID <- as.character(allSamplesWithHAZ$Subject_ID)
allSamplesWithHAZ$Time_point <- as.character(allSamplesWithHAZ$Time_point)

setClass("Subject",
         representation(subject_id = "character",
                        HAZ_map = "numeric",
                        initial_delta = "numeric",
                        midstudy_delta = "numeric"),
         prototype(subject_id = NA_character_,
                   HAZ_map = NA_real_,
                   initial_delta = NA_real_,
                   midstudy_delta = NA_real_)
)

subjectClassList <- vector(mode = "list")

for (i in 1:nrow(allSamplesWithHAZ)) {
  
  currentSubjectID <- allSamplesWithHAZ[i, "Subject_ID"]
  currentTimePoint <- allSamplesWithHAZ[i, "Time_point"]
  currentHAZ <- allSamplesWithHAZ[i, "HAZ"]
  
  if (!(currentSubjectID %in% names(subjectClassList))) {
    
    newSubject <- new("Subject")
    newSubject@subject_id <- currentSubjectID
    newSubject@HAZ_map <- structure(c(NA_real_, NA_real_, NA_real_),
                                    names = c("initial", "midstudy", "endstudy"))
    
    # get the position of the next index to append the class to w/n subjectClassList
    #   as we're not adding a new class to the list on every iteration
    currentListIndex <- length(subjectClassList) + 1
    
    subjectClassList[[currentListIndex]] <- newSubject
    names(subjectClassList)[currentListIndex] <- currentSubjectID
  }
  subjectClassList[[currentSubjectID]]
  subjectClassList[[currentSubjectID]]@HAZ_map[currentTimePoint] <- currentHAZ
  
}

# compute the delta HAZ slots
for (i in 1:length(subjectClassList)) {
  subjectClassList[[i]]@initial_delta <- subjectClassList[[i]]@HAZ_map[["midstudy"]] - subjectClassList[[i]]@HAZ_map[["initial"]]
  subjectClassList[[i]]@midstudy_delta <- subjectClassList[[i]]@HAZ_map[["endstudy"]] - subjectClassList[[i]]@HAZ_map[["midstudy"]]
}

# Now add appropriate delta HAZ to analyzed samples
analyzedSamplesWithDeltas <- analyzedSamples
analyzedSamplesWithDeltas$delta_HAZ <- NA

analyzedSamplesWithDeltas$Subject_ID <- as.character(analyzedSamplesWithDeltas$Subject_ID)

for (i in 1:nrow(analyzedSamplesWithDeltas)) {
 
  subjectID <- analyzedSamplesWithDeltas[i, "Subject_ID"]
  subjectTimePoint <- analyzedSamplesWithDeltas[i, "timepoint"]
  
  deltaHAZ <- NA
  
  if (subjectTimePoint == "initial") {
    deltaHAZ <- subjectClassList[[subjectID]]@initial_delta
  } else if (subjectTimePoint == "midstudy") {
    deltaHAZ <- subjectClassList[[subjectID]]@midstudy_delta
  }
  
  analyzedSamplesWithDeltas[i, "delta_HAZ"] <- deltaHAZ
}

# double check these values match previously calculated growth status prior
analyzedCheck <- analyzedSamplesWithDeltas
analyzedCheck <- analyzedCheck %>% 
  mutate(check = case_when(delta_HAZ > 0 ~ "Adequate",
                           delta_HAZ >= -0.30 & delta_HAZ <= 0 ~ "Moderate",
                           delta_HAZ < -0.30 ~ "Poor"))

# save
write.table(analyzedSamplesWithDeltas, file = "../documents/analyzedSamplesWithDeltas.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
