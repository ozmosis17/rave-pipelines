library(tidyverse)
library(glue)
# Helper Functions
convert_range_to_vector <- function(range_string) {
  # Remove any whitespace
  range_string <- gsub("\\s", "", range_string)

  ## If empty, return empty vector
  if(range_string==""){
    return(c())
  }

  # Split the string by comma
  ranges <- strsplit(range_string, ",")[[1]]

  # Initialize an empty vector to store individual numbers
  numbers <- c()

  # Loop through each range
  for (range in ranges) {
    # Split the range by hyphen
    range_parts <- strsplit(range, "-")[[1]]

    # If only one part exists, it's a single number
    if (length(range_parts) == 1) {
      number <- as.integer(range_parts)
      if (is.na(number)) {
        warning("Invalid number detected: ", range)
        next  # Skip to the next range
      }
      numbers <- c(numbers, number)
    } else if (length(range_parts) == 2) {
      # Convert range parts to numbers
      start <- as.integer(range_parts[1])
      end <- as.integer(range_parts[2])

      # Check for NA or NaN values
      if (is.na(start) || is.na(end)) {
        warning("Invalid range detected: ", range)
        next  # Skip to the next range
      }

      # Add numbers within the range to the vector
      numbers <- c(numbers, start:end)
    } else {
      warning("Invalid range format: ", range)
    }
  }

  # if(length(numbers)==0){
  #   stop("stop")
  # }
  # if(length(numbers)==1&&is.na(numbers)){
  #   stop("stop")
  # }
  return(numbers)
}

nanpow2db <- function(y){
  # Power to dB conversion, setting negatives and zeros to NaN
  #
  # params:
  #         y: power --required
  #
  # returns:
  #         ydB: dB (with 0s and negatives set to 0)

  # Convert negative values to 0
  y <- pmax(y, 0)

  if(length(y)==1){
    if(y==0){
      return(NaN)
    } else {
      ydB <- 10*log10(y)
    }
  } else {
    y[y==0] <- NaN
    ydB <- 10*log10(y)
  }
  return(ydB)
}

pts <- dipsaus::parse_svec("1-42,50-60")
# Read patient key data
patient_key <-
  read.csv("/Users/ozhou/Library/CloudStorage/OneDrive-TexasA&MUniversity/Karas Lab/rave-pipelines/modules/karaslab_fragility/Outcome_Classification_ML/patient_data_all_rev.csv", header = TRUE, stringsAsFactors = FALSE)|>
  mutate(
    condition = gsub("\\s.*", "", condition)
  )|>
  tibble()

patient_key <- patient_key[pts,]
patient_key$subject_code

# Set the directory path
# List of directories for each frequency range
folder <- "/Users/ozhou/Library/CloudStorage/OneDrive-TexasA&MUniversity/Karas Lab/Results_FragilityEEGDataset"

note <- "norank"
csv_files <- NULL

for (subject_code in unique(patient_key$subject_code)) {
  directory <- paste0(folder,"/",subject_code,"/",note)
  if (file.exists(directory)) {
    csv_files <- append(csv_files,list.files(directory, pattern = paste0(note,"\\.csv$"), full.names = TRUE))
  } else {
    print(paste0("patient ", subject_code, " error"))
  }
}
csv_files

frag_csv_files <- csv_files[grepl("fragility",csv_files)]
frag_csv_files

# Initialize an empty list to store dataframes
all_dataframes <- list()
all_dataframes_annotation <- list()

idx <- 1
# Loop through each list of CSV files

for (file in frag_csv_files) {
  # Extract patient code and condition from file name
  file_name <- basename(file)
  patient_code <- gsub("_.*", "", file_name)
  cnd <- gsub(".*_", "", gsub("_fragility_norank.csv", "", file_name))

  # convert "seizure" in file name to "sz" condition name
  cnd <- gsub("seizure","sz",cnd)

  ## only read the data if the patient code and condition exist in the patient key
  exist <- patient_key|>dplyr::filter(subject_code==patient_code, condition==cnd)
  if(nrow(exist)==0){
    next
  }

  # Read CSV file into a dataframe
  df <- tibble(read.csv(file))|>select(-X)

  # Add patient code, condition, and frequency as columns
  all_dataframes_annotation[[idx]] <-
  tibble(
    patient_code = patient_code,
    condition = cnd
  )
  # Store dataframe in the list
  all_dataframes[[idx]] <- df

  idx <- idx + 1
}

all_dataframes_annotation <- do.call(rbind, all_dataframes_annotation)

non_exist_patient <-
  patient_key|>
  dplyr::select(subject_code, condition)|>
  unique()|>
  anti_join(all_dataframes_annotation,
            by = c("subject_code" = "patient_code",
                   "condition" = "condition"))

if(nrow(non_exist_patient)!=0){
  warning("Patient key do not match")
  patient_key_orig <- patient_key
  patient_key <- dplyr::filter(patient_key, !(subject_code %in% non_exist_patient$subject_code & condition %in% non_exist_patient$condition))
}

all_dataframes_annotation2 <-
  all_dataframes_annotation |>
  left_join(
    patient_key,
    by = c("patient_code" = "subject_code",
           "condition" = "condition")
  )


## find out the number of readings for each electrode
times <- lapply(all_dataframes, function(x) as.numeric(gsub("^X", "", colnames(x))))
minNumReadings <- min(sapply(times, length))

# Initialize patient_data dataframe with frequency column
patient_data <- list()
for(i in seq_len(nrow(all_dataframes_annotation2))){
      patient <- all_dataframes_annotation2$patient_code[i]
      condition <- all_dataframes_annotation2$condition[i]
      ## readings are limited to the last minNumReadings readings
      idx <- rev(rev(seq_len(ncol(all_dataframes[[i]])))[1:minNumReadings])
      data_over_time_per_elec <- all_dataframes[[i]][, idx]
      colnames(data_over_time_per_elec) <- paste0("X", seq_len(ncol(data_over_time_per_elec)))

      # Extract range data for electrodes, SOZ and Resect
      electrodes <- all_dataframes_annotation2$load_electrodes[i]|>
        convert_range_to_vector()

      if(length(electrodes)==0){
        next
      }
      resect <- all_dataframes_annotation2$Resect[i]|>
        convert_range_to_vector()
      soz <- all_dataframes_annotation2$SOZ[i]|>
        convert_range_to_vector()

      # Create logical vectors indicating whether each electrode is in SOZ or Resect
      resect <- ifelse(electrodes %in% resect, 1, 0)
      soz <- ifelse(electrodes %in% soz, 1, 0)
      length(resect)
      dim(data_over_time_per_elec)

      dt <- tibble(
        patient_code = patient,
        condition = condition,
        electrode = electrodes,
        resect = resect,
        soz = soz
      )|>
      bind_cols(data_over_time_per_elec)

      patient_data <- c(patient_data, list(dt))
}

patient_data <- do.call(rbind, patient_data)


Positive <- c("subpt01", "subpt2", "subpt3", "subpt8", "subpt11", "subpt13",
              "subpt17", "subpt15", "subpt16", "subummc002", "subummc005", "subummc009",
              "subjh105", "subumf001")
Negative <- c("subpt6", "subpt7", "subpt10", "subpt12", "subpt14", "subjh101", "subjh103")

# Create outcome column
patient_data$outcome <- ifelse(patient_data$patient_code %in% Positive, 1,
                               ifelse(patient_data$patient_code %in% Negative, 0, NA))

# Convert resect to factor
patient_data$resect <- as.factor(patient_data$resect)

# Only Seizure Free
patient_data_not_sz_free <- patient_data[patient_data$outcome==0, ]
patient_data <- patient_data[patient_data$outcome==1, ]

dt_final <- patient_data

dt_final_not_sz <- patient_data_not_sz_free



