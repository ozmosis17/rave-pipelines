library(readxl)
preprocess_xls <- readxl::read_xlsx("/Users/ozhou/Downloads/Multipatient/NIHEEGDataset_pipeline_update_072224.xlsx")

pts <- 1:19
preprocess_xls$subject[pts]

# preprocess patient if not already preprocessed
print("starting preprocessing")

for (i in pts) {
  subject_code <- preprocess_xls$subject[i]
  project <- preprocess_xls$project[i]
  electrodes <- dipsaus::parse_svec(preprocess_xls$good_electrodes[i])
  sample_rate <- as.numeric(preprocess_xls$sample_rate[i])
  ictal_runs <- dipsaus::parse_svec(preprocess_xls$ictal_runs[i])
  epoch_times <- as.numeric(strsplit(preprocess_xls$epoch_times[i],",")[[1]])
  type <- preprocess_xls$type[i]
  import_format <- preprocess_xls$import_format[i]

  # Import subject from BIDS ----------------------------------------------

  pipeline_container <- raveio::pipeline_collection(tempfile(), overwrite = TRUE)

  # # you can get available inputs here:
  # raveio::pipeline("import_bids")$get_settings()

  import_bids <- pipeline_container$add_pipeline(
    "import_bids", standalone = TRUE,
    pre_hook = function(inputs, ...) {
      dipsaus::list_to_fastmap2(
        map = inputs,
        list(
          # Set Inputs Here
          overwrite = TRUE,
          # backup = FALSE,
          BIDS_dataset = project,
          BIDS_subject = paste0("sub-",subject_code),
          BIDS_runs = c(
            paste0("ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-",type,"_run-0",ictal_runs)
          ),
          BIDS_sessions = NULL
        )
      )
    }
  )

  # # you can get available inputs here:
  # raveio::pipeline("import_lfp_native")$get_settings()

  pipeline_container$add_pipeline(
    "import_lfp_native", deps = import_bids$id, standalone = TRUE,
    pre_hook = function(inputs, ...) {
      dipsaus::list_to_fastmap2(
        map = inputs,
        list(
          skip_validation = FALSE,
          import_setup__subject_code = subject_code,
          import_setup__project_name = project,
          import_channels__unit = "NA",
          import_channels__sample_rate = sample_rate,
          import_channels__electrodes = electrodes,
          # ses-presurgery_task-ictal_acq-ecog_run-01 -> presurgery_ictal_ecog_01
          import_blocks__session_block = c(
            # Now I can import just like normal RAVE subject
            paste0("presurgery_ictal_",type,"_0",ictal_runs)
          ),

          # names(raveio::IMPORT_FORMATS)[c(1:5,7)]
          import_blocks__format = import_format,
          force_import = TRUE
        )
      )
    }
  )

  # Import subject from raw_data -----------------------

  # pipeline_container <- raveio::pipeline_collection(tempfile(), overwrite = TRUE)
  #
  # pipeline_container$add_pipeline(
  #   "import_lfp_native", deps = import_bids$id, standalone = TRUE,
  #   pre_hook = function(inputs, ...) {
  #     dipsaus::list_to_fastmap2(
  #       map = inputs,
  #       list(
  #         skip_validation = FALSE,
  #         import_setup__subject_code = subject_code,
  #         import_setup__project_name = project,
  #         import_channels__unit = "NA",
  #         import_channels__sample_rate = 1000,
  #         import_channels__electrodes = electrodes,
  #         # ses-presurgery_task-ictal_acq-ecog_run-01 -> presurgery_ictal_ecog_01
  #         import_blocks__session_block = c(
  #           # Now I can import just like normal RAVE subject
  #           "presurgery_ictal_ecog_01",
  #           "presurgery_ictal_ecog_02",
  #           "presurgery_ictal_ecog_03",
  #           "presurgery_ictal_ecog_04"
  #         ),
  #
  #         # names(raveio::IMPORT_FORMATS)[c(1:5,7)]
  #         import_blocks__format = "Single BrainVision file (.vhdr+.eeg, .vhdr+.dat) per block",
  #         force_import = TRUE
  #       )
  #     )
  #   }
  # )

  # Preprocess -------------
  # build to view the order of pipelines
  pipeline_container$build_pipelines()

  # Workaround
  scheduler <- pipeline_container$get_scheduler()
  scheduler$set_settings(dry_run = FALSE, collection_root_path = "../../",
                         error_action = "error")
  scheduler$eval(scheduler$target_table$Names)
  # pipeline_container$run()

  subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))

  # run notch filter @ 60, 120, 180 Hz
  notch_pipeline <- raveio::pipeline("notch_filter")
  notch_pipeline$set_settings(
    project_name = project,
    subject_code = subject_code,
    notch_filter_lowerbound = c(59,118,178),
    notch_filter_upperbound = c(61,122,182)
  )
  notch_pipeline$run('apply_notch')

  # generate epoch
  sz_trials <- subject$blocks[ictal_runs]

  epoch_table <- data.frame(
    Block = sz_trials,
    Time = epoch_times,
    Trial = seq_len(length(sz_trials)),
    Condition = paste0(rep("sz",length(sz_trials)),ictal_runs),
    Duration = rep(NA,length(sz_trials))
  )

  raveio::safe_write_csv(
    epoch_table,
    file.path(subject$meta_path, "epoch_seizure_onset.csv")
  )

  # generate reference
  # rave::start_rave2() # alternatively, create reference.csv through rave UI

  raveio::generate_reference(paste0(project,'/',subject_code), electrodes)

  `%?<-%` <- dipsaus::`%?<-%`

  electrode_table <- subject$get_electrode_table()
  electrode_table$LabelPrefix %?<-% "All"

  # common average reference
  reference_table <- lapply(split(electrode_table, electrode_table$LabelPrefix), function(sub) {
    data.frame(
      Electrode = sub$Electrode,
      Group = sub$LabelPrefix,
      Reference = paste0('ref-',dipsaus::deparse_svec(electrodes)),
      Type = "Common Average Reference"
    )
  })

  # # bipolar reference
  # reference_table <- lapply(split(electrode_table, electrode_table$LabelPrefix), function(sub) {
  #   reference <- c(sprintf("ref_%d", sub$Electrode[-1]), "noref")
  #   data.frame(
  #     Electrode = sub$Electrode,
  #     Group = sub$LabelPrefix,
  #     Reference = reference,
  #     Type = "Bipolar Reference"
  #   )
  # })

  reference_table <- do.call("rbind", reference_table)
  rownames(reference_table) <- NULL

  # save as reference table
  raveio::safe_write_csv(
    reference_table,
    file.path(subject$meta_path, "reference_car.csv")
  )

  # label electrodes using channels.tsv file
  channels <- read.delim(paste0(raveio::rave_directories(subject_code, project, .force_format = "BIDS")$bids_subject_path,"/ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-",type,"_run-01_channels.tsv"))
  if ('name' %in% colnames(channels)){
    electrodes_table <- read.csv(paste0(subject$meta_path,'/electrodes.csv'))
    electrodes_table$Label <- channels$name[subject$electrodes]
    raveio::safe_write_csv(electrodes_table, file.path(paste0(subject$meta_path,'/electrodes.csv')), row.names = FALSE)
  } else {
    showNotification('Target file is not in the proper format!')
  }
}
