ravedash::debug_modules(module_root = rstudioapi::getActiveProject())

pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")

library(readxl)
library(stringr)

export_path <- "/Volumes/bigbrain/Fragility2024/Results_FragilityLambdaSearch"

pts <- dipsaus::parse_svec("1-22")
pipeline_xls <- read.csv("/Users/ozhou/Library/CloudStorage/OneDrive-TexasA&MUniversity/Karas Lab/patient_data_all_fragility_clean.csv")
#pipeline_xls <- readxl::read_xlsx("/Users/ozhou/Library/CloudStorage/OneDrive-TexasA&MUniversity/Karas Lab/FragilityEEGDataset_pipeline.xlsx")
pipeline_xls$subject[pts]

# not in use
# SOZ_table <- data.frame(
#   subject = stringr::str_sub(pipeline_xls$subject[pts], 4),
#   SOZ_elec = NA,
#   #SOZc_elec = NA,
#   SOZ_elec_i = NA
#   #SOZc_elec_i = NA
# )

for(i in pts){

  # subject_code <- "K009"
  # project <- "Retrostudy"
  # electrodes <- dipsaus::parse_svec("1-126")
  # display <- electrodes


  subject_code <- pipeline_xls$subject_code[i]
  project <- pipeline_xls$project_name[i]
  electrodes <- dipsaus::parse_svec(pipeline_xls$load_electrodes[i])
  display <- electrodes # display all electrodes

  epoch_name <- pipeline_xls$epoch_file_name[i]
  reference_name <- pipeline_xls$reference_name[i]
  condition <- pipeline_xls$condition[i]

  soz <- dipsaus::parse_svec(pipeline_xls$SOZ[i])
  sozc <- electrodes[!(electrodes%in%soz)]

  if(!all(c(soz,sozc) %in% electrodes)){
    warning("Not all electrodes specified in soz are loaded! Will omit unloaded electrodes from soz.")
    soz <- soz[soz %in% electrodes]
    sozc <- sozc[sozc %in% electrodes]
  }

  # display <- dipsaus::parse_svec(pipeline_xls$display_electrodes[i])
  # if(is.null(display)){
  #   display <- electrodes
  # }
  #sample_rate <- as.numeric(pipeline_xls$sample_rate[i])
  #ictal_runs <- dipsaus::parse_svec(pipeline_xls$ictal_runs[i])
  #epoch_times <- as.numeric(strsplit(pipeline_xls$epoch_times[i],",")[[1]])
  #type <- pipeline_xls$type[i]
  #import_format <- pipeline_xls$import_format[i]

  subject_check <- raveio::validate_subject(paste0(project,"/",subject_code),
                                            method = "basic", verbose = FALSE)
  subject_check$paths$data_path$valid

  #HUP dataset
  # HUP_i <- 5
  #
  # HUPxls <- readxl::read_xlsx("/Volumes/OFZ1_T7/karaslab/rave_data/bids_dir/HUPDataset/participantsHUP121923.xlsx")
  # subject_code <- stringr::str_sub(HUPxls$participant_id[HUP_i], 4)
  # project <- "HUPDataset"
  # electrodes <- dipsaus::parse_svec(HUPxls$`good electrodes`[HUP_i])
  # display <- electrodes
  #
  # ictal_runs <- seq_len(HUPxls$`ictal run`[HUP_i])
  # sample_rate <- 512
  # type <- tolower(HUPxls$implant[HUP_i])
  # import_format <- names(raveio::IMPORT_FORMATS)[3]

  # preprocess patient if not already preprocessed
  # if (!subject_check$paths$data_path$valid) {
  if (FALSE) {
    print("starting preprocessing")
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
    epoch_table <- data.frame(
      Block = subject$blocks,
      Time = epoch_times,
      Trial = seq_len(length(subject$blocks)),
      Condition = rep("sz_onset",length(subject$blocks)),
      Duration = rep(NA,length(subject$blocks))
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

  print(paste0("starting pipeline for pt: ", subject_code))
  # Fragility ----------------------------------------------
  # Add `path` to force using devel pipeline
  fragility_pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")
  #fragility_pipeline$target_table
  #raveio::pipeline_visualize(fragility_pipeline$pipeline_path)
  #fragility_pipeline$get_settings()

  # set subject object from rave
  subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))
  elec_list <- subject$get_electrode_table()

  # create export directory for this subject
  export <- file.path(export_path, subject_code)
  raveio::dir_create2(export)

  # extracts trial_num from between the parentheses in condition
  trial_num <- substr(condition, unlist(gregexpr("\\(", condition)) + 1, unlist(gregexpr("\\)", condition)) - 1)
  trial_num <- as.numeric(trial_num)

    fragility_pipeline$set_settings(
      project_name = project,
      subject_code = subject_code,
      epoch_name = epoch_name,
      epoch_time_window = c(-10,20),
      reference_name = reference_name,
      load_electrodes = electrodes,
      display_electrodes = display,
      trial_num = trial_num,
      t_window = 250,
      t_step = 125,
      sz_onset = 0,
      lambda = FALSE,
      threshold_start = 0,
      threshold_end = 20,
      threshold = 0.5,
      soz = soz,
      sozc = sozc
    )

    # # display image results in R ---------------------------------------
    # #env <- fragility_pipeline$load_shared()
    # source("./modules/karaslab_fragility/R/shared-plots.R")
    # subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))
    #
    # results <- c(fragility_pipeline$run(c("repository", "adj_frag_info","threshold_elec")))
    #
    # voltage reconstruction
    # do.call(voltage_recon_plot, c(results[1:2],
    #                               list(fragility_pipeline$get_settings("t_window"),
    #                                    fragility_pipeline$get_settings("t_step"),
    #                                    fragility_pipeline$get_settings("trial_num"),
    #                                    timepoints = 1:1000,
    #                                    elec_num = 1,
    #                                    lambda = fragility_pipeline$get_settings("lambda"))
    #                               ))
    #
    # # fragility map
    # do.call(fragility_map_plot, c(results,
    #                               list(fragility_pipeline$get_settings("display_electrodes"),
    #                                    fragility_pipeline$get_settings("sz_onset"),
    #                                    elec_list = subject$get_electrode_table(),
    #                                    'sort_fmap' = 1,
    #                                    'height' = 14)
    #                               ))

    # for plotting to pdf ---------------------------------------
    # env <- fragility_pipeline$load_shared()
    source("./modules/karaslab_fragility/R/shared-plots.R")

    tryCatch(
      error = function(e){
        if (exists(export)) {
          file.rename(export, file.path(export_path, paste0(subject_code,"_ERROR")))
        }
      },{
        results <- c(fragility_pipeline$run(c("repository", "adj_frag_info","threshold_elec")))

        # force evaluation
        #env <- c(fragility_pipeline$eval(c("repository", "adj_frag_info","threshold_elec")), shortcut = TRUE)
        #results <- list(repository = env[[1]]$repository, adj_frag_info = env[[1]]$adj_frag_info, threshold_elec = env[[1]]$threshold_elec)

        # save unranked results
        output_files(results$repository,results$adj_frag_info$frag_norank,fragility_pipeline$get_settings(),export,"norank")

        # save ranked results
        output_files(results$repository,results$adj_frag_info$frag,fragility_pipeline$get_settings(),export,"ranked")

        # print results to pdf
        pdf_path <- file.path(export, paste0(subject_code,'_seizure',trial_num/2,format(Sys.time(), "%m-%d-%Y_%H%M%S"),'.pdf'))
        grDevices::pdf(pdf_path, width = 12, height = 7)
        par(mfrow=c(2,1),mar=rep(2,4))

        # voltage reconstruction
        g <- do.call(voltage_recon_plot, c(results[1:2],
                                           list(fragility_pipeline$get_settings("t_window"),
                                                fragility_pipeline$get_settings("t_step"),
                                                fragility_pipeline$get_settings("trial_num"),
                                                timepoints = 1:1000,
                                                elec_num = 1,
                                                lambda = fragility_pipeline$get_settings("lambda"))
        ))
        print(g)

        # old fragility heatmap
        do.call(fragility_map_plot, c(results,
                                      list(fragility_pipeline$get_settings("display_electrodes"),
                                           fragility_pipeline$get_settings("sz_onset"),
                                           elec_list = elec_list,
                                           'sort_fmap' = 1,
                                           'height' = 14,
                                           threshold = fragility_pipeline$get_settings("threshold"))
        ))

        grDevices::dev.off()
    })

    # save electrodes above threshold to SOZ_table (not in use)
    # threshold_elec_i <- as.numeric(results$threshold_elec$elecnames)
    # if (!all(elec_list$Label == 'NoLabel')) {
    #   elec_i <- match(threshold_elec_i, elec_list$Electrode)
    #   threshold_elec_names <- paste0(elec_list$Label[elec_i], collapse = ", ")
    # } else {
    #   threshold_elec_names <- dipsaus::deparse_svec(threshold_elec_i)
    # }
    #
    # SOZ_table$SOZ_elec[match(paste0(subject_code),SOZ_table$subject)] <- threshold_elec_names
    # #SOZ_table$SOZc_elec[match(paste0(subject_code),SOZ_table$subject)]
    # SOZ_table$SOZ_elec_i[match(paste0(subject_code),SOZ_table$subject)] <- dipsaus::deparse_svec(threshold_elec_i)
    # #SOZ_table$SOZc_elec_i[match(paste0(subject_code),SOZ_table$subject)]

}

raveio::safe_write_csv(
  SOZ_table,
  file.path(export, paste0(project, '_fragility.csv'))
)

# export_pdf(
#   {
#     par(mfrow=c(2,1),mar=rep(2,4))
#
#     do.call(voltage_recon_plot, c(results,
#                                   list(fragility_pipeline$get_settings()$t_window,
#                                        fragility_pipeline$get_settings()$t_step,
#                                        fragility_pipeline$get_settings()$trial_num,
#                                        fragility_pipeline$get_settings()$signalScaling,
#                                        timepoints = 1:fragility_pipeline$get_settings()$t_window,
#                                        elec_num = 1)
#     ))
#
#     do.call(fragility_map_plot, c(results,
#                                   list(fragility_pipeline$get_settings()$display_electrodes,
#                                        fragility_pipeline$get_settings()$sz_onset,
#                                        elec_list = subject$get_electrode_table(),
#                                        'sort_fmap' = 1,
#                                        'height' = 14)
#     ))
#
#   },
#   path = pdf_path
# )

# DIPSAUS DEBUG START
# repository <- results$repository
# adj_frag_info <- results$adj_frag_info
# trial_num = 2
# t_window = 250
# t_step = 125
# timepoints = 1:1000
# elec_num = 1
# percentile = 0.1
# display_electrodes <- c(33,34,62:69,88:91)
