ravedash::debug_modules(module_root = rstudioapi::getActiveProject())

pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")

library(readxl)
library(stringr)

pts <- dipsaus::parse_svec("1,3,5,7-23,25-26,31,35")

pipeline_xls$subject[pts]

for(i in pts)
pipeline_xls <- readxl::read_xlsx("/Volumes/OFZ1_T7/karaslab/rave_data/bids_dir/FragilityEEGDataset/FragilityEEGDataset_pipeline.xlsx")
subject_code <- stringr::str_sub(pipeline_xls$subject[i], 5)
project <- pipeline_xls$project[i]
electrodes <- dipsaus::parse_svec(pipeline_xls$good_electrodes[i])
display <- dipsaus::parse_svec(pipeline_xls$display_electrodes[i])
if(is.null(display)){
  display <- electrodes
}
sample_rate <- as.numeric(pipeline_xls$sample_rate[i])
ictal_runs <- dipsaus::parse_svec(pipeline_xls$ictal_runs[i])
epoch_times <- as.numeric(strsplit(pipeline_xls$epoch_times[i],",")[[1]])
type <- pipeline_xls$type[i]
import_format <- pipeline_xls$import_format[i]


# subject_code <- "jh103"
# project <- "FragilityEEGDataset"
# electrodes <- c(1:4,7:12,15:23,25:33,47:63,65:66,69:71,73:110)
# #display <- electrodes
# display <- c(1:2,15:21,25:33,75:78,86,94)
# ictal_runs <- 1:3
# type <- "ecog"
# import_format <- names(raveio::IMPORT_FORMATS)[4]
# sample_rate <- 1000
#
# subject_code <- "pt01"
# project <- "FragilityEEGDataset"
# electrodes <- c(1:4,7:24,26:36,42:43,46:54,56:70,72:95)
# #display <- electrodes
# display <- c(33,34,62:69,88:91)
#
# ictal_runs <- 1:4
# sample_rate <- 1000
# type <- "ecog"
# import_format <- names(raveio::IMPORT_FORMATS)[4]

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
raveio::generate_reference(paste0(project,'/',subject_code), electrodes)

`%?<-%` <- dipsaus::`%?<-%`

electrode_table <- subject$get_electrode_table()
electrode_table$LabelPrefix %?<-% "All"

reference_table <- lapply(split(electrode_table, electrode_table$LabelPrefix), function(sub) {
  data.frame(
    Electrode = sub$Electrode,
    Group = sub$LabelPrefix,
    Reference = paste0('ref-',dipsaus::deparse_svec(electrodes)),
    Type = "Common Average Reference"
  )
})

# # bipolar by LabelPrefix
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
  #row.names = FALSE <- test if this is necessary
)

# rave::start_rave2() # need to create reference.csv through rave UI

# label electrodes using channels.tsv file
channels <- read.delim(paste0(raveio::rave_directories(subject_code, project, .force_format = "BIDS")$bids_subject_path,"/ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-",type,"_run-01_channels.tsv"))
if ('name' %in% colnames(channels)){
  electrodes_table <- read.csv(paste0(subject$meta_path,'/electrodes.csv'))
  electrodes_table$Label <- channels$name[subject$electrodes]
  raveio::safe_write_csv(electrodes_table, file.path(paste0(subject$meta_path,'/electrodes.csv')), row.names = FALSE)
} else {
  showNotification('Target file is not in the proper format!')
}

# Fragility ----------------------------------------------
# Add `path` to force using devel pipeline
fragility_pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")
#fragility_pipeline$target_table
#raveio::pipeline_visualize(fragility_pipeline$pipeline_path)
#fragility_pipeline$get_settings()

fragility_pipeline$set_settings(
  project_name = project,
  subject_code = subject_code,
  epoch_name = "seizure_onset",
  epoch_time_window = c(-10,10),
  reference_name = "car",
  load_electrodes = electrodes,
  display_electrodes = display,
  trial_num = 1,
  t_window = 250,
  t_step = 125,
  sz_onset = 0,
  lambda = 0.0001,
  threshold_start = 0,
  threshold_end = 10,
  threshold = 0.5
)

# display image results ---------------------------------------
#env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")
subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))

results <- c(fragility_pipeline$run(c("repository", "adj_frag_info","threshold_elec")))

# voltage reconstruction
do.call(voltage_recon_plot, c(results[1:2],
                              list(fragility_pipeline$get_settings("t_window"),
                                   fragility_pipeline$get_settings("t_step"),
                                   fragility_pipeline$get_settings("trial_num"),
                                   timepoints = 1:1000,
                                   elec_num = 1)
                              ))

# fragility map
do.call(fragility_map_plot, c(results,
                              list(fragility_pipeline$get_settings("display_electrodes"),
                                   fragility_pipeline$get_settings("sz_onset"),
                                   elec_list = subject$get_electrode_table(),
                                   'sort_fmap' = 1,
                                   'height' = 14)
                              ))

# for plotting to pdf ---------------------------------------
# env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")

#subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))

export_path <- file.path(subject$note_path, "karaslab_fragility")
#export_path <- file.path("/Users/ozhou/Library/CloudStorage/OneDrive-TexasA&MUniversity/Karas Lab/Fragility 2.0 Project/FragilityEEGDataset/")
raveio::dir_create2(export_path)

results <- c(fragility_pipeline$run(c("repository", "adj_frag_info","threshold_elec")))

pdf_path <- file.path(export_path, paste0(subject$subject_code,'_',format(Sys.time(), "%m-%d-%Y_%H%M%S"),'.pdf'))
grDevices::pdf(pdf_path, width = 12, height = 7)
par(mfrow=c(2,1),mar=rep(2,4))
do.call(voltage_recon_plot, c(results[1:2],
                              list(fragility_pipeline$get_settings("t_window"),
                                   fragility_pipeline$get_settings("t_step"),
                                   fragility_pipeline$get_settings("trial_num"),
                                   timepoints = 1:1000,
                                   elec_num = 1)
))
do.call(fragility_map_plot, c(results,
                              list(fragility_pipeline$get_settings("display_electrodes"),
                                   fragility_pipeline$get_settings("sz_onset"),
                                   elec_list = subject$get_electrode_table(),
                                   'sort_fmap' = 1,
                                   'height' = 14)
))
grDevices::dev.off()

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
# trial_num = 1
# t_window = 250
# t_step = 125
# timepoints = 1:1000
# elec_num = 1
# percentile = 0.1
# display_electrodes <- c(33,34,62:69,88:91)
