subject_code <- "jh103"
project <- "FragilityEEGDataset"
electrodes <- c(1:4,7:12,15:23,25:33,47:63,65:66,69:71,73:110)
#display <- electrodes
display <- c(1:2,15:21,25:33,75:78,86,94)
subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))

subject_code <- "pt01"
project <- "FragilityEEGDataset"
electrodes <- c(1:4,7:24,26:36,42:43,46:54,56:70,72:95)
#display <- electrodes
display <- c(33,34,62:69,88:91)
subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))

# Import subject_code from BIDS and preprocess ----------------------------------------------

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
        BIDS_dataset = "FragilityEEGDataset",
        BIDS_subject = paste0("sub-",subject_code),
        BIDS_runs = c(
          paste0("ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-ecog_run-01"),
          paste0("ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-ecog_run-02"),
          paste0("ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-ecog_run-03")
          #paste0("ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-ecog_run-04")
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
        import_channels__sample_rate = 1000,
        import_channels__electrodes = electrodes,
        # ses-presurgery_task-ictal_acq-ecog_run-01 -> presurgery_ictal_ecog_01
        import_blocks__session_block = c(
          # Now I can import just like normal RAVE subject
          "presurgery_ictal_ecog_01",
          "presurgery_ictal_ecog_02",
          "presurgery_ictal_ecog_03"
          # "presurgery_ictal_ecog_04"
        ),

        # names(raveio::IMPORT_FORMATS)[c(1:5,7)]
        import_blocks__format = "Single BrainVision file (.vhdr+.eeg, .vhdr+.dat) per block",
        force_import = TRUE
      )
    )
  }
)

# pipeline_container$add_pipeline(
#   "notch_filter", deps = import_lfp$id, standalone = TRUE,
#   pre_hook = function(inputs, ...) {
#     dipsaus::list_to_fastmap2(
#       map = inputs,
#       list(
#         project_name = project,
#         subject_code = subject_code,
#         notch_filter_lowerbound = c(59,118,178),
#         notch_filter_upperbound = c(61,122,182)
#       )
#     )
#   }
# )

# build to view the order of pipelines
pipeline_container$build_pipelines()

pipeline_container$run()

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
  Time = c("59.892",
           "59.892",
           "59.892"),
  Trial = seq_len(length(subject$blocks)),
  Condition = c("SZ EVENT # (PB SZ)",
                "SZ EVENT # (PB SZ)",
                "SZ EVENT # (PB SZ)"),
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
channels <- read.delim(paste0("/Volumes/OFZ1_T7/karaslab/rave_data/bids_dir/FragilityEEGDataset/sub-",subject_code,"/ses-presurgery/ieeg/sub-",subject_code,"_ses-presurgery_task-ictal_acq-ecog_run-01_channels.tsv"))
if ('name' %in% colnames(channels)){
  electrodes_table <- read.csv(paste0(subject$meta_path,'/electrodes.csv'))
  electrodes_table$Label <- channels$name[subject$electrodes]
  raveio::safe_write_csv(electrodes_table, file.path(paste0(subject$meta_path,'/electrodes.csv')), row.names = FALSE)
} else {
  showNotification('Uploaded file is not in the proper format!')
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
  nlambda = 16,
  ncores = NA
)

# for voltage reconstruction ----------------------------------------------
# env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")
vplot_data <- c(fragility_pipeline$run(c("repository", "A", "t_window", "t_step", "trial_num")),
                list('timepoints' = 1:250),
                'elec_num' = 1)
do.call(voltage_recon_plot, vplot_data)
# for plotting fragility map ---------------------------------------
# env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")

fplot_data <- c(fragility_pipeline$run(c("repository", "f_info", "displayed_elec", "sz_onset")),
                list(elec_list = subject$get_electrode_table(), 'sort_fmap' = 1, 'height' = 14))
do.call(fragility_map_plot, fplot_data)
# for plotting to pdf ---------------------------------------
# env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")

subject <- fragility_pipeline$run("subject")
#export_path <- file.path(subject$note_path, "karaslab_fragility")
export_path <- file.path("/Users/ozhou/Library/CloudStorage/OneDrive-TexasA&MUniversity/Karas Lab/Fragility 2.0 Project/FragilityEEGDataset/")
raveio::dir_create2(export_path)

fplot_data <- c(fragility_pipeline$run(c("repository", "f_info", "displayed_elec", "sz_onset")),
                list(elec_list = subject$get_electrode_table(), 'sort_fmap' = 1, 'height' = 14))

vplot_data <- c(fragility_pipeline$run(c("repository", "A", "t_window", "t_step", "trial_num")),
                list('timepoints' = 1:500),
                'elec_num' = 1)

pdf_path <- file.path(export_path, paste0(subject$subject_code,'_',format(Sys.time(), "%m-%d-%Y_%H%M%S"),'.pdf'))
grDevices::pdf(pdf_path, width = 12, height = 7)
par(mfrow=c(2,1),mar=rep(2,4))
do.call(fragility_map_plot, fplot_data)
do.call(voltage_recon_plot, vplot_data)
grDevices::dev.off()

# export_pdf(
#   {
#     par(mfrow=c(2,1),mar=rep(2,4))
#     do.call(fragility_map_plot, fplot_data)
#     do.call(voltage_recon_plot, vplot_data)
#   },
#   path = pdf_path
# )


# TODO --------------------------------------------------------------------
# clean up fragility plot function DONE
# make plot function for reconstructed voltage visualization DONE
# add function to calculate mean squared error for voltage reconstruction
