# Import subject from BIDS and preprocess ----------------------------------------------

pipeline_container <- raveio::pipeline_collection(tempfile(), overwrite = TRUE)

# # you can get available inputs here:
# raveio::pipeline("import_bids")$get_settings()

subject <- "jh103"
project <- "FragilityEEGDataset"
electrodes <- "1:4,7:12,15:23,25:33,47:63,65:66,69:71,73:110"
display <- "1:2,15:21,25:33,75:78,86,94"

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
        BIDS_subject = paste0("sub-",subject),
        BIDS_runs = c(
          paste0("ses-presurgery/ieeg/sub-",subject,"_ses-presurgery_task-ictal_acq-ecog_run-01"),
          paste0("ses-presurgery/ieeg/sub-",subject,"_ses-presurgery_task-ictal_acq-ecog_run-02"),
          paste0("ses-presurgery/ieeg/sub-",subject,"_ses-presurgery_task-ictal_acq-ecog_run-03")
          #paste0("ses-presurgery/ieeg/sub-",subject,"_ses-presurgery_task-ictal_acq-ecog_run-04")
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
        import_setup__subject_code = subject,
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
#         subject_code = subject,
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
  subject_code = subject,
  notch_filter_lowerbound = c(59,118,178),
  notch_filter_upperbound = c(61,122,182)
)
notch_pipeline$run('apply_notch')

raveio::generate_reference(paste0(project,'/',subject), electrodes)
rave::start_rave2() # generate reference table

# Fragility ----------------------------------------------
# Add `path` to force using devel pipeline
fragility_pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")
#fragility_pipeline$target_table
#raveio::pipeline_visualize(fragility_pipeline$pipeline_path)
#fragility_pipeline$get_settings()

fragility_pipeline$set_settings(
  project_name = project,
  subject_code = subject,
  load_electrodes = electrodes,
  display_electrodes = display,
  reference_name = "car",
  trial_num = 1,
  t_window = 1000,
  t_step = 500
)

# for voltage reconstruction ----------------------------------------------
# env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")
vplot_data <- c(fragility_pipeline$run(c("repository", "A", "t_window", "t_step", "trial_num")),
                list('timepoints' = 1:500),
                'elec_num' = 1)
do.call(voltage_recon_plot, vplot_data)
# for plotting fragility map ---------------------------------------
# env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")

fplot_data <- c(fragility_pipeline$run(c("repository", "f_info", "displayed_elec", "sz_onset")),
                'sort_fmap' = 1, 'height' = 14)
do.call(fragility_map_plot, fplot_data)
# for plotting to pdf ---------------------------------------
# env <- fragility_pipeline$load_shared()
source("./modules/karaslab_fragility/R/shared-plots.R")

subject <- fragility_pipeline$run("subject")
export_path <- file.path(subject$note_path, "karaslab_fragility", "exports")
raveio::dir_create2(export_path)

fplot_data <- c(fragility_pipeline$run(c("repository", "f_info", "displayed_elec", "sz_onset")),
                'sort_fmap' = 1, 'height' = 14)

vplot_data <- c(fragility_pipeline$run(c("repository", "A", "t_window", "t_step", "trial_num")),
                list('timepoints' = 1:500),
                'elec_num' = 1)

pdf_path <- file.path(export_path, paste0(Sys.time(),'.pdf'))
grDevices::pdf(pdf_path)
par(mfrow=c(2,1),mar=rep(2,4))
do.call(fragility_map_plot, fplot_data)
do.call(voltage_recon_plot, vplot_data)
grDevices::dev.off()

export_pdf(
  {
    par(mfrow=c(2,1),mar=rep(2,4))
    do.call(fragility_map_plot, fplot_data)
    do.call(voltage_recon_plot, vplot_data)
  },
  path = pdf_path
)


# TODO --------------------------------------------------------------------
# clean up fragility plot function DONE
# make plot function for reconstructed voltage visualization DONE
# add function to calculate mean squared error for voltage reconstruction
