#ravedash::debug_modules(module_root = rstudioapi::getActiveProject())

pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")

library(readxl)
library(stringr)

#pts <- dipsaus::parse_svec("1,3,5,7-23,25-26,31,35")
pts <- dipsaus::parse_svec("1")

pipeline_xls <- readxl::read_xlsx("/Volumes/bigbrain/Fragility2024/FragilityEEGDataset_pipeline_update_050224.xlsx")
#pipeline_xls <- readxl::read_xlsx("/Volumes/OFZ1_T7/karaslab/rave_data/bids_dir/FragilityEEGDataset/FragilityEEGDataset_pipeline.xlsx")
pipeline_xls$subject[pts]


for(i in pts){
#i=1
  # subject_code <- "K009"
  # project <- "Retrostudy"
  # electrodes <- dipsaus::parse_svec("1-126")
  # display <- electrodes

#  i=4

  subject_code <- pipeline_xls$subject[i]
  project <- pipeline_xls$project[i]
  electrodes <- dipsaus::parse_svec(pipeline_xls$good_electrodes[i])
  display <- electrodes # display all electrodes
  # display <- dipsaus::parse_svec(pipeline_xls$display_electrodes[i])
  # if(is.null(display)){
  #   display <- electrodes
  # }
  sample_rate <- as.numeric(pipeline_xls$sample_rate[i])
  ictal_runs <- dipsaus::parse_svec(pipeline_xls$ictal_runs[i])
  epoch_times <- as.numeric(strsplit(pipeline_xls$epoch_times[i],",")[[1]])
  type <- pipeline_xls$type[i]
  import_format <- pipeline_xls$import_format[i]

  subject_check <- raveio::validate_subject(paste0(project,"/",subject_code),
                                            method = "basic", verbose = FALSE)

  print(paste0("starting pipeline for pt: ", subject_code))
  # Fragility ----------------------------------------------
  # Add `path` to force using devel pipeline
  fragility_pipeline <- raveio::pipeline("karaslab_fragility", paths = "/Volumes/bigbrain/Fragility2024/rave_pipeline/modules/")
  #fragility_pipeline$target_table
  #raveio::pipeline_visualize(fragility_pipeline$pipeline_path)
  #fragility_pipeline$get_settings()

  # set subject object from rave
  subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))
  elec_list <- subject$get_electrode_table()

  # create export directory for this subject
  export_path <- file.path("/Volumes/bigbrain/Fragility2024/FragilityResultsACLNoNorm", subject_code)
  # export_path <- file.path("/Volumes/OFZ1_T7/karaslab/FragilityResults", subject_code)
  raveio::dir_create2(export_path)


  for (trial_num in (ictal_runs*2)) {

    fragility_pipeline$set_settings(
      project_name = project,
      subject_code = subject_code,
      epoch_name = paste0(subject_code,"_seizure"),
      epoch_time_window = c(-10,30),
      reference_name = "car",
      load_electrodes = electrodes,
      display_electrodes = display,
      trial_num = trial_num,
      t_window = 250,
      t_step = 125,
      sz_onset = 0,
      lambda = 0.001,
      threshold_start = 0,
      threshold_end = 30,
      threshold = 0.5
    )

    trial_num


    # for plotting to pdf ---------------------------------------
    # env <- fragility_pipeline$load_shared()
    source("/Volumes/bigbrain/Fragility2024/rave_pipeline/modules/karaslab_fragility/R/shared-plots.R")

    # test <- fragility_pipeline$eval("repository")
    # str(test$repository$voltage$data_list[[1]][])
    tryCatch(
      error = function(e){
        if (exists(export_path)) {
          file.rename(export_path, file.path("/Volumes/bigbrain/Fragility2024/FragilityResultsACL", paste0(subject_code,"_ERROR")))
        }
      },{
        results <- c(fragility_pipeline$run(c("repository", "adj_frag_info")))

        # save fragility matrix results to csv
        raveio::safe_write_csv(
          results$adj_frag_info$frag,
          file.path(export_path, paste0(subject_code, "_seizure", trial_num/2,"_fragilitynorank.csv"))
        )

        # env <- c(fragility_pipeline$eval(c("repository", "adj_frag_info","threshold_elec")), shortcut = TRUE)
        # results <- list(repository = env[[1]]$repository, adj_frag_info = env[[1]]$adj_frag_info, threshold_elec = env[[1]]$threshold_elec)

        # print results to pdf
        pdf_path <- file.path(export_path, paste0(subject_code,'_seizure',trial_num/2,'_',format(Sys.time(), "%m-%d-%Y_%H%M%S"),'.pdf'))
        grDevices::pdf(pdf_path, width = 12, height = 7)
        par(mfrow=c(2,1),mar=rep(2,4))

        # fragility heatmap
        do.call(fragility_map_plot, c(results,
                                      list(fragility_pipeline$get_settings("display_electrodes"),
                                           fragility_pipeline$get_settings("sz_onset"),
                                           elec_list = elec_list,
                                           'sort_fmap' = 1,
                                           'height' = 14,
                                           threshold = fragility_pipeline$get_settings("threshold"))
        ))

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
}

raveio::safe_write_csv(
  SOZ_table,
  file.path(export_path, paste0(project, '_fragility.csv'))
)
