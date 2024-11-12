ravedash::debug_modules(module_root = rstudioapi::getActiveProject())

pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")

library(readxl)
library(stringr)

export_path <- "/karaslab/Results_FragilityEEGDataset"

pts <- dipsaus::parse_svec("94-102,131-139")
pipeline_xls <- read.csv("E:/karaslab/rave-pipelines/modules/karaslab_fragility/Outcome_Classification_ML/patient_data_all_rev.csv")
pipeline_xls$subject[pts]

for(i in pts){
  # check pt metadata if needed
  subject_code <- pipeline_xls$subject_code[i]
  project <- pipeline_xls$project_name[i]
  subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))
  #print(subject$reference_names)
}

for(i in pts){

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

  subject_check <- raveio::validate_subject(paste0(project,"/",subject_code),
                                            method = "basic", verbose = FALSE)
  subject_check$paths$data_path$valid

  print(paste0("starting pipeline for pt: ", subject_code, ", ", condition))

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

  trial_num <- match(condition, subject$get_epoch(epoch_name)$table$Condition)

  fragility_pipeline$set_settings(
    project_name = project,
    subject_code = subject_code,
    epoch_name = epoch_name,
    epoch_time_window = c(-10,20),
    reference_name = reference_name,
    load_electrodes = electrodes,
    display_electrodes = display,
    condition = condition,
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

  # for saving output files ---------------------------------------
  # env <- fragility_pipeline$load_shared()
  source("./modules/karaslab_fragility/R/shared-functions.R")
  source("./modules/karaslab_fragility/R/shared-plots.R")

  tryCatch(
    error = function(e){
      if (file.exists(export)) {
        file.create(file.path(export, paste0(subject_code,"_",fragility_pipeline$get_settings("condition"),"_ERROR")))
      }
    },{
      results <- c(fragility_pipeline$run(c("repository", "adj_frag_info","quantiles")))

      # force evaluation
      #env <- c(fragility_pipeline$eval(c("repository", "adj_frag_info","threshold_elec")), shortcut = TRUE)
      #results <- list(repository = env[[1]]$repository, adj_frag_info = env[[1]]$adj_frag_info, threshold_elec = env[[1]]$threshold_elec)

      moving_average_width <- 10

      # save unranked results
      output_files(results$repository, results$adj_frag_info$frag, results$quantiles,
                   fragility_pipeline$get_settings(),export,"norank", moving_average_width)

      # save ranked results
      output_files(results$repository, results$adj_frag_info$frag_ranked, results$quantiles,
                   fragility_pipeline$get_settings(),export,"ranked", moving_average_width)

      # save R2 to csv with optimal lambdas
      raveio::safe_write_csv(
        rbind(results$adj_frag_info$R2, results$adj_frag_info$lambdas),
        file.path(export, paste0("/",subject_code,"_",fragility_pipeline$get_settings("condition"),"_R2.csv"))
      )

      # print results to pdf
      pdf_path <- file.path(export, paste0(subject_code,'_',fragility_pipeline$get_settings("condition"),"_reconstruction.pdf"))
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

      # # old fragility heatmap
      # do.call(fragility_map_plot_old, c(results,
      #                               list(fragility_pipeline$get_settings("display_electrodes"),
      #                                    fragility_pipeline$get_settings("sz_onset"),
      #                                    elec_list = elec_list,
      #                                    'sort_fmap' = 1,
      #                                    'height' = 14,
      #                                    threshold = fragility_pipeline$get_settings("threshold"))
      # ))

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
# trial_num = 3
# t_window = 250
# t_step = 125
# timepoints = 1:1000
# elec_num = 1
# percentile = 0.1
# display_electrodes <- c(33,34,62:69,88:91)
