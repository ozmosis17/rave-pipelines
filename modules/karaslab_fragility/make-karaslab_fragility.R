library(targets)
library(raveio)
source("common.R", local = TRUE, chdir = TRUE)
._._env_._. <- environment()
lapply(sort(list.files(
  "R/", ignore.case = TRUE,
  pattern = "^shared-.*\\.R", 
  full.names = TRUE
)), function(f) {
  source(f, local = ._._env_._., chdir = TRUE)
})
targets::tar_option_set(envir = ._._env_._.)
rm(._._env_._.)
...targets <- list(`__Check_settings_file` = targets::tar_target_raw("settings_path", 
    "settings.yaml", format = "file"), `__Load_settings` = targets::tar_target_raw("settings", 
    quote({
        yaml::read_yaml(settings_path)
    }), deps = "settings_path", cue = targets::tar_cue("always")), 
    input_trial_num = targets::tar_target_raw("trial_num", quote({
        settings[["trial_num"]]
    }), deps = "settings"), input_ncores = targets::tar_target_raw("ncores", 
        quote({
            settings[["ncores"]]
        }), deps = "settings"), input_load_electrodes = targets::tar_target_raw("load_electrodes", 
        quote({
            settings[["load_electrodes"]]
        }), deps = "settings"), input_t_window = targets::tar_target_raw("t_window", 
        quote({
            settings[["t_window"]]
        }), deps = "settings"), input_epoch_name = targets::tar_target_raw("epoch_name", 
        quote({
            settings[["epoch_name"]]
        }), deps = "settings"), input_t_step = targets::tar_target_raw("t_step", 
        quote({
            settings[["t_step"]]
        }), deps = "settings"), input_epoch_time_window = targets::tar_target_raw("epoch_time_window", 
        quote({
            settings[["epoch_time_window"]]
        }), deps = "settings"), input_reference_name = targets::tar_target_raw("reference_name", 
        quote({
            settings[["reference_name"]]
        }), deps = "settings"), input_subject_code = targets::tar_target_raw("subject_code", 
        quote({
            settings[["subject_code"]]
        }), deps = "settings"), input_nlambda = targets::tar_target_raw("nlambda", 
        quote({
            settings[["nlambda"]]
        }), deps = "settings"), input_project_name = targets::tar_target_raw("project_name", 
        quote({
            settings[["project_name"]]
        }), deps = "settings"), input_sz_onset = targets::tar_target_raw("sz_onset", 
        quote({
            settings[["sz_onset"]]
        }), deps = "settings"), load_subject = targets::tar_target_raw(name = "subject", 
        command = quote({
            .__target_expr__. <- quote({
                library(raveio)
                subject <- RAVESubject$new(project_name = project_name, 
                  subject_code = subject_code, strict = TRUE)
                print(subject)
                subject$epoch_names
                subject$reference_names
                subject$blocks
                subject$electrodes
                all(subject$notch_filtered)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(subject)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "subject", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "subject", target_expr = quote({
                {
                  library(raveio)
                  subject <- RAVESubject$new(project_name = project_name, 
                    subject_code = subject_code, strict = TRUE)
                  print(subject)
                  subject$epoch_names
                  subject$reference_names
                  subject$blocks
                  subject$electrodes
                  all(subject$notch_filtered)
                }
                subject
            }), target_depends = c("project_name", "subject_code"
            )), deps = c("project_name", "subject_code"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"), load_electrodes = targets::tar_target_raw(name = "loading_elec", 
        command = quote({
            .__target_expr__. <- quote({
                loading_elec <- dipsaus::parse_svec(load_electrodes)
                loading_elec <- subject$electrodes[subject$electrodes %in% 
                  loading_elec]
                if (!length(loading_elec)) {
                  stop("No valid electrode to load!")
                }
            })
            tryCatch({
                eval(.__target_expr__.)
                return(loading_elec)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "loading_elec", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "loading_elec", target_expr = quote({
                {
                  loading_elec <- dipsaus::parse_svec(load_electrodes)
                  loading_elec <- subject$electrodes[subject$electrodes %in% 
                    loading_elec]
                  if (!length(loading_elec)) {
                    stop("No valid electrode to load!")
                  }
                }
                loading_elec
            }), target_depends = c("load_electrodes", "subject"
            )), deps = c("load_electrodes", "subject"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"), load_voltage = targets::tar_target_raw(name = "repository", 
        command = quote({
            .__target_expr__. <- quote({
                repository <- raveio::prepare_subject_voltage_with_epoch(subject = subject, 
                  electrodes = loading_elec, epoch_name = epoch_name, 
                  reference_name = reference_name, time_windows = epoch_time_window)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(repository)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "repository", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "repository", target_expr = quote({
                {
                  repository <- raveio::prepare_subject_voltage_with_epoch(subject = subject, 
                    electrodes = loading_elec, epoch_name = epoch_name, 
                    reference_name = reference_name, time_windows = epoch_time_window)
                }
                repository
            }), target_depends = c("subject", "loading_elec", 
            "epoch_name", "reference_name", "epoch_time_window"
            )), deps = c("subject", "loading_elec", "epoch_name", 
        "reference_name", "epoch_time_window"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"), find_adjacency = targets::tar_target_raw(name = "A", 
        command = quote({
            .__target_expr__. <- quote({
                A <- generate_adjacency_array(repository = repository, 
                  trial_num = trial_num, t_window = t_window, 
                  t_step = t_step, nlambda = nlambda)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(A)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "A", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "A", target_expr = quote({
                {
                  A <- generate_adjacency_array(repository = repository, 
                    trial_num = trial_num, t_window = t_window, 
                    t_step = t_step, nlambda = nlambda)
                }
                A
            }), target_depends = c("repository", "trial_num", 
            "t_window", "t_step", "nlambda")), deps = c("repository", 
        "trial_num", "t_window", "t_step", "nlambda"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"), find_fragility = targets::tar_target_raw(name = "f_info", 
        command = quote({
            .__target_expr__. <- quote({
                f_info <- generate_fragility_matrix(A = A, elec = repository$electrode_list, 
                  ncores = 4)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(f_info)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "f_info", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "f_info", target_expr = quote({
                {
                  f_info <- generate_fragility_matrix(A = A, 
                    elec = repository$electrode_list, ncores = 4)
                }
                f_info
            }), target_depends = c("A", "repository")), deps = c("A", 
        "repository"), cue = targets::tar_cue("thorough"), pattern = NULL, 
        iteration = "list"), plot_fragility_map = targets::tar_target_raw(name = "plot_fragility_map_timestamp", 
        command = quote({
            .__target_expr__. <- quote({
                f_info$norm <- f_info$norm[as.character(loading_elec), 
                  ]
                elecsort <- sort(as.numeric(attr(f_info$norm, 
                  "dimnames")[[1]]))
                fsort <- as.numeric(attr(sort(f_info$avg), "names"))
                y <- rev(elecsort)
                f_info$norm <- f_info$norm[as.character(y), ]
                x <- 1:dim(f_info$norm)[2]
                m <- t(f_info$norm)
                attr(m, "xlab") = "Time (s)"
                attr(m, "ylab") = "Electrode"
                attr(m, "zlab") = "Fragility"
                tp <- repository$voltage$dimnames$Time
                yi = seq_along(y)
                if (length(y) > 10) {
                  .seq = seq(1, length(y), length.out = 10)
                  y = y[.seq]
                  yi = .seq
                }
                xtime <- round(seq(tp[1], tp[length(tp)], length.out = 9), 
                  digits = 2)
                xi <- seq(1, length(x), length.out = 9)
                secs <- seq(tp[1], tp[length(tp)])
                onset <- seq(1, length(x), length.out = length(secs))[match(sz_onset, 
                  secs)]
                plot_fragility_map_timestamp <- Sys.time()
                ravebuiltins:::draw_many_heat_maps(list(list(data = m, 
                  x = x, y = seq_along(elecsort), has_trials = TRUE, 
                  range = 0:1)), axes = c(FALSE, FALSE), PANEL.LAST = ravebuiltins:::add_decorator(function(...) {
                  abline(v = onset, lty = 2, lwd = 2)
                  mtext(y, side = 2, line = -1, at = yi, cex = (ravebuiltins:::rave_cex.lab * 
                    0.8), las = 1)
                  mtext(xtime, side = 1, line = 1, at = xi, cex = (ravebuiltins:::rave_cex.lab * 
                    0.8), las = 1)
                }, ravebuiltins:::spectrogram_heatmap_decorator()))
            })
            tryCatch({
                eval(.__target_expr__.)
                return(plot_fragility_map_timestamp)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "plot_fragility_map_timestamp", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "plot_fragility_map_timestamp", target_expr = quote({
                {
                  f_info$norm <- f_info$norm[as.character(loading_elec), 
                    ]
                  elecsort <- sort(as.numeric(attr(f_info$norm, 
                    "dimnames")[[1]]))
                  fsort <- as.numeric(attr(sort(f_info$avg), 
                    "names"))
                  y <- rev(elecsort)
                  f_info$norm <- f_info$norm[as.character(y), 
                    ]
                  x <- 1:dim(f_info$norm)[2]
                  m <- t(f_info$norm)
                  attr(m, "xlab") = "Time (s)"
                  attr(m, "ylab") = "Electrode"
                  attr(m, "zlab") = "Fragility"
                  tp <- repository$voltage$dimnames$Time
                  yi = seq_along(y)
                  if (length(y) > 10) {
                    .seq = seq(1, length(y), length.out = 10)
                    y = y[.seq]
                    yi = .seq
                  }
                  xtime <- round(seq(tp[1], tp[length(tp)], length.out = 9), 
                    digits = 2)
                  xi <- seq(1, length(x), length.out = 9)
                  secs <- seq(tp[1], tp[length(tp)])
                  onset <- seq(1, length(x), length.out = length(secs))[match(sz_onset, 
                    secs)]
                  plot_fragility_map_timestamp <- Sys.time()
                  ravebuiltins:::draw_many_heat_maps(list(list(data = m, 
                    x = x, y = seq_along(elecsort), has_trials = TRUE, 
                    range = 0:1)), axes = c(FALSE, FALSE), PANEL.LAST = ravebuiltins:::add_decorator(function(...) {
                    abline(v = onset, lty = 2, lwd = 2)
                    mtext(y, side = 2, line = -1, at = yi, cex = (ravebuiltins:::rave_cex.lab * 
                      0.8), las = 1)
                    mtext(xtime, side = 1, line = 1, at = xi, 
                      cex = (ravebuiltins:::rave_cex.lab * 0.8), 
                      las = 1)
                  }, ravebuiltins:::spectrogram_heatmap_decorator()))
                }
                plot_fragility_map_timestamp
            }), target_depends = c("f_info", "loading_elec", 
            "repository", "sz_onset")), deps = c("f_info", "loading_elec", 
        "repository", "sz_onset"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"))
