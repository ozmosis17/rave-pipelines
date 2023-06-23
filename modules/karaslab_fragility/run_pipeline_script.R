# Add `path` to force using devel pipeline
fragility_pipeline <- raveio::pipeline("karaslab_fragility", paths = "./modules/")
fragility_pipeline$target_table
raveio::pipeline_visualize(fragility_pipeline$pipeline_path)
fragility_pipeline$get_settings()

fragility_pipeline$set_settings(
  project_name = 'OnsetZone',
  subject_code = 'PT026',
  load_electrodes = "1-10"
)

plot_data <- fragility_pipeline$run(c("repository", "f_info", "loading_elec", "sz_onset"))
source("./modules/karaslab_fragility/R/shared-plots.R")
# env <- fragility_pipeline$load_shared()


subject <- fragility_pipeline$run("subject")
export_path <- file.path(subject$note_path, "karaslab_fragility", "exports")
raveio::dir_create2(export_path)

pdf_path <- file.path(export_path, "junk.pdf")

export_pdf <- function(expr, path, env = parent.frame(),
                       quoted = FALSE, width = 12, height = 7, useDingbats = FALSE, ...) {
  force(path)
  if(!quoted) {
    expr <- substitute(expr)
  }
  grDevices::pdf(path, width = width, height = height, useDingbats = useDingbats, ...)
  on.exit({
    grDevices::dev.off()
  }, add = TRUE, after = TRUE)
  eval(expr, envir = env)
}

export_pdf(
  {
    do.call(fragility_map_timestamp_plot, plot_data)
  },
  path = pdf_path
)








