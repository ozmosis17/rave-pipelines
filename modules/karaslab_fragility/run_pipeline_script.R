fragility_pipeline <- raveio::pipeline("karaslab_fragility")
fragility_pipeline$target_table
raveio::pipeline_visualize(fragility_pipeline$pipeline_path)
fragility_pipeline$get_settings()

fragility_pipeline$set_settings(
  project_name = 'OnsetZone',
  subject_code = 'PT026',
  load_electrodes = "1-10"
)

fragility_pipeline$run()
