---
title: "RAVE Fragility"
output:
  html_document: default
  pdf_document: default
---



## Introduction

A RAVE pipeline markdown is an interactive notebook that can keep your notes, code blocks, and corresponding results together, generating reports in various formats such as `PDF`, `html`, `Word`, `PowerPoint`. 

The note parts are simply `markdown`s - the same as `jupyter notebook`, or github documentations. The code blocks support `R`, `python`, `c++`. When you hit the `Knit` button in this code editor panel, the r-markdown file will be compiled, generating reproducible documentation.

With carefully designed structures, this r-markdown file will automatically generate `RAVE` pipeline scripts during the compilation. The pipeline script can be used by `RAVE` to convert your pipeline into interactive dashboard application. (This feature is currently under development)

## "RAVE" Pipeline Code Block

A `RAVE` pipeline markdown code block starts with ` ```{rave ... `. The block label following `rave` informative description of the target. After the target, the following RAVE-specific parameters configures how the block should be treated:

* `language`: specifies the programming language used; choices are: `R`, `python`
* `export`: variable name to be exported that will be available to the rest chunks
* `depends`: indicates current block depends on variables generated from other blocks; this helps `RAVE` to build non-interactive pipeline scripts internally. For blocks written in `R`, the dependence can be automatically determined.

Other parameters are available at [this `rmarkdown` book](https://bookdown.org/yihui/rmarkdown/)

## An Example

In the rest of the documentation, let's import the subject power data, baseline, and plot the collapsed mean as image.

#### Step 1: Create `RAVE` subject's instances

Noting that all the items in the `settings.yaml` are available as variables.


```r
# Load subject instance
subject <- raveio::RAVESubject$new(project_name = project_name,
                                   subject_code = subject_code,
                                   strict = TRUE)
```

With `export="subject"`, the subject variable will be registered for the following chunks to use. Be aware that all other variables created in this block will not be exposed.

#### Step 2: Initialize and load voltage data

Initialize the electrode instances and register the epoch, reference information


```r
loading_elec <- dipsaus::parse_svec(load_electrodes)
loading_elec <- subject$electrodes[subject$electrodes %in% loading_elec]
if(!length(loading_elec)) {
  stop("No valid electrode to load!")
}
```

```r
displayed_elec <- dipsaus::parse_svec(display_electrodes)
displayed_elec <- subject$electrodes[subject$electrodes %in% displayed_elec]
if(!length(displayed_elec)) {
  stop("No valid electrode to load!")
}
```

Start to load voltage. Here also create cache to the `RAVE` cache directory.


```r
repository <- raveio::prepare_subject_voltage_with_epoch(
  subject = subject,
  electrodes = loading_elec,
  epoch_name = epoch_name,
  reference_name = reference_name,
  time_windows = epoch_time_window
)
```

#### Step 3: Find Adjacency Array and Fragility Matrix


```r
adj_frag_info <- calc_adj_frag(
  repository = repository,
  trial_num = trial_num,
  t_window = t_window,
  t_step = t_step,
  lambda = lambda
)
#> Loading required package: iterators
#> Loading required package: parallel
#> Loading required package: Matrix
#> Loaded glmnet 4.1-8
```

## Build, Visualize, & Run

Please make sure the following code block is at the end of your pipeline file. This block will build the pipeline and generate a `make-karaslab_fragility.R` script with your pipeline markdown file. `RAVE` will use the generated pipeline script to execute the pipeline in the dashboard application, or in massive production mode.


```
#> Building target subject [R]
#> Building target loading_elec [R]
#> Building target displayed_elec [R]
#> Building target repository [R]
#> Building target adj_frag_info [R]
```


Once the pipeline script `make-karaslab_fragility.R` is built, you can visualize and execute the pipeline without the need of re-knit this document. Notice we use `r` block instead of `rave`. (This is because the code blocks are not part of pipeline targets.)


```{=html}
<div class="visNetwork html-widget html-fill-item" id="htmlwidget-f2a22ed0446b834d2fba" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-f2a22ed0446b834d2fba">{"x":{"nodes":{"name":["adj_frag_info","display_electrodes","displayed_elec","epoch_name","epoch_time_window","lambda","load_electrodes","loading_elec","project_name","reference_name","repository","settings","settings_path","signalScaling","subject","subject_code","sz_onset","t_step","t_window","trial_num","signalscaling"],"type":["stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem","stem"],"status":["outdated","outdated","outdated","outdated","outdated","outdated","outdated","outdated","outdated","outdated","outdated","outdated","uptodate","outdated","outdated","outdated","outdated","outdated","outdated","outdated","outdated"],"seconds":[410.199,0,0,0.001,0,0,0,0.001,0,0,1.345,0,0,0,0.01,0,0.001,0,0.001,0,null],"bytes":[4785102,109,109,66,54,55,193,224,72,56,6321845,534,900,53,66,58,49,53,53,51,null],"branches":[null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],"label":["adj_frag_info","display_electrodes","displayed_elec","epoch_name","epoch_time_window","lambda","load_electrodes","loading_elec","project_name","reference_name","repository","settings","settings_path","signalScaling","subject","subject_code","sz_onset","t_step","t_window","trial_num","signalscaling"],"color":["#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#354823","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5","#78B7C5"],"id":["adj_frag_info","display_electrodes","displayed_elec","epoch_name","epoch_time_window","lambda","load_electrodes","loading_elec","project_name","reference_name","repository","settings","settings_path","signalScaling","subject","subject_code","sz_onset","t_step","t_window","trial_num","signalscaling"],"level":[7,3,5,3,3,3,3,5,3,3,6,2,1,3,4,3,3,3,3,3,3],"shape":["dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot"]},"edges":{"from":["settings","settings","settings","settings","subject","loading_elec","epoch_name","reference_name","epoch_time_window","display_electrodes","subject","settings","settings","settings","repository","trial_num","t_window","t_step","lambda","load_electrodes","subject","settings","settings","settings","settings","settings","project_name","subject_code","settings","settings","settings_path"],"to":["epoch_name","reference_name","load_electrodes","display_electrodes","repository","repository","repository","repository","repository","displayed_elec","displayed_elec","lambda","t_step","subject_code","adj_frag_info","adj_frag_info","adj_frag_info","adj_frag_info","adj_frag_info","loading_elec","loading_elec","t_window","epoch_time_window","signalscaling","sz_onset","project_name","subject","subject","signalScaling","trial_num","settings"],"arrows":["to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to","to"]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot","physics":false},"manipulation":{"enabled":false},"edges":{"smooth":{"type":"cubicBezier","forceDirection":"horizontal"}},"physics":{"stabilization":false},"interaction":{"zoomSpeed":0.1},"layout":{"hierarchical":{"enabled":true,"direction":"LR"}}},"groups":null,"width":null,"height":null,"idselection":{"enabled":false,"style":"width: 150px; height: 26px","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":false,"hideColor":"rgba(200,200,200,0.5)","highlight":false},"main":{"text":"","style":"font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:20px;text-align:center;"},"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","highlight":{"enabled":true,"hoverNearest":false,"degree":{"from":1,"to":1},"algorithm":"hierarchical","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":true,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"},"legend":{"width":0.2,"useGroups":false,"position":"right","ncol":1,"stepX":100,"stepY":100,"zoom":true,"nodes":{"label":["Outdated","Up to date","Stem"],"color":["#78B7C5","#354823","#899DA4"],"shape":["dot","dot","dot"]},"nodesToDataframe":true},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);"},"evals":[],"jsHooks":[]}</script>
```









