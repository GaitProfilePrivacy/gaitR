library(tidyverse)
library(reticulate)
library(R.matlab)
library(gganimate)
library(gifski)
library(stringr)
library(purrr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# This is an R package named 'gaitR'

# a custom function to read and tidy the gait_profiles mat file
gait_reader <- function(link){
  # extracting the gait sublist from the read mat file
  gait_profiles = readMat(link) %>% `[[`('gait')

  # Matlab reads the cell data in a table which is read columnwise
  # R reads the data into a single list rowwise so the computations below
  # are to recreate the Matlab file structure into R
  starting_cells = 1 + (0:4)*length(gait_profiles)/5
  ending_cells = (1:5) * (length(gait_profiles)/5)

  # extracting gait time, acceleration values & experiment time from the R list
  subject_num = str_extract(link, pattern = '[:digit:]+' )

  gait_time = gait_profiles[starting_cells[1]:ending_cells[1]] %>%
    unlist(recursive = F) %>% map(.f = c)

  acc_x = gait_profiles[starting_cells[2]:ending_cells[2]] %>%
    unlist(recursive = F) %>% map(.f = c)

  acc_y = gait_profiles[starting_cells[3]:ending_cells[3]] %>%
    unlist(recursive = F) %>% map(.f = c)

  acc_z = gait_profiles[starting_cells[4]:ending_cells[4]] %>%
    unlist(recursive = F) %>% map(.f = c)

  # computed from the data but not stored since it does not add new information
  # when compared to the gait_time column
  # the computation of exp_time is only shown here to allow others to utilize
  # it instead of gait_times if needed
  exp_time = gait_profiles[starting_cells[5]:ending_cells[5]] %>%
    unlist(recursive = F) %>% map(.f = c)

  # storing the results in a tibble to match the original format
  gait_profiles_tbl = tibble(subject_num, gait_time, acc_x, acc_y, acc_z) %>%
    # creating a gait_cycle number/id column
    rownames_to_column('gait_cycle') %>%
    # moving the newly created gait_cycle column after subject_num
    relocate(gait_cycle, .after = subject_num)

  # creating a long/tidy format of the data to facilitate the calculation of the
  # resultant/magnitude of the acceleration and the visualization of the data
  gait_profiles_tbl_long = gait_profiles_tbl %>%
    unnest( cols = c(gait_time, acc_x, acc_y, acc_z) ) %>%
    mutate( acc_magnitude = sqrt(acc_x^2 + acc_y^2 + acc_z^2) )

  # saving select columns to reduce size of the returned data
  # obviously, one can comment the next two lines to return the entire data
  # if needed
  gait_profiles_tbl_long = gait_profiles_tbl_long %>%
    select(subject_num, gait_cycle, gait_time, acc_magnitude)

  return(gait_profiles_tbl_long)
}

extract_gait_profiles <- function(){
  file_names = c(paste0('Subject',
                        str_pad(1:15, width = 2, side = 'left', pad = '0'),
                        '_aZ_seg.mat'))

  # dropping subject 13 since we do not have them in our dat
  file_names = file_names[-13]

  # creating a vector of file link locations based on the github repo
  # https://github.com/saebragani/Gait-Hotelling/tree/master/Data/matFiles/Profiles
  file_links = paste0(
    'https://github.com/saebragani/Gait-Hotelling/blob/master/Data/matFiles/Profiles/',
    file_names,
    '?raw=true')

  # utilizing the custom gait_reader function to read all the data
  # map_df allows the storage of all the data into a singular data frame
  all_subj_gaits = map_df(.x = file_links, .f = gait_reader)

  # printing a summary and the structure of the all_subj_gaits object
  # glimpse(all_subj_gaits)

  return(all_subj_gaits)
}

plot_profile_samples <- function(all_subj_gaits, sample_subjects, cycle_id){
  # selecting four arbitrary subjects to facilitate the visualization of data
  # sample_subjects = c('02', '06', '10', '14')

  # creating a named character vector which will be used in renaming the panels
  sample_subject_labels = c('02' = 'Subject 02',
                            '06' = 'Subject 06',
                            '10' = 'Subject 10',
                            '14' = 'Subject 14')

  # plotting the data for select subjects and a sample of seven gait cycles
  all_subj_gaits %>%
    # selecting the sample subjects
    filter(subject_num %in% sample_subjects) %>%
    # selecting seven arbitrary gait cycles
    filter(gait_cycle %in% cycle_id) %>%
    # specifying the aesthetics
    ggplot( aes(x = gait_time, y = acc_magnitude, group = gait_cycle,
                color = gait_cycle) ) +
    facet_wrap(facets = ~ subject_num, nrow = 2, ncol = 2,
               labeller = as_labeller(sample_subject_labels)) +
    geom_line(size = 1) +
    scale_color_brewer(palette = 'Dark2') +
    scale_x_continuous(breaks = scales::pretty_breaks(11)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    theme_bw() +
    labs(x = 'Gait Cycle Time (s)',
         y = bquote('Acceleration Magnitude ' ~(m/s^2)),
         title = 'Seven Sample Gait Cycles for Four Subjects',
         subtitle = 'Two-hump patterns are somewhat different and within subject variability can be large',
         color = 'Gait Cycle') +
    theme(legend.position = 'bottom') +
    guides(color = guide_legend(nrow = 1))
}

plot_animated <- function(all_subj_gaits, subject_id){
  # filtering the data to only subject 01
  subj01_gaits = all_subj_gaits %>% filter(subject_num %in% subject_id)

  # finding the last gait cycle for subject 01
  num_gait_cycles_subj01 = subj01_gaits %>% pull(gait_cycle) %>% tail(n=1)

  p1 = subj01_gaits %>%
    # filtering every 100th gait cycle
    filter(gait_cycle %in% as.character( seq(1, num_gait_cycles_subj01, 100) ) ) %>%
    # specifying the aesthetics
    ggplot( aes(x = gait_time, y = acc_magnitude, group = gait_cycle) ) +
    # specifying the width of the line
    geom_line(size = 1) +
    # prettying up the chart
    scale_x_continuous(breaks = scales::pretty_breaks(11)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    theme_bw() +
    # setting the transition for the animation to be based on the gait cycle
    transition_manual( gait_cycle ) +
    labs(
      title = "Subject 01's Gait Profiles Sampled at Every 100th Gait Cycle",
      subtitle = 'Sampled Gait Cycle: {frame}' ,
      x = 'Gait Cycle Time (s)',
      y ='Acceleration Magnitude (m/s^2)')

  # setting the animation speed to one frame per second
  animate(p1, fps = 1)
}
