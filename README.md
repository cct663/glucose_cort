# Purpose

This repository contains the complete set of code and data required to reproduce all of the figures and analyses included in both the main text and supplementary material for Taff et al. *Individual variation in natural or manipulated corticosterone does not covary with circulating glucose in a wild bird.*

# Raw data

The raw data folder has all of the data used in the manuscript saved in tab delimited text files. There are three files and the description of column names in each file are given at the end of this readme document.

- **data_glucose_cort.txt** This is the main data file used in the manuscript that includes glucose and corticosterone measurements from adults and nestlings across all four populations.

- **adult_acth_validation.txt** These data are from the acth validation experiment on adults reported in the supplement.

- **nestling_acth_validation.txt** These data are from the acth validation experiment on nestlings reported in the supplement.

# R scripts

This folder contains all of the scripts and associated output from the manuscript. The single script in the folder reproduces all figures and analyses for both the main manuscript and the supplement and is annotated throughout to explain each section. The folder also contains figures and tables produced by the script, but these can also be re-created by running the script itself. The only figure or table not reproduced directly from code is table S2, which lists sample sizes for different groups.

# Column names

Here I describe the column names for each of the raw data files included in the analysis.

## data_glucose_cort

- *band* - unique band number for this individual
- *state* - state indicating which of the four populations this bird is from
- *site* - name of site/location within the state
- *box* - nest box number (unique within site)
- *year* - year of observation
- *date* - day of year for this sample
- *class* - adult or nestling
- *age* - age: SY = second year, ASY = after second year, AHY = after hatch year, 12_15 = 12 to 15 day old nestling
- *sex* - male or female
- *experiment* - larger experiment that this bird was a pert of
- *treatment* - treatment that this individual was associated with (for nestlings, this treatment may have been at the adult or nest)
- *cap_num* - capture number within the year of sampling
- *b_cort* - baseline corticosterone
- *s_cort* - stress-induced corticosterone
- *d_cort* - corticosterone after dexamethasone injection
- *a_cort* - corticosterone after cortrosyn injection
- *inj_type* - type of injection this individual received
- *ss_type* - type of stress series, B = baseline, S = stress-induced, D = post-dex, A## = post cortrosyn with minutes post as ##
- *b_lat* - latency from capture to sample for baseline (seconds)
- *s_lat* - latency from capture to sample for induced (minutes)
- *d_a_lat* - latency from capture to dex or cortrosyn sample (minutes)
- *b_gluc* - baseline glucose
- *s_gluc* - stress-induced glucose
- *d_gluc* - post-dexamethasone glucose
- *a_gluc* - post-cortrosyn glucose
- *stage* - breeding stage of sample
- *post_trt* - was the sample taken after some treatments had occurred
- *mass* - mass at this capture
- *bhead* - head plus bill length in mm
- *fwing* - flat wing length (mm) at this capture
- *massd15* - mass on day 15 for nestlings (day 12 and 15 are combined in one row)
- *time_start* - time in minutes after midnight that capture started
- *notes* - notes about this capture

## adult_acth_validation

- *band* - unique band number of this individual
- *treatment* - sequence of samples taken in three bleeds, B = baseline, C = control/saline, A = acth/cortrosyn
- *inj1* - type of injection for first injection (all saline)
- *inj2* - type of injection for second injection (saline or acth)
- *unit* - location where bird was captured
- *nest* - nest box (unique within location)
- *cort1* - corticosterone from first sample
- *cort2* - corticosterone from second sample
- *cort3* - corticosterone from third sample

## nestling_acth_validation

- *band* - unique band number of this individual
- *unit_box* - location/unit and nest box number of this nestling
- *unit* - location of this individual
- *nest* - nest box number (unique within unit)
- *typecort* - redundant with baseline cort, can be ignored
- *treatment* - treatment applied to this individual (saline or acth injection)
- *cort1* - corticosterone measurement for first, baseline, sample
- *cort2* - corticosterone measurement for second sample
- *cort3* - corticosteorne measurement for third sample