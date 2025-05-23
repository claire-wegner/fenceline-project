README
================

## About the project

The data and R script included in this repository were used in the
analysis for the research article by Wegner et al. (2025)
titled “Behavioural responses of dairy cows and calves to fenceline
weaning after 4 or 6 months of full cow-calf contact”. The script also
includes the creation of Fig. 4, Fig. 5 and Fig. 6. For details on study
design, please refer to the materials and methods as written in the
article.

## Description of variables

For further details on how each variable was collected and/or
calculated, please refer to the full article.

### id

The ID of the individual cow or calf.

### obs_day

The observation day in relation to the day of fenceline weaning (day 0;
not included in dataset) for day 1 - 11 after weaning. ‘BL’ refers to
the 6 days immediately prior to each weaning event.

### daily_L\_m

The daily time spent lying down in minutes.

### daily_steps

The daily number of steps taken.

### date

The exact date on which the data pertains to. The day of weaning is not
included in the dataset (weaning for 4MO = 2022-06-27; 6MO =
2022-08-24). The date shown for obs_day = ‘BL’ is the first day included
in this value; in reality, 6 days of data were used to create average
values for step count and lying time, which is thereafter referred to as
the ‘baseline’ level for either behaviour.

### trt

4_mo = fenceline weaning after 4 months of cow-calf contact; 6_mo =
fenceline weaning after 6 months of cow-calf contact

### parity

Parity of the cow.

### cow_calf

Whether the individual is a ‘cow’ or ‘calf’.

### percent_near

The percentage of time spent by a related cow-calf pair in close
proximity. Close proximity is defined as there being less than 4 m
(indoor areas) or 8 m (outdoor areas) between cow and calf.

### F_min_hourly

The total minutes per hour spent on feed-seeking behaviours; data only
collected for calves.

### voc_yes

The mean number of 5-min intervals per hour containing at least 1
vocalization.

### voc_perc_yes

The mean percentage of 5-minute intervals per hour containing at least 1
vocalization, as calculated by comparing the number of intervals
containing vocalizations to the total possible intervals per hour.
