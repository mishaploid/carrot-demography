Carrot demography and split time plots
================
Sarah Turner-Hissong

- <a href="#1-setup" id="toc-1-setup">1 Setup</a>
- <a href="#2-population-size-history"
  id="toc-2-population-size-history">2 Population size history</a>
  - <a href="#21-load-smc-estimate-results"
    id="toc-21-load-smc-estimate-results">2.1 Load SMC++ estimate
    results</a>
  - <a href="#22-plot-results-from-smc-estimate-fig-3c"
    id="toc-22-plot-results-from-smc-estimate-fig-3c">2.2 Plot results from
    SMC++ estimate (Fig. 3c)</a>
  - <a href="#23-check-bootstrap-results-for-each-population"
    id="toc-23-check-bootstrap-results-for-each-population">2.3 Check
    bootstrap results for each population</a>
  - <a href="#24-determine-timepoint-of-minimum-ne-for-each-population"
    id="toc-24-determine-timepoint-of-minimum-ne-for-each-population">2.4
    Determine timepoint of minimum Ne for each population</a>
- <a href="#3-joint-demography--split-times"
  id="toc-3-joint-demography--split-times">3 Joint demography + split
  times</a>
  - <a href="#31-load-smc-split-results"
    id="toc-31-load-smc-split-results">3.1 Load SMC++ split results</a>
  - <a href="#32-process-results-from-smc-estimate"
    id="toc-32-process-results-from-smc-estimate">3.2 Process results from
    SMC++ estimate</a>
  - <a href="#33-filter-and-order-smc-split-results"
    id="toc-33-filter-and-order-smc-split-results">3.3 Filter and order
    SMC++ split results</a>
  - <a href="#34-summarize-split-times"
    id="toc-34-summarize-split-times">3.4 Summarize split times</a>
  - <a href="#35-plot-results-from-smc-split-fig-3d"
    id="toc-35-plot-results-from-smc-split-fig-3d">3.5 Plot results from
    SMC++ split (Fig. 3d)</a>
- <a href="#4-combine-panels-for-fig-3cd"
  id="toc-4-combine-panels-for-fig-3cd">4 Combine panels for Fig 3c/d</a>

This notebook walks through plotting results from SMC++ for carrot
populations.

Reference: Coe, K., Bostan, H., Rolling, W. et al. Population genomics
identifies genetic signatures of carrot domestication and improvement
and uncovers the origin of high-carotenoid orange carrots. Nat. Plants
9, 1643–1658 (2023). <https://doi.org/10.1038/s41477-023-01526-6>

Plots of divergence time inspired by the lovely figures in this
publication: <https://www.nature.com/articles/s41588-018-0215-8>

# 1 Setup

``` r
library(tidyverse)
library(jsonlite)
library(cowplot)

# color palette for carrot populations
# update 9/27/2024 - correction for colors
# incorrect colors: Landrace A (purple), Landrace B (yellow)
# correct colors: Landrace A (yellow), Landrace B (purple)
cols <- c("Wild" = "#2F7EB8", 
          "Landrace_A" = "#FDD658", 
          "Landrace_B" = "#9950A0",
          "Early_Cultivar" = "#FF7B33", 
          "Improved_Cultivar" = "#E61321")

# colors for population codes from split times (formatted without '_')
cols2 <- c("Wild" = "#2F7EB8", 
           "LandraceA" = "#FDD658", 
           "LandraceB" = "#9950A0",
           "EarlyCultivar" = "#FF7B33", 
           "ImprovedCultivar" = "#E61321")
```

# 2 Population size history

## 2.1 Load SMC++ estimate results

This step reads in the results from `smcpp estimate`.  
Input files were generated using the `smcpp plot` command with the
`--csv` flag.

See snakemake rules for details and
<https://github.com/popgenmethods/smcpp> for more information on SMC++.

``` r
# results from smc++ estimate 
smc_estimate_results <- read_csv("results/carrot_all_pops_smc_estimate_no_timepoints.csv") %>%
  # reorder pop labels based on color palette 
  mutate(label = fct_relevel(label, names(cols))) %>% 
  # remove landrace A wild samples (these are likely feral)
  filter(!label == 'LandraceAWild')

head(smc_estimate_results)
```

    ## # A tibble: 6 × 5
    ##   label              x       y plot_type plot_num
    ##   <fct>          <dbl>   <dbl> <chr>        <dbl>
    ## 1 Early_Cultivar 0     109139. path             0
    ## 2 Early_Cultivar 0.133 109139. path             0
    ## 3 Early_Cultivar 0.154 112263. path             0
    ## 4 Early_Cultivar 0.177 121929. path             0
    ## 5 Early_Cultivar 0.205 139396. path             0
    ## 6 Early_Cultivar 0.236 167236. path             0

``` r
# results from smc++ estimate for parametric bootstrapping 
# 10 replicates, genomic data resampled in 5 Mb blocks 
smc_estimate_bootstrap <- read_csv("results/carrot_all_pops_smc_estimate_no_timepoints_bootstrap.csv") %>%
  # reorder pop labels based on color palette 
  mutate(label = fct_relevel(label, names(cols))) %>% 
  # remove landrace A wild samples (these are likely feral)
  filter(!label == 'LandraceAWild')

head(smc_estimate_bootstrap)
```

    ## # A tibble: 6 × 5
    ##   label              x       y plot_type plot_num
    ##   <fct>          <dbl>   <dbl> <chr>        <dbl>
    ## 1 Early_Cultivar 0      88193. path             0
    ## 2 Early_Cultivar 0.127  88193. path             0
    ## 3 Early_Cultivar 0.146  90653. path             0
    ## 4 Early_Cultivar 0.169  98182. path             0
    ## 5 Early_Cultivar 0.194 111583. path             0
    ## 6 Early_Cultivar 0.224 132526. path             0

## 2.2 Plot results from SMC++ estimate (Fig. 3c)

This step plots the marginal and bootstrapped results from
`smcpp estimate`.

Note: the colors of Landrace A and Landrace B are swapped compared to
the paper figure - will need to follow up on this.

``` r
estimate_plot <- 
  ggplot() + 
  # shade region of hypothesized carrot domestication window
  annotate('rect', 
           xmin = 522, 
           xmax = 1100, 
           ymin = 0, 
           ymax = Inf,
           alpha = 0.5,
           fill = 'gray') + 
  # add lines for bootstrapped estimates 
  geom_path(data = smc_estimate_bootstrap,
            aes(x = x,
                y = y,
                color = label,
                group = plot_num),
            alpha = 0.1) + 
  # add lines for marginal estimates 
  geom_path(data = smc_estimate_results, 
            aes(x = x, 
                y = y,
                color = label), 
            size = 1.5,
            lineend = 'round') +
  scale_color_manual(values = cols) + 
  # use a log10 scale for x and y axes and set number of breaks 
  # note: set x min to 80 yrs
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(80,NA),
                expand = c(0,0)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 3),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks() +
  xlab("Years ago") +
  ylab("Effective population size") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(.8,.8)) 

estimate_plot 
```

![](summarize_smc_results_files/figure-gfm/plot%20smcpp%20estimate%20results-1.png)<!-- -->

## 2.3 Check bootstrap results for each population

This code facets the bootstrapping results from SMC++ estimate to view
results by population.

``` r
estimate_bootstrap_plot <- ggplot(smc_estimate_bootstrap) + 
  # add hypothesized window of carrot domestication 
  annotate('rect', 
           xmin = 522, 
           xmax = 1100, 
           ymin = 0, 
           ymax = Inf,
           alpha = 0.5,
           fill = 'gray') + 
  geom_line(aes(x = x, 
                y = y,
                group = plot_num,
                color = label),
            size = 1,
            alpha = 0.5) +
  scale_color_manual(values = cols) + 
  # use a log10 scale for x and y axes and set number of breaks 
  # note: set x min to 80 yrs
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(80,NA)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 2),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks() +
  xlab("Generations ago") +
  ylab("Effective population size") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        legend.box.background = element_rect(colour = "black")) +
  facet_wrap(~label)

estimate_bootstrap_plot
```

![](summarize_smc_results_files/figure-gfm/plot%20smcpp%20estimate%20by%20population-1.png)<!-- -->

## 2.4 Determine timepoint of minimum Ne for each population

``` r
# minimum and maximum Ne for each population 
smc_estimate_results %>% 
  # exclude estimates within the last 100 years - not reliable 
  # also exlude wild samples - no domestication bottleneck 
  filter(x >= 100,
         !label == "Wild") %>% 
  group_by(label) %>% 
  # filter for the lowest observed Ne in each population
  filter(y == min(y)) %>% 
  select(-plot_type, -plot_num)
```

    ## # A tibble: 4 × 3
    ## # Groups:   label [4]
    ##   label                 x     y
    ##   <fct>             <dbl> <dbl>
    ## 1 Early_Cultivar    1068.  953.
    ## 2 Improved_Cultivar 1092.  895.
    ## 3 Landrace_A        1503. 1361.
    ## 4 Landrace_B        1690. 1206.

# 3 Joint demography + split times

## 3.1 Load SMC++ split results

This step reads in the results from `smcpp split` (JSON format).

``` r
# results from smc++ split 
smc_split_results <- list.files(path = "results/smc_split", 
                                pattern = "*.json",
                                recursive = TRUE,
                                full.names = TRUE) %>%
  set_names(word(., start = 3, end = 3, sep = "/")) %>% 
  map(~RJSONIO::fromJSON(., flatten = TRUE)) 

str(smc_split_results[[1]])
```

    ## List of 5
    ##  $ alpha        : num 1
    ##  $ hidden_states:List of 2
    ##   ..$ Early_Cultivar   :List of 18
    ##   .. ..$ : num 0
    ##   .. ..$ : num 2.67e-05
    ##   .. ..$ : num 3.64e-05
    ##   .. ..$ : num 4.84e-05
    ##   .. ..$ : num 6.58e-05
    ##   .. ..$ : num 9.1e-05
    ##   .. ..$ : num 0.000128
    ##   .. ..$ : num 0.000187
    ##   .. ..$ : num 0.000309
    ##   .. ..$ : num 0.265
    ##   .. ..$ : num 0.604
    ##   .. ..$ : num 0.978
    ##   .. ..$ : num 2.38
    ##   .. ..$ : num 6.07
    ##   .. ..$ : num 18.8
    ##   .. ..$ : num 36.3
    ##   .. ..$ : num 80.4
    ##   .. ..$ : NULL
    ##   ..$ Improved_Cultivar:List of 18
    ##   .. ..$ : num 0
    ##   .. ..$ : num 2.75e-05
    ##   .. ..$ : num 3.73e-05
    ##   .. ..$ : num 5.04e-05
    ##   .. ..$ : num 6.91e-05
    ##   .. ..$ : num 9.58e-05
    ##   .. ..$ : num 0.000136
    ##   .. ..$ : num 0.000202
    ##   .. ..$ : num 0.000359
    ##   .. ..$ : num 0.299
    ##   .. ..$ : num 0.582
    ##   .. ..$ : num 0.93
    ##   .. ..$ : num 2.59
    ##   .. ..$ : num 6
    ##   .. ..$ : num 19.2
    ##   .. ..$ : num 37
    ##   .. ..$ : num 81.1
    ##   .. ..$ : NULL
    ##  $ model        :List of 4
    ##   ..$ class : chr "SMCTwoPopulationModel"
    ##   ..$ model1:List of 6
    ##   .. ..$ N0          : num 1250
    ##   .. ..$ class       : chr "SMCModel"
    ##   .. ..$ knots       : num [1:8] 2.67e-05 4.84e-05 9.10e-05 1.87e-04 2.65e-01 ...
    ##   .. ..$ pid         : chr "Early_Cultivar"
    ##   .. ..$ spline_class: chr "CubicSpline"
    ##   .. ..$ y           : num [1:8] 5.173 5.638 7.118 9.914 0.509 ...
    ##   ..$ model2:List of 6
    ##   .. ..$ N0          : num 1250
    ##   .. ..$ class       : chr "SMCModel"
    ##   .. ..$ knots       : num [1:8] 2.75e-05 5.04e-05 9.58e-05 2.02e-04 2.99e-01 ...
    ##   .. ..$ pid         : chr "Improved_Cultivar"
    ##   .. ..$ spline_class: chr "CubicSpline"
    ##   .. ..$ y           : num [1:8] 2.69 3.32 5.348 9.458 0.464 ...
    ##   ..$ split : num 0.302
    ##  $ rho          : num 0.01
    ##  $ theta        : num 1e-04

``` r
# results from smc++ split for parametric bootstrapping 
# 10 replicates, genomic input resampled in 5 Mb blocks 
smc_split_bootstrap_results <- list.files("results/smc_split_bootstrap",
                                    pattern = "*.json",
                                    recursive = TRUE,
                                    full.names = TRUE) %>% 
  set_names(word(., start = 3, end = 3, sep = "/")) %>% 
  map(~RJSONIO::fromJSON(., flatten = TRUE))
```

## 3.2 Process results from SMC++ estimate

This step calculates split times from the outputs of `smcpp split` to
convert coalescent scaling to time in generations using:
$2 \cdot N0 \cdot split\ time \cdot generation\ time$

For carrot, assumed a time of two years per generation.

``` r
# create an object to store results 
split_times <- NULL

# loop through each population comparison to convert from coalescent scaling to years
# using a conversion of two years/generation 
for (i in names(smc_split_results)) {
  pop_pair <- i
  N0 <- smc_split_results[[i]]$model$model1$N0
  # split time 
  split <- smc_split_results[[i]]$model$split 
  # convert from coalescent scaling to years
  split_years <- 2*N0*split*2
  df <- data.frame(pop_pair, N0, split, split_years)
  # combine into a data frame 
  split_times <- rbind(split_times, df)
}

head(split_times)
```

    ##                         pop_pair   N0     split split_years
    ## 1 EarlyCultivar_ImprovedCultivar 1250 0.3018658    1509.329
    ## 2             EarlyCultivar_Wild 1250 3.9721659   19860.830
    ## 3          ImprovedCultivar_Wild 1250 1.7884016    8942.008
    ## 4        LandraceA_EarlyCultivar 1250 0.3853637    1926.818
    ## 5     LandraceA_ImprovedCultivar 1250 0.4332227    2166.114
    ## 6            LandraceA_LandraceB 1250 0.1628104     814.052

``` r
# repeat for bootstrapped results :) 
bootstrap_split_times <- NULL

for (i in names(smc_split_bootstrap_results)) {
  pop_pair <- i
  N0 <- smc_split_bootstrap_results[[i]]$model$model1$N0
  # split time 
  split <- smc_split_bootstrap_results[[i]]$model$split 
  # convert from coalescent scaling to years
  # two years/generation
  split_years <- 2*N0*split*2
  df <- data.frame(pop_pair, N0, split, split_years)
  bootstrap_split_times <- rbind(bootstrap_split_times, df)
}
```

## 3.3 Filter and order SMC++ split results

This step adds some project-specific filtering and reorders population
comparisons based on split time.

``` r
# add preferred factor ordering and filtering for split time estimates 
split_times_filtered <- split_times %>% 
  # order population pairs by minimum split time 
  mutate(pop_pair = fct_reorder(pop_pair, split_years, .fun = min)) %>% 
  # split population pair into two fields for plotting 
  separate(pop_pair, into = c('pop1', 'pop2'), sep = '_', remove = FALSE) %>% 
  # remove Landrace A wild samples from consideration 
  filter(!pop2 == 'LAwild') %>% 
  # focus only on population comparisons that make sense based on phylogeny
  filter(pop_pair %in% c("EarlyCultivar_ImprovedCultivar",
                         "EarlyCultivar_Wild",
                         "LandraceA_Wild",
                         "LandraceB_Wild",
                         "LandraceA_EarlyCultivar",
                         "LandraceB_EarlyCultivar"))

# same for bootstrapped results 
bootstrap_split_times_filtered <- bootstrap_split_times %>% 
  mutate(pop_pair = word(pop_pair, 1, 2, sep = "_"),
         pop_pair = fct_reorder(pop_pair, split_years, .fun = min)) %>% 
  separate(pop_pair, into = c('pop1', 'pop2'), sep = '_', remove = FALSE) %>% 
  filter(!pop2 == 'LAwild') %>% 
  # set preferred order for pop 2 (helps with plotting) 
  mutate(pop2 = fct_relevel(pop2, 
                            'Wild', 
                            'EarlyCultivar', 
                            'ImprovedCultivar')) %>% 
  # focus only on population comparisons that make sense based on phylogeny
  filter(pop_pair %in% c("EarlyCultivar_ImprovedCultivar",
                         "EarlyCultivar_Wild",
                         "LandraceA_Wild",
                         "LandraceB_Wild",
                         "LandraceA_EarlyCultivar",
                         "LandraceB_EarlyCultivar"))
```

## 3.4 Summarize split times

``` r
split_times_filtered %>% 
  group_by(pop_pair) %>% 
  distinct(split_years) %>% 
  arrange(-split_years)
```

    ## # A tibble: 6 × 2
    ## # Groups:   pop_pair [6]
    ##   pop_pair                       split_years
    ##   <fct>                                <dbl>
    ## 1 EarlyCultivar_Wild                  19861.
    ## 2 LandraceA_Wild                      13262.
    ## 3 LandraceB_Wild                       8755.
    ## 4 LandraceA_EarlyCultivar              1927.
    ## 5 EarlyCultivar_ImprovedCultivar       1509.
    ## 6 LandraceB_EarlyCultivar              1131.

``` r
bootstrap_split_times_filtered %>% 
  group_by(pop_pair) %>% 
  summarize_at('split_years', .funs = c(min = min, max = max, median = median, sd = sd)) %>% 
  arrange(-median)
```

    ## # A tibble: 6 × 5
    ##   pop_pair                         min    max median     sd
    ##   <fct>                          <dbl>  <dbl>  <dbl>  <dbl>
    ## 1 EarlyCultivar_Wild             9554. 79773. 14971. 21039.
    ## 2 LandraceA_Wild                 6037. 25010. 10845.  7307.
    ## 3 LandraceB_Wild                 6326. 51566. 10804. 17479.
    ## 4 LandraceA_EarlyCultivar        2137. 83350.  2804. 25532.
    ## 5 LandraceB_EarlyCultivar        1254. 48834.  1998. 14856.
    ## 6 EarlyCultivar_ImprovedCultivar  543. 17671.   788.  7146.

## 3.5 Plot results from SMC++ split (Fig. 3d)

This step plots results from `smcpp split` to summarize estimates of
divergence time between populations.

``` r
# pull median split time from bootstrap estimates
median_bootstrap_split_times <- bootstrap_split_times_filtered %>% 
  group_by(pop_pair) %>% 
  summarize_at(.vars = 'split_years', median) 

median_bootstrap_split_times
```

    ## # A tibble: 6 × 2
    ##   pop_pair                       split_years
    ##   <fct>                                <dbl>
    ## 1 EarlyCultivar_ImprovedCultivar        788.
    ## 2 LandraceB_EarlyCultivar              1998.
    ## 3 LandraceA_EarlyCultivar              2804.
    ## 4 LandraceA_Wild                      10845.
    ## 5 LandraceB_Wild                      10804.
    ## 6 EarlyCultivar_Wild                  14971.

``` r
# create the plot 
split_plot <- ggplot() + 
  # add hypothesized window of carrot domestication 
  annotate('rect',
         ymin = 522,
         ymax = 1100,
         xmin = -Inf,
         xmax = Inf,
         alpha = 0.5,
         fill = 'gray') +
  # add violin plots for distribution of split times from bootstrapping
  geom_violin(data = bootstrap_split_times_filtered,
              aes(pop_pair, split_years, fill = pop2),
              scale = "width",
              alpha = 0.8) +
  # add points for actual values of split times from bootstrapping (10 reps)
  geom_point(data = bootstrap_split_times_filtered,
           aes(pop_pair, split_years),
           alpha = 0.2,
           color = 'black') + 
  # add crosses for median split time from bootstrapping 
  geom_point(data = median_bootstrap_split_times,
             aes(pop_pair, split_years),
             color = 'black',
             size = 2.5,
             shape = 4) + 
  # add diamonds for marginal estimates of split time from non-bootstrapped samples 
  geom_point(data = split_times_filtered,
             aes(pop_pair, split_years),
             size = 2.5,
             shape = 23,
             fill = 'black',
             color = 'black') + 
  # use a log10 scale for x and y axes and set number of breaks 
  # note: set min to 100 years for time axis 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(100,100000)) +
  annotation_logticks(sides = "b") +
  scale_color_manual(values = cols2) + 
  scale_fill_manual(values = cols2) + 
  coord_flip() +
  xlab("") +
  ylab("Split Time (years ago)") +
  theme_classic() +
  theme(legend.position = 'none')

split_plot
```

![](summarize_smc_results_files/figure-gfm/plot%20smcpp%20split%20times-1.png)<!-- -->

# 4 Combine panels for Fig 3c/d

``` r
plot_grid(estimate_plot, split_plot,
          rel_widths = c(0.5, 0.5))
```

![](summarize_smc_results_files/figure-gfm/combine%20plots-1.png)<!-- -->

``` r
# save figure to a png
ggsave("results/figures/2023-06-08_carrot_demography_results.png",
       width = 10,
       height = 4)

# save figure to a pdf 
ggsave("results/figures/2023-06-08_carrot_demography_results.pdf",
       width = 10,
       height = 4)
```
