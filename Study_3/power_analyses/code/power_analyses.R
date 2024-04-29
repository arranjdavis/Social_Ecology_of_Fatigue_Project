"
Authors: Arran J. Davis
Emails: arran.davis@anthro.ox.ac.uk | davis.arran@gmail.com
Affiliation: Social Body Lab, Institute of Human Sciences, University of Oxford
Date: 8 April 2024
"

library(extrafont)
library(faux)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(tidyr)
library(broom.mixed)
library(purrr)
library(censReg)
library(plm)
library(raincloudplots)

#clean environment
rm(list = ls())

#set working directory
setwd(getSrcDirectory()[1])
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################################################################################

### GRAPH THEMES ###

#fonts
quartzFonts(avenir = c("Avenir Book", "Avenir Black", "Avenir Book Oblique", "Avenir Black Oblique"))

#theme
header_size = 10
axis_size = 10

#theme for plots
avenir_theme = theme(text=element_text(size=header_size,family='avenir'),
                     axis.text.x = element_text(color = 'black', size = axis_size, vjust = 1),
                     axis.text.y = element_text(color = 'black', size = axis_size),
                     axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"),
                     axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 0, l = 0), face = "bold"),
                     panel.background = element_blank(),
                     panel.grid.major.x = element_line(color = '#e7e7e7'),
                     panel.grid.major.y = element_line(color = '#e7e7e7'),
                     legend.key = element_blank(),
                     legend.title = element_text(color = 'black', size = axis_size, face = "bold"),
                     plot.title = element_text(hjust = 0.5, face = "bold", size = axis_size + 1))

#set the font
par(family = 'avenir')

################################################################################################################################################

### MODEL DATA CREATION - 2 KM TIMES BY EXPERIMENTAL CONDITION ###

#number of participants and observations per condition
participant_n = 20
item_n = 1

#model intercept
overall_intercept = 8

#calculate a 2.5% difference over an eight-minute trial (it is 12 seconds, or 20% of a minute)
t1 = 60 * 8
t1_p = t1 * 0.025
improvement = 1 - ((t1 - t1_p) / t1)

#fixed effects in minutes (around 5% faster 5 km times the social support condition)
condition_fixed_effect = -0.2
order_fixed_effect = -0.05

#random intercept SD for participants (grouping is at participant level)
random_intercept_sd = 1 * 0.1
condition_random_slope_sd = .3 * 0.05
order_random_slope_sd = .3 * 0.05

#correlation between random effects
random_effects_correlation = 0   

#error term standard deviation
error_sd = .1

#create dataframe
data = add_random(participant = participant_n, condition = 1) %>%
        add_within("participant", condition = c("control", "support")) %>%
        add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
        add_within("condition.y", session  = c("first", "second")) %>%
        add_between("participant", gender = c("female", "male"), .prob = c(.3, .7), .shuffle = TRUE) %>%
        add_between("participant", age = 15:18, .prob = c(.4, .3, .2, .1), .shuffle = TRUE) %>%
        add_recode("condition.y", "condition_numeric", control = 0, support = 1) %>%
        add_recode("session", "session_number", first = 1, second = 2) %>%
        add_ranef("participant", rand_int = random_intercept_sd, 
                  rand_slope_cond = condition_random_slope_sd,
                  rand_slope_order = order_random_slope_sd,
                  .cors = random_effects_correlation) %>%
        add_ranef(error_term = error_sd)

#remove duplicate rows (this is how the counterbalancing is created)
data = data %>% dplyr::filter((condition.y == "control" & part_odd_even != session_number | 
                                 (condition.y == "support" & part_odd_even == session_number)))

#rename outcome and transform variables
data$condition = data$condition.y
data$participant = as.factor(data$participant)

#drop unneeded columns
data = data[ , c("participant", "condition",  "condition_numeric", "session_number", "gender", "age", 
                 "rand_int", "rand_slope_cond", "rand_slope_order",  "error_term")]

#create outcome variable (assuming no effect of age or gender)
t_data = data %>% mutate(time = overall_intercept + rand_int + 
                         ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                         ((order_fixed_effect + rand_slope_order) * session_number) + error_term)

#check mean differences between conditions
t_data %>% group_by(condition) %>% summarise(mean(time))

#save the data
write.csv(t_data, "../data/social_support_effect_on_2km_times_sample_plot_data.csv", row.names = FALSE)                             

### PLOT DATA ###

#prepare data for plotting (sort by condition and participant)
t_data_sorted = t_data[with(t_data, order(condition, participant)),]
time_plot_data = data_1x1(array_1  = t_data_sorted$time[1:26],
                          array_2  = t_data_sorted$time[27:52],
                          jit_distance = .09,
                          jit_seed = 321)

#plot and save the data (citation: https://github.com/jorvlan/open-visualizations)
time_plot = raincloud_1x1_repmes(data = time_plot_data,
                                 colors = (c('#56B4E9', '#009E73')),
                                 fills = (c('#56B4E9', '#009E73')),
                                 line_color = 'gray',
                                 line_alpha = 2,
                                 size = 6,
                                 alpha = .6,
                                 align_clouds = FALSE) +
              scale_y_continuous(label = c("7:15", "7:30", "7:45", "8:00", "8:15"), breaks = c(seq(7.25, 8.25, 0.25)), limits = c(7.25, 8.25)) +
              scale_x_continuous(breaks = c(1,2), labels = c("Control", "Support")) +
              xlab("Experimental condition") +  
              ylab("Time for 2 km rowing trial (min)") +
              avenir_theme

ggsave("../plots/time_2km_by_condition.jpg", time_plot, height = 7.5, width = 7.5)

################################################################################################################################################

### DATA SIMULATIONS FOR POWER ANALYSIS OF EFFECT OF SOCIAL SUPPORT ON 2 KM TIMES ###

#create a function for power simulations (based on above data creation procedure)
power_simulation_t = function(participant_n = 40, item_n = 1,
                              overall_intercept = 8, random_intercept_sd = 1,
                              condition_fixed_effect = 1.5, condition_random_slope_sd = .3,
                              order_fixed_effect = 1.5 * .3, order_random_slope_sd = .3,
                              random_effects_correlation = 0, error_sd = .1,
                                ...) {
  
  #number of participants and observations per condition
  participant_n = participant_n
  item_n = item_n
  
  #model intercept
  overall_intercept = overall_intercept
  
  #fixed effects
  condition_fixed_effect = condition_fixed_effect
  order_fixed_effect = order_fixed_effect
  
  #random intercept SD for participants (grouping is at participant level)
  random_intercept_sd = random_intercept_sd
  condition_random_slope_sd = condition_random_slope_sd
  order_random_slope_sd = order_random_slope_sd
  
  #correlation between random effects
  random_effects_correlation = random_effects_correlation
  
  #error term standard deviation
  error_sd = error_sd
  
  #create dataframe
  data = add_random(participant = participant_n, condition = 1) %>%
                    add_within("participant", condition = c("control", "support")) %>%
                    add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
                    add_within("condition.y", session  = c("first", "second")) %>%
                    add_recode("condition.y", "condition_numeric", control = 0, support = 1) %>%
                    add_recode("session", "session_number", first = 1, second = 2) %>%
                    add_ranef("participant", rand_int = random_intercept_sd, 
                              rand_slope_cond = condition_random_slope_sd,
                              rand_slope_order = order_random_slope_sd,
                              .cors = random_effects_correlation) %>%
                    add_ranef(error_term = error_sd)
                  
  #remove duplicate rows (this is how the counterbalancing is created)
  data = data %>% dplyr::filter((condition.y == "control" & part_odd_even != session_number | 
                                (condition.y == "support" & part_odd_even == session_number)))
  
  #rename outcome and transform variables
  data$condition = data$condition.y
  data$participant = as.factor(data$participant)
  
  #drop unneeded columns
  data = data[ , c("participant", "condition",  "condition_numeric", "session_number", "rand_int", "rand_slope_cond", "rand_slope_order",  "error_term")]
  
  #create outcome variable (assuming no effect of age or gender)
  data = data %>% mutate(time = overall_intercept + rand_int + 
                        ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                        ((order_fixed_effect + rand_slope_order) * session_number) +
                        error_term)
  
  #run model
  time_mod = lmer(time ~ condition + session_number + (1 | participant), data = data)
  
  #return a dataframe of the model results
  broom.mixed::tidy(time_mod)

}

### ### ###

#calculate a 2.5% difference over an eight-minute trial (it is 12 seconds, or 20% of a minute)
t1 = 60 * 8
t1_p = t1 * 0.025
improvement = 1 - ((t1 - t1_p) / t1)

#set the social support effects (overall mean is a 2 km time of 8 minutes); effects given in minutes (2.5% 8 minutes is 12 seconds)
overall_mean = 8
improvement_percent = overall_mean * -0.025
social_support_effects = c(improvement_percent * .5, improvement_percent, improvement_percent * 1.5, improvement_percent * 2, improvement_percent * 2.5, improvement_percent * 3)

#run repeated simulations with dataset variants (variables defined above)
simulations = crossing(replications = 1:1000,
                       participant_n = c(15, 20, 25, 30, 35, 40), item_n = 1,
                       overall_intercept = 8, random_intercept_sd = 0.5,
                       condition_fixed_effect = social_support_effects, condition_random_slope_sd = 0.25,
                       order_fixed_effect = 0.5, order_random_slope_sd = 0.1,
                       random_effects_correlation = 0, error_sd = 0.1) %>%
  mutate(analysis = pmap(., power_simulation_t)) %>%
  unnest(analysis)

#save simulation data
write.csv(simulations, "../data/time_2km_social_support_power_analysis_data_sim.csv", row.names = FALSE)                             

################################################################################################################################################

### PLOT POWER ANALYSIS OF EFFECT OF SOCIAL SUPPORT ON 2 KM TIMES ###

#load the simulation data
simulations = read.csv("../data/time_2km_social_support_power_analysis_data_sim.csv")

#create dataset to plot power analysis for condition main effect
sim_results_support_t = filter(simulations, effect == "fixed", term == "conditionsupport") %>%
                          dplyr::group_by(condition_fixed_effect, participant_n) %>% 
                          dplyr::summarise(power = mean(p.value < .05), .groups = "drop")

#order the factor
unique(sim_results_support_t$condition_fixed_effect)
sim_results_support_t$condition_fixed_effect = ordered(sim_results_support_t$condition_fixed_effect, c("-0.1", "-0.2", "-0.3", "-0.4", "-0.5", "-0.6"))

#plot power analysis
pa_support_t = ggplot(aes(condition_fixed_effect, participant_n, fill = power), data = sim_results_support_t) +
                geom_tile() +
                geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
                scale_fill_viridis_c(name = "Power",
                                     limits = c(0, 1), 
                                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                scale_x_discrete(labels = c("-1.25%", "-2.5%", "-3.75%", "-5%", "-6.25%", "-7.5%")) + 
                scale_y_continuous(breaks = c(15, 20, 25, 30, 35, 40)) +
                xlab("Percentage decrease in 2 km rowing trial time in social support conditions") + 
                ylab("Participant sample size") +
                avenir_theme

#save the plot
ggsave("../plots/power_analysis_social_support_effects_on_2km_time.jpg", pa_support_t, height = 7.5, width = 10)

################################################################################################################################################

### MODEL DATA CREATION FOR EXPERIMENTAL CONDITION BY TIME INTERACTION EFFECT ON PERCEIVED FATIGUE ###

#number of participants and observations per condition
participant_n = 25
item_n = 1

#model intercept
overall_intercept = 1

#fixed effects (based off Study 1 results)
condition_fixed_effect = 0.010
order_fixed_effect = 0.004
time_fixed_effect = 0.3
cond_time_interaction_fixed_effect = -0.04

#random intercept SD for participants (grouping is at participant level)
random_intercept_sd = 0.5
condition_random_slope_sd = 0.2
order_random_slope_sd = 0.1
time_random_slope = 0.1
cond_time_interaction_sd = 0.02

#correlation between random effects
random_effects_correlation = 0   

#error term standard deviation
error_sd = .1

#number of observations (24 over 20 minutes of rowing)
n_obvs = 24

#create dataframe
data = add_random(participant = participant_n, condition = 1) %>%
       add_within("participant", condition = c("control", "support")) %>%
       add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
       add_within("condition.y", session  = c("first", "second")) %>%
       add_within("participant", time = seq(1, n_obvs, 1)) %>%
       add_recode("condition.y", "condition_numeric", control = 0, support = 1) %>%
       add_recode("session", "session_number", first = 1, second = 2) %>%
       add_ranef("participant", rand_int = random_intercept_sd, 
                 rand_slope_cond = condition_random_slope_sd,
                 rand_slope_order = order_random_slope_sd,
                 rand_slope_time = time_random_slope,
                 rand_slope_interaction = cond_time_interaction_sd,
                .cors = random_effects_correlation) %>%
       add_ranef(error_term = error_sd)

#remove duplicate rows (this is how the counterbalancing is created)
data = data %>% dplyr::filter((condition.y == "control" & part_odd_even != session_number | 
                                 (condition.y == "support" & part_odd_even == session_number)))

#rename outcome and transform variables
data$condition = data$condition.y
data$participant = as.factor(data$participant)
data$time = as.numeric(as.character(data$time))

#drop unneeded columns
data = data[ , c("participant", "condition", "condition_numeric", "session_number", "time", 
                 "rand_int", "rand_slope_cond", "rand_slope_order", "rand_slope_time", "rand_slope_interaction",  "error_term")]

#create outcome variable (assuming no effect of age or gender)
data = data %>% mutate(effort = overall_intercept + rand_int + 
                       ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                       ((order_fixed_effect + rand_slope_order) * session_number) +
                       ((time_fixed_effect + rand_slope_time) * time) +
                       ((cond_time_interaction_fixed_effect + rand_slope_interaction) * condition_numeric * time) +
                       error_term)

#round the effort data
data$effort = ifelse(data$effort > 10, 10, data$effort)
data$effort = ifelse(data$effort < 0, 0, data$effort)
data$effort = round(data$effort, 0)

#save the data
write.csv(data, "../data/social_support_and_time_interaction_on_effort_sample_plot_data.csv", row.names = FALSE)                             

#run model
time_condition_int_mod = lmer(effort ~ condition * time + session_number + (1 | participant), data = data)
summary(time_condition_int_mod)

### ### ###

#plot the model data
fatigue_plot = ggplot(data = data, aes(x = as.factor(time), y = effort, fill = condition)) +
                geom_boxplot() +
                scale_y_continuous(breaks = c(seq(0,10,1)), limits = c(0, 10)) +
                ylab("Perceived fatigue") +
                xlab("Measurement number") + 
                scale_fill_discrete(name = "Condition") +
                avenir_theme

ggsave("../plots/fatigue_by_time_condition_power_analysis_data.jpg", fatigue_plot, height = 7.5, width = 10)

### ### ###

#plot the model predictions
plot_dat = plot_model(time_condition_int_mod, type = "pred", terms = c("time", "condition"))
plot_dat = plot_dat$data
plot_dat$Condition = plot_dat$group_col
plot_dat$Condition = ifelse(plot_dat$Condition == "control", "Control",
                            ifelse(plot_dat$Condition =="support", "Support", NA))

time_cond_int_simulation = ggplot(plot_dat, aes(x = x, y = predicted, group = Condition)) + 
                            geom_line(aes(color = Condition)) +
                            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
                            scale_x_continuous(limits = c(1, 24),
                                               breaks = c(seq(1, 24, 1))) +
                            scale_y_continuous(breaks = seq(2, 10, 1)) + 
                            xlab("Measurement number") + 
                            ylab("Perceived fatigue") +
                            avenir_theme

#save the plot and the data
ggsave("../plots/fatigue_by_time_condition_interaction_model_predictions.jpg", time_cond_int_simulation, width = 10, height = 7.5)
write.csv(data, "../data/fatigue_time_condition_interaction_sim_model_data.csv", row.names = FALSE)   

################################################################################################################################################

### DATA SIMULATIONS FOR POWER ANALYSIS OF EXPERIMENTAL CONDITION BY TIME INTERACTION EFFECT ON PERCEIVED FATIGUE ###

#number of participants and observations per condition
participant_n = 25
item_n = 1

#model intercept
overall_intercept = 2

#fixed effects (take from model run on Study 1 results)
condition_fixed_effect = 0.010
order_fixed_effect = 0.004
time_fixed_effect = 0.43
cond_time_interaction_fixed_effect = -0.04

#random intercept SD for participants (grouping is at participant level)
random_intercept_sd = 0.5 * 1.5
condition_random_slope_sd = 0.2 * 1.5
order_random_slope_sd = 0.1 * 1.5
time_random_slope = 0.1 * 1.5
cond_time_interaction_sd = 0.02 * 1.5

#correlation between random effects
random_effects_correlation = 0   

#error term standard deviation
error_sd = .1

#create a function for power simulations (based on above data creation procedure)
fatigue_simulation_cond_time_int = function(participant_n = 25, item_n = 1,
                                            overall_intercept = overall_intercept, random_intercept_sd = random_intercept_sd,
                                            condition_fixed_effect = condition_fixed_effect, condition_random_slope_sd = condition_random_slope_sd,
                                            order_fixed_effect = order_fixed_effect, order_random_slope_sd = order_random_slope_sd,
                                            time_fixed_effect = time_fixed_effect, time_random_slope = 100,
                                            cond_time_interaction_fixed_effect = cond_time_interaction_fixed_effect, cond_time_interaction_sd = cond_time_interaction_sd,
                                            random_effects_correlation = random_effects_correlation, error_sd = error_sd,
                                            total_measurements = 20,
                                            ...) {
  
  #create dataframe
  data = add_random(participant = participant_n, condition = 1) %>%
    add_within("participant", condition = c("control", "support")) %>%
    add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
    add_within("condition.y", session  = c("first", "second")) %>%
    add_within("participant", time = seq(1, total_measurements, 1)) %>%
    add_recode("condition.y", "condition_numeric", control = 0, support = 1) %>%
    add_recode("session", "session_number", first = 1, second = 2) %>%
    add_ranef("participant", rand_int = random_intercept_sd, 
              rand_slope_cond = condition_random_slope_sd,
              rand_slope_order = order_random_slope_sd,
              rand_slope_time = time_random_slope,
              rand_slope_interaction = cond_time_interaction_sd,
              .cors = random_effects_correlation) %>%
    add_ranef(error_term = error_sd)
  
  #remove duplicate rows (this is how the counterbalancing is created)
  data = data %>% dplyr::filter((condition.y == "control" & part_odd_even != session_number | 
                                   (condition.y == "support" & part_odd_even == session_number)))
  
  #rename outcome and transform variables
  data$condition = data$condition.y
  data$participant = as.factor(data$participant)
  data$time = as.numeric(as.character(data$time))
  
  #drop unneeded columns
  data = data[ , c("participant", "condition", "condition_numeric", "session_number", "time", 
                   "rand_int", "rand_slope_cond", "rand_slope_order", "rand_slope_time", "rand_slope_interaction",  "error_term")]
  
  #create outcome variable (assuming no effect of age or gender)
  data = data %>% mutate(effort = overall_intercept + rand_int + 
                           ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                           ((order_fixed_effect + rand_slope_order) * session_number) +
                           ((time_fixed_effect + rand_slope_time) * time) +
                           ((cond_time_interaction_fixed_effect + rand_slope_interaction) * condition_numeric * time) +
                           error_term)
  
  #round the effort data
  data$effort = ifelse(data$effort > 10, 10, data$effort)
  data$effort = round(data$effort, 0)
  
  #run model
  time_condition_int_mod = lmer(effort ~ condition * time + session_number + (1 | participant), data = data)
  summary(time_condition_int_mod)
  
  #return a dataframe of the model results
  broom.mixed::tidy(time_condition_int_mod)
  
}

### ### ###

#set the size of the interaction effect
cond_time_interaction_fixed_effects = c(-0.02, -0.03, -0.04, -0.05, -0.06)

#run repeated simulations with dataset variants
simulations = crossing(replications = 1:1000,
                       participant_n = c(15, 20, 25, 30, 35, 40), item_n = 1,
                       overall_intercept = overall_intercept, random_intercept_sd = random_intercept_sd,
                       condition_fixed_effect = condition_fixed_effect, condition_random_slope_sd = condition_random_slope_sd,
                       order_fixed_effect = order_fixed_effect, order_random_slope_sd = order_random_slope_sd,
                       time_fixed_effect = time_fixed_effect, time_random_slope = time_random_slope,
                       cond_time_interaction_fixed_effect = cond_time_interaction_fixed_effects, cond_time_interaction_sd = cond_time_interaction_sd,
                       random_effects_correlation = random_effects_correlation, error_sd = error_sd,
                       total_measurements = 20) %>%
  mutate(analysis = pmap(., fatigue_simulation_cond_time_int)) %>%
  unnest(analysis)

#save simulation data
write.csv(simulations, "../data/fatigue_social_support_time_interaction_power_analysis_data_sim.csv", row.names = FALSE)   

################################################################################################################################################

### PLOT POWER ANALYSIS OF EXPERIMENTAL CONDITION BY TIME INTERACTION EFFECT ON PERCEIVED FATIGUE ###

#load the simulation data
simulations = read.csv("../data/fatigue_social_support_time_interaction_power_analysis_data_sim.csv")

#create dataset to plot power analysis for social support by early life adversity effect on time to exhaustion
fatigue_support_time_int_simulation = filter(simulations, effect == "fixed", term == "conditionsupport:time") %>% 
                                       dplyr::group_by(cond_time_interaction_fixed_effect, participant_n) %>% 
                                       dplyr::summarise(power = mean(p.value < .05), .groups = "drop")

#order the factor
unique(fatigue_support_time_int_simulation$cond_time_interaction_fixed_effect)
fatigue_support_time_int_simulation$cond_ela_interaction_fixed_effect = ordered(fatigue_support_time_int_simulation$cond_time_interaction_fixed_effect, 
                                                                                c("-0.02", "-0.03", "-0.04", "-0.05", "-0.06"))

#plot power analysis
x_axis_text = expression(bold(paste("Model ", bolditalic("b"), "-coefficient for experimental condition by measurement number interaction on perceived fatigue")))
fatigue_support_time = ggplot(aes(cond_ela_interaction_fixed_effect, participant_n, fill = power), data = fatigue_support_time_int_simulation) +
                        geom_tile() +
                        geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
                        scale_y_continuous(breaks = c(15, 20, 25, 30, 35, 40)) +
                        scale_fill_viridis_c(name = "Power",
                                           limits = c(0, 1), 
                                           breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                        xlab(x_axis_text) + 
                        ylab("Participant sample size") +
                        avenir_theme

#save the plot
ggsave("../plots/power_analysis_social_support_by_time_interaction_effect_on_perceived_fatigue.jpg", fatigue_support_time, height = 7.5, width = 10)
