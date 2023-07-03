"
Authors: Arran J. Davis
Emails: arran.davis@anthro.ox.ac.uk | davis.arran@gmail.com
Affiliation: Social Body Lab, Institute of Human Sciences, University of Oxford
Date: 25 August 2022
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

### MODEL DATA CREATION - TIME TO TASK FAILURE BY EXPERIMENTAL CONDITION ###

#number of participants and observations per condition
participant_n = 25
item_n = 1

#model intercept
overall_intercept = 6

#fixed effects in minutes (around 15% longer performances in the social support condition)
condition_fixed_effect = 1.5
order_fixed_effect = 1.5

#random intercept SD for participants (grouping is at participant level)
random_intercept_sd = 1 * 2.5
condition_random_slope_sd = .3 * 2.5
order_random_slope_sd = .3 * 2.5

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

#save the data
write.csv(t_data, "../data/social_support_effect_on_ttf_sample_plot_data.csv", row.names = FALSE)                             

### PLOT DATA ###

#prepare data for plotting (sort by condition and participant)
t_data_sorted = t_data[with(t_data, order(condition, participant)),]
time_plot_data = data_1x1(array_1  = t_data_sorted$time[1:26],
                          array_2  = t_data_sorted$time[27:52],
                          jit_distance = .09,
                          jit_seed = 321)

#plot and save the data (citation: https://github.com/jorvlan/open-visualizations)
ttf_plot = raincloud_1x1_repmes(data = time_plot_data,
                                colors = (c('#56B4E9', '#009E73')),
                                fills = (c('#56B4E9', '#009E73')),
                                line_color = 'gray',
                                line_alpha = 2,
                                size = 6,
                                alpha = .6,
                                align_clouds = FALSE) +
  scale_y_continuous(label = comma, breaks = c(seq(0, 20, 5)), limits = c(0, 20)) +
  scale_x_continuous(breaks = c(1,2), labels = c("Control", "Support")) +
  xlab("Condition") +  
  ylab("Time to task failure (min)") +
  avenir_theme

ggsave("../plots/ttf_by_condition.jpg", ttf_plot, height = 7.5, width = 7.5)

################################################################################################################################################

### DATA SIMULATIONS FOR POWER ANALYSIS OF EFFECT OF SOCIAL SUPPORT ON TIME TO TASK FAILURE ###

#create a function for power simulations (based on above data creation procedure)
power_simulation_t = function(participant_n = 40, item_n = 1,
                              overall_intercept = 6, random_intercept_sd = 1,
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

#set the social support effects (overall mean from Study 1)
overall_mean = (9.69 + 8.15) / 2
five_percent = overall_mean * 0.05
social_support_effects = c(five_percent * .5, five_percent, five_percent * 1.5, five_percent * 2, five_percent * 2.5, five_percent * 3)

#run repeated simulations with dataset variants (variables defined above)
simulations = crossing(replications = 1:1000,
                       participant_n = c(15, 20, 25, 30, 35, 40), item_n = 1,
                       overall_intercept = overall_intercept, random_intercept_sd = random_intercept_sd,
                       condition_fixed_effect = social_support_effects, condition_random_slope_sd = condition_random_slope_sd,
                       order_fixed_effect = order_fixed_effect, order_random_slope_sd = order_random_slope_sd,
                       random_effects_correlation = 0, error_sd = error_sd) %>%
  mutate(analysis = pmap(., power_simulation_t)) %>%
  unnest(analysis)

#save simulation data
write.csv(simulations, "../data/time_social_support_power_analysis_data_sim.csv", row.names = FALSE)                             

################################################################################################################################################

### PLOT POWER ANALYSIS OF EFFECT OF SOCIAL SUPPORT ON TIME TO TASK FAILURE ###

#load the simulation data
simulations = read.csv("../data/time_social_support_power_analysis_data_sim.csv")

#create dataset to plot power analysis for condition main effect
sim_results_support_t = filter(simulations, effect == "fixed", term == "conditionsupport") %>%
                          dplyr::group_by(condition_fixed_effect, participant_n) %>% 
                          dplyr::summarise(power = mean(p.value < .05), .groups = "drop")

#order the factor
unique(sim_results_support_t$condition_fixed_effect)
sim_results_support_t$condition_fixed_effect = ordered(sim_results_support_t$condition_fixed_effect, c("0.223", "0.446", "0.669", "0.892", "1.115", "1.338"))

#plot power analysis
pa_support_t = ggplot(aes(condition_fixed_effect, participant_n, fill = power), data = sim_results_support_t) +
                geom_tile() +
                geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
                scale_fill_viridis_c(name = "Power",
                                     limits = c(0, 1), 
                                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                scale_x_discrete(labels = c("2.5%", "5%", "7.5%", "10%", "12.5%", "15%")) + 
                scale_y_continuous(breaks = c(15, 20, 25, 30, 35, 40)) +
                xlab("Percentage increase in time to task failure in social support condition") + 
                ylab("Participant sample size") +
                avenir_theme

#save the plot
ggsave("../plots/power_analysis_social_support_effects_on_ttf.jpg", pa_support_t, height = 7.5, width = 10)

################################################################################################################################################

### MODEL DATA CREATION FOR EXPERIMENTAL CONDITION BY TIME INTERACTION EFFECT ON PERCEIVED EFFORT ###

#number of participants and observations per condition
participant_n = 25
item_n = 1

#model intercept
overall_intercept = 3.4

#fixed effects (take from model run on Study 1 results)
condition_fixed_effect = 0.010
order_fixed_effect = 0.004
time_fixed_effect = 0.43
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

#create dataframe
data = add_random(participant = participant_n, condition = 1) %>%
       add_within("participant", condition = c("control", "support")) %>%
       add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
       add_within("condition.y", session  = c("first", "second")) %>%
       add_within("participant", time = seq(1, 15, 1)) %>%
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

#save the data
write.csv(data, "../data/social_support_and_time_interaction_on_effort_sample_plot_data.csv", row.names = FALSE)                             

#run model
time_condition_int_mod = lmer(effort ~ condition * time + session_number + (1 | participant), data = data)
summary(time_condition_int_mod)

### ### ###

#plot the model data
effort_plot = ggplot(data = data, aes(x = as.factor(time), y = effort, fill = condition)) +
                geom_boxplot() +
                scale_y_continuous(breaks = c(seq(0,10,1)), limits = c(0, 10)) +
                ylab("Perceived effort") +
                xlab("Measurement number") + 
                scale_fill_discrete(name = "Condition") +
                avenir_theme

ggsave("../plots/effort_by_condition_power_analysis_data.jpg", effort_plot, height = 7.5, width = 10)

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
                            scale_x_continuous(limits = c(1, 15),
                                               breaks = c(seq(1, 15, 1))) +
                            scale_y_continuous(breaks = seq(2, 10, 1)) + 
                            xlab("Measurement number") + 
                            ylab("Perceived effort") +
                            avenir_theme

#save the plot and the data
ggsave("../plots/effort_time_condition_interaction_data_simulation.jpg", time_cond_int_simulation, width = 10, height = 7.5)
write.csv(data, "../data/effort_time_condition_interaction_sim_model_data.csv", row.names = FALSE)   

################################################################################################################################################

### DATA SIMULATIONS FOR POWER ANALYSIS OF EXPERIMENTAL CONDITION BY TIME INTERACTION EFFECT ON PERCEIVED EFFORT ###

#number of participants and observations per condition
participant_n = 25
item_n = 1

#model intercept
overall_intercept = 3.4

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
effort_simulation_cond_time_int = function(participant_n = 25, item_n = 1,
                                           overall_intercept = overall_intercept, random_intercept_sd = random_intercept_sd,
                                           condition_fixed_effect = condition_fixed_effect, condition_random_slope_sd = condition_random_slope_sd,
                                           order_fixed_effect = order_fixed_effect, order_random_slope_sd = order_random_slope_sd,
                                           time_fixed_effect = time_fixed_effect, time_random_slope = 100,
                                           cond_time_interaction_fixed_effect = cond_time_interaction_fixed_effect, cond_time_interaction_sd = cond_time_interaction_sd,
                                           random_effects_correlation = random_effects_correlation, error_sd = error_sd,
                                           ...) {
  
  #create dataframe
  data = add_random(participant = participant_n, condition = 1) %>%
    add_within("participant", condition = c("control", "support")) %>%
    add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
    add_within("condition.y", session  = c("first", "second")) %>%
    add_within("participant", time = seq(1, 15, 1)) %>%
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

#set the social support by early life adversity interaction effect
cond_time_interaction_fixed_effects = c(-0.02, -0.03, -0.04, -0.05, -0.06)

#run repeated simulations with dataset variants
simulations = crossing(replications = 1:1000,
                       participant_n = c(15, 20, 25, 30, 35, 40), item_n = 1,
                       overall_intercept = overall_intercept, random_intercept_sd = random_intercept_sd,
                       condition_fixed_effect = condition_fixed_effect, condition_random_slope_sd = condition_random_slope_sd,
                       order_fixed_effect = order_fixed_effect, order_random_slope_sd = order_random_slope_sd,
                       time_fixed_effect = time_fixed_effect, time_random_slope = time_random_slope,
                       cond_time_interaction_fixed_effect = cond_time_interaction_fixed_effects, cond_time_interaction_sd = cond_time_interaction_sd,
                       random_effects_correlation = random_effects_correlation, error_sd = error_sd) %>%
  mutate(analysis = pmap(., effort_simulation_cond_time_int)) %>%
  unnest(analysis)

#save simulation data
write.csv(simulations, "../data/effort_social_support_time_interaction_power_analysis_data_sim.csv", row.names = FALSE)   

################################################################################################################################################

### PLOT POWER ANALYSIS OF EXPERIMENTAL CONDITION BY TIME INTERACTION EFFECT ON PERCEIVED EFFORT ###

#load the simulation data
simulations = read.csv("../data/effort_social_support_time_interaction_power_analysis_data_sim.csv")

#create dataset to plot power analysis for social support by early life adversity effect on time to exhaustion
effort_support_time_int_simulation = filter(simulations, effect == "fixed", term == "conditionsupport:time") %>% 
                                      dplyr::group_by(cond_time_interaction_fixed_effect, participant_n) %>% 
                                      dplyr::summarise(power = mean(p.value < .05), .groups = "drop")

#order the factor
unique(effort_support_time_int_simulation$cond_time_interaction_fixed_effect)
effort_support_time_int_simulation$cond_ela_interaction_fixed_effect = ordered(effort_support_time_int_simulation$cond_time_interaction_fixed_effect, 
                                                                         c("-0.02", "-0.03", "-0.04", "-0.05", "-0.06"))

#plot power analysis
x_axis_text = expression(bold(paste("Model ", bolditalic("b"), "-coefficient for experimental condition by measurement number interaction on perceived effort")))
effort_support_time = ggplot(aes(cond_ela_interaction_fixed_effect, participant_n, fill = power), data = effort_support_time_int_simulation) +
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
ggsave("../plots/power_analysis_social_support_by_time_interaction_on_perceived_effort.jpg", effort_support_time, height = 7.5, width = 10)
