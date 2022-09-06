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

### MODEL DATA CREATION - CUMULATIVE WORK BY EXPERIMENTAL CONDITION ###

#number of participants and observations per condition
participant_n = 40
item_n = 1

#model intercept
overall_intercept = 37800

#fixed effects (2.5%)
condition_fixed_effect = 378 * 2.5
order_fixed_effect = 378 * 2.5

#random intercept SD for participants (grouping is at participant level)
random_intercept_sd = 2500
condition_random_slope_sd = 300
order_random_slope_sd = 300

#correlation between random effects
random_effects_correlation = 0   

#error term standard deviation
error_sd = 500

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
w_data = data %>% mutate(w = overall_intercept + rand_int + 
                           ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                           ((order_fixed_effect + rand_slope_order) * session_number) +
                           error_term)

### PLOT DATA ###

#plot the trimmed dataset
w_condition_simulation = ggplot(w_data, aes(condition, w)) + 
                          geom_violin() +
                          geom_boxplot(show.legend = FALSE) +
                          scale_x_discrete(labels = c("Control", "Support")) +
                          scale_y_continuous(labels = scales::comma,
                                             breaks = c(seq(33000, 47000, 1000))) + 
                          xlab("Experimental condition") + 
                          ylab("Cumulative work (W)") +
                          avenir_theme

ggsave("../plots/w_by_condition_data_simulation.jpg", w_condition_simulation, width = 10, height = 7.5)

################################################################################################################################################

### DATA SIMULATIONS FOR POWER ANALYSIS OF EFFECT OF SOCIAL SUPPORT ON CUMULATIVE WORK ###

#create a function for power simulations (based on above data creation procedure)
power_simulation_w = function(participant_n = 40, item_n = 1,
                              overall_intercept = 37800, random_intercept_sd = 2500,
                              condition_fixed_effect = 378 * 2.5, condition_random_slope_sd = 300,
                              order_fixed_effect = 378 * 2.5, order_random_slope_sd = 300,
                              random_effects_correlation = 0, error_sd = 500,
                                ...) {
  
  #number of participants and observations per condition
  participant_n = participant_n
  item_n = item_n
  
  #model intercept
  overall_intercept = overall_intercept
  
  #fixed effects (2.5%)
  condition_fixed_effect = condition_fixed_effect
  order_fixed_effect = order_fixed_effect
  
  #random intercept SD for participants (grouping is at participant level)
  random_intercept_sd = 2500
  condition_random_slope_sd = 300
  order_random_slope_sd = 300
  
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
  data = data %>% mutate(W = overall_intercept + rand_int + 
                        ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                        ((order_fixed_effect + rand_slope_order) * session_number) +
                        error_term)
  
  #run model
  w_mod = lmer(W ~ condition + session_number + (1 | participant), data = data)
  
  #return a dataframe of the model results
  broom.mixed::tidy(w_mod)

}

### ### ###

#set the social support effects
overall_mean = 37800
one_percent = overall_mean * 0.01
social_support_effects = c(one_percent * .5, one_percent, one_percent * 1.5, one_percent * 2, one_percent *2.5)

#run repeated simulations with dataset variants
simulations = crossing(replications = 1:1000,
                       participant_n = c(25, 30, 35, 40, 45), item_n = 1,
                       overall_intercept = 37800, random_intercept_sd = 2500,
                       condition_fixed_effect = social_support_effects, condition_random_slope_sd = 300,
                       order_fixed_effect = c(378*2.5), order_random_slope_sd = 300,
                       random_effects_correlation = 0, error_sd = 500) %>%
  mutate(analysis = pmap(., power_simulation_w)) %>%
  unnest(analysis)

#save simulation data
write.csv(simulations, "../data/w_social_support_power_analysis_data_sim.csv", row.names = FALSE)                             

################################################################################################################################################

### PLOT POWER ANALYSIS OF EFFECT OF SOCIAL SUPPORT ON CUMULATIVE WORK ###

#load the simulation data
simulations = read.csv("../data/w_social_support_power_analysis_data_sim.csv")

#create dataset to plot power analysis for condition main effect
sim_results_support_w = filter(simulations, effect == "fixed", term == "conditionsupport") %>%
                          dplyr::group_by(condition_fixed_effect, participant_n) %>% 
                          dplyr::summarise(power = mean(p.value < .05), .groups = "drop")

#order the `condition_fixed_effect` factor
sim_results_support_w$condition_fixed_effect = ordered(sim_results_support_w$condition_fixed_effect, c("189", "378", "567", "756", "945"))

#plot power analysis
pa_support_w = ggplot(aes(condition_fixed_effect, participant_n, fill = power), data = sim_results_support_w) +
                geom_tile() +
                geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
                scale_fill_viridis_c(name = "Power",
                                     limits = c(0, 1), 
                                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                scale_x_discrete(labels = c("0.5%", "1.0%", "1.5%", "2.0%", "2.5%")) + 
                scale_y_continuous(breaks = c(25, 30, 35, 40, 45)) +
                xlab("Percentage increase in cumulative work in social support condition") + 
                ylab("Participant sample size") +
                avenir_theme

#save the plot
ggsave("../plots/power_analysis_social_support_effects_on_cumulative_work.jpg", pa_support_w, height = 7.5, width = 10)

################################################################################################################################################

### MODEL DATA CREATION FOR EXPERIMENTAL CONDITION BY DEVELOPMENTAL EXPERIENCES OF SOCIAL SUPPORT INTERACTION EFFECT ON CUMULATIVE WORK ###

#number of participants and observations per condition
participant_n = 40
item_n = 1

#model intercept
overall_intercept = 37800

#fixed effects (2.5%)
condition_fixed_effect = 378 * 2.5
order_fixed_effect = 378 * 2.5
ela_fixed_effect = 0.0
cond_ela_interaction_fixed_effect = 250

#random intercept SD for participants (grouping is at participant level)
random_intercept_sd = 2500
condition_random_slope_sd = 300
order_random_slope_sd = 300
ela_random_slope = 100
cond_ela_interaction_sd = 50

#correlation between random effects
random_effects_correlation = 0   

#error term standard deviation
error_sd = 200

#create dataframe
data = add_random(participant = participant_n, condition = 1) %>%
       add_within("participant", condition = c("control", "support")) %>%
       add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
       add_within("condition.y", session  = c("first", "second")) %>%
       add_between("participant", early_life_adversity = c(rnorm(40, 0, 1))) %>%
       add_recode("condition.y", "condition_numeric", control = 0, support = 1) %>%
       add_recode("session", "session_number", first = 1, second = 2) %>%
       add_ranef("participant", rand_int = random_intercept_sd, 
                 rand_slope_cond = condition_random_slope_sd,
                 rand_slope_order = order_random_slope_sd,
                 rand_slope_ela = ela_random_slope,
                 rand_slope_interaction = cond_ela_interaction_sd,
                .cors = random_effects_correlation) %>%
       add_ranef(error_term = error_sd)

#remove duplicate rows (this is how the counterbalancing is created)
data = data %>% dplyr::filter((condition.y == "control" & part_odd_even != session_number | 
                                 (condition.y == "support" & part_odd_even == session_number)))

#rename outcome and transform variables
data$condition = data$condition.y
data$participant = as.factor(data$participant)
data$early_life_adversity = as.numeric(as.character(data$early_life_adversity))

#drop unneeded columns
data = data[ , c("participant", "condition", "condition_numeric", "session_number", "early_life_adversity", 
                 "rand_int", "rand_slope_cond", "rand_slope_order", "rand_slope_ela", "rand_slope_interaction",  "error_term")]

#create outcome variable (assuming no effect of age or gender)
data = data %>% mutate(w = overall_intercept + rand_int + 
                         ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                         ((order_fixed_effect + rand_slope_order) * session_number) +
                         ((ela_fixed_effect + rand_slope_ela) * early_life_adversity) +
                         ((cond_ela_interaction_fixed_effect + rand_slope_interaction) * condition_numeric * early_life_adversity) +
                         error_term)

#run model
w_ela_int_mod = lmer(w ~ condition * early_life_adversity + session_number + (1 | participant), data = data)
summary(w_ela_int_mod)

#get model data
plot_dat = plot_model(w_ela_int_mod, type = "pred", terms = c("early_life_adversity", "condition"))
plot_dat = plot_dat$data
plot_dat$Condition = plot_dat$group_col
plot_dat$Condition = ifelse(plot_dat$Condition == "control", "Control",
                            ifelse(plot_dat$Condition =="support", "Support", NA))

### ### ###

#load model data (so that the description of the model results in the PDF document are consistent)
plot_dat = read.csv("../data/w_social_support_desp_interaction_data_sim_model_data_for_plotting.csv")

#plot the model
w_ela_int_simulation = ggplot(plot_dat, aes(x = x, y = predicted, group = Condition)) + 
                        geom_line(aes(color = Condition)) +
                        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
                        scale_x_continuous(limits = c(-2, 2),
                                           breaks = c(seq(-2, 2, .5)), 
                                           labels = c("Low", " ", " ", " ", "", "", " ", " ", "High")) +
                        scale_y_continuous(labels = scales::comma,
                                           breaks = c(seq(37000, 43000, 1000))) + 
                        xlab("Developmental exerperiences of social support") + 
                        ylab("Cumulative work (W)") +
                        guides(colour = guide_legend(title="Experimental condition")) +
                        avenir_theme

#save the plot and the data
ggsave("../plots/w_condition_desp_interaction_data_simulation.jpg", w_ela_int_simulation, width = 10, height = 7.5)
write.csv(plot_dat, "../data/w_social_support_desp_interaction_data_sim_model_data_for_plotting.csv", row.names = FALSE)   

################################################################################################################################################

### DATA SIMULATIONS FOR POWER ANALYSIS OF EXPERIMENTAL CONDITION BY DEVELOPMENTAL EXPERIENCES OF SOCIAL SUPPORT INTERACTION EFFECT ON CUMULATIVE WORK ###

#create a function for power simulations (based on above data creation procedure)
power_simulation_int_w = function(participant_n = 40, item_n = 1,
                                  overall_intercept = 37800, random_intercept_sd = 2500,
                                  condition_fixed_effect = 945, condition_random_slope_sd = 300,
                                  order_fixed_effect = 945, order_random_slope_sd = 300,
                                  ela_fixed_effect = 0.0, ela_random_slope = 100,
                                  cond_ela_interaction_fixed_effect = -200, cond_ela_interaction_sd = 50,
                                  random_effects_correlation = 0, error_sd = 10,
                                    ...) {
  
  #create dataframe
  data = add_random(participant = participant_n, condition = 1) %>%
    add_within("participant", condition = c("control", "support")) %>%
    add_between("participant", part_odd_even = c(2, 1), .shuffle = TRUE) %>%
    add_within("condition.y", session  = c("first", "second")) %>%
    add_between("participant", early_life_adversity = c(rnorm(40, 0, 1))) %>%
    add_recode("condition.y", "condition_numeric", control = 0, support = 1) %>%
    add_recode("session", "session_number", first = 1, second = 2) %>%
    add_ranef("participant", rand_int = random_intercept_sd, 
              rand_slope_cond = condition_random_slope_sd,
              rand_slope_order = order_random_slope_sd,
              rand_slope_ela = ela_random_slope,
              rand_slope_interaction = cond_ela_interaction_sd,
              .cors = random_effects_correlation) %>%
    add_ranef(error_term = error_sd)
  
  #remove duplicate rows (this is how the counterbalancing is created)
  data = data %>% dplyr::filter((condition.y == "control" & part_odd_even != session_number | 
                                   (condition.y == "support" & part_odd_even == session_number)))
  
  #rename outcome and transform variables
  data$condition = data$condition.y
  data$participant = as.factor(data$participant)
  data$early_life_adversity = as.numeric(as.character(data$early_life_adversity))
  
  #drop unneeded columns
  data = data[ , c("participant", "condition", "condition_numeric", "session_number", "early_life_adversity", 
                   "rand_int", "rand_slope_cond", "rand_slope_order", "rand_slope_ela", "rand_slope_interaction",  "error_term")]
  
  #create outcome variable (assuming no effect of age or gender)
  data = data %>% mutate(w = overall_intercept + rand_int + 
                        ((condition_fixed_effect + rand_slope_cond) * condition_numeric) +
                        ((order_fixed_effect + rand_slope_order) * session_number) +
                        ((ela_fixed_effect + rand_slope_ela) * early_life_adversity) +
                        ((cond_ela_interaction_fixed_effect + rand_slope_interaction) * condition_numeric * early_life_adversity) +
                        error_term)
  
  #run model
  w_ela_int_mod = lmer(w ~ condition * early_life_adversity + session_number + (1 | participant), data = data)
  
  #return a dataframe of the model results
  broom.mixed::tidy(w_ela_int_mod)
  
}

### ### ###

#set the social support by early life adversity interaction effect
cond_ela_interaction_fixed_effect = c(100, 150, 200, 250, 300)

#run repeated simulations with dataset variants
simulations = crossing(replications = 1:1000,
                       participant_n = c(25, 30, 35, 40, 45), item_n = 1,
                       overall_intercept = 37800, random_intercept_sd = 2500,
                       condition_fixed_effect = 945, condition_random_slope_sd = 300,
                       order_fixed_effect = 945, order_random_slope_sd = 300,
                       ela_fixed_effect = 0, ela_random_slope = 100,
                       cond_ela_interaction_fixed_effect = cond_ela_interaction_fixed_effect, cond_ela_interaction_sd = 50,
                       random_effects_correlation = 0, error_sd = 200) %>%
  mutate(analysis = pmap(., power_simulation_int_w)) %>%
  unnest(analysis)

#save simulation data
write.csv(simulations, "../data/w_social_support_ela_interaction_power_analysis_data_sim.csv", row.names = FALSE)   

################################################################################################################################################

### PLOT POWER ANALYSIS OF EXPERIMENTAL CONDITION BY DEVELOPMENTAL EXPERIENCES OF SOCIAL SUPPORT INTERACTION EFFECT ON CUMULATIVE WORK ###

#load the simulation data
simulations = read.csv("../data/w_social_support_desp_interaction_power_analysis_data_sim.csv")

#create dataset to plot power analysis for social support by early life adversity effect on time to exhaustion
w_support_ela_int_simulation = filter(simulations, effect == "fixed", term == "conditionsupport:early_life_adversity") %>% 
                                dplyr::group_by(cond_ela_interaction_fixed_effect, participant_n) %>% 
                                dplyr::summarise(power = mean(p.value < .05), .groups = "drop")

#order the `cond_ela_interaction_fixed_effect` factor
w_support_ela_int_simulation$cond_ela_interaction_fixed_effect = ordered(w_support_ela_int_simulation$cond_ela_interaction_fixed_effect, 
                                                                         c("100", "150", "200", "250", "300"))

#plot power analysis
x_axis_text = expression(bold(paste("Model ", bolditalic("b"), "-coefficient for experimental condition by developmental experiences of social support interaction on cumulative work")))
pa_support_ela_w = ggplot(aes(cond_ela_interaction_fixed_effect, participant_n, fill = power), data = w_support_ela_int_simulation) +
                      geom_tile() +
                      geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
                      scale_y_continuous(breaks = c(25, 30, 35, 40, 45)) +
                      scale_fill_viridis_c(name = "Power",
                                           limits = c(0, 1), 
                                           breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                           labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                      xlab(x_axis_text) + 
                      ylab("Participant sample size") +
                      avenir_theme

#save the plot
ggsave("../plots/power_analysis_social_support_by_desp_interaction_effect_on_cumulative_work.jpg", pa_support_ela_w, height = 7.5, width = 10)
