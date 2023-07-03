"
Authors: Arran J. Davis
Emails: arran.davis@anthro.ox.ac.uk | davis.arran@gmail.com
Affiliation: Social Body Lab, Institute of Human Sciences, University of Oxford
Date: 01/08/2022
"

library(ggplot2)

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
