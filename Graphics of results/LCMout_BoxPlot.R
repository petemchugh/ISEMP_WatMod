# This wee bit of code is designed to make box-and-whisker plots of simulation outputs
# from Life Cycle Model simulations; in the imported dataset, columns are
# scenarios, rows are reps, and the value is a 10-yr geo mean at the end of simis

# Clear the memory to be safe
rm(list = ls(all=TRUE))

# Load necessary packages
library(ggplot2)
library(reshape)

# Interactively set the directory for writing output (tables, graphics, etc.)
this_is_the_place<-file.path(choose.dir())

# Read in the data (navigating, interactively)
data1<-read.table(file.choose(),header=TRUE,sep=",") #Data for model fitting

# Plot of scenarios
########################################################
# Make the box plot of selected scenarios (you have to tell it what those are...)
# Here is base vs. maturation of current projects plus climax state (NTP)
sub.bp1<-data.frame(data1$base,data1$NTP.c.p,data1$Rest.c.p) #subset the data
mdata <- data.frame(melt(sub.bp1)) #massage it to have the right shape

ggplot(mdata, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(alpha = .6,size = 1) + 
  scale_fill_brewer(palette = "Set1") + 
  #stat_summary(fun.y = "mean", geom = "point", shape= 25, size= 8, fill= "white") +
  ylab("Spawner Abundance\n") + 
  #theme(axis.text=element_text(size=14),axis.title=element_text(size=18))+
  theme_bw(base_size=33) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_x_discrete(limits=c("data1.base","data1.Rest.c.p","data1.NTP.c.p"),
                   labels=c("Status\nQuo", "Ongoing\nRestoration","Fully\nRestored")) +
  theme(axis.title.x=element_blank()) +
  labs(fill='sample')+theme(legend.position="none")
########################################################

# Plot of life histories
########################################################
# Make the box plot of selected scenarios (you have to tell it what those are...)
# Here is base vs. no repeat spawners plus no residents
sub.bp1<-data.frame(data1$base,data1$no.repeat,data1$no.resident) #subset the data
mdata <- data.frame(melt(sub.bp1)) #massage it to have the right shape

ggplot(mdata, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(alpha = .6,size = 1) + 
  scale_fill_brewer(palette = "Set1") + 
  #stat_summary(fun.y = "mean", geom = "point", shape= 25, size= 8, fill= "white") +
  ylab("Spawner Abundance\n") + 
  #theme(axis.text=element_text(size=14),axis.title=element_text(size=18))+
  theme_bw(base_size=33) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_x_discrete(limits=c("data1.base","data1.no.repeat","data1.no.resident"),
                   labels=c("Residents\n+ Repeats", "No\nRepeats","No\nResidents")) +
  theme(axis.title.x=element_blank()) +
  labs(fill='sample')+theme(legend.position="none")
########################################################