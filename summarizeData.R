# March 30, 2015
# Summarize processed SIPP data
# Calculate number of observations, total number switches, 
# number of EE switches, number of UE switches, and probability of switching.
# Plot results, save to ./Figures/ directory.
# Precondition: processData.R has been run.
library(dplyr)
library(zoo)
library(stats)
library(reshape2)
library(ggplot2)

# 1996 Panel --------------------------------------------------------------
processed96 <- readRDS("./Data/processed96.RData")

# Tabulate number of observations
numObs <- group_by(processed96, date) %>%       # group by date
        summarize(numObs = n()) %>%             # count number of observations
        melt(id = "date") %>%                   # reshape into long
        mutate(panel = 1996)                    # add panel id

# Tabulate weighted number of switches
numSwitches <-  group_by(processed96, date) %>%
        summarize(numSwitchesRaw = sum(wpfinwgt[switchedOcc], 
                                       na.rm = T),
                  numSwitches2d = sum(wpfinwgt[switched2d], 
                                      na.rm = T),
                  numSwitchesEE = sum(wpfinwgt[switchedOcc & EE], 
                                      na.rm = T),
                  numSwitchesUE = sum(wpfinwgt[switchedOcc & UE], 
                                      na.rm = T)) %>%
        melt(id = "date") %>%
        mutate(panel = 1996)

# Calculate probability of switching using raw code
prSwitching <-  group_by(processed96, date) %>%
        summarize(prswitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE),
                  occObs = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  occEEObs = sum(wpfinwgt[EE], na.rm = TRUE),
                  occUEObs = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 1996)

# Calculate probability of switching using SOC 2d code
prSwitching2d <-  group_by(processed96, date) %>%
        summarize(prSwitched2d = weighted.mean(switched2d[EE | UE], wpfinwgt[EE | UE],
                                               na.rm = TRUE),
                  prSwitched2dEE = weighted.mean(switched2d[EE], wpfinwgt[EE],
                                                 na.rm = TRUE),
                  prSwitched2dUE = weighted.mean(switched2d[UE], wpfinwgt[UE],
                                                 na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 1996)

# Remove 1996 data from environment
rm(processed96)

# 2001 Panel --------------------------------------------------------------
processed01 <- readRDS("./Data/processed01.RData")

# Tabulate number of observations, add to numObs
numObs <- group_by(processed01, date) %>%       # group by date
        summarize(numObs = n()) %>%             # count number of observations
        melt(id = "date") %>%                   # rehsape into long
        mutate(panel = 2001) %>%                # add panel id
        bind_rows(numObs)                       # append to numObs

# Tabulate weighted number of switches, add to numSwitches
numSwitches <-  group_by(processed01, date) %>%
        summarize(numSwitchesRaw = sum(wpfinwgt[switchedOcc], 
                                       na.rm = T),
                  numSwitches2d = sum(wpfinwgt[switched2d], 
                                      na.rm = T),
                  numSwitchesEE = sum(wpfinwgt[switchedOcc & EE], 
                                      na.rm = T),
                  numSwitchesUE = sum(wpfinwgt[switchedOcc & UE], 
                                      na.rm = T)) %>%
        melt(id = "date") %>%
        mutate(panel = 2001) %>%
        bind_rows(numSwitches)

# Calculate probability of switching, add to prSwitching
prSwitching <-  group_by(processed01, date) %>%
        summarize(prswitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 2001) %>%
        bind_rows(prSwitching)

# Calculate probability of switching using SOC 2d code
prSwitching2d <-  group_by(processed01, date) %>%
        summarize(prSwitched2d = weighted.mean(switched2d[EE | UE], wpfinwgt[EE | UE],
                                               na.rm = TRUE),
                  prSwitched2dEE = weighted.mean(switched2d[EE], wpfinwgt[EE],
                                                 na.rm = TRUE),
                  prSwitched2dUE = weighted.mean(switched2d[UE], wpfinwgt[UE],
                                                 na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 2001) %>%
        bind_rows(prSwitching2d)

# Remove 2001 data from environment
rm(processed01)

# 2004 Panel --------------------------------------------------------------
processed04 <- readRDS("./Data/processed04.RData")

# Tabulate number of observations, add to numObs
numObs <- group_by(processed04, date) %>%       # group by date
        summarize(numObs = n()) %>%             # count number of observations
        melt(id = "date") %>%                   # rehsape into long
        mutate(panel = 2004) %>%                # add panel id
        bind_rows(numObs)                       # append to numObs

# Tabulate weighted number of switches, add to numSwitches
numSwitches <-  group_by(processed04, date) %>%
        summarize(numSwitchesRaw = sum(wpfinwgt[switchedOcc], 
                                       na.rm = T),
                  numSwitches2d = sum(wpfinwgt[switched2d], 
                                      na.rm = T),
                  numSwitchesEE = sum(wpfinwgt[switchedOcc & EE], 
                                      na.rm = T),
                  numSwitchesUE = sum(wpfinwgt[switchedOcc & UE], 
                                      na.rm = T)) %>%
        melt(id = "date") %>%
        mutate(panel = 2004) %>%
        bind_rows(numSwitches)

#Calculate probability of switching, add to prSwitching
prSwitching <-  group_by(processed04, date) %>%
        summarize(prswitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 2004) %>%
        bind_rows(prSwitching)

# Calculate probability of switching using SOC 2d code
prSwitching2d <-  group_by(processed04, date) %>%
        summarize(prSwitched2d = weighted.mean(switched2d[EE | UE], wpfinwgt[EE | UE],
                                               na.rm = TRUE),
                  prSwitched2dEE = weighted.mean(switched2d[EE], wpfinwgt[EE],
                                                 na.rm = TRUE),
                  prSwitched2dUE = weighted.mean(switched2d[UE], wpfinwgt[UE],
                                                 na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 2004) %>%
        bind_rows(prSwitching2d)

# Remove 2004 data from environment
rm(processed04)

# 2008 Panel --------------------------------------------------------------
processed08 <- readRDS("./Data/processed08.RData")

# Tabulate number of observations, add to numObs
numObs <- group_by(processed08, date) %>%       # group by date
        summarize(numObs = n()) %>%             # count number of observations
        melt(id = "date") %>%                   # rehsape into long
        mutate(panel = 2008) %>%                # add panel id
        bind_rows(numObs)                       # append to numObs

# Tabulate weighted number of switches, add to numSwitches
numSwitches <-  group_by(processed08, date) %>%
        summarize(numSwitchesRaw = sum(wpfinwgt[switchedOcc], 
                                       na.rm = T),
                  numSwitches2d = sum(wpfinwgt[switched2d], 
                                      na.rm = T),
                  numSwitchesEE = sum(wpfinwgt[switchedOcc & EE],
                                      na.rm = T),
                  numSwitchesUE = sum(wpfinwgt[switchedOcc & UE], 
                                      na.rm = T)) %>%
        melt(id = "date") %>%
        mutate(panel = 2008) %>%
        bind_rows(numSwitches)


#Calculate probability of switching, add to prSwitching
prSwitching <-  group_by(processed08, date) %>%
        summarize(prswitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 2008) %>%
        bind_rows(prSwitching)

# Calculate probability of switching using SOC 2d code
prSwitching2d <-  group_by(processed08, date) %>%
        summarize(prSwitched2d = weighted.mean(switched2d[EE | UE], wpfinwgt[EE | UE],
                                               na.rm = TRUE),
                  prSwitched2dEE = weighted.mean(switched2d[EE], wpfinwgt[EE],
                                                 na.rm = TRUE),
                  prSwitched2dUE = weighted.mean(switched2d[UE], wpfinwgt[UE],
                                                 na.rm = TRUE)) %>%
        melt(id = "date") %>%
        mutate(panel = 2008) %>%
        bind_rows(prSwitching2d)

# Remove 2008 data from environment
rm(processed08)

# Plots -------------------------------------------------------------------

# Number of observations
png("./Figures/numObs.png", width = 782, height = 569)
ggplot(numObs, aes(date, value, color = factor(panel))) + 
        geom_line() + geom_point() +
        ggtitle("Number of Observations in the SIPP") + 
        xlab("Date") + ylab("Individual Observations") +
        labs(color = "Panel")
dev.off()

# Adjust factor names for better plotting
levels(numSwitches$variable) <- c("Raw", "Two-digit SOC", 
                                  "Employment-to-Employment",
                                  "Unemployment-to-Employment")

# Number of raw and converted SOC switches
numSwitchesTotal <- subset(numSwitches,
                           variable == "Raw" | variable == "Two-digit SOC")

png("./Figures/numSwitchesTotal.png", width = 782, height = 569)
ggplot(numSwitchesTotal, aes(date, value, color = factor(panel))) +
        geom_line() + geom_point() + facet_grid(.~variable) +
        ggtitle("Number of Switches by Occupation Code Type") +
        xlab("Date") + ylab("Number of Switches (weighted)") +
        labs(color = "Panel")
dev.off()

# Number of EE and UE switches
numSwitchesByFlow <- subset(numSwitches, 
                            variable == "Employment-to-Employment" | 
                                    variable == "Unemployment-to-Employment")

png("./Figures/numSwitchesByFlow.png", width = 782, height = 569)
ggplot(numSwitchesByFlow, aes(date, value, color = factor(panel))) + 
        geom_line() + geom_point() + facet_grid(.~variable) +
        ggtitle("Number of Switches by Labor Force Flow") +
        xlab("Date") + ylab("Number of Switches (weighted)") +
        labs(color = "Panel")
dev.off()

# Number of UE switches only
numSwitchesUE <- subset(numSwitches, variable == "Unemployment-to-Employment")
png("./Figures/numSwitchesUE.png", width = 782, height = 569)
ggplot(numSwitchesUE, aes(date, value, color = factor(panel))) +
        geom_line() + geom_point() +
        ggtitle("Number of Unemployment-to-Employment Switches") +
        xlab("Date") + ylab("Number of switches (weighted)") +
        labs(color = "Panel")
dev.off()

# Probability of switching
levels(prSwitching$variable) <- c("Pooled", "EE", "UE")
png("./Figures/prSwitching.png", width = 782, height = 569)
ggplot(prSwitching, aes(date, value, color = variable)) + 
        geom_point() + geom_line() +
        ggtitle("Probability of Switching Occupations") +
        xlab("Date") + ylab("Probability") +
        labs(color = "Job switch type")
dev.off()

# Probability of switching 2d
levels(prSwitching2d$variable) <- c("Pooled", "EE", "UE")
png("./Figures/prSwitching2d.png", width = 782, height = 569)
ggplot(prSwitching2d, aes(date, value, color = variable)) + 
        geom_point() + geom_line() +
        ggtitle("Probability of Switching SOC2D Occupations") +
        xlab("Date") + ylab("Probability") +
        labs(color = "Job switch type")
dev.off()
        