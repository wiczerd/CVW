# April 10, 2015
# Calculate the unemployment rate in the SIPP, compare to headline
# unemployment rate from the BLS.
library(dplyr)
library(reshape2)
library(stats)
library(ggplot2)

# ESR 1 --With job entire month, worked all weeks.
# ESR 2 --With job entire month, missed 1 or more weeks, but not because of a layoff.
# ESR 3 --With job entire month, missed 1 or more weeks because of a layoff.
# ESR 4 --With job part of month, but not because of a layoff or looking for work.
# ESR 5 --With job part of month, some time spent on layoff or looking for work.
# ESR 6 --No job in month, spent entire month on layoff or looking for work.
# ESR 7 --No job in month, spent part of month on layoff or looking for work.
# ESR 8 --No job in month, no time spent on layoff or looking for work.

haver <- read.csv("./Data/unrate.csv", skip = 1) %>%
        mutate(date = as.Date(date2, "%m/%d/%Y"),
               year = format(date, "%Y"),
               month = format(date, "%m")) %>%
        select(-starts_with("date")) %>%
        rename()

sipp96 <- readRDS("./Data/sipp96.RData")
laborStats <- sipp96 %>%
        mutate(unemp = (esr == 3 | esr == 5 | esr == 6 | esr == 7)) %>%
        group_by(date) %>%
        summarize(unRate = weighted.mean(unemp[esr != 8], wpfinwgt[esr != 8], na.rm = TRUE))
rm(sipp96)

sipp01 <- readRDS("./Data/sipp01.RData")
laborStats <- sipp01 %>%
        mutate(unemp = (esr == 3 | esr == 5 | esr == 6 | esr == 7)) %>%
        group_by(date) %>%
        summarize(unRate = weighted.mean(unemp[esr != 8], wpfinwgt[esr != 8], na.rm = TRUE)) %>%
        bind_rows(laborStats)
rm(sipp01)

sipp04 <- readRDS("./Data/sipp04.RData")
laborStats <- sipp04 %>%
        mutate(unemp = (esr == 3 | esr == 5 | esr == 6 | esr == 7)) %>%
        group_by(date) %>%
        summarize(unRate = weighted.mean(unemp[esr != 8], wpfinwgt[esr != 8], na.rm = TRUE)) %>%
        bind_rows(laborStats)
rm(sipp04)

sipp08 <- readRDS("./Data/sipp08.RData")
laborStats <- sipp08 %>%
        mutate(unemp = (esr == 3 | esr == 5 | esr == 6 | esr == 7)) %>%
        group_by(date) %>%
        summarize(unRate = weighted.mean(unemp[esr != 8], wpfinwgt[esr != 8], na.rm = TRUE)) %>%
        bind_rows(laborStats)
rm(sipp08)

merged <- laborStats %>%
        mutate(year = format(date, "%Y"),
               month = format(date, "%m")) %>%
        select(-date) %>%
        left_join(haver) %>%
        mutate(unRate = round(unRate*100, 1),
               date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        select(-month, -year) %>%
        melt(id = "date")

levels(merged$variable) <- c("SIPP", "BLS (NSA)")
png("./Figures/unrates.png", width = 782, height = 569)
ggplot(merged, aes(date, value, color = variable)) + 
        geom_line() +
        geom_point() +
        ggtitle("Comparison of BLS and SIPP Unemployment Rates") +
        xlab("Date") + ylab("Unemployment rate") +
        labs(color = "Rate")
dev.off()
        