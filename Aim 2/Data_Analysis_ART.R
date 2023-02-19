rm(list = ls())


if(!require(ARTool)){install.packages("ARTool")}
if(!require(emmeans)){install.packages("emmeans")}
if(!require(multcomp)){install.packages("multcomp")}
if(!require(rcompanion)){install.packages("rcompanion ")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(psych)){install.packages("psych")}

#Set working directory and read excel into data
work_dir <- "~/GitHub/ARES_Simulink_SR/ARES_Simulink"
setwd(work_dir)

Data <- read.csv("DataRMS.csv")
#Data <- read.csv("DataTimeIn.csv")
#Data <- read.csv("DataSmooth.csv")
#Data <- read.csv("DataLPC.csv")
#Data <- read.csv("DataCAL.csv")
#Data <- read.csv("DataAud.csv")
#Data <- read.csv("DataTact.csv")
#Data <- read.csv("DataComp.csv")


#Data = read.table(Outcome_total_headers, header = TRUE)
Data$Condition = factor(Data$Condition, levels=unique(Data$Condition))
Data$Maps = factor(Data$Maps, levels=unique(Data$Maps))
Data$Subject = factor(Data$Subject, levels=unique(Data$Subject))
###  Check the data frame

library(psych)

headTail(Data)

str(Data)

summary(Data)
### Aligned ranks anova        
library(ARTool)

model = art(Performance ~ Maps + Condition*Subject + Error(Subjects),
            data = Data)


anova(model)


### Post-hoc comparisons 
model.lm = artlm(model, "Condition")

library(emmeans)

marginal = emmeans(model.lm, 
                   ~ Condition)

pairs(marginal, adjust = "tukey")
