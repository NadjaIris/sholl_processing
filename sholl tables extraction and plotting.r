
# import libraries
rm(list=ls())
library("tidyverse")
library(dplyr)
library(stringr)
library(ggfortify)
library(ggpubr)
library(rstatix)
library(Rmisc)
library(tidyr)
library(Hmisc)
library(readxl)
library(svDialogs)
library(readr)
library(ggthemes)
library(strex)
#import the data and transform it into a readable table--------------------------- 
#adjust the working directory here: 
tracings <- dlg_dir("select the directory for your sholl tables:", default = "", gui = .GUI)
ANSWERtracings <- tracings$res

setwd(ANSWERtracings)

#give the directory for the information file containing genotype, age, etc.:
infos <- dlg_open("select your file that contains genotype, morphology, animal, cell number, etc. for each cell:", default = "", gui = .GUI)

check <- read_delim(infos$res, 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

ldf <- list() # creates a list
#adjust directory of input files here:
listcsv <- dir(ANSWERtracings,pattern = "^cell.*.csv$") # creates the list of all the csv files in the directory

infos <- dlg_input("type how many columns (number of circles) you want to save:", default = "", gui = .GUI)
radii <- as.numeric(infos$res)

infos <- dlg_input("type how big the stepping size for the radii were (e.g. 10 for 10microns):", default = "", gui = .GUI)
radiistep <- as.numeric(infos$res)

#intersections get imported as a list.
for (k in 1:length(listcsv)){
  ldf[[k]] <- read.csv(listcsv[k])
}
names(ldf) <- listcsv

#intersections and radius are put into one table:
inters <- data.frame(NULL)
for (i in 1:length(listcsv)) {
  radius <- ldf[[i]][1]
  intersections <- ldf[[i]][2]
  inters <- rbind(inters, c(radius, intersections))
}
rm(intersections)
rm(listcsv)
#add names of the cells to the table: 
namescol <- data.frame(NULL)
for (i in 1:length(ldf)) {
  namen <- data.frame(rep(names(ldf)[i],length(ldf[[i]][,2])))
  namescol <- rbind(namescol,namen)
}
rm(namen)
inters$cellname <- namescol
rm(namescol)
rm(ldf)

#merge cellname (via cellnr) with info on genotype and age: 
split1 <- str_split(inters$cellname, pattern="ll", simplify=T)
split2 <- str_split(split1, pattern="_", simplify=T)[-1,1]
inters$cellnr <- as.factor(split2)

inters <- merge(inters, check, by="cellnr")
inters$Radius <- as.numeric(inters$Radius)
rm(split1)
rm(split2)

#type the cell number(s) of cells you want to exclude: 
exclude <- dlg_input("type the cell number(s) of cells you want to exclude like this c(1, 2,...)")
exclude <- str_to_vec(exclude$res)
inters <- inters[!(inters$cellnr %in% exclude),]

#save the intersections table: 
types <- dlg_input("did you trace neurites, dendrites or axons?", default="", gui=.GUI) #specificy for the name + graph later.
yname <- types$res
filename <- paste(types$res, "intersections_table.csv", sep="_")
saving <- dlg_dir("choose the directory to save your intersections table", default = "", gui = .GUI)

write.csv(inters, paste(saving$res, filename,sep="/"))

#graphs----------------------------------------------------------------------------
# Calculate mean and standard errors for each genotype group
inters_summary <- aggregate(Inters. ~ Radius + genotype, 
                            data = inters, 
                            function(x) c(mean = mean(x), se = sd(x) / sqrt(length(!is.na(x)))))


names(inters_summary) <- c("Radius", "genotype", "Inters")

#extract subcateogries into "normal" ones. 
dfinters <- data.frame("genotype"=inters_summary$genotype,
                       "radius"=inters_summary$Radius,
                       "means"=inters_summary[,3][,1],
                       "se"=inters_summary[,3][,2])
mylabs <- c(axis.title.x = element_text("Radius"),
            axis.title.y = element_text(paste(yname, "intersections", sep=" ")))

# Plotting of traces. 
ggplot(dfinters, aes(x=radius, y=means, group=genotype))+
  geom_line()+
  geom_ribbon(aes(ymin=means-se, ymax=means+se, fill=genotype), alpha=0.5)+
  theme_classic()+
  labs(x="Radius number", y=paste(yname, "Intersections", sep=" "), 
       title="Intersections between genotypes")+
  scale_x_continuous(limits=c(0, 40))+
  scale_color_manual(name="", labels=c("Ctrl","KO"), values=c('black',"red")) +
  scale_fill_manual(name="",
                    labels=c("Ctrl","KO"), 
                    values=c('black',"red"), 
                    na.value=NA)+
  theme(legend.position=c(0.5, 0.8), 
        legend.text = element_text(size=25),
        axis.text=element_text(size=30),
        axis.title=element_text(size=30),
        title = element_text(size=30),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  stat_compare_means(method="wilcox.test",size=5, label="p.signif", hide.ns=T)

#save plot as pdf: 
plotname1 <- paste(yname, "intersections", sep="_")

pdf(paste(plotname1, ".pdf",sep=""))
ggplot(dfinters, aes(x=radius, y=means, group=genotype))+
  geom_line()+
  geom_ribbon(aes(ymin=means-se, ymax=means+se, fill=genotype), alpha=0.5)+
  theme_classic()+
  labs(x="Radius number", y=paste(yname, "Intersections", sep=" "), 
       title="Intersections between genotypes")+
  scale_x_continuous(limits=c(0, 40))+
  scale_color_manual(name="", labels=c("Ctrl","KO"), values=c('black',"red")) +
  scale_fill_manual(name="",
                    labels=c("Ctrl","KO"), 
                    values=c('black',"red"), 
                    na.value=NA)+
  theme(legend.position=c(0.5, 0.8), 
        legend.text = element_text(size=25),
        axis.text=element_text(size=30),
        axis.title=element_text(size=30),
        title = element_text(size=30),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  stat_compare_means(method="wilcox.test",size=5, label="p.signif", hide.ns=T)
dev.off()

#extract p values between genotypes: 
pvalues <- data.frame("radius"=c(1:radii))

for (i in 1:radii){
  ctrl <- inters$Inters.[inters$Radius==i & inters$genotype=="Ctrl"]
  KO <- inters$Inters.[inters$Radius==i & inters$genotype=="KO"]
  pvalues[i,2] <- wilcox.test(ctrl, KO)$p.value
}
#plot p values: 
ggplot(pvalues,aes(x=radius, y=statistic))+
  geom_line(size=1.2)+
  geom_abline(intercept = 0.05, slope=0, color="orange", size=1.1)+
  theme_classic()+
  theme(axis.title=element_text(size=30),
        title = element_text(size=30),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))+
  scale_y_continuous(breaks = c(0, 0.05,0.25, 0.5, 0.75,1),
                     labels=c(0, 0.05,0.25, 0.5, 0.75,1))

#save p value plot: 
plotsave <- dlg_input("give the p value plot a name (+.pdf):")
plotname2 <- paste(yname, "p values", sep="_")

pdf(paste(plotname2, ".pdf", sep=""), width=10, height=5)
ggplot(pvalues,aes(x=radius, y=statistic))+
  geom_line(size=1.2)+
  geom_abline(intercept = 0.05, slope=0, color="orange", size=1.1)+
  theme_classic()+
  theme(axis.title=element_text(size=30),
        title = element_text(size=30),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))+
  scale_y_continuous(breaks = c(0, 0.05,0.25, 0.5, 0.75,1),
                     labels=c(0, 0.05,0.25, 0.5, 0.75,1))

dev.off()