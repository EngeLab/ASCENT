####Code for a data-adaptive way of selecting a gamma for multipcf: 

#Prerequisites 
#Run snakemake --configfile=config/Myconfig.yaml --until find_clones with multiple gammas 
#For example 0.1,0.25,0.5,0.75,1,2,3,4,5,7.5,10,15 and 20 
#Depending on binsize this will take some time - for most data 10-40kb bins work well for this stage of the pipeline. 40kb bins will run a lot faster through mpcf. 

#When the run has finished list the mpcf result files and get how many breakpoints (lines) are in each file named like: results/{patient_id}/clones/{patient_id}-mpcf-g{gamma}{binsize}.txt.gz
#Read this file into R 

source("workflow/scripts/clone_functions_forPaper.R")

df<-read_tsv("/wrk/resources/ASCENT/file_line_counts.tsv", col_names = c("condition", "breakpoints"))

df <- df %>%
  mutate(gamma = str_extract(condition, "(?<=g)[0-9\\.]+"),
         gamma = as.numeric(gamma))
ggplot(df, aes(x=gamma, y=breakpoints))+
  geom_point()+
  theme_dntr(axis_ticks = TRUE)

#Often the lowest gammas create > 10.000 breakpoints so often it is good to look only at <1000 segments 

ggplot(df%>%filter(breakpoints<1000), aes(x=gamma, y=breakpoints))+
  geom_point()+
  theme_dntr(axis_ticks = TRUE)

#Here it is most of the time not good to choose the lowest point - as this will often represent ~45 segments which is how many segments completely diploid cells have, since each chromosome arm is a segment
#Instead choose the point where the big jumps between gamma values have stopped, this will most likely be a suitable gamma value

