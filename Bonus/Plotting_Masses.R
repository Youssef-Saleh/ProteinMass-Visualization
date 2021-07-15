# install.packages("renv",repos = "http://cran.us.r-project.org")
# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# renv::init(bare = TRUE)
# repos <- c(CRAN = "https://cloud.r-project.org")
# options(repos = repos)
# renv::install('isoband')
# renv::install("stringi", type = "win.binary")
# renv::install('ggplot2')
# renv::install('reshape2', rebuild = TRUE)
renv::restore()

library(ggplot2)
library(reshape2)

#setwd("E:/University 4/GP/bouns/") #put directory of your folder

f.data<-read.csv("Original_Mutated_Masses.csv")

dir.create("MassVisualization")
dir.create("MassChange")

# splitting the data into 50 sampled image
neW <-(split(f.data,(seq(nrow(f.data))-1)%%50))

for (i in seq_along(neW))
{
  dfm <- melt(neW[[i]][,c('Gene.Name','Mass','Original.Mass')],id.vars = 1)
  
  #Plotting Mass difference using line
  plot<-ggplot(neW[[i]],aes(y=MassDifference,x=Gene.Name, group = 1))+ 
    geom_line(color="blue",size=2) + geom_point(size=2)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(file=paste0("MassVisualization/MassVisualization",i,".png"))
  
  # Visualize the difference between mutated Mass and original mass using bar graph
  Plot2<-ggplot(data = dfm, aes(x = Gene.Name, y = value)) + 
    geom_bar(aes(fill = variable),stat = "identity",position = "dodge")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(file=paste0("MassChange/MassChange",i,".png"))
}


