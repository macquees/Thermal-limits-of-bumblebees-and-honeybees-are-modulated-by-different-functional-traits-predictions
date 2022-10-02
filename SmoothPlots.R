## this file creates the plots for figures 2 and 3

setwd("~/SUSPOLL/Code")
#install.packages("ggplot2")
library(ggplot2)

#colorblind friendly palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# black, goldenrod, light blue, green, yellow, dark blue, orange, mauve
# To use for fills, add
# scale_fill_manual(values=cbbPalette)
# To use for line and point colors, add
# scale_colour_manual(values=cbbPalette)


### Read in data set
MaxAirTemp_HB = read.csv('MaxAirTemp_HB.csv',header=F)
MaxAirTemp_BB = read.csv('MaxAirTemp_BB.csv',header=F)

masses = seq(from=0.01, by=0.0005, to=0.16)   #full range of mass values

### Do the smoothing 
Bees = data.frame(mass = masses, slowBB = MaxAirTemp_BB[,1],
                 defaultBB = MaxAirTemp_BB[,2],
                 fastBB = MaxAirTemp_BB[,3],
                 slowHB = MaxAirTemp_HB[,1],
                 defaultHB = MaxAirTemp_HB[,2],
                 fastHB = MaxAirTemp_HB[,3])
span.lo=0.55
slowBB.lo = loess(slowBB~masses,Bees,span=span.lo)
midBB.lo = loess(defaultBB~masses,Bees,span=span.lo)
fastBB.lo = loess(fastBB~masses,Bees,span=span.lo)

slowBB_smooth = predict(slowBB.lo,data.frame(mass=seq(0.01,0.16,0.0005)))
midBB_smooth = predict(midBB.lo,data.frame(mass=seq(0.01,0.16,0.0005)))
fastBB_smooth = predict(fastBB.lo,data.frame(mass=seq(0.01,0.16,0.0005)))

slowHB.lo = loess(slowHB~masses,Bees,span=span.lo)
midHB.lo = loess(defaultHB~masses,Bees,span=span.lo)
fastHB.lo = loess(fastHB~masses,Bees,span=span.lo)

slowHB_smooth = predict(slowHB.lo,data.frame(mass=seq(0.01,0.16,0.0005)))
midHB_smooth = predict(midHB.lo,data.frame(mass=seq(0.01,0.16,0.0005)))
fastHB_smooth = predict(fastHB.lo,data.frame(mass=seq(0.01,0.16,0.0005)))


### Create the data frame for ggplot
Slow = c(slowBB_smooth,slowHB_smooth)
Mid = c(midBB_smooth,midHB_smooth)
Fast = c(fastBB_smooth,fastHB_smooth)
Bee = c(rep("Bumblebee",length(masses)),rep("Honeybee",length(masses)))

Bees_smooth = data.frame(Bee = factor(Bee), Slow_Temp = Slow, Mid_Temp = Mid, Fast_Temp = Fast, Mass = c(masses,masses))

### Make a custom color palette
myColors <- c(cbbPalette[3],cbbPalette[7])
names(myColors) <- c("Bumblebee","Honeybee")
custom_colors <- scale_colour_manual(name = "Species Names", values = myColors)


### Establish the plot 
p = ggplot(Bees_smooth,aes(x=Mass,y=Mid_Temp,group=Bee,color=Bee)) + 
  geom_ribbon(aes(x = Mass, ymin=Slow_Temp,ymax=Fast_Temp,fill=Bee,color=Bee),  #shaded region
              alpha=0.3,colour = NA) +
  geom_line(size=1) +  #line for the middle speed
  geom_segment(aes(x = 0.0349, y = 40.3, xend = 0.0465, yend = 40.3), #HB thorax mass range
               colour = cbbPalette[7],size=2) +  
  geom_segment(aes(x = 0.014, y = 40.7, xend = 0.132, yend = 40.7),  #BB thorax mass range
               colour = cbbPalette[3],size=2) + 
  labs(x='mass of thorax (g)', 
       y='Maximum air temperature for sustained flight (C)') + 
  theme_bw()+
  theme(legend.title=element_blank()) +
  scale_colour_manual(values=c("Bumblebee"=cbbPalette[3],"Honeybee"=cbbPalette[7]))+
  scale_fill_manual(values=c("Bumblebee"=cbbPalette[3],"Honeybee"=cbbPalette[7])) +
 annotate("text", x = c(0.145,0.145,0.145,0.145,0.145,0.145), y = c(20,30,33,22,30,40), #add text labels
           label = c("1 m/s","4.5 m/s","5.5 m/s","3.6 m/s","4.5 m/s","9 m/s"),size=6) + 
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 15))+   #fix text size
  theme(legend.text = element_text(size = 20))

### Plot it with text, adjust axis text size
p + annotate("text", x = c(0.145,0.145,0.145,0.145,0.145,0.145), y = c(20,30,33,22,30,40), 
             label = c("1 m/s","4.5 m/s","5.5 m/s","3.6 m/s","4.5 m/s","9 m/s"),size=6) + 
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 15))+
  theme(legend.text = element_text(size = 20))

### Save the plot as a png  
#ggsave(filename="~/SUSPOLL/Manuscript/MaxAirTempPlot_ggplot2.png")
png("~/SUSPOLL/Manuscript/MaxAirTempPlot_ggplot2.png",width = 680, height = 480)
print(p)
dev.off()


#######################################################################################
#### Lookup Table Plots ############
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Read in data for BB
Flying_BB = read.csv('Flying_Bumblebee_LookupTable.csv',header=F)
colnames(Flying_BB) = c("Air_Temp","Thorax_Temp")
Resting_BB = read.csv('Resting_Bumblebee_LookupTable.csv',header=F)
colnames(Resting_BB)= c("Air_Temp", "Cooling_Off", "Cooling_On")
LethalTemp_BB = 45  #lethal thorax/air temp in C
CoolingTemp_BB = 42  #thorax temp where cooling begins
FlyingTemp_BB = 30  #thorax temp where flight can begin


### Read in data for HB
Flying_HB = read.csv('Flying_Honeybee_LookupTable.csv',header=F)
colnames(Flying_HB) = c("Air_Temp","Thorax_Temp")
Resting_HB = read.csv('Resting_Honeybee_LookupTable.csv',header=F)
colnames(Resting_HB)= c("Air_Temp", "Cooling_Off", "Cooling_On")
LethalTemp_HB = 52  #lethal thorax/air temp in C
CoolingTemp_HB = 47.9  #thorax temp where cooling begins
FlyingTemp_HB = 35  #thorax temp where flight can begin



#plot for BB
plot(Flying_BB[1:13,],lty=1,type="l",xlim=c(0,50),ylim=c(25,52),col=cbbPalette[8],
     main = "Bumblebee",ylab="Equilibrium thorax temperature (C)",xlab = "Air temperature (C)",
     cex.axis=1.5,cex.lab=1.5,lwd=2)
points(Flying_BB[14:17,],lty=2,type="l",col=cbbPalette[8],lwd=2)
points(Resting_BB$Air_Temp,Resting_BB$Cooling_Off,col=cbbPalette[4],lty=1,type="l",lwd=2)
points(Resting_BB$Air_Temp,Resting_BB$Cooling_On,col=cbbPalette[4],lty=2,type="l",lwd=2)
# abline(h=FlyingTemp_BB,col=cbbPalette[3],lwd=2)
# abline(h=CoolingTemp_BB,col=cbbPalette[3])
# abline(h=LethalTemp_BB,col=cbbPalette[3],lwd=2)
polygon(x=c(-5,55,55,-5),y=c(FlyingTemp_BB,FlyingTemp_BB,LethalTemp_BB,LethalTemp_BB)
       ,col=rgb(86/255,180/255,233/255,alpha=0.25),border=NA)


### plot for HB
plot(Flying_HB[1:14,],lty=1,type="l",xlim=c(0,50),ylim=c(25,52),col=cbbPalette[8],
     main = "Honeybee",ylab="Equilibrium thorax temperature (C)",xlab = "Air temperature (C)",
     cex.axis=1.5,cex.lab=1.5,lwd=2)
points(Flying_HB[15:18,],lty=2,type="l",col=cbbPalette[8],lwd=2)
points(Resting_HB$Air_Temp,Resting_HB$Cooling_Off,col=cbbPalette[4],lty=1,type="l",lwd=2)
points(Resting_HB$Air_Temp,Resting_HB$Cooling_On,col=cbbPalette[4],lty=2,type="l",lwd=2)
# abline(h=FlyingTemp_HB,col=cbbPalette[7],lwd=2)
# abline(h=CoolingTemp_HB,col=cbbPalette[7])
# abline(h=LethalTemp_HB,col=cbbPalette[7],lwd=2)
polygon(x=c(-5,55,55,-5),y=c(FlyingTemp_HB,FlyingTemp_HB,LethalTemp_HB,LethalTemp_HB)
        ,col=rgb(213/255,94/255,0/255,alpha=0.3),border=NA)


### plot the legend
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("topleft",legend=c("flying","resting","with cooling \n behaviour","without cooling \n behaviour",
                          "bumblebee \n limits for flight","honeybee \n limits for flight"),
       lty=c(1,1,2,1,1,1),col=c(cbbPalette[8],cbbPalette[4],1,1,cbbPalette[3],cbbPalette[7]))









### plot the legend
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("topleft",legend=c("flying","resting","with cooling behaviour","without cooling behaviour",
                          "bumblebee limits for flight","honeybee limits for flight"),
       lty=c(1,1,2,1,1,1),col=c(cbbPalette[8],cbbPalette[4],1,1,cbbPalette[3],cbbPalette[7]))








