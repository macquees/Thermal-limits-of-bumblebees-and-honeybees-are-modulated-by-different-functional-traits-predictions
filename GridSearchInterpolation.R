## this file furns the ABC fitting for E and i_0

setwd("~/SUSPOLL/Code")
library(akima) #for interpolation
library(MASS)  #for fancy histogram plot
library(fields) #used of image.plot to get legends

Bumblebee=T   #only one can be true
Honeybee=F 

Euclid_Cutoff = 10

if(Bumblebee==TRUE){
  # read in the model results
  Fly_P1_P2 = read.csv('Fly_P1_P2_Bumblebee.csv',header=FALSE)
  Cool_P1_P2 = read.csv('Cool_P1_P2_Bumblebee.csv',header=FALSE)
  # target for ability to fly is 10\dC
  # target for cooling is 20-25\dC,
  fly_min = 9
  fly_max = 11
  cool_min = 20
  cool_max = 25
  
  ## #these are the parameter values used in the simulations
  i0_min = 0.001349728  #shouldn't go lower than resting i0...
  i0_max = 0.15  #Heinrich max number (rounded up?)
  i0_step = 0.0003
  i0_axis = seq(from=i0_min,by=i0_step,to=i0_max)
  E_min = 0
  E_max = 0.7
  E_step = 0.002
  E_axis = seq(from=E_min,by=E_step,to=E_max)
  
  shiver_target = 10
  diverge_target = 22.5
  
}

if(Honeybee==TRUE){
  # read in the model results
  Fly_P1_P2 = read.csv('Fly_P1_P2_Honeybee.csv',header=FALSE)
  Cool_P1_P2 = read.csv('Cool_P1_P2_Honeybee.csv',header=FALSE)
  # target for ability to fly is 10\dC
  # target for cooling is 35\dC,
  fly_min = 9
  fly_max = 11
  cool_min = 34
  cool_max = 36
  
  i0_default = 3.20*10^-2;     #I_flying
  i0_min = 4.52*10^(-4);  #shouldn't go lower than resting P1... 
  i0_max = 1.5*i0_default;  #some maximum
  i0_step = 0.0001;
  i0_axis = seq(from=i0_min,by=i0_step,to=i0_max)

  E_min = 0
  E_max = 0.7
  E_step = 0.002
  E_axis = seq(from=E_min,by=E_step,to=E_max)
  # E_min = 0.63
  # E_max = 0.63
  # E_step = 0.002
  # E_axis = 0.63
  
  shiver_target = 10
  diverge_target = 35
  
  
}

#plot the heatmaps
par(mfrow=c(2,1))
imagePlot(i0_axis,E_axis,as.matrix(Fly_P1_P2),main='air temp able to fly')
imagePlot(i0_axis,E_axis,as.matrix(Cool_P1_P2),main='air temp cooling starts')
# points(i0_default,E_default,col='red',pch=20)

imagePlot(i0_axis,E_axis,abs(as.matrix(Fly_P1_P2)-shiver_target), 
           main='able to fly (difference from target)')
imagePlot(i0_axis,E_axis,abs(as.matrix(Cool_P1_P2)-diverge_target),
           main='cooling starts (difference from target)')


i0_values_list = NULL
E_values_list = NULL
Dist_values_list = NULL

FlyTemp_list = NULL
CoolTemp_list = NULL


i0_values_interp_list = NULL
E_values_interp_list = NULL
Dist_values_interp_list = NULL
FlyTemp_interp_list = NULL
CoolTemp_interp_list = NULL

for(i in 1:1000){
  # sample the target temperatures
  # target for ability to fly is 10\dC
  # target for cooling is 20-25\dC, so use midpoint
  Fly_target = runif(n=1,min=fly_min,max=fly_max)
  Cool_target = runif(n=1,min=cool_min,max=cool_max)
  
  # calculate the euclidean distance for each model result
  Fly_difference = Fly_P1_P2-Fly_target
  Cool_difference = Cool_P1_P2-Cool_target
  Euclid_dist = sqrt(Fly_difference^2+Cool_difference^2)

  
  ######################################################
  # ###### First just use  euclidian distance cutoff ####
  # 
  # 
  # # # find the minimum distance
  # cutoff_points = which(Euclid_dist<=Euclid_Cutoff,arr.ind=TRUE)
  # i0_values_cutoff = i0_axis[cutoff_points[,1]]
  # E_values_cutoff = E_axis[cutoff_points[,2]]
  # Dist_values_cutoff = Euclid_dist[cutoff_points]
  # 
  # FlyTemp_cutoff = Fly_P1_P2[cutoff_points]
  # CoolTemp_cutoff = Cool_P1_P2[cutoff_points]
  # 
  # # save it and the param values
  # i0_values_list = c(i0_values_list,i0_values_min)
  # E_values_list = c(E_values_list,E_values_min)
  # Dist_values_list = c(Dist_values_list,Dist_values_min)
  # FlyTemp_list = c(FlyTemp_list,FlyTemp_cutoff)
  # CoolTemp_list = c(CoolTemp_list,CoolTemp_cutoff)
  
  
  ###### Euclidian distance minimum, single point ####
  
  # # find the minimum distance
  min_point_indices = which(Euclid_dist==min(Euclid_dist),arr.ind=TRUE)
  min_i0_values = i0_axis[min_point_indices[,1]]
  min_E_values = E_axis[min_point_indices[,2]]
  min_E_values = E_axis[min_point_indices[,2]]
  min_Dist_values = Euclid_dist[min_point_indices]
  min_Fly_Temp = Fly_P1_P2[min_point_indices]
  min_Cool_Temp = Cool_P1_P2[min_point_indices]
  
  # save it and the param values
  i0_values_list = c(i0_values_list,min_i0_values)
  E_values_list = c(E_values_list,min_E_values)
  Dist_values_list = c(Dist_values_list,min_Dist_values)
  FlyTemp_list = c(FlyTemp_list,min_Fly_Temp)
  CoolTemp_list = c(CoolTemp_list,min_Cool_Temp)
  
  
  # 
  
  # #########################################################
  # ###### Then use interpolation ###########################
  # ## select the points to interpolate with
  # Nr30 = sort(as.matrix(Euclid_dist))[300]
  # candidates = which(Euclid_dist<=Nr30,arr.ind=TRUE)
  # #rows for i0 and colums for f
  # i0_values = i0_axis[candidates[,1]]
  # E_values = E_axis[candidates[,2]]
  # Dist_values = Euclid_dist[candidates]
  # FlyTemp = Fly_P1_P2[candidates]
  # CoolTemp = Cool_P1_P2[candidates]
  # 
  # # Nr30 = sort(as.matrix(Fly_difference))[30]
  # # candidates = which(Fly_difference<=Nr30,arr.ind=TRUE)
  # # #rows for i0 and colums for f
  # # i0_values = i0_axis[candidates[,1]]
  # # E_values = E_axis[candidates[,2]]
  # # Dist_values = Fly_difference[candidates]
  # 
  # interpolation_results = interpp(x=i0_values, y=E_values, z=Dist_values,
  #                                 xo=seq(min(i0_values), max(i0_values), length = length(candidates)+50),
  #                                 yo=seq(min(E_values), max(E_values), length = length(candidates)+50),
  #                                        linear=FALSE, extrap=TRUE,jitter.random=FALSE)
  # #maybe use https://astrostatistics.psu.edu/su07/R/html/mgcv/html/te.html instead
  # 
  # Best_Dist = min(interpolation_results$z)
  # Best_indexs = which(interpolation_results$z==Best_Dist)
  # Best_i0s = interpolation_results$x[Best_indexs]
  # Best_Es = interpolation_results$y[Best_indexs]
  # 
  # Best_i0= mean(Best_i0s)
  # Best_E = mean(Best_Es)
  # #actually save the mean or median only, not all of them
  # 
  # i0_values_interp_list = c(i0_values_interp_list,Best_i0)
  # E_values_interp_list = c(E_values_interp_list,Best_E)
  # Dist_values_interp_list = c(Dist_values_interp_list,Best_Dist)
  # 
  
}

# par(mfrow=c(2,2))
par(mfrow=c(1,1))
plot(i0_values_list,E_values_list, main = 'candidates at minimum distance')
par(mfrow=c(3,1))
hist(i0_values_list,main = 'i_0 values')
hist(E_values_list,main = 'E values')
hist(Dist_values_list, main ='Euclidean distance values')
par(mfrow=c(2,1))
hist(FlyTemp_list, main ='Flying temperature values')
hist(CoolTemp_list, main ='Cooling temperature values')

# plot(i0_values_interp_list,E_values_interp_list, main = 'candidates based on interpolation')
# hist(i0_values_interp_list,main = 'i_0 values')
# hist(E_values_interp_list,main = 'E values')
# hist(Dist_values_interp_list, main ='Euclidean distance values')





# from  https://www.r-bloggers.com/2014/09/5-ways-to-do-2d-histograms-in-r/
##### Addendum: 2D Histogram + 1D on sides (from Computational ActSci w R) #######
#http://books.google.ca/books?id=YWcLBAAAQBAJ&pg=PA60&lpg=PA60&dq=kde2d+log&source=bl&ots=7AB-RAoMqY&sig=gFaHSoQCoGMXrR9BTaLOdCs198U&hl=en&sa=X&ei=8mQDVPqtMsi4ggSRnILQDw&redir_esc=y#v=onepage&q=kde2d%20log&f=false
n_breaks_i0 = 100
n_breaks_E = 100
n_breaks_histogram = 15
h1_minvalue <- hist(i0_values_list, n=n_breaks_histogram, plot=F)
h2_minvalue <- hist(E_values_list, n = n_breaks_histogram, plot=F)
top_minvalue <- max(h1_minvalue$counts, h2_minvalue$counts)
k_minvalue <- kde2d(i0_values_list, E_values_list, n=c(n_breaks_i0,n_breaks_E))  #choose n to match histograms
# margins
oldpar <- par()
par(mar=c(5,5,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k_minvalue,xlab="i_0 value",ylab="E value") #plot the image
# imagePlot(k_minvalue,xlab="i_0 value",ylab="E value",col = hcl.colors(12, "YlOrRd", rev = TRUE))
par(mar=c(0,4,1,0.5))
barplot(h1_minvalue$counts, axes=F, ylim=c(0, top_minvalue), space=0, col='red')
par(mar=c(4,0,0.5,1))
barplot(h2_minvalue$counts, axes=F, xlim=c(0, top_minvalue), space=0, col='red', horiz=T)

#find the maximum density and the corresponding i0,E values
max_density_minvalue = which(k_minvalue$z==max(k_minvalue$z),arr.ind=TRUE)
i0_guess_minvalue = k_minvalue$x[max_density_minvalue[1]]
E_guess_minvalue = k_minvalue$y[max_density_minvalue[2]]
i0_guess_minvalue
E_guess_minvalue


if(Bumblebee==T){
  write.csv(i0_values_list,file="Fit_i0_values_list_BB.csv")
  write.csv(E_values_list,file="Fit_E_values_list_BB.csv")
  write.csv(Dist_values_list,file="Fit_Dist_values_list_BB.csv")
  write.csv(FlyTemp_list,file="Fit_FlyTemp_list_BB.csv")
  write.csv(CoolTemp_list,file="Fit_CoolTemp_list_BB.csv")
}


if(Honeybee==T){
  write.csv(i0_values_list,file="Fit_i0_values_list_HB.csv")
  write.csv(E_values_list,file="Fit_E_values_list_HB.csv")
  write.csv(Dist_values_list,file="Fit_Dist_values_list_HB.csv")
  write.csv(FlyTemp_list,file="Fit_FlyTemp_values_list_HB.csv")
  write.csv(CoolTemp_list,file="Fit_CoolTemp_values_list_HB.csv")
}




# # from  https://www.r-bloggers.com/2014/09/5-ways-to-do-2d-histograms-in-r/
# ##### Addendum: 2D Histogram + 1D on sides (from Computational ActSci w R) #######
# #http://books.google.ca/books?id=YWcLBAAAQBAJ&pg=PA60&lpg=PA60&dq=kde2d+log&source=bl&ots=7AB-RAoMqY&sig=gFaHSoQCoGMXrR9BTaLOdCs198U&hl=en&sa=X&ei=8mQDVPqtMsi4ggSRnILQDw&redir_esc=y#v=onepage&q=kde2d%20log&f=false
# h1_interp <- hist(i0_values_interp_list, breaks=25, plot=F)
# h2_interp <- hist(E_values_interp_list, breaks=25, plot=F)
# top_interp <- max(h1_interp$counts, h2_interp$counts)
# k_interp <- kde2d(i0_values_interp_list, E_values_interp_list, n=25)
# # margins
# oldpar <- par()
# par(mar=c(3,3,1,1))
# layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
# image(k_interp, main='from interpolation') #plot the image
# par(mar=c(0,2,1,0))
# barplot(h1_interp$counts, axes=F, ylim=c(0, top_interp), space=0, col='red')
# par(mar=c(2,0,0.5,1))
# barplot(h2_interp$counts, axes=F, xlim=c(0, top_interp), space=0, col='red', horiz=T)
# 
# #find the maximum density and the corresponding i0,f values
# max_density_interp = which(k_interp$z==max(k_interp$z),arr.ind=TRUE)
# i0_guess_interp = k_interp$x[max_density_interp[1]]
# E_guess_interp = k_interp$y[max_density_interp[2]]
# 
# 
# 
# 
# 
# 
# 
# # # select the 10-30 points to interpolate with
# # candidates = which(Euclid_dist<=2,arr.ind=TRUE)
# # #rows for i0 and colums for f
# # i0_values = i0_axis[candidates[,1]]
# # E_values = E_axis[candidates[,2]]
# # Dist_values = Euclid_dist[candidates]
# # 
# # interpolation_results = interpp(x=i0_values, y=E_values, z=Dist_values,
# #         xo=seq(min(i0_values), max(i0_values), length = length(candidates)+20),
# #         yo=seq(min(E_values), max(E_values), length = length(candidates)+20))
# # 
# # Best = min(interpolation_results$z)
# # Best_index = which(interpolation_results$z==Best)
# # Best_i0 = interpolation_results$x[Best_index]
# # Best_E = interpolation_results$y[Best_index]
# # 
# # plot(interpolation_results$x,interpolation_results$y, main = 'candidates within a distance of 1 from the criteria')
# # plot(interpolation_results$z)
