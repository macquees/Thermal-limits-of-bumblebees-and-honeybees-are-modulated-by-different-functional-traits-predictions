%%%% This file is used to make figure 3
%% setup 
%BB thorax mass [0.014,0.132]
%HB thorax mass [0.0349,0.0465]
varymasses = 0.01:0.0005:0.16;   %full range
% varymasses = 0.001:0.05:0.2;    %testing the code limited range
pcolors=['b' 'y' 'r'];
MaxAirTemp_BB = zeros(length(varymasses),3);  %row for each mass, column for each flight speed
MaxAirTemp_HB = zeros(length(varymasses),3);  %row for each mass, column for each flight speed
BB = [1,0];  %do BB first
HB = [0,1];  %do HB second

CoolingOn = true;   %only for cooling behaviour 
CoolingOff = false;  
CoolingSwitch = false;  

Flying = true;      %only for flying bee
Resting = false;      
Shivering = false; 

Fitting = false;     %only true if running fitting procedure
VaryMass = true;     %here we do vary the mass


%% loop through everything 
for Bee = 1:2
%%%% Set model switches
Bumblebee = BB(Bee);    %only one of Bumblebee and Honeybee can be true
Honeybee = HB(Bee);

%%%%%%%%%%% Storage vectors %%%%%%%%%%%%%%%%

for m = 1:length(varymasses) 
    m_thorax = varymasses(m);  %set the thorax mass (will be given correct name in Parameters_Script
    
    for v_bee = 1:3  %vary the flight speed
    %% Set parameters 
    Parameters_Script;  

    Temps_median = zeros(51,2); %rows for temperature, store equilibrium temp
    Temps_median(:,2) = (0:50).'; %2nd column is the temperature
   
    %% Calculate when bee goes into thermal danger zone even with cooling
    for i = 1:51  %go through temps 0 to 50 
        j = 3; %flying bee only 
   
        %%% Reset the Environmental Parameters that depend on air temp 
        T_aC = i-1;    %varying air temp
        T_aK = T_aC+273.15;     %air temp in K
        %nu = 2.791*10^(-7)*T_aK^(0.7355)/Pr;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
        mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
        rho = Pr/(R_specific*T_aK);  %in humid air air at sea level 
        nu = mu/rho; %kinematic viscosity for humid air %The Shock Absorber Handbook
        kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));


        %%% Solve the model 
        plotting = false;    %if true, will show plots of each solution curve
        logging=  false;    %log scale is for calculating rates
        tspan = 0:2000;     %specify how long to solve for
        maxy = 70+273.15;   %cutoff to reduce computation time

        SolveModel_Script;

        %%% Calculate mean/median values 
        y = calc_length(tspan,y,maxy);  %check if bee maxed out temp and fill in if necessary
        if max(y) == maxy   %if temp reached maximum, use that
            Temps_median(i,1) = maxy-273.15;
        else  %otherwise, take the mean/median of the end
            Temps_median(i,1) = median(y(500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
        end



    end
%     ["Bee" Bee "mass" m "speed" v] %print out to show where we are

    %% Calculate the max air temp for flight
    MaxAirTemp_indices = find(Temps_median(:,1) <= CoolingTemp);  %get the array indices where flying bee with cooling thorax temp >= lethal
    length_of_indices = length(MaxAirTemp_indices);
    
    if Honeybee==true 
        if isempty(MaxAirTemp_indices)
            MaxAirTemp_HB(m,v_bee) = 0;  %if it always reaches lethal point, set to min air temp
        elseif length_of_indices == 51  %all of the air temps are less than lethal thorax point
            MaxAirTemp_HB(m,v_bee) = maxy;  %if it never reaches lethal point, set to max air temp
        else
            MaxAirTemp_HB(m,v_bee) = Temps_median(MaxAirTemp_indices(length_of_indices),2);  %if it does, get the air temp for the last index 
        end
    end
    
    if Bumblebee==true
        if isempty(MaxAirTemp_indices)
            MaxAirTemp_BB(m,v_bee) = 0;  %if it always reaches lethal point, set to min air temp
        elseif length_of_indices == 51  %all of the air temps are less than lethal thorax point
            MaxAirTemp_BB(m,v_bee) = maxy;  %if it never reaches lethal point, set to max air temp
        else
            MaxAirTemp_BB(m,v_bee) = Temps_median(MaxAirTemp_indices(length_of_indices),2);  %if it does, get the air temp for the last index 
        end
    end
    
    end
end
end

%%save the data for plotting in R
writematrix(MaxAirTemp_BB,'MaxAirTemp_BB.csv');
writematrix(MaxAirTemp_HB,'MaxAirTemp_HB.csv');

%%% The following can be used to plot curves in matlab, but they won't be
%%% as smooth as in the R script
% %% plot BB curves, slightly smothed
% plot(varymasses,smooth(MaxAirTemp_BB(:,1),0.1,'loess'),'r')
% hold on
% plot(varymasses,smooth(MaxAirTemp_BB(:,2),0.1,'loess'),'r')
% plot(varymasses,smooth(MaxAirTemp_BB(:,3),0.1,'loess'),'r')
% 
% %plot HB curves, smothed
% plot(varymasses,smooth(MaxAirTemp_HB(:,1),0.1,'loess'),'b')
% plot(varymasses,smooth(MaxAirTemp_HB(:,2),0.1,'loess'),'b')
% plot(varymasses,smooth(MaxAirTemp_HB(:,3),0.1,'loess'),'b')
% 
% %Fix window and labels
% xlim([0.01 0.16])
% xlabel('mass of thorax (g)')
% ylabel('Maximum air temperature for sustained flying')
% 
% 
% %add ranges to top
% x = [0.052];  %midpoint of BB range
% y = [44];    %suitable y-value
% err = [0.0286];  %half width of BB range
% errorbar(x,y,err,'horizontal', 'color', 'r')
% plot(0.057,44,'color','b','Marker','.', 'MarkerSize', 12)  %default BB m_th
% 
% x = [0.0411];  %midpoint of HB range
% y = [43.5];    %suitable y-value
% err = [0.0075];  %half width of HB range
% errorbar(x,y,err,'horizontal', 'color', 'b')
% plot(0.0407,43.5,'color','r','Marker','.', 'MarkerSize', 12)  %default HB m_th
% 
