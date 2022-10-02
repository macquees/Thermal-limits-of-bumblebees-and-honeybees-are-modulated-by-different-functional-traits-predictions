%%%%%% Used to produce figure 2, and can be modified to produce 
%%%%%% individual solution curves


%%%%%%%%%%% Model Switches %%%%%%%%%%%%%%%%
Bumblebee = true;    %only one of Bumblebee and Honeybee can be true
Honeybee = false;

Cooling_Off = [1 0 0];  %first cooling off 
Cooling_Switch = [0 1 0];  %second switching (but actually skip it)
Cooling_On = [0 0 1];  %third cooling on
RestingBehaviour = [1 0 0];  %first resting
ShiveringBehaviour = [0 1 0];  %second do shivering (but skip for now)
FlyingBehaviour = [0 0 1];  %third do flying

Fitting = false;     %only true when running fitting procedure
VaryMass = false;    %only true when running mass/flight speed procedure

%%%%%%%%%%% Storage vectors %%%%%%%%%%%%%%%%
Temps = zeros(51,9); %rows for temperature, colums for resting/shivering/flying
Temps_median = zeros(51,9); %rows for temperature, colums for resting/shivering/flying
Temps_mean = zeros(51,9); %rows for temperature, colums for resting/shivering/flying

pcolors=['b' 'g' 'r'];

for cool_behaviour = [1 3]  %off/switch/on
for metabolic_state = [1 3]  %resting/shivering/flying

    %%%%%%%%%%% Model Switches %%%%%%%%%%%%%%%%
    CoolingOff = Cooling_Off(cool_behaviour);   %only one of cooling can be true here
    CoolingSwitch = false;  %not using switching
    CoolingOn = Cooling_On(cool_behaviour);  

    Resting = RestingBehaviour(metabolic_state);      %only one of metabolic states can be true here
    Shivering = ShiveringBehaviour(metabolic_state); 
    Flying = FlyingBehaviour(metabolic_state);

    j = nonzeros([Resting Shivering Flying].*[1 2 3])';  %whichever one it is
    
    if CoolingOff == true %physiological model with abdomen off    
        saving_index = j;   %index for storage vector is 1, 2, or 3
    elseif CoolingSwitch == true %behavioural model with abdomen switching
        saving_index = j+3;   %index for storage vector
    elseif CoolingOn == true %behavioural model with abdomen cooling always on    
        saving_index = j+6;   %index for storage vector
    end


    %%%%%%%%%%%% load all the default parameters %%%%%%%%%%%%%%%%%%%%
    Parameters_Script;
    
    %% go through temps 0 to 50 
     for i = 1:51  
%      for i = 20   %use this to run a single temperature value

    %% Reset the Environmental Parameters that depend on air temp 
    T_aC = i-1;    %varying air temp
    T_aK = T_aC+273.15;     %air temp in K
    %nu = 2.791*10^(-7)*T_aK^(0.7355)/Pr;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
    mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
    rho = Pr/(R_specific*T_aK);  %in humid air air at sea level 
    nu = mu/rho; %kinematic viscosity for humid air %The Shock Absorber Handbook
    kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));


  
    %% Solve the model 
    plotting = true;    %if true, will show plots of each solution curve
    logging=  false;    %log scale is for calculating rates
    tspan = 0:2000;     %specify how long to solve for
    maxy = 70+273.15;   %cutoff to reduce computation time

    SolveModel_Script;
    
    %% Calculate mean/median values 
    y = calc_length(tspan,y,maxy);  %check if bee maxed out temp and fill in if necessary
    if max(y) == maxy   %if temp reached maximum, use that
        Temps_median(i,saving_index) = maxy-273.15;
        Temps_mean(i,saving_index) = maxy-273.15;
    else  %otherwise, take the mean/median of the end
        Temps_median(i,saving_index) = median(y(500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
        Temps_mean(i,saving_index) = mean(y(500:length(tspan)))-273.15;     %but save both to compare just in case
    end
    [i j saving_index] %print out to show where we are


    end
end
end


%% Make data for plotting the combined lookup table in R

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %identify when flying bee goes above flight possible temp
    fly = find(Temps_median(:,3) >= FlyingTemp);  %get the array indices where flying (no cooling model) bee temp>=30
    if isempty(fly) %if the bee never gets above 30, set temp to  max
        flyindex= maxy-273.15;
    else
        flyindex = fly(1)-1;   %otherwise, get the air temp for the first index in shiver30
    end

    %identify when bee thorax reaches critical cooling point
    coolstart = find(Temps_median(:,3) >= CoolingTemp);  %get the array indices where flying (no cooling model) bee temp>=42
    if isempty(coolstart)
        coolstartindex = maxy-273.15;  %if it never reaches 42, set to max
    else
        coolstartindex = coolstart(1)-1;  %if they do, get the air temp for the first index in paDiverge
    end
    
    %identify when bee thorax goes above critical cooling point for cooling
    %on model
    coolincrease = find(Temps_median(:,9) >= CoolingTemp);  %get the array indices where flying (cooling model) bee temp>=42
    if isempty(coolincrease)
        coolincreaseindex = maxy-273.15;  %if it never reaches 42, set to max
    else
        coolincreaseindex = coolincrease(1);  %if they do, get the air temp for the first index in paDiverge
    end
    
    
    %identify when bee thorax reaches thermal max
    lethal = find(Temps_median(:,9) >= LethalTemp);  %get the array indices where flying (cooling model) bee temp>=42
    if isempty(lethal)
        lethalindex = maxy-273.15;  %if it never reaches 45, set to max
    else
        lethalindex = lethal(1);  %if they do, get the air temp for the first index in paDiverge
    end



Flying_Temps = [transpose(Temps_median(flyindex:coolstartindex,3)),Temps_median(coolstartindex,3),... 
                Temps_median(coolstartindex,3), Temps_median(coolstartindex,3), Temps_median(51,9)];
Flying_Air = [air(flyindex:coolstartindex), air(coolstartindex), air(coolincreaseindex),air(coolincreaseindex), air(51)];
Flying_Data = transpose([Flying_Air;Flying_Temps]);

Resting_Data = [transpose(air),Temps_median(:,1),Temps_median(:,7)];

if Honeybee == true
writematrix(Flying_Data,'Flying_Honeybee_LookupTable.csv');
writematrix(Resting_Data,'Resting_Honeybee_LookupTable.csv');
end

if Bumblebee == true
writematrix(Flying_Data,'Flying_Bumblebee_LookupTable.csv');
writematrix(Resting_Data,'Resting_Bumblebee_LookupTable.csv');
end

%% use the following to generate plots in matlab

% %plot thorax temp as a function of air temp
% figure(4) %create a new plot window for the lookup table plot
% hold on
% air=0:50;
% thorax=air;
% %plines = ['-' '-' '-' '--' '--' '--'];
%  plot(air,Temps_median(:,1),'b-')  %cooling off, resting
%  plot(air,Temps_median(:,2),'g')  %cooling off, shivering
%  plot(air,Temps_median(:,3),'r')  %cooling off, flying
%  plot(air,Temps_median(:,4),'b-.')  %physiological (cooling switching), resting (-. for dash-dotted line)
%  plot(air,Temps_median(:,5),'g-.')  %physiological (cooling switching), shivering
%  plot(air,Temps_median(:,6),'r-.')  %physiological (cooling switching), flying
%  plot(air,Temps_median(:,7),'b--')  %cooling on, resting (-. for dash-dotted line)
%  plot(air,Temps_median(:,8),'g--')  %cooling on, shivering
%  plot(air,Temps_median(:,9),'r--')  %cooling on, flying
% ylim([20,LethalTemp])
% %plot(air,thorax,'k:')
% % title('Thorax temp - median')
% xlabel('Air Temperature (C)')
% ylabel('Equilibrium Thorax Temperature (C)')
% 
% yline(FlyingTemp);  %30 is the min thorax temp for flight, Heinrich1983
% yline(LethalTemp);   %lethal thorax temp 
% yline(CoolingTemp);   %lethal thorax temp 
% %yline(43,'k:');
% II=area([0,50],[CoolingTemp CoolingTemp],LethalTemp,'EdgeColor', 'none', 'FaceColor', 'r');
% alpha(0.1)
% fill(CoolingTemp,LethalTemp, [0.9 0.9 0.9])
% %hold off;
% %legend('resting','flying','environmental','behavioural','T_{th} = T_{air}','Location','southeast');
% %legend('resting, physiological','flying, physiological','resting, behavioural','flying, behavioural','T_{th} = T_{air}','Location','southeast');
% %legend('resting','shivering','flying','Location','southeast');
% %legend('cooling always off','cooling always on','Location','southeast');
% 
% legend('resting','flying','Location','northeast')
% 
% %plot thorax temp as a function of air temp
% figure(5) %create a new plot window for the lookup table plot
% hold on
% air=0:50;
% thorax=air;
% plot(air(flyindex:coolstartindex),Temps_median(flyindex:coolstartindex,3),'r')  %cooling off, flying
% %plot(air(coolincreaseindex:51),Temps_median(coolincreaseindex:51,9),'b--')  %cooling on, flying
% plot(air([coolstartindex coolincreaseindex]),[Temps_median(coolstartindex,3) Temps_median(coolstartindex,3)],'r--')  %cooling on, flying
% plot([air(coolincreaseindex) air(51)],[Temps_median(coolstartindex,3) Temps_median(51,9)],'r--')  %cooling on, flying
% ylim([20,53])
% xlabel('Air Temperature (C)')
% ylabel('Equilibrium Thorax Temperature (C)')
% 
% II=area([0,50],[CoolingTemp CoolingTemp],LethalTemp,'EdgeColor', 'none', 'FaceColor', 'r');
% alpha(0.1)
% fill(CoolingTemp,LethalTemp, [0.9 0.9 0.9])
% yline(FlyingTemp,'r');  %30 is the min thorax temp for flight, Heinrich1983
% yline(LethalTemp,'r');   %lethal thorax temp 
% yline(CoolingTemp,'r');   %lethal thorax temp 
% 
% plot(air,Temps_median(:,1),'k')  %cooling off, resting
% plot(air,Temps_median(:,7),'k--')  %cooling on, resting (-. for dash-dotted line)
% ylim([25,LethalTemp])
