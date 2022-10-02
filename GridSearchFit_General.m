%%% This script carries evaluates the model at each point on the grid of 
%%% i_0 and E values for the ABC fitting in R

%%Goal of script: determine the air temperature at which the flying bee's thorax
%%temp gets above flying min and cooling critical point 

%%%%%%%%%%% Model Switches %%%%%%%%%%%%%%%%
Bumblebee = true;    %only one of Bumblebee and Honeybee can be true
Honeybee = false;
CoolingOff = true;   
%no cooling because we're trying to fit using the temperature at which the
%cooling should begin




%grid for P1 and P2

if Bumblebee==true
%P1_default = 0.06229515;     %I_flying, Kammer1974, converted to W
P1_min = 0.001349728;  %shouldn't go lower than resting P1... 
P1_max = 0.15;  %Heinrich number
P1_step = 0.0003;

%P1_max = 0.04;   %smaller range
%P1_step = 0.0001;

P1_axis = [P1_min:P1_step:P1_max];
P1_grid_size = length(P1_axis);
end



%P2_default = 0.63;  %E 
P2_min = 0;
P2_step = 0.002;
P2_max = 0.7;
P2_axis = [P2_min:P2_step:P2_max];
P2_grid_size = length(P2_axis);

%Storage matrices
%rows for P1, columns for r
Fly = zeros(P1_grid_size,P2_grid_size);  
Cool = zeros(P1_grid_size,P2_grid_size);
%ThermalDanger = zeros(P1_grid_size,P2_grid_size);

%tic
%ticBytes(gcp);
parfor i = 1:P1_grid_size
    for j = 1:P2_grid_size
        P1_guess = P1_axis(i);
        P2_guess = P2_axis(j);
        KeyTemps = RunModelABC_v2(P1_guess,P2_guess,Bumblebee,Honeybee,CoolingOff);
        Fly(i,j) = KeyTemps(1); %how far from warming goal
        Cool(i,j) = KeyTemps(2);  %how far from divergence goal
%         ThermalDanger(i,j) = KeyTemps(3);  %temp where bee goes above 42C
    end
end
%tocBytes(gcp);
%toc

if Bumblebee==true
    writematrix(Fly,'Fly_P1_P2_Bumblebee.csv');
    writematrix(Cool,'Cool_P1_P2_Bumblebee.csv');
%     writematrix(ThermalDanger,'ThermalDanger_P1_P2_Bumblebee.csv');
end

if Honeybee==true
    writematrix(Fly,'Fly_P1_P2_Honeybee.csv');
    writematrix(Cool,'Cool_P1_P2_Honeybee.csv');
%     writematrix(ThermalDanger,'ThermalDanger_P1_P2_Honeybee.csv');
end



clear  %clear workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Now Do Honeybee %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% pause(1800)

%%%%%%%%%%% Model Switches %%%%%%%%%%%%%%%%
Bumblebee = false;    %only one of Bumblebee and Honeybee can be true
Honeybee = true;
CoolingOff = true;   
%no cooling because we're trying to fit using the temperature at which the
%cooling should begin
Flying = true;



%grid for P1 and P2

if Honeybee==true
P1_default = 3.20*10^-2;     %I_flying
P1_min = 4.52*10^-4;  %shouldn't go lower than resting P1... 
P1_max = 1.5*P1_default;  %
P1_step = 0.0001;

% P1_max = 0.02;   %smaller range
% P1_step = 0.00005;

P1_axis = [P1_min:P1_step:P1_max];
P1_grid_size = length(P1_axis);
end

P2_default = 0.63;  %E 
P2_min = 0;
P2_step = 0.002;
P2_max = 0.7;
P2_axis = [P2_min:P2_step:P2_max];
P2_grid_size = length(P2_axis);
% P2_default = 0.63;  %E 
% P2_axis = 0.63;
% P2_grid_size = 1;

%Storage matrices
%rows for P1, columns for r
Fly = zeros(P1_grid_size,P2_grid_size);  
Cool = zeros(P1_grid_size,P2_grid_size);
% ThermalDanger = zeros(P1_grid_size,P2_grid_size);

%tic
%ticBytes(gcp);
parfor i = 1:P1_grid_size
    for j = 1:P2_grid_size
        P1_guess = P1_axis(i);
        P2_guess = P2_axis(j);
        KeyTemps = RunModelABC_v2(P1_guess,P2_guess,Bumblebee,Honeybee,CoolingOff);
        Fly(i,j) = KeyTemps(1); %how far from warming goal
        Cool(i,j) = KeyTemps(2);  %how far from divergence goal
%         ThermalDanger(i,j) = KeyTemps(3);  %temp where bee goes above 42C
    end
end
%tocBytes(gcp);
%toc

if Bumblebee==true
    writematrix(Fly,'Fly_P1_P2_Bumblebee.csv');
    writematrix(Cool,'Cool_P1_P2_Bumblebee.csv');
%     writematrix(ThermalDanger,'ThermalDanger_P1_P2_Bumblebee.csv');
end

if Honeybee==true
    writematrix(Fly,'Fly_P1_P2_Honeybee.csv');
    writematrix(Cool,'Cool_P1_P2_Honeybee.csv');
%     writematrix(ThermalDanger,'ThermalDanger_P1_P2_Honeybee.csv');
end