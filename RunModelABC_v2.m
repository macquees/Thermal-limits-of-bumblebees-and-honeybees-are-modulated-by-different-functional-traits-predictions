function KeyTemps = RunModelABC_v2(P1_guess, P2_guess,Bumblebee,Honeybee,CoolingOff)
% this function runs the model specifically for the fitting procedure 

Temps_median = zeros(51,4); %rows for temperature, colums for resting/shivering/flying
Temps_median(:,4) = (0:50).'; %fourth column is the temperature


%%%%%% Constant values %%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant
%h_fg = 2.33*10^6;  %latent heat of vaporization of water, J/kg 
A = 9.1496*10^10;  %clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3;  %clausius-clapyron constant B for water in K
rh =  0.6908;   %mean relative humidity from Arrian's data
%D_A = 2.06*10^-5;  %diffusion coefficient of air into itself in m^2/s - %calculated version is below
MM_air = 0.0289652;   %molar mass of dry air in g/mol
MM_vapor = 0.018016;  %molar mass of water vapor in kg/mol

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%

P = 332.3878; %mean solar irradiance from all of Arrian's data
T_gC = 17.1;                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
Pr = 1.013*10^5; %atmospheric pressure in N/m^2 (value used in Sidebotham)
R_specific = 287.058;  %J/kg/K for dry air
%kappa = 0.02534;   %15C dry air at sea level    %thermal conductivity of air (fill in an equation )
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%%%%%%%% Bee Parameters same for BB and HB %%%%%%%%%%%%%%%%
%s = 0.9;  %fraction of internal temp at surface - default value calculated from Church1960 data in C
%s = 0.9141875;  %fraction of internal temp at surface - BB fitted value
%s = 1; %temporarily make surface and internal temp equal 
s = 0.9965;  %ratio calculated from Church1960 data converted to K
c = 3.349;  %specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
C_l = 2.429809*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = 2.9;
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
maxy=52+273.15;    %don't solve above 52C because it's not biologically relevant
tspan = 0:2000; 
%E_opts = [0.63 0.63 0.63];    %Brown2004 activation energy
%E_opts = [0.63 0 0]; 
%r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
%r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
%r = 0.03;    %playing around with value to get desired behaviour
%R_0 = 0.000616/2;  %radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue
   

if Bumblebee==true
   %%%%%%%%%%%%%%% Bee Parameters %%%%%%%%%%%%%%%

    A_th = 9.3896*10^(-5)  ; %thorax surface area in m^2, from Church1960
    A_h = 3.61375*10^(-5)  ; %head surface area in m^2, from Cooper1985 - will need to update this to BB
    M_b = 0.149;   %mass of the bee in g, Joos1991, default
    %M_b = 0.035;   %mass of the bee in g, Joos1991, minimum
    %M_b = 0.351;   %mass of bee in g, Joos1991, maximum
    M_th = 0.057;  %mass of thorax in g, Joos1991
    %M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
    l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0 0.1 4.1];   %make the bee be out of wind when resting/shivering, default
    %v_options = [0 0 5.5];   %make the bee be out of wind when resting/shivering, maximum (Osborne2013)
    %v_options = [0 0 1];   %make the bee be out of wind when resting/shivering, minimum (Osborne2013)
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 42+273.15;     %median temp for  abdomen cooling
    %I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
    %I_flying = 0.06229515;     %Kammer1974, converted to W
    %I_flying = 0.0018375;     %fitted value
    %I_flying = 0.2097035;     %Heinrich1975, converted to W
    %I_flying = 0.03;   %experiment with the value
    %masses = [M_b M_b M_b];   %reference weight for Kammer only data is just M_b for now
    %masses = [M_b (0.25+0.60)/2 (0.25+0.60)/2];   %reference weight for Heinrich (flying)
    masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
    %RefTemps = [25+273.15, 19.55556+273.15, 19.55556+273.15];   %Reference temp is 25C for Kammer (resting), 35-44C for Heinrich (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
    %norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3

    y0 = 30+273.15;   %initial temperature of the bee's thorax in K 
    LethalTemp = 45;  %lethal thorax/air temp in C
    CoolingTemp = 42;   %thorax temp where cooling begins
    FlyingTemp = 30;      %thorax temp where flight can begin
end

if Honeybee==true
    A_th = 4.5*10^(-5)  ; %thorax surface area in m^2, from ???
    A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 
    M_b = 0.100;   %mass of the bee in g, Joos1991
    M_th = 0.0407;  %mass of thorax in g, Joos1991
    l_th = 0.004;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = 0.91;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0.1 0.1 5.6];   %make the bee be out of wind when resting/shivering, default
    % v_options = [0.1 0.1 10];   %make the bee be out of wind when resting/shivering, maximum (fill in ref)
    %v_options = [0.1 0.1 1];   %make the bee be out of wind when resting/shivering, minimum (fill in ref)
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 47.9+273.15;     %median temp for  evaporative cooling

    %I_resting = 5.65*(80/1000)*(1/1000); %Rothe1989, mW/g -> W, 80mg reference mass
    %I_flying = 0.4*(80/1000); %Nachtigall1989, W/g -> W, 80mg reference mass%    
    %I_resting = 0.001349728;     %%experiment with value
    %I_flying = 0.004;     %experiment with value
    masses = [0.08 0.08 0.08];   %reference weight for Rothe/Nachtigal (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp 
    %norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3

    y0 = 39+273.15;   %initial temperature of the bee's thorax in K (fill in ref)
    LethalTemp = 52;  %lethal thorax/air temp in C
    CoolingTemp = 47.9;   %thorax temp where cooling begins
    FlyingTemp = 35;      %thorax temp where flight can begin
end


for i = 1:51  %go through temps 0 to 50 
    j = 3;  %flying only


    %%%%%%%%%%% Environmental Parameters that depend on air temp %%%%%%%%%%%%%%%%
    T_aC = i-1;    %varying air temp
    T_aK = T_aC+273.15;     %air temp in K
    %nu = 2.791*10^(-7)*T_aK^(0.7355)/Pr;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
    mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
    rho = Pr/(R_specific*T_aK);  %in humid air air at sea level 
    nu = mu/rho; %kinematic viscosity for humid air %The Shock Absorber Handbook
    kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));

    %%%%%%%%% Bee parameters that depend on metabolic state j
    v = v_options(j);
    %D_A = v*(2*R_0);    %diffusion coefficient equivalent to advection
    M_ref = masses(j);
    T_ref = RefTemps(j);
    %i_0 = norm_constants(j);   %being fitted
    
    if Bumblebee==true  %these are the parameters being fit for BB
        i_0 = P1_guess;
        E = P2_guess;
    end

    if Honeybee==true  %these are the parameters being fit for HB
        i_0 = P1_guess;
        E = P2_guess;
    end
    
    
    %%%%%% Solar Radiation %%%%%%%
    % Does not depend on T_h
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S = (alpha_si*epsilon_a*A_th*P)/(M_th*c) + (alpha_so*epsilon_a*A_th*a*P)/(M_th*c);   %solar radiation thorax
    S_h = (alpha_si*epsilon_a*A_h*P)/(M_th*c) + (alpha_so*epsilon_a*A_h*a*P)/(M_th*c);   %solar radiation head

    %%%%%% Thermal Radiation %%%%%%%
    % Term 2 does depend on T_h
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    R1 = (alpha_th*epsilon_a*A_th*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee thorax (from sky + earth)
    R2 = (epsilon_e*A_th*sigma)/(M_th*c);   %thermal radiation out of bee thorax
    R1_h = (alpha_th*epsilon_a*A_h*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee head (from sky + earth)
    R2_h = (epsilon_e*A_h*sigma)/(M_th*c);   %thermal radiation out of bee head


    %%%%%%%% Convection %%%%%%%%%
    %1st term does depend on T_h
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h  = (C_l*kappa/l_th)*(v*l_th/nu)^n;
    C1 = (h*A_th*s)/(M_th*c);           %thorax, will be multiplied by T_th, *s is for thorax surface temperature
    C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
    C1_h = (h*A_h)/(M_th*c);           %head
    C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



    %%%%%%%%%%% Metabolic %%%%%%%%%
    %Does  depend on T_th 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i_0=I, depends on T_th, relative to ref temp/mass



%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatfluxhd.m, which may need to be
%renamed
%Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));   %using this will make ode23 stop solving at T_th=maxy
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy),'RelTol',1e-12,'AbsTol',1e-14);

%solve and plot solution curves separately for passive or active model 

if CoolingOff == true %physiological model with cooling off    
    saving_index = j;   %index for storage vector is 1, 2, or 3
    Ab = 0;
    Ev1 = 0;
    Ev2 = 0;   %no abdomen or evaporative cooling
    CoolingSwitch_indicator = 0;  %no switching
    [t,y] = ode23(@(t,y) heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2), tspan, y0, Opt); %for use with constant environmental conditions
    y = calc_length(tspan,y,maxy);  %check if bee maxed out temp and fill in if necessary
    if max(y) == maxy   %if temp reached maximum, use that
        Temps_median(i,saving_index) = maxy-273.15;
    else  %otherwise, take the mean/median of the end
        Temps_median(i,saving_index) = median(y(500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
    end
end

end

if Bumblebee==true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %identify when shivering bee goes above 30C
    fly30 = find(Temps_median(:,3) >= FlyingTemp);  %get the array indices where flying (active model) bee temp>=30
    if isempty(fly30) %if the bee never gets above 30, set temp to  max
        fly30Temp = maxy;
    else
        fly30Temp = Temps_median(fly30(1),4);   %otherwise, get the air temp for the first index in fly30
    end

    %identify when bee thorax reaches critical cooling point
    cool42 = find(Temps_median(:,3) >= CoolingTemp);  %get the array indices where flying (active model) bee temp>=42
    if isempty(cool42)
        cool42Temp = maxy;  %if it never reaches 42, set to max
    else
        cool42Temp = Temps_median(cool42(1),4);  %if they do, get the air temp for the first index in paDiverge
    end

    % %identify when active flying bee goes above 42C
    % fly42 = find(Temps_median(:,4) >= 42);  %get the array indices where shivering bee temp>=42
    % if isempty(fly42)
    %     fly42Temp = 50;  %if it never goes above 42C, set temp to 50 as max
    % else
    %     fly42Temp = Temps_median(fly42(1),1);   %otherwise, get the air temp for the first index in fly42
    % end
    KeyTemps = [fly30Temp cool42Temp];
end

if Honeybee==true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %identify when shivering bee goes above 35C
    fly35 = find(Temps_median(:,3) >= FlyingTemp);  %get the array indices where flying (active model) bee temp>=35
    if isempty(fly35) %if the bee never gets above 35, set temp to  max
        fly35Temp = maxy;
    else
        fly35Temp = Temps_median(fly35(1),4);   %otherwise, get the air temp for the first index in fly35
    end

    %identify when bee thorax reaches critical cooling point
    cool47p9 = find(Temps_median(:,3) >= CoolingTemp);  %get the array indices where flying (active model) bee temp>=47.9
    if isempty(cool47p9)
        cool47p9Temp = maxy;  %if it never reaches 47.9, set to max
    else
        cool47p9Temp = Temps_median(cool47p9(1),4);  %if they do, get the air temp for the first index in paDiverge
    end

    % %identify when active flying bee goes above 42C
    % fly42 = find(Temps_median(:,4) >= 42);  %get the array indices where shivering bee temp>=42
    % if isempty(fly42)
    %     fly42Temp = 50;  %if it never goes above 42C, set temp to 50 as max
    % else
    %     fly42Temp = Temps_median(fly42(1),1);   %otherwise, get the air temp for the first index in fly42
    % end
    KeyTemps = [fly35Temp cool47p9Temp];
end




