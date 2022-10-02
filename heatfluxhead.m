function dydt = heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2)
%this function calculates dT_th/dt


if Honeybee == true 
    X_sfc = (rh*A*exp(B/(y-delta_T_h)))/Pr;
    Y_sfc = 1/( 1 + ((1-X_sfc)/X_sfc)*(MM_air/MM_vapor) );
else 
    Y_sfc = 0;
end

if CoolingSwitch_indicator == 1
    f = (1./(1+exp(-3.*(y - T_mK))));
else
    f = 1;
end

% S_h = 0;
% R1_h = 0;
% R2_h = 0;
% C1_h = 0;
% C2_h = 0;

dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))...   %thorax baseline
    + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h)... %head baseline
    - Ab*(y-T_aK)*f...   %abdomen 
    - (Ev1+Ev2*log( 1 - Y_sfc ))*f;   %evaporative



% (S + R1 - R2*T_th^4 - C1*T_th - C2 + I*exp(-E/(k*T_th)))...   %thorax baseline
%     + (S_h + R1_h - R2_h*(T_th-delta_T_h)^4 - C1_h*(T_th-delta_T_h) - C2_h)... %head baseline
%     - Ab*(T_th-T_aK)*f   %abdomen 