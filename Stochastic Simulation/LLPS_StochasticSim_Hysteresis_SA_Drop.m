clc
clear
close all


Para 
  
Cellnum=1000;
Lara0=([0:.25:5])*0.001;
Lara=Lara0(1);
[~,ytau00] = ode23s(@LLPS_ODE2_CellGrowth,[0 5000],[.01 .01 0.01]);
y00=ceil(ytau00(end,1:2)*Vt*Omega)
Tmax=20; % check this part
%%
for i=12:length(Lara0)

    Lara=Lara0(i);

    for cell=1:Cellnum  

        [i cell] 
        [ttau0,ytau0, N_dilute0, N_droplet0,concentration_dilute0,concentration_droplet0,concentration_ave0] = ExactStoch_BatchGrowth_EndPoint(@(t,y) LLPS_SDE2_CellGrowth(t,y),[0 Tmax],y00);

%         hold on 
%         plot(ttau0, concentration_ave0(:,2)) 
%         pause(0.01)

        y_off2on(i,cell)=ytau0(end,end);
        concentration_ave_off2on(i,cell)=concentration_ave0(end,end);
        concentration_dilute_off2on(i,cell)=concentration_dilute0(end);
        concentration_droplet_off2on(i,cell)=concentration_droplet0(end);

    end

save LLPS_StochasticSim_Hysteresis_SA_Drop
end
%% 
Lara0=([0:.25:5])*0.001;
Lara=Lara0(end);
[~,ytau00] = ode23s(@LLPS_ODE2_CellGrowth,[0 5000],[.1 .1 0.01]);
y00=ceil(ytau00(end,1:2)*Vt*Omega)
Tmax=20;

for i=1:length(Lara0)

    Lara=Lara0(i);

    for cell=1:Cellnum  

%         [i cell] 
        [ttau0,ytau0, N_dilute0, N_droplet0,concentration_dilute0,concentration_droplet0,concentration_ave0] = ExactStoch_BatchGrowth_EndPoint(@(t,y) LLPS_SDE2_CellGrowth(t,y),[0 Tmax],y00);

%         hold on 
%         plot(ttau0, concentration_ave0(:,2)) 
%         pause(0.01)

        y_on2off(i,cell)=ytau0(end,end);
        concentration_ave_on2off(i,cell)=concentration_ave0(end,end);
        concentration_dilute_on2off(i,cell)=concentration_dilute0(end);
        concentration_droplet_on2off(i,cell)=concentration_droplet0(end);

    end

save LLPS_StochasticSim_Hysteresis_SA_Drop
end

 
