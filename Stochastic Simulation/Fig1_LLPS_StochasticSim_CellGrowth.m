clc
clear
% close all 


muFold=[0; 1];
for i=1:2 % [no condensate, condensate] 

    Para 
 
    mu=muFold(i)*mu0;

    Lara0=([0:.25:5])*0.001;
    Lara=Lara0(end);
    [~,ytau00] = ode23s(@LLPS_ODE2_CellGrowth,[0 5000],[.1 .1 0.01]);
    
    Tmax=20; 
    Lara=0.001;

    [ttau0,ytau0, N_dilute0, N_droplet0,concentration_dilute0,concentration_droplet0,concentration_ave0,vt_rec] = ExactStoch_BatchGrowth(@(t,y) LLPS_SDE2_CellGrowth(t,y),[0 Tmax],ceil(ytau00(end,1:2)*Vt*Omega));
 
    Sa=(Cmin+(Cmax-Cmin)*Lara^3/(Lara^3+Km^3));     
    ProductionRate_dilute0=(km*Sa*concentration_dilute0.^2./(Sa*concentration_dilute0.^2+1)+km0);
    ProductionRate_droplet0=(km*Sa*concentration_droplet0.^2./(Sa*concentration_droplet0.^2+1)+km0); 
 
  
    subplot(3,3,i) 
    hold on 
    plot(ttau0, concentration_ave0(:,2))
    ylabel('Protein Concentration in Dilute Phase')
    ylim([0 3]) 
    box on 
    xlabel('Time')

    subplot(3,3,i+3) 
    hold on
    box on 
    plot(ttau0, N_droplet0) 
    xlabel('Time') 


    subplot(3,3,i+6) 
    hold on 
    plot(ttau0, concentration_droplet0)
    plot(ttau0, max(ProductionRate_droplet0,ProductionRate_dilute0))
    ylabel('ProductionRate')
    ylim([0 10])
 
    box on 
    xlabel('Time')
end
 





