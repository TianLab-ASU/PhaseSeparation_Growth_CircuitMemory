clc
clear
% close all 

 
for cell=1:100 % [no condensate, condensate] 
    cell
    Para  

    Lara0=([0:.25:5])*0.001;
    Lara=Lara0(end);
    [~,ytau00] = ode23s(@LLPS_ODE2_CellGrowth,[0 5000],[.1 .1 0.01]);
    
    Tmax=20; 
    Lara=0.001;

    [ttau0,ytau0, N_dilute0, N_droplet0,concentration_dilute0,concentration_droplet0,concentration_ave0,vt_rec] = ExactStoch_BatchGrowth(@(t,y) LLPS_SDE2_CellGrowth(t,y),[0 Tmax],ceil(ytau00(end,1:2)*Vt*Omega));
 
    Sa=(Cmin+(Cmax-Cmin)*Lara^3/(Lara^3+Km^3));     
    ProductionRate_dilute0=(km*Sa*concentration_dilute0.^2./(Sa*concentration_dilute0.^2+1)+km0);
    ProductionRate_droplet0=(km*Sa*concentration_droplet0.^2./(Sa*concentration_droplet0.^2+1)+km0); 
  
    P_ave(cell,:)=interp1(ttau0,concentration_ave0(:,end),0:0.001:10);
    Drop_size(cell,:)=interp1(ttau0,N_droplet0,0:0.001:10); 
end
 

%%

subplot(2,3,3)


imagesc(0:0.001:10,1:100,P_ave)
xlabel('Time (h)')
ylabel ('cell index')


subplot(2,3,6)

imagesc(0:0.001:10,1:100,Drop_size)
xlabel('Time (h)')
ylabel ('cell index')

