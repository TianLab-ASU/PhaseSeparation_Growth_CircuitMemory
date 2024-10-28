function dy = LLPS_ODE2_CellGrowth(t,y)
  
global Vt V_TF NA
% global mu gamma kBT D 
global Omega Drop_Factor
global Lara Cmin Cmax Km
global km0 km dm kp dp 
global lagphase Nmax1 growthrate  


M=y(1); %*(Vt-V_drop);
P=y(2); %*(Vt-V_drop); 
Nt=ceil(P*Vt*Omega);  
CellNum=y(3);  

N_droplet=0:Nt;
N_dilute=Nt-N_droplet;
[F, ~, ~, ~, ~]=FreeEnergy(N_dilute,N_droplet);

[~, index]=min(F);
N_dilute_SS=N_dilute(index);   
N_droplet_SS=N_droplet(index);     
V_droplet=N_droplet_SS*V_TF;
V_dilute=Vt-V_droplet;  
    

Sa=(Cmin+(Cmax-Cmin)*Lara^3/(Lara^3+Km^3)); 

if N_droplet_SS==0 % no condensate
%    concentration_droplet=0;  
%    concentration_dilute=N_dilute_SS/V_dilute/Omega;  
   concentration_dilute=P;
   vm=km*Sa*concentration_dilute^2/(Sa*concentration_dilute^2+1)+km0;  
else
%    concentration_droplet=N_droplet_SS/V_droplet/Omega; 
   concentration_droplet=1/V_TF/Omega*Drop_Factor; 
   vm=km*Sa*concentration_droplet^2/(Sa*concentration_droplet^2+1)+km0; 
end
% GrowthRATE=growthrate/(AarC/J2+1)*(1-CellNum/Nmax1)*(t>lagphase
GrowthRATE=growthrate*(1-CellNum/Nmax1)*(t>lagphase);

dy=[vm-dm*M-GrowthRATE*M;
    kp*M-dp*P-GrowthRATE*P;
    GrowthRATE*CellNum; % growth with circuit 
    ];
 