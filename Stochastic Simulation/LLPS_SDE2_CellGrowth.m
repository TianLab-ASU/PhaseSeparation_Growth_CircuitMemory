function [S, P, K, N_dilute_SS, N_droplet_SS] = LLPS_SDE2_CellGrowth(t,y)
  
global Vt V_TF  
global Omega Drop_Factor
global Lara Cmin Cmax Km
global km0 km dm kp dp 
% global mu gamma kBT D
% global DoublingTime
 
  
M=y(1); %*(Vt-V_drop);
Nt=y(2); %*(Vt-V_drop); 

% delta_F=FreeEnergy(N_dilute+1,N_droplet-1)-FreeEnergy(N_dilute,N_droplet);
% kin=6*D*N_dilute/(V_dilute^(2/3));  %6*D*(Vt-V_drop)^(1/3)/v*phi_plus;
% kout=kin*exp(-delta_F/kBT); 

N_droplet=0:Nt;
N_dilute=Nt-N_droplet;
[F, F_dilute, F_droplet, F_SurfaceTension, Vol_Frac_dilute1]=FreeEnergy(N_dilute,N_droplet);

[Fmin, index]=min(F);
N_dilute_SS=N_dilute(index);   
N_droplet_SS=N_droplet(index);   % this could be used if the droplet promote transcription 
% Vol_Frac_dilute=Vol_Frac_dilute1(index);  
V_droplet=N_droplet_SS*V_TF;
V_dilute=Vt-V_droplet;   


Sa=(Cmin+(Cmax-Cmin)*Lara^3/(Lara^3+Km^3));  
if N_droplet_SS==0 % no condensate
%    concentration_droplet=0;  
   concentration_dilute=N_dilute_SS/V_dilute/Omega;   
%    concentration_dilute=P;
   vm=(km*Sa*concentration_dilute^2/(Sa*concentration_dilute^2+1)+km0)*Omega*Vt;  
else
%    concentration_droplet=N_droplet_SS/V_droplet/Omega; 
   concentration_droplet=1/V_TF/Omega*Drop_Factor;   
   vm=(km*Sa*concentration_droplet^2/(Sa*concentration_droplet^2+1)+km0)*Omega*Vt; 
end 
  

S = [  % How many of each chemical is used as substractes 
    0  0  % 0-->M
    1  0  % 0-->P
    1  0  % M-->0 
    0  1  % P-->0 
    ];
    
P = [  % How many does each chemical species is made as products
    1  0 % 0-->M
    1  1 % 0-->P
    0  0 % M-->0 
    0  0 % P-->0 
    ];
 
K  = [     
    vm
    kp
    dm
    dp 
    ];  

K=real(K);  