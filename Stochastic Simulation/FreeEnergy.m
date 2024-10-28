function [F, F_dilute, F_droplet, F_SurfaceTension, Vol_Frac_dilute]=FreeEnergy(N_dilute,N_droplet)

    global Vt V_TF
    global mu gamma kBT

    V_droplet=N_droplet*V_TF;
    V_dilute=Vt-V_droplet; 


    Vol_Frac_dilute=N_dilute*V_TF./V_dilute; 

    f_droplet=-mu/V_TF;

    f_dilute=kBT/V_TF.*(Vol_Frac_dilute.*log(Vol_Frac_dilute)-Vol_Frac_dilute);

    A=4*pi*(3/4/pi*V_TF*N_droplet).^(2/3);   
   
    F_dilute=V_dilute.*f_dilute;
    F_droplet=V_droplet.*f_droplet;
    F_SurfaceTension=gamma*A;
    F=F_dilute+F_droplet+F_SurfaceTension;
end