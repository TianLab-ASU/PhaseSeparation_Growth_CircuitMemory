

%condenstate parameters
global Vt Vt0 V_TF NA
global mu gamma kBT D  
global Omega Drop_Factor


Vt=4*10^(-19);% E.coli volume ≈ 0.4 µm^3=4*10^(-19)^3 m^3
Vt0=Vt;
V_TF=10^(-25); % TF molecular volume m^3
NA=6.0*10^23;% Avogadro constant
Omega=NA*10^(-3); % converting ode variable coentration (micromolar or μM) to number N: N=C*V*Omega 
% Drop_Factor=0.00029;
Drop_Factor=0.0005;

gamma=10^(-6); kBT=1.38*10^(-23)*305;  D=0.01;
mu0=4*10^(-20); 
mu=mu0;

global lagphase Nmax1 growthrate

lagphase=0;Nmax1=1;

%Circuit parameters
global Lara Cmin Cmax Km
global km0 km dm kp dp 

Lara=2*0.001;

km=4; km0=km*0.04; dm=8; kp=10; dp=2; 
Cmin=0.6; Cmax=3; Km=3*0.001;

growthrate=1.3;