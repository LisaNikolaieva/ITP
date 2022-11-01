global DO_PARALLEL;

DO_PARALLEL = 0;
%%

NN0=2.5e4;%5e3;
s0=1;  
s1=-1;
s2=-1;

%% Grid parameters
L = 150;    Lx=L; Ly=L; 

N = 256;    Nx=N; Ny=N; 


%% Physical parameters
a_s=5.3112e-9;                  % m     Rb87 s-wave scattering length
omega_r = 2*pi*110;             % Hz
omega_z = 2*pi*245;             % Hz
kappa=omega_z/omega_r;
mass=87/6.02*1e-26;             % kg    Rb87 mass
hbar=1.06e-34;

l_r=sqrt(hbar/(mass*omega_r));  % m 
l_z=sqrt(hbar/(mass*omega_z));  % m

g3D_ph=4*pi*a_s*hbar^2/mass;    % m^5*kg/s^2
g2D_ph=g3D_ph/(sqrt(2*pi)*l_z); % m^4*kg/s^2

R_i_ph=20e-6;          % m
R_o_ph=55e-6;         % m

a_ph=10e-6;         % half-height diameter of barrier


%% dimension 2D
k_b=1.380649e-23;                                           % J/K
t_mult=1/omega_r;
w_mult=omega_r;
r_mult=l_r;                                                 % m
r_mult_microm=r_mult*1e6;                                   % microm
V_mult=hbar*omega_r;                                        % J
mu_mult=hbar*omega_r;                                       % J
mu_mult_nK=mu_mult/k_b*1e9;                                 % nK
Psi_mult=1/l_r;                        




%% Dimensionless parameters


g=sqrt(8*pi)*a_s/l_z;              % 2D dimensionless

R_i=R_i_ph/r_mult;                 % dimensionless
R_o=R_o_ph/r_mult;                 % dimensionless
a=a_ph/r_mult;                     % dimensionless


%% ITP parameters
niter = 20000; 
dt_itp = 0.01; 
MU_DIF = 1.e-9; 


%% Initialize grid
rx = linspace(-Lx/2,Lx/2,Nx);
ry = linspace(-Ly/2,Ly/2,Ny);

hx = rx(2)-rx(1);
hy = ry(2)-ry(1);

kx = [ (0:Nx/2-1)*2*pi/Lx -(Nx/2:-1:1)*2*pi/Lx];
ky = [ (0:Ny/2-1)*2*pi/Ly -(Ny/2:-1:1)*2*pi/Ly];

dV = hx * hy;


%% 
[X,Y] = meshgrid(rx,ry);
[KX,KY] = meshgrid(kx,ky);
if DO_PARALLEL
    kk = gpuArray((KX.^2+KY.^2));
else
    kk = ((KX.^2+KY.^2));
end




%% Potential

V=0.5*min((sqrt(X.^2+Y.^2)-R_i).^2,(sqrt(X.^2+Y.^2)-R_o).^2);

if DO_PARALLEL
    V = gpuArray(V);
end

%% Operators

eD_hp_itp = exp(0.5*dt_itp*(-0.5)*kk);




%% Clear
clear k KX KY;

h1=figure;