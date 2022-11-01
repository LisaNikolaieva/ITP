%Config;
%%

k_b=1.380649e-23;                                           % J/K
t_mult=1/omega_r;
r_mult=l_r;                                                 % m
r_mult_microm=r_mult*1e6;                                   % microm
V_mult=hbar*omega_r;                                        % J
mu_mult=hbar*omega_r;                                       % J
mu_mult_nK=mu_mult/k_b*1e9;                                 % nK
Psi_mult=1/l_r;                                             % m^-1
rho_mult_3D_z0=(Psi_mult/sqrt(sqrt(pi)*l_z))^2*mass;       % kg/m^3





Psi_3D_z0=Psi*Psi_mult/sqrt(sqrt(pi)*l_z);



f1=figure('visible', 'on',  'Position', [50 -50 2000 1000]);

figure(f1);
ax = gca;
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.5;
ax.TickLabelInterpreter='latex'; 
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

    
%% Density       
ax1=subplot(1,2,1);

    hold on    
    pc1=surf(rx*r_mult_microm, ry*r_mult_microm, 1e-6*abs(Psi_3D_z0).^2);
    pc1.EdgeColor='none';
    view(2);
    axis tight
    daspect ([1 1 1])
    box on
    cb1=colorbar;
    cb1.Location='eastoutside';
    cb1.Label.String='$\rho$, $1/$cm$^3$';
    cb1.Label.Interpreter='latex';
    cb1.Label.FontSize=16;
    cb1.Label.Position=[1 -3.1 0];
    cb1.Label.Rotation=0;
    cb1.TickLabelInterpreter='latex';
    cb1.FontSize=16;
    cm1=colormap(ax1,parula);
    xlabel('x, $\mu$m', 'interpreter','latex','FontSize', 45);
    ylabel('y, $\mu$m','FontSize', 45, 'interpreter','latex');
    ax = gca;
    ax.FontSize=16;
    ax.LabelFontSizeMultiplier = 1.5;
    ax.TickLabelInterpreter='latex'; 

%% Phase
ax2=subplot(1,2,2);

    su2=surf(rx*r_mult_microm, ry*r_mult_microm, angle(Psi_3D_z0));
    su2.EdgeColor='none';
    view(2);
    axis tight
    daspect ([1 1 1])
    box on
    cb2=colorbar;
    cb2.Location='eastoutside';
    cb2.Label.String='$\arg \Psi$';
    cb2.Label.Interpreter='latex';
    cb2.Label.FontSize=16;
    cb2.Label.Position=[1 -3.1 0];
    cb2.Label.Rotation=0;
    cb2.TickLabelInterpreter='latex';
    cb2.FontSize=16;
    cm2=colormap(ax2,hsv);
    xlabel('x, $\mu$m', 'interpreter','latex','FontSize', 45);
    ylabel('y, $\mu$m','FontSize', 45, 'interpreter','latex');
    ax = gca;
    ax.FontSize=16;
    ax.LabelFontSizeMultiplier = 1.5;
    ax.TickLabelInterpreter='latex'; 

    
    
