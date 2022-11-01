Config;
niter=200;

%% Initializing
mu_TF=16;
rho_TF=(mu_TF-V)/g.*(mu_TF>=V);

Psi=sqrt(rho_TF).*exp(1i*s0*atan2(Y,X)).*exp(1i*s1*atan2(Y,(X-0.5*(R_i+R_o)))).*exp(1i*s2*atan2(Y,(X+0.5*(R_i+R_o))));
Psi = Psi.*sqrt(NN0/((sum(sum(abs(Psi).^2)))*dV));   

if DO_PARALLEL
    Psi = gpuArray(Psi);
end
%% ITP
tic

    PsiL=Psi;
    PsiR=Psi;
    khi=0.5;
    dt_itp=0.1;
    
for i=1:niter
    I(i)=i;
    
%     [Psi, MU(i), dE(i)]=SSFM(eD_hp_itp,dt_itp,V,g,NN0,kk,dV,Psi);
    
    [PsiL, MUL(i), dEL(i)]=SSFM(eD_hp_itp,dt_itp,V,g,NN0,kk,dV,Psi);
    [PsiS, MUS(i), dES(i)]=SSFM(eD_hp_itp,khi*dt_itp,V,g,NN0,kk,dV,Psi);
    
    if (MUL(i)<MUS(i))||(dEL(i)<dES(i))
        Psi=PsiL;
        khi=sqrt(khi);
        MU(i)=MUL(i);
    else
        Psi=PsiS;
        dt_itp=khi*dt_itp;
        eD_hp_itp = exp(0.5*dt_itp*(-0.5)*kk);
        khi=khi^2;
        MU(i)=MUS(i);
    end
    
    

     

    
    
    %% Figure 
    figure(h1);
    h1.Color='k';
    set(gca,'Color','k');
    set(gca,'xcolor','[0.55 0.55 0.55]') 
    set(gca,'ycolor','[0.55 0.55 0.55]') 
    set(gcf, 'InvertHardCopy', 'off');
    ax = gca;
    ax.FontSize=20;
    ax.LabelFontSizeMultiplier = 2.5;
    ax.TickLabelInterpreter='latex'; 
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.Box='on';
    xlabel('i','FontSize',40,'Interpreter','latex');
    ylabel('$\mu$','FontSize',40,'Interpreter','latex');
    hold on
    p=plot(I,MU,'c.-');
    p.MarkerSize=15;
    
    
    
    
    
    
    
    
    
end
toc
MU(end)


%% Residual
muF=sum(sum(HPsi(kk,V,g,Psi).*conj(Psi)))/sum(sum((abs(Psi.^2))))
DeltaPsi=ifftn(0.5*kk.*fftn(Psi))+(V+g*(abs(Psi).^2)-muF).*Psi;


mu=real(muF);
disp(mu);

% clear MU

Pic_maker;


