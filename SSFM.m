function [Psi, mu, dE] = SSFM(eD_hp_itp,dt_itp,V,g,NN0,kk,dV,Psi)

    Psi = ifftn(eD_hp_itp.*fftn(Psi));
    Psi = exp(-(V+g*(abs(Psi)).^2)*dt_itp).*Psi;
    Psi = ifftn(eD_hp_itp.*fftn(Psi));
    
    exp_mu = sqrt(NN0/(sum(sum((abs(Psi).^2)))*dV));
    Psi = Psi*exp_mu;

    mu = log(exp_mu)/dt_itp;
    H_Psi=HPsi(kk,V,g,Psi);
    dE=dV*sum(sum(abs(H_Psi).^2))-(dV*sum(sum(conj(Psi).*H_Psi))).^2;


end