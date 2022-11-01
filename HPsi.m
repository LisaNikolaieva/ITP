function [HPsi] = HPsi(kk,V,g,Psi)

HPsi=ifftn(0.5*kk.*fftn(Psi))+(V+g*(abs(Psi).^2)).*Psi;

end