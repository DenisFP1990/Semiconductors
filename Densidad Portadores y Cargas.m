clear; clc;
% Constantes Físicas
me = 9.11e-31;      % Masa Electrón en espacio libre [Kg]
kb_eV = 8.617e-5;   % Constante de Boltzman [ev/K]
kb_J = 1.38e-23;    % Constante de Boltzman [J/K]
hbar_J = 1.054e-34; % Constante de planck reducida [J*s]
hbar_eV = 6.582e-16;% Constante de planck reducida [eV*s]
h = 6.626e-34;      % Constante de planck [J*s]
% Masas efectivas para cálculo de densidad de estados
me_eff = me*1.12;   % Masa efectiva de electrones en Silicio
mh_eff = me*1.08;   % Masa efectiva de huecos en Silicio
% Parámetros del semiconductor
Eg = 1.12; Ec = Eg; Ev = 0; % Energía de band gap, CB y VB [eV]
T = 300;            % Temperatura de estudio [Kelvin]
Ei = (Ec-Ev)/2 + 3/4*kb_eV*T*log(me_eff/mh_eff); % Energía intrínseca [eV]
Ef = Ei+24*kb_eV*T;  % Energía de Fermi caso Extrínseco N [eV]
% Cálculo de Densidades de estados
Ecb = 1.12 : 0.01 : 2; % Variable energética para DoS de CB [eV]
Evb = -1 : 0.01 : 0;   % Variable energética para DoS de VB [eV]
Etot = -1 : 0.01 : 2;  % Variable energética para función Fermi-Dirac [eV]
gcE = 0.01^3*((2*me_eff)^1.5/(2*pi^2*hbar_J^3))*sqrt(Ecb-Ec); % [1/eV*cm^3]
gvE = 0.01^3*((2*mh_eff)^1.5/(2*pi^2*hbar_J^3))*sqrt(Ev-Evb); % [1/eV*cm^3]
fE = 1./(1+exp((Etot-Ef)/(kb_eV*T)));       % Función Fermi-Dirac
fE_CB = 1./(1+exp((Ecb-Ef)/(kb_eV*T)));     % Función Fermi-Dirac en CB
fE_VB = 1-1./(1+exp((Evb-Ef)/(kb_eV*T)));   % Función Fermi-Dirac en VB
eta_e = (Ef-Ec)/(kb_eV*T); % Parámetro eta para electrones
eta_h = (Ev-Ef)/(kb_eV*T); % Parámetro eta para huecos
fermiInt_e = fermi(1/2,eta_e); % Cálculo de integral de Fermi orden 1/2 para e
fermiInt_h = fermi(1/2,eta_h); % Cálculo de integral de Fermi orden 1/2 para h
% Cálculo de densidad de estados efectivos y densidad de portadores
Nc = (0.01^3)*2*(2*pi*me_eff*kb_J*T/h^2)^1.5; % en CB [1/cm^3]
Nv = (0.01^3)*2*(2*pi*me_eff*kb_J*T/h^2)^1.5; % en VB [1/cm^3]
n = Nc*fermiInt_e;      % Cálculo de densidad de electrones
p = Nv*fermiInt_h;      % Cálculo de densidad de huecos
% Gráficas y Resultados
subplot(1,3,1); plot(fE,Etot); grid on; hold on; plot(1-fE,Etot); hold off;
line([0 1],[Ec Ec],'Color','black','LineStyle','-.');
line([0 1],[Ev Ev],'Color','black','LineStyle','-.');
line([0 1],[Ef Ef],'Color','black','LineStyle','-.');
text(0.5,Ef+0.1,strcat('EF=',num2str(Ef,3)));
xlabel('f(E) Fermi-Dirac'); ylabel('Energía [eV]');
legend('f(E)','1-f(E)'); title('Factor Ocupación');
subplot(1,3,2); plot(gcE,Ecb); grid on; hold on; plot(gvE,Evb); hold off;
line([0 15e49],[Ec Ec],'Color','black','LineStyle','-.');
line([0 15e49],[Ev Ev],'Color','black','LineStyle','-.');
line([0 15e49],[Ef Ef],'Color','black','LineStyle','-.');
text(7e49,Ef+0.1,strcat('EF=',num2str(Ef,3)));
xlabel('DoS(E)'); ylabel('Energía [eV]');
legend('gc(E)','gv(E)'); title('Densidad Estados');
subplot(1,3,3);
plot(fE_CB.*gcE,Ecb); grid on; hold on; plot(fE_VB.*gvE,Evb); hold off;
line([0 2e49],[Ec Ec],'Color','black','LineStyle','-.');
line([0 2e49],[Ev Ev],'Color','black','LineStyle','-.');
line([0 2e49],[Ef Ef],'Color','black','LineStyle','-.');
text(1e49,Ef+0.1,strcat('EF=',num2str(Ef,3)));
xlabel('DoS*f(E)'); ylabel('Energía [eV]');
legend('f(E)*gc(E)','[1-f(E)]*gv(E)'); title('Distribución Portadores');
text(3e39,1.5,strcat('n=',num2str(n,3)));
text(3e39,-0.5,strcat('p=',num2str(p,3))); saveas(gcf,'simulacion.png');