% Gr�fica de relaci�n de bandas prohibidas y permitidas para
% Modelo Krunig-Penney
clear; clc;
% Constantes Generales
hbar = 1.054e-34;    % Constante Reducida de Planck [J*s]
me = 9.11e-31;       % Masa del electr�n libre [Kg]
meff = 1.08;         % Masa efectiva de un electr�n en el Silicio
me_eff = me*meff;    % Masa del electr�n en el Silicio
ev2J = 1.6e-19;      % Conversi�n de eV a Joules
hbar_eV = hbar/ev2J; % Constante Planck Reducida en [eV*s]
% Par�metros del programa
a = 15;              % ancho del pozo [nm]
b = 1;               % ancho de la barrera [nm]
Uo = 500e-3;         % Energ�a de la barrera de potencial [eV]
E = 0 : Uo/100 : Uo; % Valores de evaluaci�n de la energ�a
% Evaluamos la funci�n de la derecha
syms k;
fr = cos(k*(a+b));
% T�rminos de la ecuaci�n de la izquierda
alpha = sqrt(2*me_eff*E)/hbar_eV;
beta = sqrt(2*me_eff*(Uo-E))/hbar_eV;
% Graficamos para 3 casos de anchos de pozo y barrera
% 1) a=15 y b=1  [nm]
subplot(3,1,1);
% Evaluamos la funci�n de la izquierda
fl = -(beta.^2+alpha.^2)/(2*alpha.*beta)*sin(alpha*a).*sin(beta*b)+cos(alpha*a).*cos(beta*b)
fplot(fr,'--b'); hold on; plot(E,fl,'r'); grid on;
plot(E,ones(length(E)),'g'); plot(E,-ones(length(E)),'g');
axis([0 Uo -1.5 1.5]);
xlabel('Energ�a [eV]'); ylabel('f(E)')
title('a=15nm, b=1nm y Uo=500meV');
legend('cos(k(a+b)','f_i_z_q_u_i_e_r_d_a');
% 2) a=20 y b=3  [nm]
subplot(3,1,2);
a = 20; b = 3;
% Evaluamos la funci�n de la izquierda
fl = -(beta.^2+alpha.^2)/(2*alpha.*beta)*sin(alpha*a).*sin(beta*b)+cos(alpha*a).*cos(beta*b);
fplot(fr,'--b'); hold on; plot(E,fl,'r'); grid on;
plot(E,ones(length(E)),'g'); plot(E,-ones(length(E)),'g');
axis([0 Uo -1.5 1.5]);
xlabel('Energ�a [eV]'); ylabel('f(E)')
title('a=20nm, b=3nm y Uo=500meV');
legend('cos(k(a+b)','f_i_z_q_u_i_e_r_d_a');
% 2) a=50 y b=15  [nm]
subplot(3,1,3);
a = 50; b = 15;
% Evaluamos la funci�n de la izquierda
fl = -(beta.^2+alpha.^2)/(2*alpha.*beta)*sin(alpha*a).*sin(beta*b)+cos(alpha*a).*cos(beta*b);
fplot(fr,'--b'); hold on; plot(E,fl,'r'); grid on;
plot(E,ones(length(E)),'g'); plot(E,-ones(length(E)),'g');
axis([0 Uo -1.5 1.5]);
xlabel('Energ�a [eV]'); ylabel('f(E)')
title('a=50nm, b=15nm y Uo=500meV');
legend('cos(k(a+b)','f_i_z_q_u_i_e_r_d_a');
saveas(gcf,'relacionfE.png')