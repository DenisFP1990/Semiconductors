% Gr�fica de relaci�n E-k Para Pozo de Potencial finito
clear; clc;
% Constantes Generales
hbar = 1.054e-34;    % Constante Reducida de Planck [J*s]
me = 9.11e-31;       % Masa del electr�n libre [Kg]
meff = 1.08;         % Masa efectiva de un electr�n en el Silicio
me_eff = me*meff;    % Masa del electr�n en el Silicio
ev2J = 1.6e-19;      % Conversi�n de eV a Joules
hbar_eV = hbar/ev2J; % Constante Planck Reducida en [eV*s]
% Par�metros del programa
a = [0.5 1 1.5 2];   % 4 valores que tomara el ancho del pozo [nm]
Uo = 100e-3;         % Energ�a de la barrera de potencial [eV]
% Al modificar Uo se puede encontrar los diferentes casos
E = 0 : Uo/100 : Uo; % Valores de evaluaci�n de la energ�a
% T�rminos de la ecuaci�n
k = sqrt(2*me_eff*E)/hbar_eV;
alpha = sqrt(2*me_eff*(Uo-E))/hbar_eV;
% Evaluamos la funci�n de la izquierda para los 4 casos de a
for i = 1 : length(a)
  aux = tan(k*a(i));
  fl(:,i) = aux;
end
% Evaluamos la funci�n de la derecha independiente de a
fr = (2*alpha.*k)./(k.^2-alpha.^2);
% Evaluaci�n Num�rica utilizando una libreria de terceros (cint)
[E1,fE1] = cint(E,fr,E,fl(:,1));
[E2,fE2] = cint(E,fr,E,fl(:,2));
[E3,fE3] = cint(E,fr,E,fl(:,3));
[E4,fE4] = cint(E,fr,E,fl(:,4));
% Gr�ficas
subplot(2,2,1);
plot(E,fr,'r'); hold on; grid on;
plot(E,fl(:,1),'b'); plot(E1,fE1,'*');
axis([0 Uo -20 20])
xlabel('Energ�a [eV]'); ylabel('f(E)')
title('f(E)-E para a=0.5 nm y Uo=100meV');
legend('(2af*k)/(k^2-af^2)','tan(ka)');
subplot(2,2,2);
plot(E,fr,'r'); hold on; grid on;
plot(E,fl(:,2),'b'); plot(E2,fE2,'*');
axis([0 Uo -20 20])
xlabel('Energ�a [eV]'); ylabel('f(E)')
title('f(E)-E para a=1 nm y Uo=100meV');
legend('(2af*k)/(k^2-af^2)','tan(ka)');
subplot(2,2,3);
plot(E,fr,'r'); hold on; grid on;
plot(E,fl(:,3),'b'); plot(E3,fE3,'*');
axis([0 Uo -20 20])
xlabel('Energ�a [eV]'); ylabel('f(E)')
title('f(E)-E para a=1.5 nm y Uo=100meV');
legend('(2af*k)/(k^2-af^2)','tan(ka)');
subplot(2,2,4);
plot(E,fr,'r'); hold on; grid on;
plot(E,fl(:,4),'b'); plot(E4,fE4,'*');
axis([0 Uo -20 20])
xlabel('Energ�a [eV]'); ylabel('f(E)')
title('f(E)-E para a=2 nm y Uo=100meV');
legend('(2af*k)/(k^2-af^2)','tan(ka)');
saveas(gcf,'relacionfE.png')