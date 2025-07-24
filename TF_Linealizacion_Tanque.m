% clear all; close all; clc
% 
% %Parámetros
% F = 1;      % m^3/h
% V = 1;      % m^3
% R = 1.9859; % Kcal/(kmol*K)
% H = -5960;  % Kcal/kmol
% E = 11843;  % Kcal/kmol
% A = 34930800; % 1/h
% pCp = 500;  % Kcal/(m^3*K)
% UA = 150;   % Kcal/(K*h)
% 
% % Punto de operación
% CA_ss = 8.5698;    % kmol/m^3
% T_ss  = 311.2639;  % K
% Ti_ss = 300;       % K
% Tc_ss = 292;       % K
% CAi_ss = 10;       % kmol/m^3
% 
% T = linspace(300, 400, 200);  
% CA = linspace(8, 12, 200); 
% 
% k = A .* exp(-E ./ (R .* T));
% 
% % Ecuación no lineal
% nl = (F/V)*(Ti_ss - T) + (-H/pCp) .* k .* CA + (UA/(V*pCp))*(Tc_ss - T);
% 
% %linealizada
% lineal = 1.989203226 * (CA - CA_ss) + 1.0 * (Ti_ss - 300) + 0.3 * (Tc_ss - 292) - 2.349305992*(T - T_ss);
% 
% 
% figure;
% plot(T, nl, 'b', 'LineWidth', 2); hold on; 
% plot(T, lineal, 'r', 'LineWidth', 2);
% 
% plot(T_ss, 0, 'ro', 'MarkerSize', 6, 'LineWidth', 2)
% plot([T_ss T_ss], ylim, 'k--')
% plot(xlim, [0 0], 'k--')
% 
% xlabel('T [K]')
% ylabel('dT/dt')
% legend('No lineal', 'Linealizado')
% grid on;




%Funciones de transferencia

syms s
A = [-1.166879465  -0.08802902616; 1.989203226 -2.349305992]; 
B = [0; 0.3]; 
C = [0 1]; 
E = [1 0; 0 1];
D = 0;

%Entrada
I = eye(size(A));
Gp = C*inv(s*I - A)*B; % Proceso interno que hace matlab, aunque con tf está normalizado el resultado, o sea lo dividen por 10^8

sys1 = ss(A, B, C, D);
TF_Tc = tf(sys1); %Función de tranferencia para la entrada Tc


%Perturbaciones
Gd = C*inv(s*I - A)*E; %aunque con tf está normalizado el resultado, o sea lo dividen por 10^8

sys2 = ss(A, E, C, D);
TF_per_CAi = tf(sys2(1)); %Para la perturbación de CAi
TF_per_Ti = tf(sys2(2)); %Para la perturbación de Ti

Co = ctrb(A, B); %Matriz de controlabilidad
%rank(Co); %Calcula cuantas columnas L.I tiene

%{ 
Estabilidad interna
La estabilidad interna se refiere a si las variables internas 
(estados del sistema) se mantienen limitadas o no crecen indefinidamente en el tiempo
con pequeñas entradas o perturbaciones
%}
disp(" Valores propios de la matriz A son:" ) %polos
disp(eig(A))


%{ 
Estabilidad BIBO
Un sistema es BIBO estable si, para cualquier entrada acotada (es decir, que no crezca sin límites), 
la salida también está acotada.
%}
disp("Polo del sistema :")
disp(pole(tf(sys1)))
