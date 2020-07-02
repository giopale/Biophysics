Cm = 1e-2; %Membrane capacity F/m^2
F = 1.602e-19*6.022e23; %Faraday's constant
v  = 5e-15; %Cell volume, m^3
A = 1.4140e-9; %cell area, m^2
T=300; %Temperature, Kelvin
R=8.3144;   %Gas constant J/mol/K
Pk=10e-6; %K permeability, m/s
Pna=Pk/50; %Na permeability, m/s

Kout = 4;   %External ion concentration, M/m^3
Naout=145;
V0 = -65; %Resting potential, mV
K0=155;
Na0=12;

%Risultato: Y = [V, K, Na]

f1=@(t,y)[(Pk*F^2*y(1))/Cm/R/T * (y(2)*exp(F*y(1)/R/T) - Kout)/(1-exp(F*y(1)/R/T)) + (Pna*F^2*y(1)/R/Cm/T)*(y(3)*exp(F*y(1)/R/T)-Naout)/(1-exp(F*y(1)/R/T)); ... 
    (Pk*F*y(1)*A/v/R/T)*(y(2)*exp(F*y(1)/R/T)-Kout)/(1-exp(F*y(1)/R/T)); ...
    (Pna*F*y(1)*A/v/R/T)*(y(3)*exp(F*y(1)/R/T) - Naout)/(1-exp(F*y(1)/R*T))];
[T,Y]=ode15s(f1,[0,1e-2],[V0;K0;Na0]);


figure(123);
clf
tiledlayout(3,1);
nexttile
plot(T,Y(:,1));
xlabel('Time [s]');
ylabel('V');
dim = [.6 .6 .3 .3];
str = ['V finale = ',num2str(Y(end,1))];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
nexttile
plot(T,Y(:,2));
dim = [.6 .24 .3 .3];
str = ['[K] finale = ',num2str(Y(end,2))];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
ylabel('[K]');
xlabel('Time [s]');
nexttile
plot(T,Y(:,3));
dim = [.6 .001 .3 .3];
str = ['[Na] finale = ',num2str(Y(end,3))];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
ylabel('[Na]');
xlabel('Time [s]');