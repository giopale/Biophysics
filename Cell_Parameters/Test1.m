clear all
disp("Puo' prendere la Rolls, ma me la riporti col pieno.")
% Physics constants
kb=1.38e-23;
% Simulation parameters
sr= 1e5; % Hz, sampling rate of the virtual oscilloscope
dur= 7e-3; % seconds, duration of the simulation
tstart=1e-3;
tstop=5e-3;
Vin_amplitude = 10e-3;
tstep=1/sr;
time=(0:tstep:dur)';
Nstep=size(time,1)-1;

% Cell parameters
Ra= 1e7; %Ohm
Rm=1e8; %Ohm
Cm=3e-11; %F
tau=Cm*(1/Rm + 1/Ra)^(-1);


Vin=zeros(Nstep+1,1);
Vm=zeros(Nstep+1,1);
I0=zeros(Nstep+1,1);
I1=I0;

Vin=step_fun(Vin,time,tstart,tstop,Vin_amplitude);

s90 = (4*pi^2*Cm^2*4*kb*300*Ra*sr^2)/(1+4*pi*sr^2*10^12*Cm^2);
amp90= sqrt(s90*sr/2);

figure(10)
clf
pl1=plot(time,Vin);
prop = {"Input voltage", "Time [s]", "Voltage [mV]","m","y"};
SetPlot(get(gcf), prop)

for i=1:1:Nstep
	dVm=(-Vm(i)+Vin(i)*Rm/(Rm+Ra))/tau;
	I0(i) = Vm(i)/Rm+Cm*dVm;
	Vm(i+1) = Vm(i)+dVm*tstep;
	I1(i)=I0(i)+amp90*randn;
end

figure(20)
plot(time*10^3,I0*10^12, time*10^3, I1*10^12)
prop = {"Current", "Time [s]"," [mV]","m","y"};
SetPlot(get(gcf), prop)












disp("... Racheeeeeeel!")




function [out_vec] = step_fun(in_vec, time, tstart, tstop, amp)
	in_vec=in_vec-in_vec;
	in_vec(time>=tstart)= amp; %Step function: Vin_amplitude*(theta(tstart)-theta(-tstop))
	in_vec(time>=tstop)=0;
	out_vec=in_vec;
end

function [] = SetPlot(fig, prop )
	font="CMU Sans Serif"
	fontbold= "CMU Sans Serif Bold"
	ax=fig.CurrentAxes;
	ax.Title.String=prop(1);
	ax.Title.FontName=fontbold;
	ax.XLabel.String=prop(2);
	ax.XLabel.FontName=font;
	ax.YLabel.String=prop(3);
	ax.YLabel.FontName=font;
	
	switch prop{4}
	case "s"
		TitFS=18
		LabFS=12
	case "m"
		TitFS=22
		LabFS=18
	case "l"
		TitFS=26
		LabFS=20
	otherwise
		TitFS=20
		LabFS=15
	end
	ax.Title.FontSize=TitFS
	ax.XLabel.FontSize=LabFS
	ax.YLabel.FontSize=LabFS

	if(size(prop,2)>=5e-3)
		switch prop{5}
		case "y"
			ax.YLim=enlarge(ax.YLim,0.1)
		case "x"
			ax.XLim=enlarge(ax.XLim,0.1)
		case "xy"
			ax.YLim=enlarge(ax.YLim,0.1)
			ax.XLim=enlarge(ax.XLim,0.1)
		otherwise
			return
		end
	end
end

function [outint]=enlarge(inint,hwmuch)
	interv= abs(inint(2)-inint(1))
	piece= 0.5*(hwmuch*interv)
	if inint(1)<1e-12
		outint= [inint(1), inint(2)+piece]
	else
		outint= [inint(1)-piece, inint(2)+piece]
	end
end
