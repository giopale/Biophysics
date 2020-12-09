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
disp('Computing signal... ')
Ra= 1e7; %Ohm
Rm=1e8; %Ohm
Cm=3e-11; %F
tau=Cm*(1/Rm + 1/Ra)^(-1);


Vin=zeros(Nstep+1,1);
Vm=zeros(Nstep+1,1);
I0=zeros(Nstep+1,1); % Current without noise
I1=I0; % Current with noise 

Vin=step_fun(Vin,time,tstart,tstop,Vin_amplitude);

noise_amp= 1.5e-10;

f10=figure(10);
f10.Visible='off';
clf
pl1=plot(time,Vin);
prop = {"Input voltage", "Time [s]", "Voltage [mV]","m","y"};
SetPlot(get(gcf), prop)

for i=1:1:Nstep
	dVm=(-Vm(i)+Vin(i)*Rm/(Rm+Ra))/tau;
	I0(i) = Vm(i)/Rm+Cm*dVm;
	Vm(i+1) = Vm(i)+dVm*tstep;
	I1(i)=I0(i)+noise_amp*randn;
end
f20 = figure(20);
f20.Visible='off';
clf
pl2=plot(time,I1); %, time*10^3, I1*10^12)
prop = {"Current", "Time [ms]"," Current [pA]","m","y"};
SetPlot(get(gcf), prop)
hold on
disp('Saving plot... ')

%Ra, Rm estimation:
peak1=PeakEstim(time,I1,tstart,5); 
Ra1=Vin_amplitude/peak1;
I1_infty_start = ceil((tstart+5*tau)/tstep);	%I1_infty estimation from data
I1_infty_end = ceil(tstop/tstep)-10;
I1_infty=mean(I1(I1_infty_start:I1_infty_end));
Rm1=Vin_amplitude/I1_infty - Ra1;

%Cm estimation via curve fitting
x= time(time>tstart & time<tstart+13*tau);	%Select right time interval
y=I1(time>tstart & time<tstart+13*tau) - I1_infty; %subtract baseline
fitres=fit(x,y,'exp1','startpoint',[1e-9,-1/0.2545e-3]); %fit specifying initial parameters
tau1=-1/fitres.b;
Cm1=tau1*(1/Ra1 +1/Rm1);
y1=fitres(x)+I1_infty;
pl3=plot(x,y1,'r-',x,y+I1_infty,'b.');
hold off

% Errors
sigpt=noise_amp;
peak1_sig=sigpt/sqrt(5);
Ra1_sig=Ra*peak1_sig/peak1;
I1_infty_sig=sigpt/sqrt(size(I1(I1_infty_start:I1_infty_end),1));
Rm1_sig=sqrt(Rm^2*(I1_infty_sig/I1_infty)^2 + Ra1_sig^2);
bounds=confint(fitres);
bounds=-(1./bounds);
bounds(:,1)=[];
tau1_sig=(bounds(2)-bounds(1))/2;
R=(1/Ra1 + 1/Rm1);
R_sig=sqrt((Ra1_sig/Ra^2)^2 + (Rm1_sig/Rm1^2)^2);
Cm1_sig=Cm1*sqrt((tau1_sig/tau1)^2 + (R_sig/R)^2);







% --------------- Resting potential --------------------
disp('Computing signal with resting potential... ')
Vrest=-0.08; %mV, resting potential
I0r=I0-I0;
I1r=I0r;
Vmr=zeros(Nstep+1,1);
Vmr(1)=Vrest*Ra/(Rm+Ra);
for i=1:1:Nstep+1
	Vc=(Rm*Vin(i)+Ra*Vrest)/(Ra+Rm);
	dVmr=-1/tau*(Vmr(i)-Vc);
	I0r(i) = ((Vmr(i)-Vrest)/Rm + Cm*dVmr);
	I1r(i)=I0r(i)+noise_amp*randn;
	Vmr(i+1) = Vmr(i)+dVmr*tstep;
end

%Ra, Rm estimation:
peak1=PeakEstim(time,I1r,tstart,5); 
I1r_infty_start = ceil((tstart+7*tau)/tstep);	%I1_infty estimation from data
I1r_infty_end = ceil(tstop/tstep)-10;
I1r_infty=mean(I1r(I1r_infty_start:I1r_infty_end));
Ra1r=Vin_amplitude/(peak1-I1r_infty);
Rm1r=Vin_amplitude-Vrest/I1r_infty;


tau1r=0;
Cm1r=0;

f10=figure(30); %plot input voltage
f10.Visible='off';
clf
pl1=plot(time,Vin);
prop = {"Input voltage, resting pot.", "Time [s]", "Voltage [mV]","m","y"};
SetPlot(get(gcf), prop)

f40=figure(40); %plot solution for current
f40.Visible='off';
clf
pl4=plot(time,I1r);
prop = {"Current w. Resting Potential", "Time [s]", "Current [A]","m","y"};
SetPlot(get(gcf), prop)
hold on
grid on

%Cm estimation via curve fitting
xr= time(time>tstart & time<tstart+13*tau);	%Select right time interval
yr=I1r(time>tstart & time<tstart+13*tau) - I1r_infty; %subtract baseline
fitres1=fit(xr,yr,'exp1','startpoint',[1e-9,-1/0.2545e-3]); %fit specifying initial parameters
tau1r=-1/fitres.b;
Cm1r=tau1r*(1/Ra1r +1/Rm1r);
y1r=fitres(xr)+I1r_infty;
pl5=plot(xr,y1r,'r-',xr,yr+I1r_infty,'b.');
disp('Saving plot... ')

% Errors
sigpt=noise_amp;
peak1_sig=sigpt/sqrt(5);
Ra1r_sig=Ra1r*peak1_sig/peak1;
I1r_infty_sig=sigpt/sqrt(size(I1r(I1r_infty_start:I1r_infty_end),1));
Rm1r_sig=sqrt(Rm1r^2*(I1r_infty_sig/I1r_infty)^2 + Ra1r_sig^2);
bounds=confint(fitres1);
bounds=-(1./bounds);
bounds(:,1)=[];
tau1r_sig=(bounds(2)-bounds(1))/2;
R=(1/Ra1r + 1/Rm1r);
R_sig=sqrt((Ra1r_sig/Ra1r^2)^2 + (Rm1r_sig/Rm1r^2)^2);
Cm1r_sig=Cm1r*sqrt((tau1r_sig/tau1r)^2 + (R_sig/R)^2);




% Display some results on screen
% disp('Caso base')
% varnames={'I1_infty [pA]','Ra1 [MOhm]','Rm1 [MOhm]','tau1 [ms]','Cm1 [pf]'};
% T=table([I1_infty,I1_infty_sig]*1e12, [Ra1,Ra1_sig]*1e-6, [Rm1,Rm1_sig]*1e-6,[tau1,tau1_sig]*1e3,[Cm1,Cm1_sig]*1e12,'VariableNames',varnames);
% disp(T);

% disp('Caso resting potential')
% varnames={'I1r_infty [pA]','Ra1r [MOhm]','Rm1r [MOhm]','tau1r [ms]','Cm1r [pf]'};
% Tr=table([I1r_infty,I1r_infty_sig]*1e12, [Ra1r,Ra1r_sig]*1e-6, [Rm1r,Rm1_sig]*1e-6,[tau1r,tau1r_sig]*1e3,[Cm1r,Cm1r_sig]*1e12,'VariableNames',varnames);
% disp(Tr);


% ------------------ Bessel filtering ------------------
% Initializing structure fields
npoints=5;

signal(1).time=time;
signal(1).I1=I1;
signal(1).If=I1-I1;
signal(1).base=0;
signal(1).fcutoff=0;
signal(1).fig=figure('Visible',false);

fcutoff=logspace(0,1,npoints)*1e3;
disp('Filtering signals:')
for i=1:1:npoints
	signal(i).time=time;
	signal(i).I1=I1;
	signal(i).If=I1-I1;
	signal(i).fcutoff=fcutoff(i);
	%filtering
	[b, a] = besself(4,fcutoff(i)*2*pi);
	[bz, az] = impinvar(b,a,sr);
	signal(i).If=filter(bz,az,signal(1).I1);
	%fitting
	time_peak=time(signal(i).If == max(signal(i).If));
	tstart_fit=time_peak+0.1e-3;
	signal(i).base=baseline(signal(i).If,tstart, tstop, tstep,tau);
	[tauf,x, y, fitf]=fitthis(signal(i).time,signal(i).If,tstart_fit, tau, signal(i).base);

	signal(i).fig=figure('Visible',false);
	tit='Filtered signal: '+ string(round(fcutoff(i)/1e3,1)) +' kHz';
	prop = {tit, "Time [s]", "Current [A]","m","y"};
	clf
	plot(time,signal(i).If);
	hold on
	plot(x,y,'b.',x,fitf,'r-');
	SetPlot(get(gcf), prop);
	grid on
	filename='Bess'+string(ceil(fcutoff(i)))+'Hz';
	printpdf(signal(i).fig,filename);
	disp('plot '+string(i)+' over '+string(npoints) +' saved')
end

%Cm estimation via curve fitting
% i=1;
% tstart_fit=tstart+0.7e-3;
% signal(i).base=baseline(signal(i).If,tstart, tstop, tstep,tau);
% [tauf,x, y]=fitthis(signal(i).time,signal(i).If,tstart_fit, tau, signal(i).base);











disp("... Racheeeeeeel!")



function [out_vec] = step_fun(in_vec, time, tstart, tstop, amp)
	in_vec=in_vec-in_vec;
	in_vec(time>=tstart)= amp; %Step function: Vin_amplitude*(theta(tstart)-theta(-tstop))
	in_vec(time>=tstop)=0;
	out_vec=in_vec;
end

function [] = SetPlot(fig, prop )
	font="CMU Sans Serif";
	fontbold= "CMU Sans Serif Bold";
	ax=fig.CurrentAxes;
	ax.Title.String=prop(1);
	ax.Title.FontName=fontbold;
	ax.XLabel.String=prop(2);
	ax.XLabel.FontName=font;
	ax.YLabel.String=prop(3);
	ax.YLabel.FontName=font;
	
	switch prop{4}
	case "s"
		TitFS=18;
		LabFS=12;
	case "m"
		TitFS=22;
		LabFS=18;
	case "l"
		TitFS=26;
		LabFS=20;
	otherwise
		TitFS=20;
		LabFS=15;
	end
	ax.Title.FontSize=TitFS;
	ax.XLabel.FontSize=LabFS;
	ax.YLabel.FontSize=LabFS;

	if(size(prop,2)>=5e-3)
		switch prop{5}
		case "y"
			ax.YLim=enlarge(ax.YLim,0.1);
		case "x"
			ax.XLim=enlarge(ax.XLim,0.1);
		case "xy"
			ax.YLim=enlarge(ax.YLim,0.1);
			ax.XLim=enlarge(ax.XLim,0.1);
		otherwise
			return
		end
	end
end

function [outint]=enlarge(inint,hwmuch)
	interv= abs(inint(2)-inint(1));
	piece= 0.5*(hwmuch*interv);
	if inint(1)<1e-12
		outint= [inint(1), inint(2)+piece];
	else
		outint= [inint(1)-piece, inint(2)+piece];
	end
end

function [peak_val]=PeakEstim(time,signal,peak_time,mean_int)
	%This function performs a mean of the values of
	%the vector signal from time_loc to time_loc+mean_int
	%
	timest=time(2)-time(1);
	peak_idx=ceil(peak_time/timest)+1;
	peak_val=mean(signal(peak_idx:peak_idx+mean_int));
end

function []=printpdf(fig_handle,filename)

	set(fig_handle,'Units','Inches');
	pos = get(fig_handle,'Position');
	set(fig_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print(fig_handle,filename,'-dpdf','-bestfit');

end

function [I1_infty]=baseline(I1, tstart, tstop, tstep, tau)
	I1_infty_start = ceil((tstart+5*tau)/tstep);	%I1_infty estimation from data
	I1_infty_end = ceil(tstop/tstep)-10;
	I1_infty = mean(I1(I1_infty_start:I1_infty_end));
end


function [tau1,x, y, fitfun]=fitthis(xdata, ydata, tstart_fit, tau, baseline)
	x= xdata(xdata>tstart_fit & xdata<tstart_fit+4*tau);	%Select right time interval
	y= ydata(xdata>tstart_fit & xdata<tstart_fit+4*tau) - baseline; %subtract baseline
	fitres=fit(x,y,'exp1','startpoint',[1e-9,-1/0.2545e-3]); %fit specifying initial parameters
	tau1=-1/fitres.b;
	% Cm1=tau1*(1/Ra1 +1/Rm1);
	y=y+baseline;
	fitfun=fitres(x)+baseline;
	% figure(666)
	% pl3=plot(x,y1,'r-',x,y,'b.');

end
