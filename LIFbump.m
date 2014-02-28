%Leaky integrate-and-fire bump model, edited 02-28-14 by Shashaank
%Vattikuti

clear
%define some figure properties (more at end of code)
figure(1)
clf %comment out to keep prior figure
% hold on %comment in(out) to draw on(replace) figure 
fs=18;  %font size
ms=10;  %marker size
mc=[0 0 0]; %marker (face) color (RGB vector or string)

%simulation parameters
dt=single(0.05); %timestep for integration
T=2; %total simulated time in seconds
iterations=single(T*1000/dt); %number of iterations for integration

%test parameters (begin section)
N=50; %number of neurons, equal number of e and i neurons
        %if g_GABA fixed at 0.3 and g_NMDA at 0.48, then adjusting N
        %between 50 and 500 produces the density effect

%-----maximum conductance for postsynaptic neuron by channel-type

g_NMDA=.48;    %if g_GABA fixed at 0.3 and N=50, then can reduce to 50% before failure
                %around 10 x this is the max given the other parameters
g_GABA=.3; %if g_NMDA fixed at 0.48 and N=50, then can reduce to 5% (sometimes 3%) of 0.3 before failure
                %g_GABA=0.3 about max given other parameters

%---afferent current parameters 

ta1_i=0;   %(msec), start time for first input 
ta1_f=200; %(msec), stop time for first input 
mua1=0.5;
siga1=.1;
maga1=200;

ta2_i=1000; %(msec), second input
ta2_f=1600; 
mua2=0.75;
siga2=.1;
maga2=10; 

%1. Normalize physical distance to maximum length of 1.
% As such, any change in the number of neurons will not affect the synaptic
% footprints defined later. 

x=[1:N]/N; %normalized locations of neurons
dX(1,:)=single(abs(x-mua1));    %calculate first distance
dX(2,:)=single(abs(dX(1,:)-1));  %calculate second distance
w1=exp(-dX(1,:).^2./(2*siga1^2));
w2=exp(-dX(2,:).^2./(2*siga1^2));
tmp=single(w1+w2);
% tmp=exp(-((x-mua1).^2)/(2*(siga1^2))); %generate input vector for first input
Ia1(:,1)=maga1.*tmp/sum(tmp); %normalize vector to sum to one, then scale by magnitude

dX(1,:)=single(abs(x-mua2));    %calculate first distance
dX(2,:)=single(abs(dX(1,:)-1));  %calculate second distance
w1=exp(-dX(1,:).^2./(2*siga2^2));
w2=exp(-dX(2,:).^2./(2*siga2^2));
tmp=single(w1+w2);
% tmp=exp(-((x-mua2).^2)/(2*(siga2^2)));
Ia2(:,1)=maga2.*tmp/sum(tmp);


%test parameters (end section)

%construct network 


%2. Define synpatic footprint
%nomenclature according to row, column name of weight matrices (see "Define weight matrices")
sigmaEE=single(.4/3); %e to e, standard deviation is relative to max length of one
sigmaIE=single(.9/3); %e to i
sigmaEI=single(.9/3);  %i to e
%II is uniform

%3. Define weight matrices
%These are multiplied by 's' column vectors with the synapse state of each
%presynaptic neuron. As such, the rows of each weight matrix are the
%presynaptic weights into each postsynaptic neuron. 
%N(postsyn w)xN(presyn w) times NX1 (presyn mag)  = Nx1 (postsyn w mag).
%Also normalize the rows to sum to one. 

%E(post) rows x E(pre) columns, e to e
for i=1:N
    dX(1,:)=single(abs(x-x(i)));    %calculate first distance for periodic (ring) connections
    dX(2,:)=single(abs(dX(1,:)-1));  %calculate second distance
    w1=exp(-dX(1,:).^2./(2*sigmaEE^2));
    w2=exp(-dX(2,:).^2./(2*sigmaEE^2));
    EE(:,i)=single(w1+w2);
end
%normalize rows
for j=1:N
    EE(j,:)=EE(j,:)/sum(EE(j,:));
end

%I(post)xE(pre), e to i
for i=1:N
    dX(1,:)=single(abs(x-x(i)));    %calculate first distance
    dX(2,:)=single(abs(dX(1,:)-1));  %calculate second distance
    w1=exp(-dX(1,:).^2./(2*sigmaIE^2));
    w2=exp(-dX(2,:).^2./(2*sigmaIE^2));
    IE(:,i)=single(w1+w2);
end

for j=1:N
    IE(j,:)=IE(j,:)/sum(IE(j,:));
end

%E(post)xI(pre), i to e
for i=1:N
    dX(1,:)=single(abs(x-x(i)));    %calculate first distance
    dX(2,:)=single(abs(dX(1,:)-1));  %calculate second distance
    w1=exp(-dX(1,:).^2./(2*sigmaEI^2));
    w2=exp(-dX(2,:).^2./(2*sigmaEI^2));
    EI(:,i)=single(w1+w2);
end

for j=1:N
    EI(j,:)=EI(j,:)/sum(EI(j,:));
end

%I(post)xI(pre), i to i
II=1/N*ones(N);

%---SET ELECTROPHYSIOLOGIC PARAMETERS

%---------time constants in msec units----------------
tau_rise_NMDA=2;   %rise time constant for NMDA channels
tau_decay_NMDA=100;%NMDA decay time constant
tau_decay_GABA=10; %GABA decay time constant



%capacitance is in nF; the base units for a Farad (unit for capacitance) is (s^4*amp^2)/(m^2*kg)
%conductance is in uS; the base units for a Siemens (unit for channel conductance) is (s^3*amp^2)/(m^2*kg)

g_LE=.025; %leak conductance
g_LI=.02;

CapE=.5;
CapI=.2;

%---inverse capacitance
invCapE=1/CapE;
invCapI=1/CapI;


%voltage parameters (mV), base units for volts is (kg*m^2)/(s^3*amp)

V_L=-70;
V_rest=-60;
V_Th=-50;
V_reset=-55;
V_syn_NMDA=0;
V_syn_GABA=-70;

ve=V_rest*ones(N,1);
s1NMDA=zeros(N,1);         
s2NMDA=zeros(N,1);      

vi=V_rest*ones(N, 1);
sGABA=zeros(N,1);

M_e=zeros(N, 1);   
M_i=zeros(N, 1);   

spikes_e=zeros(N, 1); 
spikes_i=zeros(N, 1);



spikes_E=int8(zeros(N, iterations));

for i=2:iterations


    time=single(i*dt); %do not comment out
    
    timer=int32(time) %report sim time to user, can be commented out 
    
    %---update M and synapses---------
    
    M_e=(1./(1+exp(-0.062.*(ve))./3.57)); 
    
    M_i=(1./(1+exp(-0.062.*(vi))./3.57)); 
    
    
    s1NMDA=s1NMDA+spikes_e;
    s1NMDA=s1NMDA*(1-dt*1/tau_rise_NMDA);
    alpha=0.5;
    s2NMDA=s2NMDA+dt*(alpha*s1NMDA.*(1-s2NMDA)-s2NMDA.*1/tau_decay_NMDA);
    
    
    
    sGABA=sGABA+spikes_i;
    sGABA=sGABA.*(1-dt*1/tau_decay_GABA);
    
    
    
    %---update currents (nA)-------------
    
    I_LE=(g_LE*(ve-V_L));
    I_NMDAE=(g_NMDA*EE*s2NMDA.*M_e.*((ve-V_syn_NMDA)));
    I_GABAE=(g_GABA*EI*sGABA.*((ve-V_syn_GABA)));
    
    I_LI=(g_LI*(vi-V_L));
    I_NMDAI=(g_NMDA*IE*s2NMDA.*M_i.*((vi-V_syn_NMDA)));
    I_GABAI=(g_GABA*II*sGABA.*((vi-V_syn_GABA)));
    
    %---afferent current parameters (2)
    if time>=ta1_i & time<=ta1_f
        I_E=Ia1-I_NMDAE-I_GABAE-I_LE;
    elseif time>=ta2_i & time<=ta2_f
        I_E=Ia2-I_NMDAE-I_GABAE-I_LE;  
    else
        I_E=-I_NMDAE-I_GABAE-I_LE;
    end
    
    I_I=-I_NMDAI-I_GABAI-I_LI;
    
    
    
    %---update potentials----------
    ve2=invCapE.*dt.*I_E;
    vi2=invCapI.*dt.*I_I;
    
    ve=ve+ve2;
    
    
    vi=vi+vi2;
    
    %---check and adjust for spikes-----------
    
    spikes_e=ve>V_Th; %generates binary vector (1=spike, 0=no spike)
    ve=ve.*(1-spikes_e)+spikes_e*V_reset; %resets any neuron with v>V_th
    spikes_i=vi>V_Th;
    vi=vi.*(1-spikes_i)+spikes_i*V_reset;
    
    spikes_E(:,i)=int8(spikes_e); %store E neuron spikes for later analysis (i is current iteration)

end

f=calc_freq(spikes_E,dt); %calculate firing rates
pv=popvec(f); %calculate center of activity in degrees (population vector) 

figure(1);
set(gcf,'defaulttextinterpreter','tex');
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextFontname', 'Helvetica')
set(0,'DefaultAxesFontName', 'Helvetica')

%note plot metrics based on 100 msec bins used for frequency calculation
plot(pv,'k','linestyle','none','marker','sq','markerfacecolor',mc,'markersize',ms)
axis([0 T*10 0 360])
set(gca,'ytick',0:90:360)
set(gca,'xtick',5:5:T*10)
set(gca,'xticklabel',500:500:T*1000)
box off
title('population vector','fontsize',fs)
xlabel('milliseconds','fontsize',fs)
ylabel('\theta','fontsize',fs)
w=(ta1_f-ta1_i)/100;
rectangle('Position',[.01*ta1_i,mua1*360-.1*360,w,.2*360])
w=(ta2_f-ta2_i)/100;
rectangle('Position',[.01*ta2_i,mua2*360-.1*360,w,.2*360])
text(0.5,10,'box \approx location and dimensions of input','fontsize',fs*.8)
