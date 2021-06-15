%% NPZ isotope model, April 28, 2020

clear all;
close all;

%% Parameters

r=0.0036765; % 15N/14N of air N2

vmax=0.4; % max growth rate of phyto, /day
gmax=0.7; % max grazing rate, /day. Change to 0.6 to reach steady state
%gmax=0.7;

k1=0.5; % half saturation of assimilation
k2=0.3; % half saturation of grazing

assim_eff_zoo=0.1; % zoo assimilation efficiency 
poop_rate=0.0; % fraction of the gut that goes to poop, 0.06
sink_rate=0.0; % fraction of the dead zoo that sinks, while the rest goes to NO3, 0.9

upwelling_rate=0.000; % constant rate supply of NO3, 0.0005
up15=upwelling_rate*r;
up14=upwelling_rate*(1-r);

dt=0.05; % time intervel
time=1500; % total time in days
%time = 1500;
ntstep=time/dt;

low=10e-10;

%% Initial conditions
[NO3,phyto,zoo,F_NO3,F_phyto,d15NO3,d15phyto,d15zoo,gut_15,gut_14,gut]=deal(zeros(1,ntstep));
[NO3_15,NO3_14,phyto_15,phyto_14,zoo_15,zoo_14]=deal(zeros(1,ntstep));

NO3(1)=5;
NO3_15(1)=NO3(1)*r;
NO3_14(1)=NO3(1)*(1-r);

phyto(1)=1;
phyto_15(1)=phyto(1)*r;
phyto_14(1)=phyto(1)*(1-r);

zoo(1)=0.5;
zoo_15(1)=zoo(1)*r;
zoo_14(1)=zoo(1)*(1-r);

gut(1)=0.00001;
gut_15(1)=zoo(1)*r;
gut_14(1)=zoo(1)*(1-r);


%% Loop
for i=1:ntstep
    
    %   Calculate isotop values in delta notation
    d15NO3(i)=(NO3_15(i)/NO3_14(i)/r-1)*1000;
    d15phyto(i)=(phyto_15(i)/phyto_14(i)/r-1)*1000;
    d15zoo(i)=(zoo_15(i)/zoo_14(i)/r-1)*1000;
    
    %   The fraction of 15N
    F_NO3(i)=NO3_15(i)/(NO3_15(i)+NO3_14(i));
    F_phyto(i)=phyto_15(i)/(phyto_15(i)+phyto_14(i));
    
    %   Michaelis Menten for NO3 uptake and grazing
    MMNO3=NO3(i)/(NO3(i)+k1);
    MMphyto=phyto(i)/(phyto(i)+k2);
    
    %   Total fluxes, DO NOT forget dt
    uptake=vmax*phyto(i)*MMNO3*dt; % assimilation
    grazing=gmax*zoo(i)*MMphyto*dt; % grazing
    death=0.1*zoo(i)/(zoo(i)+1)*dt; % zoo death
   
    %   Fluxes of 15N and 14N
    uptake15=uptake*F_NO3(i)/1.005;
    uptake14=uptake*(1-F_NO3(i));
    grazing15=grazing*F_phyto(i); % this grazing is directly to the pool of gut
    grazing14=grazing*(1-F_phyto(i));
        %   Here, the Transient gut pool = old gut pool + grazing at time step(i)
    pee15=(grazing15+gut_15(i))*(1 - assim_eff_zoo - poop_rate)/1.0035;
    pee14=(grazing14+gut_14(i))*(1 - assim_eff_zoo - poop_rate);
    zoo_assim15=(grazing15+gut_15(i))*assim_eff_zoo;
    zoo_assim14=(grazing14+gut_14(i))*assim_eff_zoo;
    
    NO3_15(i+1)=NO3_15(i) - uptake15 + pee15 + zoo_15(i)*death*(1-sink_rate) + up15*cos(2*pi/10000*i)*dt;
    NO3_14(i+1)=NO3_14(i) - uptake14 + pee14 + zoo_14(i)*death*(1-sink_rate) + up14*cos(2*pi/10000*i)*dt;
    
    phyto_15(i+1)=phyto_15(i) + uptake15 - grazing15;
    phyto_14(i+1)=phyto_14(i) + uptake14 - grazing14;
    
    zoo_15(i+1)=zoo_15(i)*(1-death) + zoo_assim15;
    zoo_14(i+1)=zoo_14(i)*(1-death) + zoo_assim14;
    
    %   Gut_15 pool is changing, while gut_14 pool is not, due to the isotope
    % fractionation of peeing
    gut_15(i+1)=gut_15(i) + grazing15 - pee15 - zoo_assim15;
    gut_14(i+1)=gut_14(i) + grazing14 - pee14 - zoo_assim14;
    
    %   New total pools
    NO3(i+1)=NO3_15(i+1)+NO3_14(i+1);
    phyto(i+1)=phyto_15(i+1)+phyto_14(i+1);
    zoo(i+1)=zoo_15(i+1)+zoo_14(i+1);
    
%     if NO3(i+1) < 0
%         NO3(i+1) = low;
%     else
%         NO3(i+1) = NO3(i+1);
%     end
%     if phyto(i+1) < 0
%         phyto(i+1) = low;
%     else
%         phyto(i+1) = phyto(i+1);
%     end
%     if zoo(i+1) < 0
%         zoo(i+1) = low;
%     else
%         zoo(i+1) = zoo(i+1);
%     end
%     
end
    %   The last time step isotope value in delta notation
    d15NO3(i+1)=(NO3_15(i+1)/NO3_14(i+1)/r-1)*1000;
    d15phyto(i+1)=(phyto_15(i+1)/phyto_14(i+1)/r-1)*1000;
    d15zoo(i+1)=(zoo_15(i+1)/zoo_14(i+1)/r-1)*1000;

%% Plots
t=(0:1:ntstep).*dt;
figure1 = figure;
plot(t,NO3,t,phyto,'--',t,zoo,':','LineWidth',2);
title('Population','fontsize',16);
legend('NO3','Phyto','Zoo','location','best','fontsize',12);
xlabel('Time (days)','fontsize',16);
ylabel('Concentration (µM)','fontsize',16)
saveas(figure1, 'population_p.png')

figure2 = figure;
plot(t,d15NO3,t,d15phyto,'--',t,d15zoo,':','LineWidth',2);
title('\delta^{15}N','fontsize',16);
legend('\delta^{15}N_{NO3}','\delta^{15}N_{Phyto}','\delta^{15}N_{Zoo}','location','best','fontsize',12);
xlabel('Time (days)','fontsize',16);
ylabel(['\delta^{15}N, ', char(8240) ' vs. air'],'fontsize',16);
saveas(figure2, 'd15N_p.png')




    
    
