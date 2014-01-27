%Chiller Model Examples
%Example 5
%Load chiller state
chiller(3);

%Define water-side boundary conditions
TEWI = 16;
TCWI = 30;
TEWO_SET = 10;
MEWAT = 13.2;
MCWAT = 16.7;

%Fault settings
tstart  = 10;   %time at which fault begins
tend    = 180;  %time to reach fully developed state
fullfault  = 0.4;  %severity when fully developed, 0.0 being nominal

%Set-up plotting
FIG = figure;
set(FIG,'Position',[231 132 1128 908]);
subplot(311); axis([0 300 0 1500]);
grid on; hold on;
xlabel('s'); ylabel('kPa');
subplot(312); axis([0 300 0 50]);
grid on; hold on;
xlabel('s'); ylabel('deg C');
subplot(313); axis([0 300 0 100]);
grid on; hold on;
xlabel('s'); ylabel('kW');

%Execute chiller for 300s
%
%Update water side conditions every 2s
t = 2;
%Initialize counter and output-storage
i = 1;
output = [];
%Begin loop...

Tewi    = TEWI;
Tcwi    = TCWI;
Tewo_set= TEWO_SET;
mewat   = MEWAT;
mcwat   = MCWAT;
while(i<300)
    
    %Fault introduction
    if(i>=tstart)
        if(i<=tstart+tend)
            instfault = fullfault*(i-tstart)/tend;  %current fault severity
            mcwat = (1-instfault)*MCWAT;
        end
    end
    %Updated boundary condition
    i
    u = [Tewi;Tcwi;Tewo_set;mewat;mcwat]
    j = 1;
    while(j<=t)
        y = chiller(1,u);
        output = [output;y'];        
        j = j + 1;
    end
    figure(FIG);
    subplot(311);
    plot(i,y(2),'b.',i,y(3),'r.');
    figure(FIG);    
    subplot(312);
    plot(i,y(12),'b.',i,y(13),'r.');
    figure(FIG);    
    subplot(313);
    plot(i,y(8),'r.');
    save output;
    i = i+t;
    pause(0.1);
end
%Save state at the end
chiller(2);