%Chiller Model Examples
%Example 3
%Load chiller state
chiller(3);

%Define water-side boundary conditions
TEWI = 16;
TCWI = 30;
TEWO_SET = 10;
MEWAT = 13.2;
MCWAT = 16.7;

%Set-up plotting
FIG = figure;
set(FIG,'Position',[231 132 1128 908]);
subplot(311); axis([0 200 0 1500]);
grid on; hold on;
xlabel('s'); ylabel('kPa');
subplot(312); axis([0 200 0 50]);
grid on; hold on;
xlabel('s'); ylabel('deg C');
subplot(313); axis([0 200 0 100]);
grid on; hold on;
xlabel('s'); ylabel('kW');

%Execute chiller for 300s
%
%Update water side conditions every 10s
t = 2;
%Initialize counter and output-storage
i = 1;
output = [];
%Begin loop...
while(i<200)
    Tewi    = TEWI;
    Tcwi    = TCWI;
    Tewo_set= TEWO_SET;
    mewat   = MEWAT;
    mcwat   = MCWAT;
    if(i>=50)
        Tewi = TEWI - 2;
    end
    if(i>=75)
        Tcwi = TCWI + 2;
    end
    u = [Tewi;Tcwi;Tewo_set;mewat;mcwat];
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