%Chiller Model Examples
%Example 7
%
%Ensure that Initial_MINIMAL.txt has entries as follows:
%18.0 21.0 124.7
chiller(0);

%Set-up plotting
FIG = figure;
set(FIG,'Position',[231 132 1128 908]);
subplot(311); axis([0 1000 0 1500]);
grid on; hold on;
xlabel('s'); ylabel('kPa');
subplot(312); axis([0 1000 0 50]);
grid on; hold on;
xlabel('s'); ylabel('deg C');
subplot(313); axis([0 1000 0 100]);
grid on; hold on;
xlabel('s'); ylabel('kW');

%Execute chiller
%
%Update water side conditions every 10s
t = 10;
%Initialize counter and output-storage
i = 0;
output = [];
%Open input stream for reading BCs
fid = fopen('FFBCs.txt');
%Begin loop...
while(i<5110)
    [A,count] = fscanf(fid,'%5f');
    Tewi = A(1);
    Tcwi = A(2);
    Tewo_set = A(3);
    mewat = A(4);
    mcwat = A(5);
    u = [Tewi;Tcwi;Tewo_set;mewat;mcwat];
    j = 1;
    while(j<=t)
        y = chiller(1,u);
        output = [output;y'];        
        j = j + 1;
    end
    subplot(311);
    plot(i,y(2),'b.',i,y(3),'r.');
    subplot(312);
    plot(i,y(12),'b.',i,y(13),'r.');
    subplot(313);
    plot(i,y(8),'r.');
    save output;
    i = i+10;
    pause(0.1);
end
%Save state at the end
chiller(2);