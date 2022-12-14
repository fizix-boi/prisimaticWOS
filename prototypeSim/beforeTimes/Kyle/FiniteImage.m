%%Constants

% Physical Conastants
i = sqrt(-1);
epHe = 1.056;
ep0 = 8.854e-12;
rho = 145;
sigma = 3.78e-4;
g = 9.81;
k = 0.6911e-6;
e = 1.602e-19;
um = 1e6;
nm = 1e9;
pcm = 1e-2;

% Experimental constants
d = 13e-9;
r0 = 0e-9;
ePos = [0,0,d+r0];

el = Electron(ePos,-e);

Gates = {};
Gates{1} = Geometry([0 0 200],[400 400 400],0);
%Gates{2} = Geometry([0 0 0e-6],[20e-6 20e-6 20e-6],0);

% Simulation Constants
epsilon = 1e-15;
totalWalks = 1e3;
walkMax = inf;

%% Debug Codes

CalcV = Gates{1}.V+(1/4/pi/ep0)*e/(2*el.home(3));
%hits = zeros([2 totalWalks]);

figure(10)
subplot(1,1,1)
hold off

figure(20)
subplot(1,1,1)
hold off

tic

%% Calculations
Voltage = 0;

for walkNum = 1:totalWalks

    %Walk Start
    el.reset;
    rStep = inf;
    for i = 1:length(Gates)
        if Gates{i}.distance(el.position) < rStep
            rStep = Gates{i}.distance(el.position);
        end
    end

    %Walk To Plate
    counter = 0;
    while rStep > epsilon && counter < walkMax
        el.walk(rStep);
        rStep = inf;
        for i = 1:length(Gates)
            if Gates{i}.distance(el.position) < rStep
                rStep = Gates{i}.distance(el.position);
            end   
        end
        counter = counter+1;
    end

    %Determine Voltage
    rStep = inf;
    gateNum = 0;
    for i = 1:length(Gates)
        if Gates{i}.distance(el.position) < rStep
            rStep = Gates{i}.distance(el.position);
            gateNum = i;
        end
    end
    
    Voltage = Voltage + (Gates{i}.V - (1/4/pi/ep0)*el.charge/el.distance());
    VoltGuess = Voltage/walkNum;
    Error = abs(VoltGuess-CalcV)/abs(CalcV)*100;
    figure(10)
    subplot(2,1,1)
    plot(walkNum,Error,'O','Color','#0072BD')
    hold on
    subplot(2,1,2)
    plot(walkNum,VoltGuess,'O','Color','#D95319')
    hold on
    if mod(walkNum,totalWalks/10) == 0
        fprintf('%d : %f\n',walkNum,toc)
    end

    %hits(1,walkNum) = el.position(1);
    %hits(2,walkNum) = el.position(2);
    
end
Voltage = Voltage/totalWalks;
fprintf('Simulation : %1.2e\nTheory : %1.2e\nError : %2.2f%%\n',[Voltage,CalcV,abs(Voltage-CalcV)/abs(CalcV)*100])

% figure(2)
% subplot(1,1,1)
% hold off
% for i = 1:totalWalks
%     plot(hits(1,i),hits(2,i),'O')
%     hold on
% end
% lim = el.home(3);
% xlim([-lim lim])
% ylim([-lim lim])


    














