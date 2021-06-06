% Azimuth off of array axis.
% Elevation between 45deg to 135deg.
clear
clf
tic
%% Radar Processor
total_distBins = {};
total_speedBins = {};
finalDistBins = {};
finalSpeedBins = {};
targets = {};
PosSpeeds = {};
NegSpeeds = {};
cFlag = 0;
exitLoop = 0;

CPI = 5e-3;
scan_angle = deg2rad(45:5:135);
scan_length = length(scan_angle);
El_mode = 2;

%% Initialization
radar_vars = Init_ADC();
radar_vars.PRI_counter = 0;
noise_variance = radar_vars.noise_variance;

Rmin = 3e8*radar_vars.tau/2;
N_max = round(30e3/Rmin);

lambda = 3e8/(radar_vars.f0);
beta = 2*pi/lambda;
d = lambda/2;
NumberOfElements = 20;

%% Initialize Graphs
for p = 1:3
    axList{2*p-1} = subplot(3,2,2*p-1, polaraxes);
    axList{2*p} = subplot(3,2,2*p, polaraxes);
    plotCell{2*p-1} = polarscatter(axList{2*p-1}, 0, 0);
    plotCell{2*p} = polarscatter(axList{2*p}, 0, 0);
    legend({'Positive', 'Negative'}, 'location', 'best')
end
sgt = sgtitle("Time: " + num2str(radar_vars.sim_time));
%% PRF Scheduler
% delF = 250;
delF = 350;
Fmin = delF*3e8/2/radar_vars.f0;
NF_max = ceil(3000/3.6/Fmin);
% PRF = [5.2e3 7.8e3 9e3 19e3];
% PRF = [17e3 15e3 19e3];
% PRF = [11e3 10e3 9e3 19e3];
PRF = [17e3 19e3];
np_max = CPI*PRF;
np = floor(PRF./delF);
CPI_met = all(np./np_max <= 1);
Pfa = 1e-8;
Vt = sqrt(-2.*noise_variance.*np*log(Pfa));

M = 2;
N = 3;

Refresh_Time = sum(np./PRF*N)*length(scan_angle)*3;
%% Ping the ADC
scanNo = 1;
while(1)
    if exitLoop == 1
        break
    end
    for p = 1:3;
        if exitLoop == 1
            break
        end
        radar_vars.EL = p;
        for z = 1:scan_length
            if exitLoop == 1
                break
            end
            radar_vars.PSI = rad2deg(-pi*cos(scan_angle(z)));
            for k = 1:length(PRF)
                if exitLoop == 1
                    break
                end
                radar_vars.PRF = PRF(k);
                for j = 1:N
                    % Accumulation of np pulses
                    for i = 1:np(k)
                        % Check if Time is Up
                        if radar_vars.sim_time+Refresh_Time > 30
                            exitLoop = 1;
                            break
                        end
                        [radar_vars, My_I_and_Q_bins] = Read_ADC_X2(radar_vars);
                        my_Complex_values= My_I_and_Q_bins{1} + 1i*My_I_and_Q_bins{2};
                        if i == 1
                            total_Complex_values = my_Complex_values;
                        else
                            total_Complex_values = [total_Complex_values; my_Complex_values];
                        end
                    end
                     if exitLoop == 1
                        break
                    end
                    radar_vars.PRI_counter = 0;
                    total_Abs_values = abs(total_Complex_values);
                    
                    % Doppler Processing
                    F_norm = -0.5: (1./np(k)) : (0.5 - 1./np(k));
                    Fourier_series = fft(total_Complex_values);
                    Fourier_series = fftshift(Fourier_series,1);
                    
                    % Erase Ground Return
                    if radar_vars.EL == 1
                        ZeroPoint = size(Fourier_series,1)/2 + 1;
                        Fourier_series(ZeroPoint,:) = 0;
                    end
                    % Threshold Detection
                    total_distBins{j,1} = any(abs(Fourier_series) > Vt(k));
                    total_speedBins{j,1} = any(rot90(abs(Fourier_series), -1) > Vt(k));
                    
                end
                
                total_distBins = cell2mat(total_distBins);
                total_speedBins = cell2mat(total_speedBins);
                
                % M of N
                if N ~= 1
                    finalDistBins{z,k} = sum(total_distBins) >= M;
                    finalSpeedBins{z,k} = sum(total_speedBins) >= M;
                else
                    finalDistBins{z,k} = sum(total_distBins);
                    finalSpeedBins{z,k} = sum(total_speedBins);
                end
                % Back to Empty Cells
                total_distBins = {};
                total_speedBins = {};
            end
        end
        %% Disambiguate
        for i = 1:scan_length
            targets{i,1} = coincidence(finalDistBins(i,:), N_max)*Rmin;
            [PosSpeeds{i,1}, NegSpeeds{i,1}] = coincidenceDoppler(finalSpeedBins(i,:), NF_max);
            PosSpeeds{i,1} = PosSpeeds{i,1}*Fmin;
            NegSpeeds{i,1} = NegSpeeds{i,1}*Fmin;
        end
        %% Make Graphable Values
        DisOut = [0 0];
        for i = 1:scan_length
            for j = 1:length(targets{i})
                Dnew = [scan_angle(i) targets{i}(j)];
                DisOut = [DisOut;Dnew];
            end
        end
        DisOut(1,:) = [];
        
        PosOut = [0 0];
        for i = 1:scan_length
            for j = 1:length(PosSpeeds{i})
                Dnew = [scan_angle(i) PosSpeeds{i}(j)];
                PosOut = [PosOut;Dnew];
            end
        end
        PosOut(1,:) = [];
        
        NegOut = [0 0];
        for i = 1:scan_length
            for j = 1:length(NegSpeeds{i})
                Dnew = [scan_angle(i) NegSpeeds{i}(j)];
                NegOut = [NegOut;Dnew];
            end
        end
        NegOut(1,:) = [];
        %% Log Results
        FinalOut{scanNo, 1, p} = radar_vars.sim_time;
        FinalOut{scanNo, 2, p} = DisOut;
        FinalOut{scanNo, 3, p} = PosOut;
        FinalOut{scanNo, 4, p} = NegOut;
        
        %% Graph PPI
        a = polarscatter(axList{2*p-1}, DisOut(:,1), DisOut(:,2));
        b = polarscatter(axList{2*p}, PosOut(:,1), PosOut(:,2),'b');
        hold on
        b = polarscatter(axList{2*p}, NegOut(:,1), NegOut(:,2),'r');
        legend({'Opening', 'Closing'}, 'location', 'best')
        hold off
        axList{2*p-1}.RLim = [0 30e3];
        axList{2*p-1}.Title.String = "Distance, El Mode: " + num2str(p);
%         axList{2*p}.RLim = [0 800];
        axList{2*p}.Title.String = "Velocity, El Mode: " + num2str(p);
        drawnow                         % Update elevations as they happen
    
    end
    % Update Time
    figno = gcf;
    figno.Name = "Time: " + num2str(radar_vars.sim_time);
    sgt.String = "Sim Time: " + num2str(radar_vars.sim_time) + ...
        "s    Real Time: " + num2str(toc) + "s";
%     drawnow       % Update all elevations at once
    if radar_vars.EL == 3
        cFlag = 1;
    end
    scanNo = scanNo + 1;
end
% Last Update
drawnow
sgt.String = "Sim Time: " + num2str(radar_vars.sim_time) + ...
        "s    Real Time: " + num2str(toc) + "s";
% Separate Final Data Into Elevations
El1Out = FinalOut(:,:,1);
EL2Out = FinalOut(:,:,2);
El3Out = FinalOut(:,:,3);
%% Sim Time
Total_Time = sum(np./PRF*N)*length(scan_angle) % Old but still valid for
% 1 elvation sweep
SimTime = radar_vars.sim_time
