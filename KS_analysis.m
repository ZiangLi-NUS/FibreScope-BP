% Liza's BP filter
clc;clear; close all;

%% Change accordingly
file_count_selection = 1;
arm_girth = 32.2; % unit: cm

% above this amount will be considered as a KS peak
blood_pulse_threshold = 0.1;
% larger than this time(s) interval will be considered as a KS peak
min_peak_distance = 0.6;
% smaller than this time(s) interval will be considered as a KS peak
max_peak_distance = 1.7;
% number of consequent data to check
pulse_check_length = 5;
% Write audio file or not, 1 for write
audio_write_flag = 1;
% Save figure, 1 for write
figure_save_flag = 1;

%% Initialization
title_font_size = 40;
label_font_size = 35;
xytick_font_size = 30;
File = dir(fullfile('*.csv'));
FileNames = {File.name}'; 
file_name = FileNames{file_count_selection};
file_name = erase(file_name,'.csv');
disp(file_name);
filename_char = convertStringsToChars(file_name);
csv_to_mat(file_name);
mat_file = file_name + ".mat";
load(mat_file);
disp('[INFO] Successfully open raw data file');
Sample1 = voltage_array;
fs = fs_data;
AP_fit_f = [0.2986,-71.1601]; % convert analog pressure sensor voltage to pressure
if(mean(pressure_array)<2)
    pressure_fit = polyval(AP_fit_f,pressure_array*1000);
else
    pressure_fit = polyval(AP_fit_f,pressure_array);
end
fix_coefficient_f = [0.00417,-0.34167,7.4];
fix_coefficient = polyval(fix_coefficient_f,arm_girth);
AB_fit = [0.3313, -17.8113]; % convert airep ressure to arm feeling pressure
AB_fit(1) = AB_fit(1)/fix_coefficient;
pressure_fit = polyval(AB_fit,pressure_fit);
test_point = pressure_fit;

disp('[PROCESS] Smoothing Pressure Data');
pressure_fit = smooth(round(0.01*fs),pressure_fit);
offset = mean(pressure_fit(0.5*fs:1*fs));
pressure_fit = pressure_fit - offset;
[~, max_pos] = max(pressure_fit(1:60*fs));
front_cuttoff_count = max_pos+0.5*fs;
end_cuttoff_pos = find(pressure_fit<40);
end_cuttoff_time = find(end_cuttoff_pos>50*fs,1);
end_cuttoff_count = end_cuttoff_pos(end_cuttoff_time);
Sample1 = Sample1(front_cuttoff_count:end_cuttoff_count);
pressure_data = pressure_fit(front_cuttoff_count:end_cuttoff_count);

run_time = length(Sample1)/fs; % acquire recorded data length
sample_count = 1:length(Sample1);
sample_time = sample_count/fs;


raw_data_plot = figure(1);
set(raw_data_plot,'Position',[50,50,1000,950]);
subplot(2,1,1)
raw_data_x = 1:length(voltage_array);
raw_data_x = raw_data_x/fs;
raw_data_run_time = length(voltage_array)/fs;
plot(raw_data_x, voltage_array,'b');
set(gca,"box","on",'linewidth',2,'FontSize',xytick_font_size)
xlim([0 raw_data_run_time])
title("Raw Data", 'FontSize', title_font_size)
xlabel("Time (s)", 'FontSize', label_font_size)
ylabel("Voltage (V)", 'FontSize', label_font_size)
grid on;

subplot(2,1,2);
plot(raw_data_x, pressure_fit,'r','LineWidth',2);
set(gca,"box","on",'linewidth',2,'FontSize',xytick_font_size)
xlim([0 raw_data_run_time])
title("Pressure Data", 'FontSize', title_font_size);
xlabel('Time (s)', 'FontSize', label_font_size);
ylabel('Pressure (mmHg)', 'FontSize', label_font_size);
grid on;

%% filter
disp('[INFO] Filtering KS')
disp('[PROCESS] Running high/low pass filter');
Sample1 = high_pass(Sample1, fs, 40);
Sample1 = low_pass(Sample1, fs, 100);
disp('[PROCESS] Running notch filter');
for i = 1:3
    multicancel = 2;
    for j = 1:multicancel
        L = length(Sample1);
        sample_fft = fft(Sample1);
        P2 = abs(sample_fft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = fs*(0:(L/2))/L;
        Sample1 = notch(50*i-40,50*i+1,fs,f,P1,Sample1,1);
    end
end
disp('[PROCESS] Normalize signal')
Sample1 = normalize(Sample1);

%% Plot
filtered_data_plot = figure(2);
set(filtered_data_plot,'Position',[50,50,1000,950]);
subplot(2,1,1)
plot(sample_time,Sample1,'b');
set(gca,"box","on",'linewidth',2,'FontSize',xytick_font_size)
xlim([0 run_time])
title("Filtered korotkoff sounds",'FontSize',title_font_size)
xlabel("Time (s)",'FontSize',label_font_size)
ylabel("Amplitude",'FontSize',label_font_size)

sample_fft = fft(Sample1);
L = length(Sample1); 
P2 = abs(sample_fft/L);
P1 = P2(1:round(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:round(L/2))/L;
subplot(2,1,2)
plot(f,P1,'r');
set(gca,"box","on",'linewidth',2,'FontSize',xytick_font_size)
xlim([2 200])
title("FFT analysis",'FontSize',title_font_size)
xlabel("Frequency/Hz",'FontSize',label_font_size)
ylabel("Amplitude",'FontSize',label_font_size)

%% analysis pulse time stamp

[peak_amp,peak_loc] = findpeaks(abs(Sample1),'MinPeakProminence', blood_pulse_threshold,...
    'MinPeakDistance', min_peak_distance*fs);
peak_loc = peak_loc/fs;
if length(peak_loc) > pulse_check_length
    pulse_start_time = find_pulse_start(peak_loc,pulse_check_length ,max_peak_distance);
    pulse_end_time = find_pulse_end(peak_loc,pulse_check_length ,max_peak_distance);
    if pulse_start_time < 0.1*run_time
        pulse_start_time = 0;
    end
else
    pulse_start_time = 0;
    pulse_end_time = 0;
end


if pulse_start_time ~= 0 && pulse_end_time ~= 0
    disp(['[RESULT] Pulse start at ',num2str(pulse_start_time),...
        's, Pulse end at ',num2str(pulse_end_time),'s'])
    pulse_rate = pulse_rate_cal(peak_loc,pulse_start_time,pulse_end_time);
    disp(['[RESULT] Pulse Rate: ',num2str(pulse_rate),'bpm']);
    systolic_index = round(pulse_start_time*fs);
    diastolic_index = round(pulse_end_time*fs);
    mean_period = round(0.01*fs);
    systolic_pressure = mean(pressure_data(systolic_index-mean_period:systolic_index+mean_period));
    systolic_pressure = roundn(systolic_pressure,-1);
    diastolic_pressure = mean(pressure_data(diastolic_index-mean_period:diastolic_index+mean_period));
    diastolic_pressure = roundn(diastolic_pressure,-1);
    fprintf('[RESULT] SBP: %d mmHg, DBP: %d mmHg\n',...
        systolic_pressure,diastolic_pressure);
else
    disp('[INFO] Unable to find Blood Pressure');
end

%% Plot Filtered Data
expantion_time = 1.2; %unit:s
plot_KS = Sample1((pulse_start_time-expantion_time)*fs:(pulse_end_time+expantion_time)*fs);
amp_rate = 1/max(abs(plot_KS));
plot_KS = amp_rate * plot_KS;
pressure_plot = pressure_data(1:length(Sample1));
plot_P = pressure_plot((pulse_start_time-expantion_time)*fs:(pulse_end_time+expantion_time)*fs);
plot_t = 1:length(plot_KS);
plot_t = plot_t/fs;
filtered_plot = figure(3);
set(filtered_plot,'Position',[50,50,1200,950])
subplot(2,1,1);
plot(plot_t,plot_KS,'b');
set(gca,'FontSize',xytick_font_size, 'box', 'on', 'linewidth',2);
xlim([0 round(plot_t(end),1)]);
xticks([0 round(plot_t(end),1)])
xticklabels(["",""]);
ylim([-1 1]);
yticks([-1 0 1]);
title("Korotkoff Sound","FontSize",title_font_size);
ylabel("Amplitude","FontSize",label_font_size)
hold on;
subplot(2,1,2);
plot(plot_t,plot_P,'r',"LineWidth",2);
set(gca,'FontSize',xytick_font_size, 'box', 'on', 'linewidth',2);
xlim([0 round(plot_t(end),1)]);
xlabel("Time (s)","FontSize",label_font_size);
ylabel("Pressure (mmhg)","FontSize",label_font_size);
title("Pressure","FontSize",title_font_size);
ylim([80 150]);
hold on;
txt = ['SBP ' num2str(systolic_pressure) ' mmHg'];
text(plot_t(end)*0.55,135,txt,"FontSize",label_font_size)
txt = ['DBP ' num2str(diastolic_pressure) ' mmHg'];
text(plot_t(end)*0.55,125,txt,"FontSize",label_font_size)
txt = ['Pulse ' num2str(pulse_rate) ' bpm'];
text(plot_t(end)*0.55,115,txt,"FontSize",label_font_size)
figure_name = file_name +"_ArmGirth"+ num2str(arm_girth)+...
    "_"+num2str(pulse_start_time)+"-"+num2str(pulse_end_time)+"s";
if figure_save_flag == 1
    saveas(filtered_plot, figure_name+".png");
end

% combined plot
result_combined_plot = figure(4);
set(result_combined_plot,'Position',[50,50,1200,950])
colororder(["b", "r"])
yyaxis left;
plot(plot_t,plot_KS,'b');
set(gca,'FontSize',xytick_font_size);
ylim([-1.5 1.5]);
yticks([-1 0 1]);
ylabel("Amplitude","FontSize",label_font_size)
yyaxis right;
plot(plot_t,plot_P,"r","LineWidth",2);
set(gca,'FontSize',xytick_font_size,'linewidth',2);
xlim([0 round(plot_t(end),1)]);
xlabel("Time (s)","FontSize",label_font_size);
ylabel("Pressure (mmhg)","FontSize",label_font_size);
hold on;
if figure_save_flag == 1
    saveas(result_combined_plot, "combined_" + figure_name+".png");
end




%% Save the audio signal as a WAV file
if audio_write_flag == 1
    disp("[PROCESS] Writing filtered signal audio");
    amp_rate = 1/max(abs(plot_KS));
    plot_KS = amp_rate*plot_KS;
    wav_file = file_name + ".wav";
    audiowrite(wav_file, plot_KS, fs,"BitsPerSample",24);
end
%% Functions
function filtered = notch(low_limit, high_limit,fs,freq_array,freq_amp,raw_data,varargin)
    f = freq_array;
    P1 = freq_amp;
    Sample1 = raw_data;
    peak_range = find(f>low_limit & f<high_limit);
    [~, peak_index] = max(P1(peak_range));
    peak_index = peak_index+peak_range(1)-1;
    peak_freq = round(f(peak_index),2);
    % Notch filter 
    notch_freq = peak_freq; 
    if nargin == 6
        bw = 5;
    else
        bw = varargin{1}; %  Bandwidth for the notch filter
    end
    [bn, an] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
    Sample1 = filter(bn, an, Sample1); 
    filtered = Sample1;
    % disp(['noise peak:',num2str(peak_freq)])
end


function filtered = low_pass(raw_data, fs, freq)
    %Low-pass filter
    Fc_lp = freq; 
    [b_lp, a_lp] = butter(4, Fc_lp/(fs/2)); 
    filtered = filter(b_lp, a_lp, raw_data);
end

function filtered = high_pass(raw_data, fs, freq)
    % High-pass filter
    Fc_hp = freq;
    [b_hp, a_hp] = butter(4, Fc_hp/(fs/2), 'high');
    filtered = filter(b_hp, a_hp,raw_data);
end

function filtered = smooth(windowSize,Sample1)
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    filtered = filter(b,a,Sample1);
end

function filtered = normalize(Raw_data)
    Sample1 = Raw_data;
    norm_factor = 1/max(abs(Sample1));
    filtered = Sample1*norm_factor;
end

function start_time = find_pulse_start(peak_loc,check_length,pulse_interval)
% peak_loc: time stamp of the peak(peaks' location)
% check_length: number of consiquent data to check
% pulse_interval: the interval smaller than this can be consider as a pulse
    peak_arrary_interval = zeros(1,length(peak_loc)-1);
    for i=1:length(peak_loc)-1
        peak_arrary_interval(i) = peak_loc(i+1)-peak_loc(i);
    end
    pulse_start_flag = 0;
    for i=1:length(peak_arrary_interval)-check_length
        if peak_arrary_interval(i)<pulse_interval
            pulse_start_flag = 1;
            for j=i:i+check_length
                if peak_arrary_interval(j)>pulse_interval
                    pulse_start_flag = 0;
                end
            end
        end
        if pulse_start_flag == 1
            pulse_start_point = peak_loc(i);
            break
        else
            pulse_start_point = 0;
            continue
        end
    end
    start_time = pulse_start_point;
end

function end_time = find_pulse_end(peak_loc,check_length,pulse_interval)
% peak_loc: time stamp of the peak(peaks' location)
% check_length: number of consiquent data to check
% pulse_interval: the interval smaller than this can be consider as a pulse
    peak_arrary_interval = zeros(1,length(peak_loc)-1);
    for i=1:length(peak_loc)-1
        peak_arrary_interval(i) = peak_loc(length(peak_loc)-i+1)-peak_loc(length(peak_loc)-i);
    end
    pulse_start_flag = 0;
    for i=1:length(peak_arrary_interval)-check_length
        if peak_arrary_interval(i)<pulse_interval
            pulse_start_flag = 1;
            for j=i:i+check_length
                if peak_arrary_interval(j)>pulse_interval
                    pulse_start_flag = 0;
                end
            end
        end
        if pulse_start_flag == 1
            pulse_start_point = peak_loc(length(peak_loc)-i+1);
            break
        else
            pulse_start_point = 0;
            continue
        end
    end
    end_time = pulse_start_point;
end


function csv_to_mat(file_name)
    if isfile(file_name + ".mat") == 0
        disp("[INFO] Converting csv file to mat file");
        csv_file = file_name;
        csv_name = csv_file + ".csv";
        csv_read = readtable(csv_name);
        voltage_array = csv_read{:,2}; % KS data
        pressure_array = csv_read{:,3}; % pressure data
        data_time = csv_read{:,1};
        fs_data = round(length(data_time)/data_time(end));
        save(csv_file,"voltage_array","pressure_array","fs_data");
        pause(1);
    else
        disp("[INFO] mat file already exsist, no need for conversion");
    end
end


function p_rate = pulse_rate_cal(peak_loc,pulse_start_time,pulse_end_time)
    pulse_rate_start = find(peak_loc>=pulse_start_time);
    pulse_rate_start = pulse_rate_start(1);
    pulse_rate_end = find(peak_loc>=pulse_end_time-5);
    pulse_rate_end = pulse_rate_end(1);
    pulse_total_amount = pulse_rate_end-pulse_rate_start;
    pulse_rate_start_time = peak_loc(pulse_rate_start);
    pulse_rate_end_time = peak_loc(pulse_rate_end-1);
    p_rate = pulse_total_amount/(pulse_rate_end_time-pulse_rate_start_time)*60;
    p_rate = round(p_rate,1);
end
