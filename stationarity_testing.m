clc;clear all;close all;

%% IMPORT 

export_names={'SO01_RO1.mat'};


% Import preprocessing data 
i=1;
path_import='E:\Nastava riteh pula 23_24\Mobilnost Koper\GC\Kod\';
name_import_1=char(export_names(i))
path_all=strcat(path_import,name_import_1);
load(path_all)

ch1=3; 
ch2=33; 

S1=EEG.data(ch1,:);
S2=EEG.data(ch2,:);


% MATLAB code for testing stationarity using Augmented Dickey-Fuller (ADF) test
% Input data: two sets of signals, signal1 and signal2 (as arrays)

% Example signals (replace with actual data)
signal1=S1;
signal2=S2;

% ADF test for signal 1
[h1, pValue1] = adftest(signal1);
fprintf('Signal 1: Stationary=%d, p-value=%.4f\n', h1, pValue1);

% ADF test for signal 2
[h2, pValue2] = adftest(signal2);
fprintf('Signal 2: Stationary=%d, p-value=%.4f\n', h2, pValue2);

% Interpretation
if h1 == 1
    disp('Signal 1 is stationary.');
else
    disp('Signal 1 is not stationary.');
end

if h2 == 1
    disp('Signal 2 is stationary.');
else
    disp('Signal 2 is not stationary.');
end