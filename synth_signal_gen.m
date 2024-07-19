% Date: 19/07/2024
% Author: Zoran Sverko, zoran.sverko@riteh.uniri.hr
% ---------------------------------
% ---------------------------------
% ---------------------------------
% Article:
% Generation and Analysis of Generation Parameters of Synthetic EEG Signals for Testing Dynamic Directed Brain Connectivity Estimation Methods
% Authors: Sverko Zoran; Rogelj Peter; Vlahinic Sasa
% ---------------------------------
% ---------------------------------
% ---------------------------------

clc;clear all;close all;
%% IMPORT
% Enter the path from your computer for data import (the data used is described in the article text).
export_names={'SO01_RO1.mat'};

i=1;
path_import='E:\Nastava riteh pula 23_24\Mobilnost Koper\GC\Kod\';
name_import_1=char(export_names(i))
path_all=strcat(path_import,name_import_1);
load(path_all)


%% connectivity matrix calculation - Granger causality (GC) 

order = 19; % regression order
nrEl = size(EEG.data,1); % number of electrodes
conG2=zeros(nrEl);
tic
for c1=1:nrEl
    for c2= c1+1:nrEl
        GC = GCmodel(EEG.data([c1 c2],:), order);
        GC = max(GC,[0 0]);
        conG2(c1,c2) = GC(1);
        conG2(c2,c1) = GC(2);
    end
end
toc % Elapsed time is 35.265318 seconds.

figure(22);
imagesc(conG2); colorbar; title(['Granger causality, order=' num2str(order) ', ' EEG.setname ]);
xlabel('influencing electrode');
ylabel('influenced electrode');
axis equal;axis tight;
xticks([1:nrEl]);xticklabels({EEG.chanlocs.labels}); xtickangle(90)
yticks([1:nrEl]);yticklabels({EEG.chanlocs.labels});


%% Logic of Thinking
% We took two electrodes that are connected
% From these electrodes, we took intervals that are not the same in time but are of the same length
% We will use electrodes Fc1 (3rd channel) and F1 (33rd channel)
% The sampling rate is Fs = 160 Hz
% We will take three intervals that are not connected, each lasting 14.4 seconds
% Each interval is 768/160 = 4.8 seconds long
ch1=3; %channel 3 - FC1
ch2=33; %channel 33 - F1

S1=EEG.data(ch1,:);
S2=EEG.data(ch2,:);
data=[S1; S2];

S1_int=S1(1,1:2304); % In the article, it is T' from FC1.
S2_int=S2(1,4609:6912); % In the article, it is S' from F1.
data_int=[S1_int; S2_int];


%% Calculate GC S1_int in S2_int, we expect it to be around 0
order=19;
GC_int=zeros(2);

GC=GCmodel(double(data_int),order)
disp("GC of unrelated parts of signals:");disp(GC);

%% We calculate the GC model for the selected electrodes to obtain parameters for signal reconstruction
data=[EEG.data(ch1,:);EEG.data(ch2,:)];
order2=19;
[GC, omega1, omega2, omega12, e1, e2, e12] = GCmodel(double(data),order2);
disp("GC of original signals:");disp(GC);

%omega - from the initially considered electrodes Fc1 and F1
omegaT=omega12(1,1:12);
omegaS=omega12(1,13:24);

% signals to reconstruct from:
T1 = data(1,1:2304);
S1 = data(2,4609:6912);
eT12 = e12(1,1+order2:2304+order2);% error in the interval from 20:2324, we discarded the first 19 samples

% weighting signal
K=ones(1,2304);
K(1,769:1536)=0;

% reconstruction
R = reconstructSig2([omegaT; omegaS], [T1; S1], K , eT12); % ); %






%% Calculating GC on precisely defined intervals 1:768, 769:1536, 1537:2304 - reference functional connectivity values RFC
synth_data_xx = R;

testData=[synth_data_xx; S2_int];
disp("section 1 - k=1:")
GC_INT1=GCmodel(double(testData(:,1:768)),order)

disp("section 2 - k=0:")
GC_INT2=GCmodel(double(testData(:,769:1536)),order)

disp("section 3 - k=1:")
GC_INT3=GCmodel(double(testData(:,1537:2304)),order)



%% Synthetic and source signals
figure('Name','Signals','Position',[255.4 343 793.6 420]);
plot(synth_data_xx,'k');

hold on;
plot(S2_int,'g');

title('Synthetic and source signals','FontSize',12);
ylabel('Amplitude','FontSize',14)
xlim([0, size(synth_data_xx,2)])
% ylim([-2,3]);
xticks([0 640 1280 1920 2304])
xticklabels({'0','4','8','12','14.4'})
yticks([-300 -150 0 150 300])
% yticklabels({'-300','-150','0','150','300'})
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',14)
% b = get(gca,'YTickLabel');
% set(gca,'YTickLabel',b,'fontsize',14)
% new_b=strrep(b(:),'.',',');
% set(gca,'YTickLabel',new_b)
xline(768,'k','LineWidth',2,'HandleVisibility','off')
xline(1536,'k','LineWidth',2,'HandleVisibility','off')
xlabel({'\fontsize{14}Time [s]'},'fontsize',14);
% xlabel({'\fontsize{12}Time';'\fontsize{12}[s]';'\fontsize{12}\bf(a)'},'fontsize',12);
a=strcat('synthetic-','\it{R}',"'")
b=strcat('','\it{S}',"'")
legend(a,b,'fontsize',12)
ax=gca;
ax.XAxis.FontSize=14;
ax.YAxis.FontSize=14;

% Create textbox
annotation('textbox',...
    [0.716913978494623 0.741904761904763 0.149013440860215 0.162425965411041],...
    'String',['Interval 3, K=1',sprintf('\n'),'\it{GC}=0.7723'],...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.438654569892473 0.728571428571429 0.15207123655914 0.177003079838901],...
    'String',['Interval 2, K=0',sprintf('\n'),'\it{GC}=0.0319'],...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.182505376344086 0.724761904761905 0.156204301075268 0.180398009950253],...
    'String',{'Interval 1, K=1','\it{GC}=0.7745'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% Dynamic GC analysis with constant window size (sliding window analysis) 

f_fil_n=fun_constant_window_v1_GC(testData,64,19);
f_fil_w=fun_constant_window_v1_GC(testData,320,19);

c=strcat('Dynamic {\it GC} ',{' '},a,'->',b);
figure('Name','Dynamic GC - constant window width','Position',[53.800000000000004 343 1403.2 420])
subplot(2,1,1)
plot(f_fil_w,'k');
hold on;
TFC_GC_values=ones(1,size(f_fil_w,2));

TFC_GC_values(1,1:768)=TFC_GC_values(1,1:768).*GC_INT1(1);% RFC value on 1st interval
TFC_GC_values(1,769:1536)=TFC_GC_values(1,769:1536).*GC_INT2(1);% RFC value on 2nd interval
TFC_GC_values(1,1537:end)=TFC_GC_values(1,1537:end).*GC_INT3(1);% RFC value on 3rd interval


title(c,'FontSize',18);
ylabel('{\it GC} value','FontSize',16)
xlim([0, size(synth_data_xx,2)])
% ylim([-2,3]);
xticks([0 640 1280 1920 2304])
xticklabels({'0','4','8','12','14.4'})
xline(768,'k','LineWidth',2,'HandleVisibility','off')
xline(1536,'k','LineWidth',2,'HandleVisibility','off')
plot(TFC_GC_values,'g','LineWidth',2)

% xlabel({'\fontsize{18}Time';'\fontsize{16}[seconds]';'\fontsize{20}\bf(a)'},'fontsize',16);
xlabel({'\fontsize{20}\bf(a)'},'fontsize',16);

legend('Estimated \it{GC}','\it{RFC GC}','fontsize',16)
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;



subplot(2,1,2)
plot(f_fil_n,'k');
hold on;
TFC_GC_values=ones(1,size(f_fil_n,2));

TFC_GC_values(1,1:768)=TFC_GC_values(1,1:768).*GC_INT1(1);
TFC_GC_values(1,769:1536)=TFC_GC_values(1,769:1536).*GC_INT2(1);
TFC_GC_values(1,1537:end)=TFC_GC_values(1,1537:end).*GC_INT3(1);


% title(c,'FontSize',18);
ylabel('{\it GC} value','FontSize',16)
xlim([0, size(synth_data_xx,2)])
% ylim([-2,3]);
xticks([0 640 1280 1920 2304])
xticklabels({'0','4','8','12','14.4'})
xline(768,'k','LineWidth',2,'HandleVisibility','off')
xline(1536,'k','LineWidth',2,'HandleVisibility','off')
plot(TFC_GC_values,'g','LineWidth',2)

% xlabel({'\fontsize{18}Time';'\fontsize{16}[seconds]';'\fontsize{20}\bf(b)'},'fontsize',16);
xlabel({'\fontsize{18}Time [s]';'\fontsize{20}\bf(b)'},'fontsize',16);

% legend('Estimated \it{GC}','\it{TFC GC}','fontsize',16)
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;
