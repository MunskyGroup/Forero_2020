function [mRNAFluc_Rave,Ser5phFluc_Rave, UnphFluc_Rave, vector_minima_mRNA,lmin_mRNA] = minimaLocation(UnphFluc,Ser5phFluc,mRNAFluc,detectionThreshold)
% 
% UnphFluc=xlsread('Data_with_Oscillations.xlsx','Unph');
% Ser5phFluc=xlsread('Data_with_Oscillations.xlsx','Ser5ph');
% mRNAFluc=xlsread('Data_with_Oscillations.xlsx','mRNA');


timeFluc=[0:199]';
%% Locating minimal values
for i=1:size(UnphFluc,2)
    % To reduce noise. Smoothing data with a window size of 3 seconds.
    UnphFluc_Rave(:,i)=movmean(UnphFluc(:,i),3);
    Ser5phFluc_Rave(:,i)=movmean(Ser5phFluc(:,i),3);
    mRNAFluc_Rave(:,i)=movmean(mRNAFluc(:,i),3);
    % Locating minimal values in the threshold "MinProminence"
    
    traj = mRNAFluc_Rave(:,i);
    pure_detect = (mRNAFluc_Rave(:,i) < 1-detectionThreshold);
    pure_detect2 = ones(200,1);
    pure_detect2(pure_detect) = traj(pure_detect);

    
    lmin_mRNA(:,i)= islocalmin(pure_detect2,'MinProminence', 1-detectionThreshold,'MinSeparation',15, 'FlatSelection', 'first');%find local min mRNA

    %lmin_mRNA(:,i)= islocalmin(mRNAFluc_Rave(:,i),'MinProminence', 1-detectionThreshold,'MinSeparation',15, 'FlatSelection', 'first');%find local min mRNA
    %lmin_Unph(:,i)= islocalmin(UnphFluc_Rave(:,i),'MinProminence', detectionThreshold,'MinSeparation',15);%find local min Unph
    %lmin_Ser5ph(:,i)= islocalmin(Ser5phFluc_Rave(:,i),'MinProminence', detectionThreshold,'MinSeparation',15);%find local min Unph
    % to avoid conflicts in the code, we are ignoring the first 10 seconds
    % and the las 10 seconds
    lmin_mRNA(1:10,i)=0; lmin_mRNA(190:200,i)=0;
    %lmin_Unph(1:10,i)=0; lmin_Unph(190:200,i)=0;
    %lmin_Ser5ph(1:10,i)=0; lmin_Ser5ph(190:200,i)=0;
    % Saving time points with the location of minimal values for mRNA and UnphFluc in a cell.
    vector_minima_mRNA{i} = timeFluc(lmin_mRNA(:,i))+1;
    %vector_minima_Unph{i} = timeFluc(lmin_Unph(:,i))+1;
    %vector_minima_Ser5ph{i} = timeFluc(lmin_Ser5ph(:,i))+1;
end


end