

function [POL2_I,SER5_I,MRNA_I] = get_raw_intensities()
[POL2_I] = xlsread('Raw_Intensities Analysis-BLC_W_BG_WO_RunAve_May18','Sheet1','C2:C4001');
[SER5_I] = xlsread('Raw_Intensities Analysis-BLC_W_BG_WO_RunAve_May18','Sheet1','D2:D4001');
[MRNA_I] = xlsread('Raw_Intensities Analysis-BLC_W_BG_WO_RunAve_May18','Sheet1','E2:E4001');

POL2_I = reshape(POL2_I,200,20);
SER5_I = reshape(SER5_I,200,20);
MRNA_I = reshape(MRNA_I,200,20);



end




function [] = plot_intensity_histograms()
[POL2_I,SER5_I,MRNA_I] = get_raw_intensities();
figure;
clf;
histogram(POL2_I,30);
hold on;
histogram(SER5_I, 30);
histogram(MRNA_I,30);

end
