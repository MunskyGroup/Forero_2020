function [POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_simulated_intensities(percent,POL2_I,SER5_I,MRNA_I)
sz = size(POL2_I);
POL2_I_Norm = zeros(sz);
SER5_I_Norm  = zeros(sz);
MRNA_I_Norm = zeros(sz);


for i = 1:length(POL2_I(1,:))
    
    top5_POL2 = quantile(POL2_I(:,i),percent);
    top5_SER5 = quantile(SER5_I(:,i),percent);
    top5_MRNA = quantile(MRNA_I(:,i),percent);

    POL2_I_Norm(:,i) = (POL2_I(:,i))./(top5_POL2);
    SER5_I_Norm(:,i) = (SER5_I(:,i))./(top5_SER5);
    MRNA_I_Norm(:,i) = (MRNA_I(:,i))./(top5_MRNA);
    
    
end



POL2_I_Norm = min(POL2_I_Norm,1.5);
SER5_I_Norm = min(SER5_I_Norm,1.5);
MRNA_I_Norm = min(MRNA_I_Norm,1.5);

% for i=1:length(POL2_I(1,:))
%     POL2_I_Norm(:,i) = (POL2_I_Norm(:,i)-mean(POL2_I_Norm(:,i)))/std(POL2_I_Norm(:,i));
%     SER5_I_Norm(:,i) = (SER5_I_Norm(:,i)-mean(SER5_I_Norm(:,i)))/std(SER5_I_Norm(:,i));
%     MRNA_I_Norm(:,i) = (MRNA_I_Norm(:,i)-mean(MRNA_I_Norm(:,i)))/std(MRNA_I_Norm(:,i));
% end

end




