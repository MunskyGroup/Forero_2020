function [POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_raw_intensities(percent)
[POL2_I,SER5_I,MRNA_I] = get_raw_intensities();
[POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_simulated_intensities(percent,POL2_I,SER5_I,MRNA_I);

% for i = 1:length(POL2_I(1,:))
%     
%     top5_POL2 = quantile(POL2_I(:,i),percent);
% 
%     top5_SER5 = quantile(SER5_I(:,i),percent);
% 
%     top5_MRNA = quantile(MRNA_I(:,i),percent);
% 
% 
%     POL2_I_Norm(:,i) = (POL2_I(:,i))./(top5_POL2);
%     SER5_I_Norm(:,i) = (SER5_I(:,i))./(top5_SER5);
%     MRNA_I_Norm(:,i) = (MRNA_I(:,i))./(top5_MRNA);
% 
%     POL2_I_Norm = min(POL2_I_Norm,1.5);
%     SER5_I_Norm = min(SER5_I_Norm,1.5);
%     MRNA_I_Norm = min(MRNA_I_Norm,1.5);
%  
%     
% end
% % for i=1:length(POL2_I(1,:))
% %     POL2_I_Norm(:,i) = (POL2_I_Norm(:,i)-mean(POL2_I_Norm(:,i)))/std(POL2_I_Norm(:,i));
% %     SER5_I_Norm(:,i) = (SER5_I_Norm(:,i)-mean(SER5_I_Norm(:,i)))/std(SER5_I_Norm(:,i));
% %     MRNA_I_Norm(:,i) = (MRNA_I_Norm(:,i)-mean(MRNA_I_Norm(:,i)))/std(MRNA_I_Norm(:,i));
% % end

plottype = false;

if plottype == true
    figure(20);clf;
    subplot(3,3,1)
    histogram(POL2_I_Norm(:)); hold on;
    subplot(3,3,5)
    histogram(SER5_I_Norm(:)); hold on;
    subplot(3,3,9)
    histogram(MRNA_I_Norm(:)); hold on;
    subplot(3,3,2)
    scatter(SER5_I_Norm(:),POL2_I_Norm(:)); hold on;
    subplot(3,3,3)
    scatter(MRNA_I_Norm(:),POL2_I_Norm(:)); hold on;
    subplot(3,3,6)
    scatter(MRNA_I_Norm(:),SER5_I_Norm(:)); hold on;
    subplot(3,3,4)
    ksdensity([SER5_I_Norm(:),POL2_I_Norm(:)],'PlotFcn','contour');
    set(gca,'xlim',[-2 2],'ylim',[-2 2])
    subplot(3,3,7)
    ksdensity([MRNA_I_Norm(:),POL2_I_Norm(:)],'PlotFcn','contour');
    set(gca,'xlim',[-2 2],'ylim',[-2 2])
    subplot(3,3,8)
    ksdensity([MRNA_I_Norm(:),SER5_I_Norm(:)],'PlotFcn','contour');
    set(gca,'xlim',[-2 2],'ylim',[-2 2])
    
end

end
