function [mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,G1_rnap,G1_ser5,G1_mrna] = load_corr_data_figure(rezero_acc, acc_type)

[POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_raw_intensities(.95);


[data]=xlsread('RawCorrelationsData.xlsx','mRNA_RawAC');
lnt = floor((length(data(:,1))+1)/2);
mrna.tv = data(lnt:end,1);


[data]=xlsread('RawCorrelationsData.xlsx','Ser5ph_RawAC');
lnt = floor((length(data(:,1))+1)/2);
ser5.tv = data(lnt:end,1);


[data]=xlsread('RawCorrelationsData.xlsx','RNAP2_RawAC');
lnt = floor((length(data(:,1))+1)/2);
rnap.tv = data(lnt:end,1);


[data]=xlsread('RawCorrelationsData.xlsx','mRNA-RNAP2_RawCC');
mrna_rnap.tv = data(lnt-10:lnt+10,1);


[data]=xlsread('RawCorrelationsData.xlsx','mRNA-Ser5ph_RawCC');
mrna_ser5.tv = data(lnt-10:lnt+10,1);


[data]=xlsread('RawCorrelationsData.xlsx','Ser5ph-RNAP2_RawCC');
ser5_rnap.tv = data(lnt-10:lnt+10,1);

Nt = 200;
pol2acc100 = [];
ser5acc100 = [];
mrnaacc100 = [];

mrna_rnap100 = [];
ser5_rnap100 = [];
mrna_ser5100 = [];
G1s_rnap = zeros(1,20);
G1s_ser5 = zeros(1,20);
G1s_mrna = zeros(1,20);
for j = 1:20
    V = [POL2_I_Norm(:,j),SER5_I_Norm(:,j),MRNA_I_Norm(:,j)];
    V = (V - repmat(mean(V),200,1));%./std(V,0,1);
    [C,lags] = xcorr(V);
    C(round(length(C)/2):end,:) = C(round(length(C)/2):end,:)./(Nt:-1:1)'; %Xcorr scale fix
    C(1:round(length(C)/2)-1,:) = C(1:round(length(C)/2)-1,:)./(1:1:Nt-1)'; 

        i = 1;
        mod_rnap_100 = C(:,(i-1)*3+i);
        mod_rnap_100 = mod_rnap_100(((length(mod_rnap_100)+1)/2):((length(mod_rnap_100)+1)/2)+31);
        G1s_rnap(j) = mod_rnap_100(2);
        
        i = 2;
        mod_ser5_100 = C(:,(i-1)*3+i);
        mod_ser5_100 = mod_ser5_100(((length(mod_ser5_100)+1)/2):((length(mod_ser5_100)+1)/2)+31);
        G1s_ser5(j) = mod_ser5_100(2);

        i = 3;
        mod_mrna_100 = C(:,(i-1)*3+i);
        mod_mrna_100 = mod_mrna_100(((length(mod_mrna_100)+1)/2):((length(mod_mrna_100)+1)/2)+31);
        G1s_mrna(j) = mod_mrna_100(2);
        
        [G1s_rnap(j),G1s_ser5(j),G1s_mrna(j)] = get_g0(mod_rnap_100,mod_ser5_100,mod_mrna_100,acc_type);
        
end
    


    
G1_rnap = mean(G1s_rnap);
G1_ser5 = mean(G1s_ser5);
G1_mrna = mean(G1s_mrna);

G1s_rnap = mean(G1s_rnap)*ones(size(G1s_rnap));
G1s_ser5 = mean(G1s_ser5)*ones(size(G1s_ser5));
G1s_mrna = mean(G1s_mrna)*ones(size(G1s_mrna));




for j = 1:20
    V = [POL2_I_Norm(:,j),SER5_I_Norm(:,j),MRNA_I_Norm(:,j)];
    V = (V - repmat(mean(V),200,1));
    
    [C,lags] = xcorr(V);
    C(round(length(C)/2):end,:) = C(round(length(C)/2):end,:)./(Nt:-1:1)'; %Xcorr scale fix
    C(1:round(length(C)/2)-1,:) = C(1:round(length(C)/2)-1,:)./(1:1:Nt-1)'; 

        i = 1;
        mod_rnap_100 = C(:,(i-1)*3+i);
        mod_rnap_100 = mod_rnap_100(((length(mod_rnap_100)+1)/2):((length(mod_rnap_100)+1)/2)+31);
         mod_rnap_100 = mod_rnap_100./G1_rnap;
        %mod_rnap_100 = mod_rnap_100./G1s_rnap(j);
        
       % mod_rnap_100 = mod_rnap_100
        
        
        if j == 1
            pol2acc100 = mod_rnap_100;
            
        else
            pol2acc100 = horzcat(pol2acc100,mod_rnap_100);
        end
        
        i = 2;
        mod_ser5_100 = C(:,(i-1)*3+i);
        mod_ser5_100 = mod_ser5_100(((length(mod_ser5_100)+1)/2):((length(mod_ser5_100)+1)/2)+31);
        mod_ser5_100 = mod_ser5_100./G1_ser5;
        %mod_ser5_100 = mod_ser5_100./G1s_ser5(j);
        
        
       % mod_ser5_100 = mod_ser5_100/mod_ser5_100(2);

        if j == 1
            ser5acc100 = mod_ser5_100;
            
        else
            ser5acc100 = horzcat(ser5acc100,mod_ser5_100);
        end
        
        i = 3;
        mod_mrna_100 = C(:,(i-1)*3+i);
        mod_mrna_100 = mod_mrna_100(((length(mod_mrna_100)+1)/2):((length(mod_mrna_100)+1)/2)+31);
        mod_mrna_100 = mod_mrna_100./G1_mrna;
        %mod_mrna_100 = mod_mrna_100./G1s_mrna(j);
        %mod_mrna_100 = mod_mrna_100/mod_mrna_100(2);
        if j == 1
            mrnaacc100 = mod_mrna_100;
            
        else
            mrnaacc100 = horzcat(mrnaacc100,mod_mrna_100);
        end

        [~,I] = max(abs(C)); %lag times
        lag_offsets = lags(I);

        mod_mrna_rnap_100 = C(:,7);
        mod_mrna_rnap_100 = mod_mrna_rnap_100(((length(mod_mrna_rnap_100)+1)/2-10):((length(mod_mrna_rnap_100)+1)/2)+10);
        %mod_mrna_rnap_100 = mod_mrna_rnap_100./(sqrt(G1_mrna)*sqrt(G1_rnap));
        %mod_mrna_rnap_100 = mod_mrna_rnap_100/(sqrt(G1s_mrna(j))*sqrt(G1s_rnap(j)));
       
        %mod_mrna_rnap_100 = mod_mrna_rnap_100/mod_mrna_rnap_100(8);
       
        if j == 1
            mrna_rnap100 = mod_mrna_rnap_100;
            
        else
            mrna_rnap100 = horzcat(mrna_rnap100,mod_mrna_rnap_100);
        end
        
        mod_ser5_rnap_100 = C(:,4);
        mod_ser5_rnap_100 = mod_ser5_rnap_100(((length(mod_ser5_rnap_100)+1)/2-10):((length(mod_ser5_rnap_100)+1)/2)+10);
        %mod_ser5_rnap_100 = mod_ser5_rnap_100./(sqrt(G1_ser5)*sqrt(G1_rnap));
        %mod_ser5_rnap_100 = mod_ser5_rnap_100/(sqrt(G1s_ser5(j))*sqrt(G1s_rnap(j)));
        %mod_ser5_rnap_100 = mod_ser5_rnap_100/mod_ser5_rnap_100(8);
        if j == 1
            ser5_rnap100 = mod_ser5_rnap_100;
            
        else
            ser5_rnap100 = horzcat(ser5_rnap100,mod_ser5_rnap_100);
        end
        mod_mrna_ser5_100 = C(:,8);
        mod_mrna_ser5_100 = mod_mrna_ser5_100(((length(mod_mrna_ser5_100)+1)/2-10):((length(mod_mrna_ser5_100)+1)/2)+10);
        
        %mod_mrna_ser5_100 = mod_mrna_ser5_100./(sqrt(G1_mrna)*sqrt(G1_ser5));
        %mod_mrna_ser5_100 = mod_mrna_ser5_100/(sqrt(G1s_mrna(j))*sqrt(G1s_ser5(j)));
        %mod_mrna_ser5_100 = mod_mrna_ser5_100/mod_mrna_ser5_100(8);
        if j == 1
            mrna_ser5100 = mod_mrna_ser5_100;
            
        else
            mrna_ser5100 = horzcat(mrna_ser5100,mod_mrna_ser5_100);
        end
end



 mrna.mn_ac = mean(mrnaacc100,2);
 ser5.mn_ac = mean(ser5acc100,2);
 rnap.mn_ac = mean(pol2acc100,2);
 

%[g_mr, g_ms, g_sr] =  get_g0_cc(mean(mrna_rnap100,2)',mean(mrna_ser5100,2)',mean(ser5_rnap100,2)',cc_type);

 mrna.sem_ac = std(mrnaacc100,0,2)./sqrt(20);
 ser5.sem_ac = std(ser5acc100,0,2)./sqrt(20);
 rnap.sem_ac = std(pol2acc100,0,2)./sqrt(20);
 

 
 mrna_rnap.mn_cc = mean(mrna_rnap100,2);
 mrna_ser5.mn_cc = mean(mrna_ser5100,2);
 ser5_rnap.mn_cc = mean(ser5_rnap100,2);
 
 if rezero_acc == true
     rezero_mrna = mean(mrna.mn_ac(end-9:end));
     rezero_ser5 = mean(ser5.mn_ac(end-9:end));
     rezero_rnap = mean(rnap.mn_ac(end-9:end));
     
%      mrnaacc100 = mrnaacc100-rezero_mrna;
%      ser5acc100 = ser5acc100-rezero_ser5;
%      pol2acc100 = pol2acc100-rezero_rnap;
     
%      mrna.mn_ac = mean(mrnaacc100,2);
%      ser5.mn_ac = mean(ser5acc100,2);
%      rnap.mn_ac = mean(pol2acc100,2);

     mrna_scl= mrna.mn_ac(1)/(mrna.mn_ac(1)-rezero_mrna);
     ser5_scl = ser5.mn_ac(1)/(ser5.mn_ac(1)-rezero_ser5);
     rnap_scl = rnap.mn_ac(1)/(rnap.mn_ac(1)-rezero_rnap);

     mrna.mn_ac = (mrna.mn_ac-rezero_mrna)*mrna_scl;
     ser5.mn_ac = (ser5.mn_ac-rezero_ser5)*ser5_scl;
     rnap.mn_ac = (rnap.mn_ac-rezero_rnap)*rnap_scl;
     
     mrna.sem_ac = mrna.sem_ac*mrna_scl;
     ser5.sem_ac = ser5.sem_ac*ser5_scl;
     rnap.sem_ac = rnap.sem_ac*rnap_scl;
          
     
     
 end
 
   
 g0 = mrna_rnap.mn_cc(11);
 mrna_rnap.mn_cc = mrna_rnap.mn_cc./g0;
 mrna_rnap.sem_cc = std(mrna_rnap100,0,2)./(g0*sqrt(20));
 
 
 g0 = mrna_ser5.mn_cc(11);
  mrna_ser5.mn_cc = mrna_ser5.mn_cc./g0;
 mrna_ser5.sem_cc = std(mrna_ser5100,0,2)./(g0*sqrt(20));
 
 
 g0 = ser5_rnap.mn_cc(11);
  ser5_rnap.mn_cc = ser5_rnap.mn_cc./g0;
 ser5_rnap.sem_cc = std(ser5_rnap100,0,2)./(g0*sqrt(20));
 
 
end
