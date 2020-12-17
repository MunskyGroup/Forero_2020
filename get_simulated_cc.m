function [mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,G1_rnap,G1_ser5,G1_mrna] = get_simulated_cc(t, POL2_I_Norm,SER5_I_Norm,MRNA_I_Norm, rezero_acc)
%[POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_raw_intensities(.95);



lnt = floor((length(POL2_I_Norm(:,1))+1)/2);

Nt = length(t);
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
    [C,~] = xcorr(V);
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


    [G1s_rnap(j),G1s_ser5(j),G1s_mrna(j)] = get_g0(mod_rnap_100,mod_ser5_100,mod_mrna_100,'G0_intp');

    
    
G1_rnap = mean(G1s_rnap);
G1_ser5 = mean(G1s_ser5);
G1_mrna = mean(G1s_mrna);

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
  
%  figure;
%  plot(mrna_rnap100)
%  hold on
%  plot(mean(mrna_rnap100,2),'o')


 mrna.mn_ac = mean(mrnaacc100,2);
 ser5.mn_ac = mean(ser5acc100,2);
 rnap.mn_ac = mean(pol2acc100,2);
 
 
 
 mrna.sem_ac = std(mrnaacc100,0,2)./sqrt(20);
 ser5.sem_ac = std(ser5acc100,0,2)./sqrt(20);
 rnap.sem_ac = std(pol2acc100,0,2)./sqrt(20);
 
 
 mrna_rnap.mn_cc = mean(mrna_rnap100,2);
 mrna_ser5.mn_cc = mean(mrna_ser5100,2);
 ser5_rnap.mn_cc = mean(ser5_rnap100,2);
 
 
 if rezero_acc == true
     rezero_mrna = mean(mrna.mn_ac(end-11:end-6));
     rezero_ser5 = mean(ser5.mn_ac(19:25));
     rezero_rnap = mean(rnap.mn_ac(end-5:end));
     
     mrnaacc100 = (mrnaacc100-rezero_mrna)./(1-rezero_mrna);
     ser5acc100 = (ser5acc100-rezero_ser5)./(1-rezero_ser5);
     mrna.mn_ac = mean(mrnaacc100,2);
     ser5.mn_ac = mean(ser5acc100,2);

     pol2acc100 = (pol2acc100-rezero_rnap)./(1-rezero_rnap);
     rnap.mn_ac = mean(pol2acc100,2);


     mrna.sem_ac = std(mrnaacc100,0,2)./sqrt(20);   
     ser5.sem_ac = std(ser5acc100,0,2)./sqrt(20);
     rnap.sem_ac = std(pol2acc100,0,2)./sqrt(20);
     
     
     rezero_mrna_rnap = mean(mrna_rnap.mn_cc([1:3, end-3:end]));
     rezero_mrna_ser5 = mean(mrna_ser5.mn_cc([1:3, end-3:end]));
     rezero_ser5_rnap = mean(ser5_rnap.mn_cc([1:3, end-3:end]));
     
     
     mrna_rnap100 = (mrna_rnap100-rezero_mrna_rnap)./(1-rezero_mrna_rnap);
     mrna_ser5100 = (mrna_ser5100-rezero_mrna_ser5)./(1-rezero_mrna_ser5);
     ser5_rnap100 = (ser5_rnap100-rezero_ser5_rnap)./(1-rezero_ser5_rnap);
     
     
     mrna_rnap.mn_cc = mean(mrna_rnap100,2);
     mrna_ser5.mn_cc = mean(mrna_ser5100,2);
     ser5_rnap.mn_cc = mean(ser5_rnap100,2);
     
     mrna_rnap.sem_ac = std(mrna_rnap100,0,2)./sqrt(20);   
     mrna_ser5.sem_ac = std(mrna_ser5100,0,2)./sqrt(20);
     ser5_rnap.sem_ac = std(ser5_rnap100,0,2)./sqrt(20);
          
     
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




    

    