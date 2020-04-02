function [C,lags,t_array,POL2,SER5,TS,POL2i,SER5i,TSi] = get_stoch_corrs_Nov25(Xfree,~,delta_t,model_options,traj_length)


trypt_scale = model_options.trypt_scale;
trypt_mu = model_options.trypt_mu;
trypt_var = model_options.trypt_var;

% Initialize default parameters at zero.
parameters = zeros(10,1);
% X(1:10,1) = 0;

% Set relevant parameters as given by input.
parameters(model_options.freepars) = Xfree;


kin = parameters(1);
kout = parameters(2);
kinit = parameters(3);
kabort = parameters(4);
kesc = parameters(5);
kdphos = parameters(6);
kelong = parameters(7);
kproc = parameters(10);

if model_options.photobleach == 1
    blch_pol2 = parameters(15);
    blch_ser5 = parameters(16);
    blch_ts = parameters(17);
else
blch_pol2 = 0.0;
blch_ser5 = 0.0;
blch_ts = 0.0;
end

% Fixed geometry of the gene.
LMS2_start = 0.7;
LMS2_end = 2.9+LMS2_start;
Lgene_end = 1.6+LMS2_end;

% Approximate time to reach steady state.
burnin = (1/kin+min(1/kinit,1/kout)+min(1/kesc+Lgene_end/kelong))*10;

% Simulation time to consider after burnin
sim_time = traj_length;

% Add bursting parameters into model.
if model_options.burst==1
    kon = parameters(8);
    koff = parameters(9);
    % Adjust SS time to encompass many switches.
    burnin = burnin+10*(1/kon+1/koff);
end

% Set number of polymerases to simulate.
Nsim = max(10000,ceil((burnin + sim_time)*kin*1.1));

Nsim = min(6000000, Nsim);

% Compute times between polymerases enter 'bubble' events.

IN_delt_times = -log(rand(Nsim,1))/kin;
% Convert to times that polymerases enter 'bubble' events.
times_IN = cumsum(IN_delt_times);  
% Truncate Poly entry events to just those that occur before end of
% simulation time.
Nsim = max(1000,sum(times_IN<(sim_time+burnin)));



if model_options.trypt(1) == 1   %BLOCKING K_IN
    %times_IN = times_IN(times_IN<model_options.trypt_time);
    times_IN = times_IN(times_IN<model_options.trypt_time+trypt_mu+sqrt(trypt_var)*randn());
    Nsim = sum(times_IN<(sim_time+burnin));
end

% if ismember(model_options.trypt, [7,2,5,4])    %BLOCKING K_IN
%     %times_IN = times_IN(times_IN<model_options.trypt_time);
%     times_IN = times_IN(times_IN<model_options.trypt_time+10+sqrt(10)*randn());
%     Nsim = sum(times_IN<(sim_time+burnin));
% end

% if ismember(model_options.trypt, [7,3,6,5])    %BLOCKING K_IN
%     %times_IN = times_IN(times_IN<model_options.trypt_time);
%     times_IN = times_IN(times_IN<model_options.trypt_time+10+sqrt(10)*randn());
%     Nsim = sum(times_IN<(sim_time+burnin));
% end


t_on_off = [];
% Truncate Poly entry events to just those that occur before end of
% simulation time.
times_IN = times_IN(1:Nsim);

% Adjust model for bursting kinetics.
if model_options.burst==1
    min_number_of_bursts = 10;
    
    n_bursts = ceil(max(min_number_of_bursts,1.1*max(times_IN)/(1/kon+1/koff)));
    t_turn_on = -log(rand(n_bursts,1))/kon;
    t_turn_off = -log(rand(n_bursts,1))/koff;

    
    
    % Generate times of burst 'on' and burst 'off'
    t_on_off(1:2:n_bursts*2-1) = t_turn_on;
    t_on_off(2:2:n_bursts*2) = t_turn_off;
    t_on_off = cumsum(t_on_off);
    t_on_off_table = reshape(t_on_off,2,n_bursts)';
    
    if model_options.trypt(5) ==1  %blocking k_on
        inhib_t = model_options.trypt_time + trypt_mu + sqrt(trypt_var)*randn();
       t_on_off_table = t_on_off_table(( t_on_off_table(:,1) < inhib_t),:);
        
        
    end
    

    
    % Find indices of poly 2 entries that occur during 'off' states
    J = false(size(times_IN));
    for i = 1:length(t_on_off_table)
        J(times_IN>=t_on_off_table(i,1)&times_IN<t_on_off_table(i,2))=1;
    end
    
    % Truncate to only those that enter when in 'on' state.
    times_IN = times_IN(J);
    Nsim = size(times_IN,1);
    
    
    
    
    
end

model_options.trypt_scale

% % Potential time for POL2 to leave the 'bubble'.
times_OUT = times_IN - log(rand(Nsim,1))/(kout); 

ag = 1;
tmax = max(times_IN);
all_times_IN =[];
all_times_OUT =[];
all_times_ser5  =[];
all_times_ser5_abort =[];
all_times_ser5_dephos =[];
all_times_terminate =[];
all_times_escape =[];
all_times_entering_proc = [];

new_in = times_IN;
new_outs = times_OUT;

t_inhibit = model_options.trypt_time+trypt_mu+sqrt(trypt_var)*randn();

while ~isempty(new_in)

    all_times_IN =[all_times_IN;new_in];
    [times_OUT,times_ser5,times_ser5_abort,times_ser5_dephos,times_terminate,times_escape,new_in,new_outs,escapes,time_enter_proc] = ...
        gen_new_pol2s_new(new_in,new_outs,parameters,model_options,t_inhibit);
    
    

    all_times_OUT =[all_times_OUT;times_OUT];
    all_times_ser5  =[all_times_ser5;times_ser5];
    all_times_ser5_abort =[all_times_ser5_abort;times_ser5_abort];
    all_times_ser5_dephos =[all_times_ser5_dephos;times_ser5_dephos];
    all_times_terminate =[all_times_terminate;times_terminate];
    all_times_escape =[all_times_escape;times_escape];
    all_times_entering_proc = [all_times_entering_proc; time_enter_proc];
end

times_IN =all_times_IN;
times_OUT =all_times_OUT;
times_ser5  =all_times_ser5;
times_ser5_abort =all_times_ser5_abort;
times_ser5_dephos =all_times_ser5_dephos;
times_terminate =all_times_terminate;
times_escape =all_times_escape;

times_entering_proc = all_times_entering_proc;

Nsim = length(times_ser5);







% if model_options.trypt(7) == 1
%     times_processing_inhibited = - log(rand(Nsim,1))/(kproc*trypt_scale);
%     affected_procs = (times_entering_proc >= t_inhibit);
%     times_OUT(affected_procs) = times_escape(affected_procs) + Lgene_end/kelong + times_processing_inhibited(affected_procs);
%     times_terminate(affected_procs) = times_escape(affected_procs) + Lgene_end/kelong + times_processing_inhibited(affected_procs);
%     
% end


% % 
% % % Potential time for SER5 Phosphorylation
% % times_ser5 = times_IN - log(rand(Nsim,1))/(kinit); 
% 
% % if model_options.trypt(2) ==1   %BLOCKING K_INIT
% %     %times_IN = times_IN(times_IN<model_options.trypt_time);
% %     ser5inhibit_t = model_options.trypt_time+trypt_scale+sqrt(trypt_scale)*randn();
% %     times_ser5(times_ser5>ser5inhibit_t) = inf;
% %   
% % end
% 
% % tappend = [1];
% % cnt = 1
% % times_escape = [];
% 
% % ser5_phos = (times_ser5<=times_OUT );
% 
% % % these POL2 do not get phosphorylated and simply go away.
% % times_ser5(~ser5_phos) = inf;
% % % active_ser5 = times_ser5(ser5_phos);
% % times_ser5_abort = times_ser5;
% % times_escape = inf*ones(size(times_ser5));
% 
%tmax = max(times_IN);
% while length(tappend > 0)
%     length(tappend);
%     cnt
%     if cnt ==1
%         
%         r = rand(size(times_ser5)); 
%      
%         returns_2_pol2 = r<(kabort/(kabort+kesc))  & times_ser5 ~= inf;  %ones that abort back to POL2
%         escapes        = r>=(kabort/(kabort+kesc)) & times_ser5 ~= inf;  %ones that abort back to POL2
% %         indexes = find(cond);
% %         reverse_indexes = find(~cond);
%         
%         tappend = times_ser5(returns_2_pol2) ;
%         %tappend = tappend(tappend ~= inf);
%         
%         dt = -log(rand(length(returns_2_pol2),1))/(kabort+kesc);
%         
%         
% %         dt = -log(rand(length(tappend),1))/(kabort+kesc);
%         tappend = tappend + dt(returns_2_pol2);
% 
%         new_outs = times_OUT(returns_2_pol2);
%         times_OUT(returns_2_pol2) = tappend;
%         
%         %ones that continue to escape
% %         dt_e = -log(rand(sum(escapes),1))/(kabort+kesc);
%        
%         times_escape(escapes) = times_ser5(escapes) + dt(escapes);
%         %times_escape = [times_escape; times_ser5(~cond) + dt_e];
%         
%         times_ser5_abort(returns_2_pol2) = tappend;
%         
%     else
%         
%         r = rand(size(tappend));  %only active ser5ph pols's
%         
%         returns_2_pol2 = r<(kabort/(kabort+kesc)) ;  %ones that abort back to POL2
%         escapes = r<(kabort/(kabort+kesc)) ;  %ones that abort back to POL2
%         reverse_indexes = indexes(~returns_2_pol2);
%         indexes = indexes(returns_2_pol2);
%         tappend = tappend(returns_2_pol2);
%         
%         dt = -log(rand(sum(returns_2_pol2),1))/(kabort+kesc);
%         
%         
%         tappend = tappend + dt;
%         new_outs = times_OUT(indexes);
%         times_OUT(indexes) = tappend;
%         %new_outs = tappend - log(rand(length(tappend),1))/(kout);  %move old out here
%         
%         %ones that continue to escape
% 
%         dt_e = -log(rand(sum(~returns_2_pol2),1))/(kabort+kesc);
%         
%         times_escape(reverse_indexes) = times_ser5(reverse_indexes) + dt_e;
%         
%         times_ser5_abort(indexes) = tappend;
%         
%     end
%     %new_times_ser5 = tappend - log(rand(tappend,1))/(kinit); 
%     
%     cnt = cnt + 1;
%     size(tappend)
%     times_IN = [times_IN; tappend];
%     times_OUT = [times_OUT; new_outs];
%     times_ser5 = [times_ser5; inf*ones(size(tappend))];
%     
%     times_ser5_abort = [times_ser5_abort; inf*ones(size(tappend))];
%     times_escape = [times_escape; inf*ones(size(tappend))];
%     
%     
%     if cnt > 1e7
%         return
%     end
%     
%     
% end
% Nsim = length(times_ser5);
% times_ser5_dephos = times_ser5 - log(rand(Nsim,1))/(kdphos);
% 
% % Does this POL2 escape to start transcription?
% transcr =  (times_escape ~= inf);
% 
% % If not, it will not transcribe mRNA.
% times_escape(~transcr) = inf;
% 
% if model_options.trypt(3) == 1   %BLOCKING K_escape
%     
%     times_escape(times_escape>model_options.trypt_time+trypt_scale+sqrt(trypt_scale)*randn()) = inf;
% end
% 
% 
% 
% % Add processing if required by model.
% if model_options.processing==1
%     kproc = parameters(10);
%     times_processing = - log(rand(Nsim,1))/(kproc);
% %     times_processing = max(0,1/kproc + 1/sqrt(3)*1/kproc*randn(Nsim,1));
% else
%     times_processing = zeros(Nsim,1);
% end
% 
% % if model_options.elongation_randomness == 1
% %     random_elong = parameters(11);
% %     times_processing = max(-Lgene_end/kelong,times_processing+random_elong*randn(Nsim,1));
% % end
%   
% 
% % Updates times at which Pol2 (i) leaves the bubble or (ii) leaves due to
% % Ser5 dephosphorylation or (iii) leaves due to completion of
% % transcription.
% 
% % This case is redundant, but is left here and commented for completeness.
% % times_OUT(~ser5_phos) = times_OUT(~ser5_phos); 
% 
% % Note -- The following line from old code makes POL2 leave at same rate
% % even after Ser5 is phosphorylated.
% %      times_OUT(ser5_phos&~transcr) = min(times_OUT(ser5_phos&~transcr),times_ser5_abort(ser5_phos&~transcr));
% % I (Brian) changed this to now leave only when Ser5 aborts.
% times_OUT2 = times_OUT;
% %times_OUT(ser5_phos&~transcr) = times_ser5_abort(ser5_phos&~transcr);
% 
% 
% 
% if model_options.trypt(2) == 1    %BLOCKING K_INIT
%     %times_IN = times_IN(times_IN<model_options.trypt_time);
%     
%     times_OUT(~transcr&times_IN > ser5inhibit_t) = times_OUT2(~transcr&times_IN > ser5inhibit_t);
%     %times_OUT(times_OUT>ser5inhibit_t ) = times_OUT2(times_OUT2>ser5inhibit_t);
%    
%     %transcr = ser5_phos&(times_escape<=ser5inhibit_t & times_escape ~= inf);
%      %transcr = ser5_phos&(times_escape<=times_ser5_abort & times_escape ~= inf);
% end
% 
% if model_options.trypt(3) == 1     %BLOCKING K_escape
%     
%     times_escape(times_escape>model_options.trypt_time+trypt_scale+sqrt(trypt_scale)*randn()) = inf;
%     
%     transcr = ser5_phos&(times_escape<=times_ser5_abort & times_escape ~= inf);
%     
% end
% 
% 
% % Note -- The following line from old code makes POL2 leave at same rate
% % even after transcription begins.
% %       times_OUT(transcr) = min(times_OUT(transcr),times_escape(transcr)+Lgene_end/kelong+times_processing(transcr));
% % I (Brian) changed this to now leave only when transcription completes.
% %times_OUT(transcr) = times_escape(transcr &times_escape ~=inf)+Lgene_end/kelong+times_processing(transcr);
% % POL2_int = []
% % for it = 1:Nt    %put this into coder of its own function
% %     t = t_array(it);
% %     if it == 100
% %         x=1
% %         
% %     end
% % %     POL2_int(it) = sum(times_IN<t&times_OUT>t);
% %     POL2_int(it) = sum((t-times_IN((times_IN<t&times_OUT>t))));
% % end
% % 
% % plot(POL2_int)
% times_OUT(transcr&times_escape ~=inf) = times_escape(transcr &times_escape ~=inf)+Lgene_end/kelong+times_processing(transcr&times_escape ~=inf);
% 
% 
% if model_options.trypt(3) == 1      %BLOCKING K_escape
%     
%     times_OUT(ser5_phos&~transcr) = times_ser5_abort(ser5_phos&~transcr);
%     
% end
% 
% % initialize termination times.
% times_terminate = 0*times_OUT;
% % calculate time of termination for each mRNA
% times_terminate(transcr) = times_escape(transcr)+Lgene_end/kelong+times_processing(transcr);
% times_terminate(~transcr) = -inf;
% 
% % Note -- The following line from old code makes Ser5 leave at same rate
% % even after transcription begins.
% %      times_ser5_abort(transcr) = min(times_ser5_abort(transcr),times_terminate(transcr));
% % I (Brian) changed this to now leave only when transcription completes.
% times_ser5_abort(transcr) = times_terminate(transcr);
% % Dephosphorylation of Ser5 could happen prior to completion.
% times_ser5_dephos(transcr) = min(times_ser5_dephos(transcr),times_terminate(transcr));
% 
% % Simulated tryptolide

if sum(model_options.trypt) == 0
    %tmax = max(times_IN);
    tmin = max(0,floor(tmax)-traj_length);
else
    tmax = model_options.trypt_time+150;
    tmin = model_options.trypt_time-100;

    delta_t = 0.1;        
end 



t_array = [tmin:delta_t:tmax];
Nt = length(t_array);

POL2_int = zeros(Nt,1);
SER5_int = zeros(Nt,1);
TS_int = zeros(Nt,1);

% G = max([times_OUT,times_ser5_abort,times_escape+LMS2_start/kelong],[],2);
G = times_OUT;
J = (G>=tmin);

times_IN = times_IN(J);
times_OUT = times_OUT(J);
times_ser5 = times_ser5(J);


times_ser5_abort = times_ser5_abort(J);
times_ser5_dephos = times_ser5_dephos(J);
times_terminate = times_terminate(J);
times_escape = times_escape(J);

size(times_IN);
chip_ts_vec = [];
chip_pol2_vec = [];
chip_ser5_vec = [];
chip =0 ;
preallocated_tescape =(times_escape+LMS2_start/kelong);
fact = 1;
% Convert event times to molecule number at all times.
for it = 1:Nt    %put this into coder of its own function
    t = t_array(it);

%     POL2_int(it) = sum(times_IN<t&times_OUT>t);
    POL2_int(it) = sum(exp(-blch_pol2*(t-times_IN((times_IN<t&times_OUT>t)))));
    %     SER5_int(it) = sum(times_ser5<t&times_ser5_abort>t&times_ser5_dephos>t);
    SER5_int(it) = sum(exp(-blch_ser5*(t-times_ser5((times_ser5<t&  times_ser5_dephos>t & times_OUT > t)))));

    
    
%     TS_int(it) = sum((times_terminate>t).*min(1,max(0,(t-(times_escape+LMS2_start/kelong))/((LMS2_end-LMS2_start)/kelong))));

    if model_options.photobleach == 1  
        fact =  exp(-blch_ts*(t-(times_escape+LMS2_end/kelong)));
        fact = min(fact,1);
        fact(times_terminate<=t)=0;
        
    end 
    TS_int(it) = sum((times_terminate>t).*...
        min(1,max(0,(t-(times_escape+LMS2_start/kelong))/((LMS2_end-LMS2_start)/kelong))).*fact);

    
if chip == 1
    if mod(it,40) == 0

   ser_5_pre = times_ser5<t&times_ser5_abort>t&times_ser5_dephos>t;
     
        
    ts_pre = (times_terminate>t).*min(1,max(0,(t-preallocated_tescape)/((LMS2_end-LMS2_start)/kelong))) ;
    ts_pre2 = ts_pre;
    ts_pre2(ts_pre==0) = 5;
    
    ser5_elong = (ts_pre2 == ser_5_pre);
    
        time_elong = Lgene_end/kelong;
        chip_ts = max(0,(t-preallocated_tescape) ).*(times_terminate>t);   %./(times_terminate-times_escape).*Lgene_end;

        chip_ts_processing = chip_ts>time_elong;
        chip_ts_processing2 = (chip_ts~=0);
        chip_ts_processing2(chip_ts_processing2 == 0) = 3;

        chip_ser5_processing = chip_ts(ser5_elong==chip_ts_processing2)~=0.*Lgene_end;
        chip_ser5_processing = chip_ser5_processing(chip_ser5_processing~=0).*Lgene_end;
        chip_ts_processing = (chip_ts(chip_ts>time_elong)~=0).*Lgene_end;
        chip_ts(chip_ts==0) = Inf;

        chip_ts_elong_pre = double(chip_ts<=time_elong);
        chip_ts_elong_pre(chip_ts_elong_pre == 0) = 5;


        chip_ser5_elong = chip_ts(  ser5_elong == chip_ts_elong_pre)  /time_elong.*Lgene_end;

        chip_ts_elong = chip_ts(chip_ts<=time_elong)/time_elong.*Lgene_end;
        num_ser5 = sum(ser_5_pre);
        num_at_front_ser5 = num_ser5 - length(chip_ser5_elong) - length(chip_ser5_processing);
        chip_ser5 = [    chip_ser5_elong', chip_ser5_processing', zeros(num_at_front_ser5,1)'];



    chip_pol2 = zeros(POL2_int(it) - length(chip_ts),1);
    
    chip_ts_vec = [chip_ts_vec, chip_ts_processing', chip_ts_elong'];
    chip_pol2_vec =  [chip_pol2_vec, chip_pol2', chip_ts_processing', chip_ts_elong', chip_ser5];
    chip_ser5_vec =  [chip_ser5_vec, chip_ser5];
    end
    
end




end





% Rescaling (not used)
% if model_options.calibration == 1
%     POL2 = POL2_int*parameters(15);
%     SER5 = SER5_int*parameters(16);
%     TS = TS_int*parameters(17);
% else
    POL2 = POL2_int;
    SER5 = SER5_int;
    TS = TS_int;    
% end
    
% Compute auto- and cross-correlations (before adding noise);
% [POL2n, SER5n, TSn] = Normalize_simulated_intensities(0.95,POL2,SER5,TS);
% V = [POL2n,SER5n,TSn];
% V = (V - repmat(mean(V),Nt,1) );%./std(V,0,1);
% [C,lags] = xcorr(V);
% C(round(length(C)/2):end,:) = C(round(length(C)/2):end,:)./(Nt:-1:1)'; %Xcorr scale fix
% C(1:round(length(C)/2)-1,:) = C(1:round(length(C)/2)-1,:)./(1:1:Nt-1)'; 
% close all;
% plot([POL2_int(1:200), SER5_int(1:200), TS_int(1:200)])
% 
% figure
% plot([POL2_int, SER5_int, TS_int])
% Add shot noise.

if model_options.shot_noise == 1
   
    
%     POL2 = POL2 + parameters(12)*randn(size(POL2));
%     SER5 = SER5 + parameters(13)*randn(size(SER5));
%     TS = TS + parameters(14)*randn(size(TS));
    
    POL2 = POL2 + std(POL2)*parameters(12)*randn(size(POL2));
    SER5 = SER5 + std(SER5)*parameters(13)*randn(size(SER5));
    TS = TS + std(TS)*parameters(14)*randn(size(TS));
    
end   

% 
% if model_options.shot_noise == 1
%     POL2 =poissrnd( POL2, length(POL2));
%     SER5 =poissrnd( SER5, length(SER5));
%     TS = poissrnd( TS, length(TS));
% end    
POL2i = POL2_int;
SER5i = SER5_int;
TSi = TS_int;
try
[POL2, SER5, TS] = Normalize_simulated_intensities(0.95,POL2,SER5,TS);

catch
    xx=1
end
% Exit for tryptolide analysis.
if sum(model_options.trypt) ~= 0
    C = [];
    lags =[];
return
end


% 
% figure; 
% subplot(3,3,1)
% histogram(POL2);
% subplot(3,3,5)
% histogram(SER5);
% subplot(3,3,9)
% histogram(TS);
% subplot(3,3,2)
% scatter(SER5,POL2)
% subplot(3,3,3)
% scatter(TS,POL2)
% subplot(3,3,6)
% scatter(TS,SER5)
% subplot(3,3,4)
% 
% ksdensity([SER5(:),POL2(:)],'PlotFcn','contour');
% set(gca,'xlim',[-2 2],'ylim',[-2 2])
% subplot(3,3,7)
% ksdensity([TS(:),POL2(:)],'PlotFcn','contour');
% set(gca,'xlim',[-2 2],'ylim',[-2 2])
% subplot(3,3,8)
% ksdensity([TS(:),SER5(:)],'PlotFcn','contour');
% set(gca,'xlim',[-2 2],'ylim',[-2 2])

% % Compute auto- and cross-correlations.
V = [POL2,SER5,TS];
V = (V - repmat(mean(V),Nt,1) );%./std(V,0,1);
[C,lags] = xcorr(V);
C(round(length(C)/2):end,:) = C(round(length(C)/2):end,:)./(Nt:-1:1)'; %Xcorr scale fix
C(1:round(length(C)/2)-1,:) = C(1:round(length(C)/2)-1,:)./(1:1:Nt-1)'; 

% shift initial time to t0 = 0;
t_array = t_array-min(t_array);

end
