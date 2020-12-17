classdef SimpleModel

 properties
     %Defult properties of the model
     Nstates = 3;
     parameters = [0.433596648003337,1,0.666624850282066,0.198727194675144,15.4054402060124,0.777851806420281,1,1.98448561566709,1.42073542241570,0.408137598525255,0.924805163389989,1.22915843769510,1.14124029524139]; %(1/min) Rate of concentration reduction
    
     kon = .433596648003337;
     koff = 1000;
     kesc = 0.666624850282066;
     kproc = 0.198727194675144;
     kin =  1000*15.4054402060124;
     kout = 0.777851806420281;
     frac = 1;
     eta_rnap = 1.98448561566709;
     eta_ser5 = 1.42073542241570;
     eta_ts = 0.408137598525255;
    
     sc1 = 0.924805163389989;
     sc2 = 1.22915843769510;
     sc3 = 1.14124029524139;

         S =     [[ 1    -1     0     0     0     0];
         [0     0     1    -1    -1     0];
         [0     0     0     0     1    -1]];   
  
         
     
 end
 properties (Dependent)
   
     A
     EX
     means
     Q
     SIG
     c
     b
     W1
     W0
     stats
     
 end
 
 methods
    
     function stats = get.stats(self)
            p_esc = self.kesc/(self.kesc+ self.kout);
            p2_arrival = self.kin/1000  *  self.kon;
            m_cluster = (p2_arrival / (self.kesc + self.kout));
            av_mrna = m_cluster*(self.kesc/self.kproc); 
            
    stats = {'POLII arrival rate',p2_arrival, 'min^-1';
                 'Prob POLII escape', p_esc, '---';
                 'Average mRNA',av_mrna, 'molecules';
                 'Average time solo CTD',1/(self.kesc + self.kout), 'min';
                 'Average POLII',av_mrna+m_cluster, 'molecules';
                 'Average mRNA production',m_cluster*self.kesc,'molecules/min';
                 'mRNA completion time', 1/self.kproc,'min';
                 'Average solo CTD',m_cluster,'molecules';
                 
                 };
     end
     
     function c = get.c(self)
         c = [[0,1,1]; [0,self.frac,1];[0,0,1] ];
     end
     
     function b = get.b(self)
         b = [self.kon; 0;0;];
     end   
     
     function A = get.A(self)
         A =  self.S*self.W1;
     end
     
     function W1 = get.W1(self)
         W1 =[[-self.kon         0         0];
            [self.koff         0         0];
            [self.kin         0         0];
             [    0    self.kout         0];
             [    0    self.kesc*self.frac         0];
             [    0         0    self.kproc];];         
     end
     
     function W0 = get.W0(self)
         W0 = [[self.kon; 0; 0; 0; 0; 0;]];
     end
     
     function EX = get.EX(self)
       EX =  -self.A\self.b;
     end
     
     function means = get.means(self)
       means =  self.c*self.EX;
     end
     
     function Q = get.Q(self)
         Q = self.S*diag(self.W1*self.EX+self.W0)*self.S';
     end
     
     function SIG = get.SIG(self)
         SIG = lyap(self.A,self.Q);
     end
     
     
     
     function sigred = get_sigred(self, TT)
         
        x0 = zeros(length(self.A)^2,1);
        x0(:) = self.SIG;
        SIGc = self.c*self.SIG*self.c';
        n = length(self.A);
        PHI = spalloc(n^2,n^2,n^3);
        for i=1:n
            PHI((i-1)*n+1:i*n,(i-1)*n+1:i*n)=self.A;
        end
        fun = @(t,x)PHI*x;
        options = odeset('jacobian',PHI);
        [TT,YY] = ode23s(fun,TT,x0,options);      
        
        SIGt = 0*self.SIG;
        sigred = zeros(size(self.c,1)^2,length(TT));
        for i=1:length(TT)
            SIGt(:) = YY(i,:);
            tmp = self.c*SIGt*self.c';
            sigred(:,i) = tmp(:);
        end

        sgr = sigred;

        sigred(1,:) = sgr(1,:)/sgr(1,1);
        sigred(2,:) = self.sc1*sgr(2,:)/sqrt(sgr(1,1)*sgr(5,1));
        sigred(3,:) = self.sc2*sgr(3,:)/sqrt(sgr(1,1)*sgr(9,1));
        sigred(4,:) = self.sc1*sgr(4,:)/sqrt(sgr(5,1)*sgr(1,1));
        sigred(5,:) = sgr(5,:)/sgr(5,1);
        sigred(6,:) = self.sc3*sgr(6,:)/sqrt(sgr(5,1)*sgr(9,1));
        sigred(7,:) = self.sc2*sgr(7,:)/sqrt(sgr(9,1)*sgr(1,1));
        sigred(8,:) = self.sc3*sgr(8,:)/sqrt(sgr(9,1)*sgr(5,1));
        sigred(9,:) = sgr(9,:)/sgr(9,1);
        
        sigred([1 5 9],1) = sigred([1 5 9],1)+([self.eta_rnap,self.eta_ser5,self.eta_ts].^2)';

     end
     
     function [ana_means, ctd_ode, ts_ode] = ana_inhibs(self, tode, inhibs)


  
        
        
        W1after = [[-self.kon*inhibs(1)         0         0];
            [self.koff*inhibs(2)         0         0];
            [self.kin*inhibs(5)         0         0];
             [    0    self.kout*inhibs(6)         0];
             [    0    self.kesc*self.frac*inhibs(3)*inhibs(7)          0];
             [    0         0    self.kproc*inhibs(4)];]; 
         
         W0after = [[self.kon*inhibs(1); 0; 0; 0; 0; 0;]];
         
    
       Anew = self.S*W1after;
       
        on = [];
        ctd_ode = [];
        ts_ode = [];
            
            
        for j = 1:length(tode)

            %P = expm(A*tode(j)+ W0)*b;
        %      if j == 1
        %          b(2) = means(2);
        %          b(3) = means(3);
        %      else
        %          b(2) = 0;
        %          b(3) = 0;         
        %      end

             %P =  A^-1*(-eye(3) + expm(A*tode(j)))*b - x0;
             %P = A^-1*(-eye(3) + expm(A*tode(j)))*b ;
            xx = zeros(self.Nstates,1);
            xx(1) = self.kon;
            xx(2) = 0;
            xx(3) = 0;
            P= Anew^-1*(-eye(3) + expm(Anew*tode(j)))*xx + expm(Anew*tode(j))*[0;self.EX(2);self.EX(3)];
            on = [on, P(1)];
            ctd_ode = [ctd_ode,P(2)];
            ts_ode = [ts_ode,P(3)];
           
        end
        
        ana_means = [on; ctd_ode; ts_ode];
            
     end
     

     function [simulated_ssa, simulated_ssa_no_noise,simulated_ssa_norm] = ssa_traj_inhibs(self, time, inhibs, inhib_time)
         
        x0 = zeros(self.Nstates,1);
        
        time_var = 0;
        signal_update_rate = 0;

        W = @(x) self.W1*x + self.W0;
        
        
        W1after = [[-self.kon*inhibs(1)         0         0];
            [self.koff*inhibs(2)         0         0];
            [self.kin*inhibs(5)         0         0];
             [    0    self.kout*inhibs(6)         0];
             [    0    self.kesc*self.frac*inhibs(3)*inhibs(7)          0];
             [    0         0    self.kproc*inhibs(4)];]; 
         
         W0after = [[self.kon*inhibs(1); 0; 0; 0; 0; 0;]];
         
         Wafter = @(x) W1after*x + W0after;
         
         
        
        sol = run_single_SSA_generic_inhib(x0,self.S,W, Wafter,time,time_var,signal_update_rate, inhib_time);  

       
        [pol2_ssa_w_shot,ser5_ssa_w_shot,ts_ssa_w_shot,ts_ssa,pol2_ssa,ser5_ssa] = self.convert_ssa(sol);
     
        simulated_ssa = [pol2_ssa_w_shot'; ser5_ssa_w_shot';ts_ssa_w_shot'];
        simulated_ssa_no_noise = [pol2_ssa'; ser5_ssa';ts_ssa'];
        
        
        [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa_w_shot,ser5_ssa_w_shot,ts_ssa_w_shot);
        simulated_ssa_norm = [pol2norm'; ser5norm';tsnorm'];
     
 
        
     end
     
     
     function [simulated_ssa, simulated_ssa_no_noise,simulated_ssa_norm] = ssa_traj(self, time)
         
        x0 = zeros(self.Nstates,1);
        
        time_var = 0;
        signal_update_rate = 0;

        W = @(x) self.W1*x + self.W0;
        sol = run_single_SSA_linda(x0,self.S,W,time,time_var,signal_update_rate);  

       
        [pol2_ssa_w_shot,ser5_ssa_w_shot,ts_ssa_w_shot,ts_ssa,pol2_ssa,ser5_ssa] = self.convert_ssa(sol);
     
        simulated_ssa = [pol2_ssa_w_shot'; ser5_ssa_w_shot';ts_ssa_w_shot'];
        simulated_ssa_no_noise = [pol2_ssa'; ser5_ssa';ts_ssa'];
        
        
        [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa_w_shot,ser5_ssa_w_shot,ts_ssa_w_shot);
        simulated_ssa_norm = [pol2norm'; ser5norm';tsnorm'];
     
     end
     
     
     
     
    function [pol2_ssa_w_shot,ser5_ssa_w_shot,ts_ssa_w_shot,ts_ssa, pol2_ssa,ser5_ssa] = ...
    convert_ssa(self,sol)

        % Get integer numbers from simulations
        pol2_ssa = sol(2,:)' + sol(3,:)';
        ser5_ssa = pol2_ssa;
        ts_ssa = sol(3,:)';

        pol2_ssa_w_shot = self.add_shot_noise(pol2_ssa,self.eta_rnap,30/200);
        ser5_ssa_w_shot = self.add_shot_noise(ser5_ssa,self.eta_ser5,23/200);
        ts_ssa_w_shot = self.add_shot_noise(ts_ssa,self.eta_ts,5/200);
        % The third inputs are the fraction of data points that were assigned as
        % negative values.

        end

    function ssa_w_shot = add_shot_noise(self,ssa,eta,tr)
        std_mod = std(ssa);
        shot_mod = std_mod*eta;
        ssa_w_shot = ssa + randn(size(ssa))*shot_mod;

        % Define fraction tr as zero
        tmp = sort(ssa_w_shot);
        thresh = tmp(ceil(length(tmp)*tr));
        ssa_w_shot = ssa_w_shot-thresh;
     end


 end
end
