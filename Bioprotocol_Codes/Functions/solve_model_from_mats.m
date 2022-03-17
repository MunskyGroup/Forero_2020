function [ModelCorrelations, ModelMeans, ModelVariances] = solve_model_from_mats(S, W1, W0, c, noise_parameters, time_vec, parameters)
    
    % solve the models analytical correlations, means, and variances

    % x' = S*W1(x) + W0

    % S - stoichiometry matrix 
    % W1 - Linear Propensity matrix(x)
    % W0 - Independent Propensity matrix [ kon 0 0 0 0 0 0 ] 1xN reactions
    % c - intensity transformation matrix, ie molecule to intensity matrix
    % noise_parameters - last Nchannel^2
    % time_vec - time vector of each sample
    % parameters - parameters to solve for

    % See Forero2020 Methods for full derivation
    b = @(p)(S*W0(p));
    TT = time_vec;
    A =  S*W1(parameters);
    EX =  -A\b(parameters);
    Q = S*diag(W1(parameters)*EX+W0(parameters))*S';
    SIG = lyap(A,Q);

    x0 = zeros(length(A)^2,1);
    x0(:) = SIG;
    SIGc = c(parameters)*SIG*c(parameters)';
    n = length(A);
    PHI = spalloc(n^2,n^2,n^3);
    for i=1:n
        PHI((i-1)*n+1:i*n,(i-1)*n+1:i*n)=A;
    end
    fun = @(t,x)PHI*x;
    options = odeset('jacobian',PHI);
    [TT,YY] = ode23s(fun,TT,x0,options);      
    
    SIGt = 0*SIG;
    sigred = zeros(size(c(parameters),1)^2,length(TT));
    for i=1:length(TT)
        SIGt(:) = YY(i,:);
        tmp = c(parameters)*SIGt*c(parameters)';
        sigred(:,i) = tmp(:);
    end

    sgr = sigred;
    G0s = [];
    ind_array = [];
    k = 1;
    for i = 1:length(A)
        for j = 1:length(A)
            if i == j
                G0s(end+1) = sgr(k,1);
                ind_array(end+1) = k;
            end
            k = k + 1;
        end
    end

    scs = noise_parameters(parameters);
    scs = scs(length(A)+1:end);

    scs_inds = [];
    k = 1;
    for i = 1:length(A)
        for j = 1:length(A)
            scs_inds(end+1) = sum([i,j]);
        end
        k = k + 1;
    end

    scs_inds = scs_inds - min(scs_inds);
    
    k = 1;
    for i = 1:length(A)
        for j = 1:length(A)
            if i == j
                sigred(k,:) = sgr(k,:)/G0s(i);
            else
                sigred(k,:) = scs(scs_inds(k))*sgr(k,:)/(sqrt(G0s(i))*sqrt(G0s(j)));
            end
            k = k + 1;
        end
    end
    etas = noise_parameters(parameters);
    etas = etas(1:length(A));
    sigred(ind_array,1) = sigred(ind_array,1)+(etas.^2)';

    ModelCorrelations = sigred; %analytical correlations
    ModelMeans = c(parameters)*EX; %analytical means
    ModelVariances = SIGc; %analytical variance