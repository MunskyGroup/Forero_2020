function [ZOut,b,Z_Array,A_fill] = Run_3_State_FSP_BEM(S,w,t_Array,f,b,eps_max)
%% Initialize states.
M = size(S,2);   % Number of reactions.
N_bd = length(b); % Number of boundary functions

err=1;                   % Initialize the error in the FSP.
while err > eps_max
    %% Compute Inf. Gen. Matrix.
    Nmax = ceil(b(4:6));                    % Find the maximum number of each species.
    
    Nst = prod(Nmax+1);                     % Total number of states.
    Y = zeros(Nst,1);X = zeros(Nst,1);Z = zeros(Nst,1);
    X(:) = repmat(repmat([0:Nmax(1)]',1,Nmax(2)+1),Nmax(3)+1,1);% x-value of all states
    Y(:) = repmat(repmat([0:Nmax(2)],Nmax(1)+1,1),1,Nmax(3)+1,1); % y-value of all states
    Z(:) = bsxfun(@(x,y)x*y,ones((Nmax(1)+1)*(Nmax(2)+1),1),[0:Nmax(3)]);
    
    inds = [1:Nst];                         % indices of the FSP states
    props = w(X,Y,Z);                       % propensities of all reactions at all states
    A_Square = spalloc(Nst,Nst,Nst*(M+1));  % allocate Inf.Gen.Mat.
    C = zeros(N_bd,Nst);                    % allocate output matrix
    for mu = 1:M
        stoich_1D = S(1,mu)+S(2,mu)*(Nmax(1)+1)+S(3,mu)*(Nmax(1)+1)*(Nmax(2)+1); % transform mu^th stoichiometry vector into 1D coordinates
        indsNew = inds + stoich_1D;              % find state after mu^th reaction
        XNew = X+S(1,mu); YNew = Y+S(2,mu); ZNew = Z+S(3,mu);
        % find state after mu^th reaction
        fNew = f(XNew,YNew,ZNew);                     
        % compute boundary functions at new state.
                
        A_Square = A_Square-spdiags([props(:,mu)],0,Nst,Nst);
        % Add diagonal to the
        % reduced inf.gen.mat.
        
        Valid_constraints = (fNew<=repmat(b,size(fNew,1),1));      % determine which reactions STAY in X_J
        props_keep = props(:,mu);
        props_keep(min(Valid_constraints,[],2)~=1) = 0;
        A_Square = A_Square+spdiags(props_keep,-stoich_1D,Nst,Nst);
        % Enter the props that stay in
        % X_J into the inf.gen.mat.
        
        for j=1:N_bd                        %specify the output matrix.
            props_reject = props(:,mu)./abs(sum(Valid_constraints-1,2));
            props_reject(Valid_constraints(:,j)==1) = 0;
            C(j,:) = C(j,:)+(props_reject)';
        end
    end
    
    Valid_constraints = min((f(X,Y,Z)<=repmat(b,size(X,1),1)),[],2);
    A_Red = A_Square(Valid_constraints,Valid_constraints);  % Find A_J
    B_Red = C(:,Valid_constraints);                         % Find A_{Ji,J}
    
    A_fill = [[A_Red;B_Red],spalloc(size(A_Red,1)+N_bd,N_bd,0)]; % Combine to find FSP Inf.Gen.Mat.
    P0 = sparse(zeros(size(A_fill,1),1)); P0(1)=1;          % Specify initial condition.
    
    m=30;
    tryagain=1;
    while tryagain==1
        [P, ~, ~, ~, P_array, sinks,tryagain] = mexpv_modified(max(t_Array), A_fill, P0, 1e-6, m, 2, t_Array,eps_max,[length(A_fill)-N_bd+1:length(A_fill)]);
        m=m+10;
    end
    % Use Krylov method to integrate the ODE dP/dt = A*P, assuming that A
    % is an inf.gen.mat.
    err = sum(sinks)            % Sum the prob. in the sinks.
%     sinks
%     pause
    if err>eps_max               % If necessary, expand the projection space.
        b(sinks>eps_max/10) = b(sinks>eps_max/10)*1.2
    end
end
Valid_constraints = min((f(X,Y,Z)<=repmat(b,size(X,1),1)),[],2);
ZOut = zeros(Nmax+1);
Z_Array = zeros([size(P_array,2),Nmax+1]);
for i = 2:size(P_array,2)
    ZOut(Valid_constraints) = P(1:end-N_bd);
    ZOut = reshape(ZOut,Nmax+1);
    Z_Array(i,:,:,:) =ZOut;
end



