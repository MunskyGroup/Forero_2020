function fitresult = fitGaussian (x, y, sem)

nRepForFit =100;

%% this code is intended to apply the function fit to the experimental data. 
% the idea is to fit a Gaussian plot.
% input: x is the independent variable. y is the dependent variable. sem is
% the stander error of the mean.
% Output: fitresult returns the matlab object 'fitobject' that contains the
% parametrized fit function.

%%For Unph with Unph as min
[x, y, sem] = prepareCurveData( x, y, sem);
% Set up fittype and options.
ft = fittype( '(a*exp(-((x-b)/c)^2))+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxIter = 1e5;
opts.Weights = sem;
for k =1:nRepForFit
    % Fit model to data.
    opts.StartPoint = rand*[0.0461713906311539 -0.09 0.823457828327293 0.694828622975817];
    [fitresult, gof] = fit( x, y, ft, opts );
    selParam (k,:) = coeffvalues(fitresult);
    rmse(k) = gof.rmse;
end
[~,indexBest] = min(rmse);
opts.StartPoint = selParam(indexBest,:);
[fitresult, gof] = fit( x, y, ft, opts );

end