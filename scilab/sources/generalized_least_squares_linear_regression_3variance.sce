// Example of fitting: comparison of variances between ordinary/generalized least squares

// Options of the algorithm
options = struct( ...
    'randseed'   , 0.5        , ...// Initialization value for random generation
    'sigmay1'    , 0.2        , ...// Standard deviation of noise on Y
    'granddist1' , 'nor'      , ...// Distrib. used to generate Y noise (see doc of scilab)
    'grandp11'   , 0          , ...// First parameter of grand fct (see doc of scilab)
    'grandp12'   , 0.2        , ...// Second parameter of grand fct (see doc of scilab)
    'sigmay2'    , 0.5        , ...// Standard deviation of noise on Y
    'granddist2' , 'unf'      , ...// Distrib. used to generate Y noise (see doc of scilab)
    'grandp21'   , -0.5*(sqrt(12)*0.5),...// First parameter of grand fct (see doc of scilab)
    'grandp22'   ,  0.5*(sqrt(12)*0.5),...// Second parameter of grand fct (see doc of scilab)
    'samplesnbr' , 2000       , ...// Number of executions of the fitting
    'R'          , [1 ; 2]    ,... // Reference model used to compute Y from X
    'dX'         , 0.0025      ... // Increment of X to go from 0 to 1.
);
// Parameters of the model used to generate date
R = options.R;
// Number of simulations
N = options.samplesnbr;
// X values (inputs of the model)
X = (0:options.dX:1)';
// Number of samples
Nx = size(X,1);
// Showing the options of the script
disp(options, "Options of the script: ");
// Initialization of random generation
rand("seed", options.randseed);
// Initializing the vector to create the weight matrix
V = zeros(1,Nx);
// Odd indexes
iodd = (1:2:Nx)';
// Elements of the weight matrix
V(iodd) = options.sigmay1;
// Even indexes
ieven = (2:2:Nx)';
// Elements of the weight matrix
V(ieven) = options.sigmay2;
// Creating the weight matrix
W = diag(1 ./ V.^2);
// Matrix to compute a solution for the linear model: A * S = Y
A = [ones(Nx,1) X];
// Precomputing once for all the matrix multipyling Y for OLS
Mols = (A' * A)^-1 * A';
// Precomputing once for all the matrix multipyling Y for OLS
Mgls = (A' * W * A)^-1 * A' * W;
// Initialization of the output vector with the right size
Y = zeros(Nx,1);
// Computing the reference output (without noise)
Zr = R(1) + R(2) * X;
// Parameters errors for each simulation
Seols = zeros(2, N);
// Estimation errors
Yeols = zeros(Nx, N);
// Parameters errors for each simulation
Segls = zeros(2, N);
// Estimation errors
Yegls = zeros(Nx, N);
// Running the simulations
for i = 1:1:N
    // Elements with odd index have noise 1
    Y(iodd) = R(1) + R(2) * X(iodd) + grand(size(iodd,1), 1, options.granddist1, options.grandp11, options.grandp12);
    // Elements with even index have noise 2
    Y(ieven) = R(1) + R(2) * X(ieven) + grand(size(ieven,1), 1, options.granddist2, options.grandp21, options.grandp22);
    // Computing OLS squares using the classical formula and the precomputed matrix M
    S = Mols * Y;
    // Computing parameters errors
    Seols(:,i) = S - R;
    // Computing the output from the model
    Z = S(1) + S(2) * X;
    // Computing output error with respect to the reference model
    Yeols(:,i) = Z - Zr;
    // Classical solution using the inversion of (A' * A)
    S = Mgls * Y;
    // Computing parameters errors
    Segls(:,i) = S - R;
    // Computing the output from the model
    Z = S(1) + S(2) * X;
    // Computing output error with respect to the reference model
    Yegls(:,i) = Z - Zr;
end
// Showing diagrams
// Theoretical covariance matrix of the OLS
Qols = (A' * A)^-1 * A' * diag(V.^2) * A * (A' * A)^-1;
disp(Qols, "Theoretical covariance matrix for OLS");
// Theretical covariance matrix of the GLS
Qgls = (A' * W * A)^-1; 
disp(Qgls, "Theoretical covariance matrix for GLS");
// Observed covariance matrices for OLS (from simulations)
Qeols = cov(Seols(1,:)',Seols(2,:)')
disp(Qeols, "Observed covariance matrix for OLS");
// Observed covariance matrices for GLS (from simulations)
Qegls = cov(Segls(1,:)',Segls(2,:)')
disp(Qegls, "Observed covariance matrix for GLS");
// Theoretical standard deviation of output for OLS
sigmazols = sqrt(diag(A * Qols * A'));
// Observed mean of the output error per X-value for OLS
meanzeols = mean(Yeols, 2);
disp(max(meanzeols), 'Max bias for OLS');
// Observed the standard deviation of the output error per X-value for OLS
sigmazeols = stdev(Yeols, 2);
// Theoretical standard deviation of output for GLS
sigmazgls = sqrt(diag(A * Qgls * A'));
// Observed mean of the output error per X-value for GLS
meanzegls = mean(Yegls, 2);
disp(max(meanzegls), 'Max bias for GLS');
// Observed the standard deviation of the output error per X-value for GLS
sigmazegls = stdev(Yegls, 2);
// Showing the result
h = figure();
ha = newaxes();
title("Standard deviations for OLS and GLS", "fontsize", 4);
xlabel("X: input values", "fontsize", 2);
ylabel("Y: put values", "fontsize", 2);
// Plotting the reference model
plot(X, sigmazols, 'b');
plot(X, sigmazeols, 'k');
plot(X, sigmazgls, 'g');
plot(X, sigmazegls, 'r');
legend(['Theoretical OLS', 'Observed OLS', 'Theoretical GLS', 'Observed GLS'], 4)
// Plotting only one point of the theoretical error for the legend
//plot(X(1), Zr(1) - sigmazols(1), 'color', 'grey', 'marker', '.', 'markersize', 3);
// Plotting only one point from output error for the legend
//plot(X(1), Zr(1)+ meanzeols(1), 'color', 'red', 'marker', '.', 'markersize', 3);
// Plotting the theoretical standard deviation
//plot([X'; X'], [(Zr - sigmaz)'; (Zr + sigmaz)'], 'color', 'black', ...
 //    'marker', '.', 'markersize', 3, 'markerfac', 'black', 'markeredg', 'black');
// Plotting the standard deviation of the output error
//plot([(X+options.dX/10)'; (X+options.dX/10)'], [(Zr + meanze - sigmaze)'; (Zr + meanze + sigmaze)'], 'color', 'red', ...
//     'marker', '.', 'markersize', 3, 'markerfac', 'red', 'markeredg', 'red');
// Plotting the mean deviation of the output error
//plot(X, Zr + meanze, 'LineStyle', 'none', 'marker', '.', 'markersize', 3, 'markerfac', 'red');
//legend(['reference model', 'theoretical std dev', 'output error std dev'], 4)
// New figure for the mean of the output errors


// Difference between the model output and Y
// Yd(1,i) = mean(Z - Y);
// Computing the std. dev. taking the number of parameters into account
// Yd(2,i) = sqrt(sum(((Z - Y) - Yd(1,i)).^2) / (Nx - size(S,1)));
