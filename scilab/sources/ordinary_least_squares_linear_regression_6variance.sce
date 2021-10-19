// This script illustrates the computation of the variance of the result
// provided by the least squares estimation

// Options of the algorithm
options = struct( ...
    'randseed'   , 0.5        , ...// Initialization value for random generation
    'sigmay'     , 0.5        , ...// Standard deviation of noise on Y
    'granddist'  , 'unf'      , ...// Distrib. used to generate Y noise (see doc of scilab)
    'grandp1'    , -sqrt(12)/4, ...// First parameter of grand fct (see doc of scilab)
    'grandp2'    , sqrt(12)/4 , ...// Second parameter of grand fct (see doc of scilab)
    'samplesnbr' , 1000       , ...// Number of executions of the fitting
    'R'          , [1 ; 2]    ,... // Reference model used to compute Y from X
    'dX'         , 0.0025      ... // Increment of X to go from 0 to 1.
);
// Initialization of random generation
rand("seed", options.randseed);
// Number of tests (using N as a short notation)
N = options.samplesnbr;
// Parameters of the model used to generate date
R = options.R;
// X values (inputs of the model)
X = (0:options.dX:1)';
// Number of samples for each simulation
Nx = size(X,1)
// Parameters errors for each simulation
Se = zeros(2, N);
// Estimation errors
Ye = zeros(Nx, N);
// Mean and std. dev. of the difference between Z (recomputed output) and Y (output)
Yd = zeros(2, N);
// Building the matrix A (column 1 = 1, column 2 = X) which is the same for all simulations
A = [ones(Nx,1) X];
// Precomputing once for all the matrix multipyling Y
M = (A' * A)^-1 * A';
// Computing the covariance matrix once for all
Q = options.sigmay^2 * (A' * A)^-1;
// Computing the reference output (without noise)
Zr = R(1) + R(2) * X;
// Running the simulations
for i = 1:1:N
    // Computation of the output from the reference model + noise
    Y = Zr + grand(Nx, 1, options.granddist, options.grandp1, options.grandp2);
    // Computing least squares using the classical formula and the precomputed matrix M
    S = M * Y;
    // Computing parameters errors
    Se(:,i) = S - R;
    // Computing the output from the model
    Z = S(1) + S(2) * X;
    // Computing output error with respect to the reference model
    Ye(:,i) = Z - Zr;
    // Difference between the model output and Y
    Yd(1,i) = mean(Z - Y);
    // Computing the std. dev. taking the number of parameters into account
    Yd(2,i) = sqrt(sum(((Z - Y) - Yd(1,i)).^2) / (Nx - size(S,1)));
end
// Theoretical covariance matrix computed from model
disp(Q, "Theoretical covariance matrix")
// Covariance matrix computed from sample
Qe = cov(Se(1,:)',Se(2,:)')
disp(Qe, "Covariance matrix from samples")
// Computing the standard deviation using the model and the covariance matrix
// Equivalent operation: sigmaz = sqrt(diag(A * Q * A'))
sigmaz = sqrt(Q(1,1) + Q(2,2) * X.^2 + 2 * X * Q(1,2));
// Computing the mean of the output error per X-value
meanze = mean(Ye, 2);
// Computing the standard deviation of the output error per X-value
sigmaze = stdev(Ye, 2);
// Showing the result
h = figure();
ha = newaxes();
title("Standard deviation of ordinary least squares fitting", "fontsize", 4);
xlabel("X: input values", "fontsize", 2);
ylabel("Y: put values", "fontsize", 2);
// Plotting the reference model
plot(X, Zr, 'b');
// Plotting only one point of the theoretical error for the legend
plot(X(1), Zr(1) - sigmaz(1), 'color', 'black', 'marker', '.', 'markersize', 3);
// Plotting only one point from output error for the legend
plot(X(1), Zr(1)+ meanze(1), 'color', 'red', 'marker', '.', 'markersize', 3);
// Plotting the theoretical standard deviation
plot([X'; X'], [(Zr - sigmaz)'; (Zr + sigmaz)'], 'color', 'black', ...
     'marker', '.', 'markersize', 3, 'markerfac', 'black', 'markeredg', 'black');
// Plotting the standard deviation of the output error
plot([(X+options.dX/10)'; (X+options.dX/10)'], [(Zr + meanze - sigmaze)'; (Zr + meanze + sigmaze)'], 'color', 'red', ...
     'marker', '.', 'markersize', 3, 'markerfac', 'red', 'markeredg', 'red');
// Plotting the mean deviation of the output error
plot(X, Zr + meanze, 'LineStyle', 'none', 'marker', '.', 'markersize', 3, 'markerfac', 'red');
legend(['reference model', 'theoretical std dev', 'output error std dev'], 4)
// New figure for the mean of the output errors
h = figure();
ha = newaxes();
title("Mean of the output errors", "fontsize", 4);
xlabel("X: input values", "fontsize", 2);
ylabel("Error between the reference and the recomputed model", "fontsize", 2);
plot(X, meanze, 'r');
// New figure for the standard deviations
h = figure();
ha = newaxes();
title("Standard deviation of the errors", "fontsize", 4);
xlabel("X: input values", "fontsize", 2);
ylabel("Error between the reference and the recomputed model", "fontsize", 2);
plot(X, sigmaz, 'k');
plot(X, sigmaze, 'r');
legend(['theoretical', 'from samples']);
// New figure for the mean and the standard deviation recomputed from samples for each simulation
h = figure();
ha = newaxes();
title("Standard deviation of the errors", "fontsize", 4);
xlabel("Iterations", "fontsize", 2);
ylabel("Error between the recomputed and the initial outputs", "fontsize", 2);
plot(Yd(1,:), 'k');
plot(Yd(2,:), 'r');
legend(['mean', 'standard deviation']);
