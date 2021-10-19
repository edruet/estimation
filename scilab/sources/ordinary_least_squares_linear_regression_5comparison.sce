// Example of fitting: ordinary least squares (OLS-) used for a linear model
// This script checks the results provided by different methods:
// - Direct analytical formula for linear model
// - Solution using the inverse of (A' A)
// - Solution using the resolution of the linear system rather than the matrix inversion
// - The gradient descent

// Command to disable scilab warnings due to functions reload
funcprot(0);

// Options of the algorithm
options = struct( ...
    'randseed'   , 0       ,...// Initialization value for random generation
    'sigmay'     , 0.2     ,...// Standard deviation of noise on Y
    'R'          , [1 ; 2] ,...// Reference model used to compute Y from X
    'dX'         , 0.0025  ,...// Increment of X to go from 0 to 1.
    'graditermax', 1000    ,...// Max nb of iterations for gradient descent
    'alpha'      , 0.5     ,...// Learning rate for gradient descent
    'dJthreshold', 0.0000001  ...// Threshold to stop the gradient descent 
);
// Parameters of the model used to generate date
R = options.R;
// X values (inputs of the model)
X = (0:options.dX:1)';
// Initialization of random generation
rand("seed", options.randseed);
// Computation of the output
Y = R(1) + R(2) * X + options.sigmay * rand(size(X,1),1,'normal');
// Matrix to compute a solution for the linear model: A * S = Y
A = [ones(size(X,1),1) X];
// Computing the output values using the reference model
Zref = A * R;
// Direct solution using the linear model (exploiting structure of matrix A)
N = size(Y,1);
sumx = sum(X);
sumy = sum(Y);
sumx2 = sum(X.^2);
sumxy = sum(X .* Y);
Sdir = [(sumx2 * sumy - sumx * sumxy) ; (-sumx * sumy + N * sumxy)] / (N * sumx2 - sumx^2);
Zdir = A * Sdir;
// Solution using the classical analytical solution
Sinv = (A' * A)^-1 * A' * Y;
Zinv = A * Sinv;
// Solution using the resolution of a linear system (linsolve resolves A * x + b = 0) 
Slin = linsolve(A' * A, A' * -Y);
Zlin = A * Slin;
// Solution using gradient descent (S = S - alpha * (d/dSj) J(S))
Sgrad = [0; 0];
for i = 1:1:options.graditermax
    dJ = A' * ((A * Sgrad) - Y);
    Sgrad = Sgrad - (options.alpha / N) * dJ;
    if abs(dJ) < options.dJthreshold
        break;
    end 
end
Zgrad = A * Sgrad;
// Show the result
h = gcf();
ha = gca();
title("Ordinary least squares fitting", "fontsize", 4);
xlabel("X: input values", "fontsize", 2);
ylabel("Y: put values", "fontsize", 2);
// Plotting the data X-Y in blue
plot(X, Y, 'b.');
// Plotting the reference model in green
plot(X, A * R, 'g')
// Plotting the direct solution in cyan
plot(X, A * Sdir, 'cx', 'MarkerSize', 10);
// Plotting the solution using inversion in magenta
plot(X, A * Sinv, 'm+', 'MarkerSize', 10);
// Plotting the solution using linear system resolution in red
plot(X, A * Slin, 'ro', 'MarkerSize', 6);
// Plotting the solution computed with gradient descent in black
plot(X, A * Sgrad, 'k.', 'MarkerSize', 4);

legend(['reference', 'data', 'OLS-direct', 'OLS-inversion', 'OLS-linsolve', 'OLS-gradient'])

// Summary of options
disp(options, "Options of the script: ")
// Checking the mean and standard deviation of the recomputed output
disp(R', "Reference model");
disp(Sdir', "Direct solution");
disp(Sinv', "Classical solution");
disp(Slin', "Solution with linear resolution");
disp(Sgrad', "Solution using gradient descent");
disp(i, "Number of iterations for gradient descent");
disp(mean(Y - Zdir), "Mean of the output error");
disp(stdev(Y - Zdir), "Standard deviation of the output error");
