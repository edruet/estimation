// Example of fitting: ordinary least squares used for a linear model

// Options of the algorithm
options = struct( ...
    'randseed'   , 0       ,...// Initialization value for random generation
    'sigmay'     , 0.2     ,...// Standard deviation of noise on Y
    'R'          , [1 ; 2] ,...// Reference model used to compute Y from X
    'dX'         , 0.0025   ...// Increment of X to go from 0 to 1.
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
// Direct solution using the linear model (exploiting structure of matrix A)
// The matrix (A' * A) is a 2x2 matrix that can be inverted with the sample
// formula: [a11 a12; a21 a22]^-1 = [a11 -a12; -a21 a22] / (a11 * a22 - a12 * a21)
N = size(Y,1);
sumx = sum(X);
sumy = sum(Y);
sumx2 = sum(X.^2);
sumxy = sum(X .* Y);
// Computing of the parameters of the model with the detailed 
S = [(sumx2 * sumy - sumx * sumxy) ; (-sumx * sumy + N * sumxy)] / (N * sumx2 - sumx^2);
// Computing the estimates for Y-values based on solution S
Z = A * S;

// Summary of options
disp(S - R, "OLS parameters errors", S, "OLS solution", options, "Options of the script: ");

// Showing diagrams
// Plotting the data X-Y in blue
plot(X, Y, 'b.');
// Plotting the reference model in green
plot(X, A * R, 'g', 'LineWidth', 2)
// Plotting the direct solution in red
plot(X, A * S, 'r', 'LineWidth', 2);
title("Ordinary least squares fitting", "fontsize", 4);
xlabel("X: input values", "fontsize", 2);
ylabel("Y: put values", "fontsize", 2);
legend(['Samples', 'Reference model', 'Computed model'], 4);
