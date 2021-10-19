// Example of fitting: generalized least squares used for a linear model

// Options of the algorithm
options = struct( ...
    'randseed'   , 0       ,...// Initialization value for random generation
    'sigmay1'    , 0.2     ,...// Standard deviation 1 of noise on Y
    'sigmay2'    , 1       ,...// Standard deviation 2 of noise on Y
    'R'          , [1 ; 2] ,...// Reference model used to compute Y from X
    'dX'         , 0.0025   ...// Increment of X to go from 0 to 1.
);
// Parameters of the model used to generate date
R = options.R;
// X values (inputs of the model)
X = (0:options.dX:1)';
// Number of samples
N = size(X,1);
// Initialization of random generation
rand("seed", options.randseed);
// Initializing the vector to create the weight matrix
V = zeros(1,N);
// Initialization of the output vector with the right size
Y = zeros(N,1);
// Computation of the output with alternate noises
i = (1:2:N)';
Y(i) = R(1) + R(2) * X(i) + options.sigmay1 * rand(size(i,1),1,'normal');
V(i) = options.sigmay1;
i = (2:2:N)';
Y(i) = R(1) + R(2) * X(i) + options.sigmay2 * rand(size(i,1),1,'normal');
V(i) = options.sigmay2;
// Creating the weight matrix
W = diag(1 ./ (V.^2));
// Matrix to compute a solution for the linear model: W * A * S = W * Y
A = [ones(N,1) X];
// Classical solution using the inversion of (A' * A)
S = (A' * W * A)^-1 * A' * W * Y;
// Computing the estimates for Y-values based on solution S
Z = A * S;
// Show the result
// Summary of options
disp(S - R, "GLS parameters errors", S, "GLS solution", options, "Options of the script: ");
// Showing diagrams
// Plotting the data X-Y in blue
plot(X, Y, 'b.');
// Plotting the reference model in green
plot(X, A * R, 'g', 'LineWidth', 2)
// Plotting the direct solution in red
plot(X, A * S, 'r', 'LineWidth', 2);
title("Generalized least squares fitting", "fontsize", 4);
xlabel("X: input values", "fontsize", 2);
ylabel("Y: put values", "fontsize", 2);
legend(['Samples', 'Reference model', 'Computed model'], 4);
// Cost function
J = ((A * S) - Y)' * W^-1 * ((A * S) - Y);
disp(J, "Cost function");
