// Example of fitting: generalized least squares used for a linear model

// Options of the algorithm
options = struct( ...
    'randseed'   , 0       ,...// Initialization value for random generation
    'sigmay1'    , 0.2     ,...// Standard deviation 1 of noise on Y
    'sigmay2'    , 1       ,...// Standard deviation 2 of noise on Y
    'R'          , [1 ; 2] ,...// Reference model used to compute Y from X
    'dX'         , 0.0025  ,...// Increment of X to go from 0 to 1.
    'graditermax', 10000   ,...// Max nb of iterations for gradient descent
    'alpha'      , 0.5     ,...// Learning rate for gradient descent
    'dJthreshold', 0.0000001  ...// Threshold to stop the gradient descent 
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
// Creating the inverse of the weight matrix for the determination of dJ
W_1 = diag(V.^2);
// Matrix to compute a solution for the linear model: A * S = Y
A = [ones(N,1) X];
// Resolution using gradient descent
// Initial values: [0, 0]
S = [0; 0];
// Iterate until the variation of the cost function is below the specified threshold
for i = 1:1:options.graditermax
    // Computing the delta of cost function for each parameter: A[ixj] * (A*S-Y)[1xj]
    dJ = A' * W_1 * ((A * S) - Y);
    // Updating the parameters of the model
    S = S - (options.alpha / N) * dJ;
    // Checking the norm of the update vector for cost function
    if norm(dJ) < options.dJthreshold
        break;
    end 
end
// Computing the estimates for Y-values based on solution S
Z = A * S;
// Showing the results
// Summary of options
disp(i, "Number of iterations", S - R, "GLS parameters errors", ...
     S, "GLS solution", options, "Options of the script: ");
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
// Cost function computed from the solution using gradient descent
J = ((A * S) - Y)' * W_1 * ((A * S) - Y);
disp(J, "Cost function of solution using gradient descent");
// Cost function computed from the analytical formula for GLS
W = diag(1 ./ V.^2);
S = (A' * W * A)^-1 * A' * W * Y;
J = ((A * S) - Y)' * W_1 * ((A * S) - Y);
disp(J, "Cost function of solution using analytical formula");
