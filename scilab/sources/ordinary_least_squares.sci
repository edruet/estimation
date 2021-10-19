// Function that executes the ordinary least squares estimation
// Inputs parameters:
//    X: input values of the model
//    Y: output of the model
//    options: structure with some additional parameters
// Outputs parameters:
//    S: solution of the OLS equation A * S = y
//    Z: estimate for the output based on the computed model
//    Q: covariance matrix of the parameteres of the model
function [S, Z, Q] = fit_ols(X, Y, options)
    // Building the matrix A (column 1 = 1, column 2 = X)
    A = [ones(size(X,1),1) X];
    // Computing least squares using inversion of A
    S = (A' * A)^-1 * A' * Y;
    // Computing the estimates from X and from the solution S
    Z = S(1) + S(2) * X;
    // Computing the covariance matrix
    Q = options.sigmay^2 * (A' * A)^-1;
endfunction
