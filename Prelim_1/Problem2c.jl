using DifferentialEquations
using Plots

#variables
function equation(dA,A,S,t)
    a_x = 1.5
    b_x = 5.0
    z_x = 0.4
    n_zx = 2.7
    n_xz = 2.7
    del_z = 1.0
    x_z = 1.5

    dA[1] = (a_x + b_x * S)/(1 + S + (A[2]/z_x)^n_zx) - A[1]
    dA[2] = 1/(1 + (A[1]/x_z)^n_xz) - (del_z * A[2])
end

#initialization
tspan = (0.0,100.0)
u0 = [0.0; 0.0]
S = range(0.001, length = 10^4, stop = 100)
X = zeros(length(S))
Z = zeros(length(S))

#looping to get the values of X and Z for different values of S
for i in 1:length(S)
    prob = ODEProblem(equation,u0,tspan,S[i])
    solution = solve(prob)
    X[i] = solution[1,end]
    Z[i] = solution[2,end]
end

plot(S,X, xaxis = ("S", :log), yaxis = ("X", :log), label =("X vs S)"))
