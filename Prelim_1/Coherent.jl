using DifferentialEquations
using Plots
gr(size=(500,500), show = true)

#Variables
a_x = 3.9 * 10^(-2)
a_y = 4.3 * 10^(-3)
b_x = 6.1
b_y = 5.7
del_y = 1.05
del_z = 1.04
z_x = 1.3 * 10^(-5)
y_z = 11 * 10^(-3)
x_z = 12 * 10^(-2)
x_y = 7.9 * 10^(-4)
nzx = 2.32
nxz = 2.0
nxy = 2.0
nyz = 2.0

function equation(dA, A, p, t)
    S = 100
    dA[1] = ((a_x + b_x * S)/(1 + S + (A[3]/z_x)^nzx)) - A[1]            #dX/dt
    dA[2] = ((a_y + b_y * S)/(1 + S + (A[1]/x_y)^nxy)) - (del_y * A[2])  #dY/dt
    dA[3] = 1/(1+((A[1]/x_z)^nxz)+((A[2]/y_z)^nyz)) - (del_z * A[3])     #dZ/dt
end

tspan = (0.0, 100.0)
# steady state value of X, Y, Z
u0 = [5.734114738, 0.005147128, 0.000420886]
prob1 = ODEProblem(equation,u0,tspan)
solution1 = solve(prob1)
# 25% above steady state value of X, Y, Z
u0 = [7.167643422,	0.006433911,	0.000526108]
prob2 = ODEProblem(equation,u0,tspan)
solution2 = solve(prob2)

#25% below steady state value of X, Y, Z
u0 = [4.300586053,	0.003860346,	0.000315665]
prob3 = ODEProblem(equation,u0,tspan)
solution3 = solve(prob1)

 plot(solution1, vars=(0,3), xlabel="Time", ylabel="Concentration (Z)", label=["Cell 1"])
 plot!(solution2, vars=(0,3), xlabel="Time", ylabel="Concentration (Z)", label=["Cell 2"])
 plot!(solution3, vars=(0,3), xlabel="Time", ylabel="Concentration (Z)", label=["Cell 3"])
