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

function equation1(dA, A, p, t)
    S = 0.02
    dA[1] = ((a_x + b_x * S)/(1 + S + (A[3]/z_x)^nzx)) - A[1]            #dx/dt
    dA[2] = ((a_y + b_y * S)/(1 + S + (A[1]/x_y)^nxy)) - (del_y * A[2])  #dy/dt
    dA[3] = 1/(1+((A[1]/x_z)^nxz)+((A[2]/y_z)^nyz)) - (del_z * A[3])     #dz/dt
end

function equation2(dA, A, p, t)
    S = 10
    dA[1] = ((a_x + b_x * S)/(1 + S + (A[3]/z_x)^nzx)) - A[1]            #dx/dt
    dA[2] = ((a_y + b_y * S)/(1 + S + (A[1]/x_y)^nxy)) - (del_y * A[2])  #dy/dt
    dA[3] = 1/(1+((A[1]/x_z)^nxz)+((A[2]/y_z)^nyz)) - (del_z * A[3])     #dz/dt
end

function equation3(dA, A, p, t)
    S = 10^5
    dA[1] = ((a_x + b_x * S)/(1 + S + (A[3]/z_x)^nzx)) - A[1]            #dx/dt
    dA[2] = ((a_y + b_y * S)/(1 + S + (A[1]/x_y)^nxy)) - (del_y * A[2])  #dy/dt
    dA[3] = 1/(1+((A[1]/x_z)^nxz)+((A[2]/y_z)^nyz)) - (del_z * A[3])     #dz/dt
end

u0 = [0.0, 0.0, 0.0]
tspan = (0.0, 100.0)
prob1 = ODEProblem(equation1,u0,tspan)
solution1 = solve(prob1)
prob2 = ODEProblem(equation2,u0,tspan)
solution2 = solve(prob2)
prob3 = ODEProblem(equation3,u0,tspan)
solution3 = solve(prob3)

 plot(solution1, vars=(0,1), xlabel="Time", ylabel="X", label=["S=0.02"])
 plot!(solution2, vars=(0,1), xlabel="Time", ylabel="X", label=["S=10"])
 plot!(solution3, vars=(0,1), xlabel="Time", ylabel="X", label=["S=10^5"])
