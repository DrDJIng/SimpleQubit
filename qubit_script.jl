# Note: to see proper speedup, must import these packages first, or run the file twice.
using DifferentialEquations
using LinearAlgebra
using Plots
using LaTeXStrings

#Including the file containing the function I want to use
#Essentially a script inside a function, to utilise Julia's strong typing system
include("doubleQubit.jl")

#Defining the solution as the output of the function defined in above import
sol = @time doubleQubit();

#Plotting the results, with axis labels, legend labels
@time plot(sol[3], real(sol[1]), xaxis = "Time (t)", yaxis = L"\left< \sigma^z \right>", label = "Qubit 1")
#plot! plots on already existing axis
@time plot!(sol[3], real(sol[2]), xaxis = "Time (t)", yaxis = L"\left< \sigma^z \right>", label = "Qubit 2")
#Inclusion of the relaxation operator, and other expectation values, is left as an exercise for the reader.
