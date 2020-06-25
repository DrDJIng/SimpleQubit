function doubleQubit()
    # Spin operators
    sX = [0im 1.0; 1.0 0]
    sY = [0 -1.0im; 1.0im 0]
    sZ = [1.0 0im; 0 -1.0]
    # Creation and annihilation operators
    sUp = [0im 1.0; 0 0]
    sDown = Adjoint(sUp)
    # Qubit energies, in arbitrary energy units
    q1Energy = 10.0
    q2Energy = 10.0
    # Coupling energy
    g = 1.0
    # Define identity matrix
    identMatrix = Matrix{Complex{Float64}}(I,2,2)
    # Hamiltonian, with XX YY coupling, linking single excitation states together
    H = q1Energy .* kron(sZ, identMatrix) + q2Energy .* kron(identMatrix, sZ) + g .* (kron(sX, identMatrix) * kron(identMatrix, sX) + kron(sY, identMatrix) * kron(identMatrix, sY))
    # Initial density matrix
    rhoInit = kron(sUp * sDown, sDown * sUp) #Initialised in energetic state.
    # Lindblad operator
    L = kron(sZ, identMatrix)
    # Dephasing rate
    depRate = 1.0
    # Setting up the function to be parsed
    f(rho, p, t) = -1.0im .* (H * rho - rho * H) + depRate .* (L * rho * Adjoint(L) - 0.5 .* (Adjoint(L)*L*rho + rho*Adjoint(L)*L))
    # Time span for solving
    tspan = (0.0, 5)
    # Defining the problem as an ODE problem (Julia is great like this)
    prob = ODEProblem(f, rhoInit, tspan)
    # Solving the problem with solve, with manual solution spacings to maintain matrix size equivalence with matlab
    sol = solve(prob, Tsit5(), saveat = 0.01)
    # Defining complex64 zero sets for expectation value calculations
    expZ1 = zeros(ComplexF64, length(sol.t),)
    exp1Z = zeros(ComplexF64, length(sol.t),)
    # Calculating the Sz expectation values
    for ic = 1:length(sol.t)
      expZ1[ic] = tr(kron(sZ, identMatrix) * sol[:,:,ic])
      exp1Z[ic] = tr(kron(identMatrix, sZ) * sol[:,:,ic])
    end
    # Constructing output structure
    expZZ = (expZ1, exp1Z, sol.t)
    # Returning expectation values
    return expZZ
end
