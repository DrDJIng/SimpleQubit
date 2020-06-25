function runner
tic
sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1.0 0; 0 -1.0];
sUp = [0 1.0; 0 0];
sDown = sUp';
q1Energy = 10.0;
q2Energy = 10.0;
g = 1.0;
identMatrix = eye(2);
H = q1Energy .* kron(sZ, identMatrix) + q2Energy .* kron(identMatrix, sZ) + g .* (kron(sX, identMatrix) * kron(identMatrix, sX) + kron(sY, identMatrix) * kron(identMatrix, sY));
rhoInit = kron(sUp*sDown, sDown*sUp);
L = kron(sZ,identMatrix);
relRate = 1.0;
tspan = [0.0, 5.0];
[t,y] = ode45(@(t, rho) lindblad(rho, H, relRate, L), tspan, rhoInit);
expZ1 = zeros(length(t),1);
expZ2 = zeros(length(t),1);
for ic = 1:length(t)
    expZ1(ic) = trace(kron(sZ, identMatrix) * reshape(y(ic,:),4,4));
    expZ2(ic) = trace(kron(identMatrix, sZ) * reshape(y(ic,:),4,4));
end

plot(t, real(expZ1), t, real(expZ2))
legend("sz1","sz2")
xlabel("Time (t)");
ylabel("<\sigma^z>");
toc
end