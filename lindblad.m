function f = lindblad(rho, H, relRate, L)
rho = reshape(rho, 4, 4);
f = -1i .* (H * rho - rho * H) + relRate .* (L * rho * (L') - 0.5 .* ((L')*L*rho + rho*(L')*L));
f = reshape(f, 16, 1);
end