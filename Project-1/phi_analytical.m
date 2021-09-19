function y=phi_analytical(Pe,L,phi_o,phi_L, X)
y = (exp(X*Pe/L) - 1)/(exp(Pe) - 1);
y = y*(phi_L - phi_o);
y=y + phi_o;
end
