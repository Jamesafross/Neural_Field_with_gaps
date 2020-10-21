

function J(u1,u2,u3,kappa,v_syn,beta,Alpha)


    return [-u3.-kappa.+2 .*u2 2 .*u1      0     -u1      ;
            -2*pi^2*u1     -u3+2*u2  0      v_syn   ;
            0              0        -Alpha  Alpha   ;
            Alpha*beta     0         0     -Alpha   ]
end
