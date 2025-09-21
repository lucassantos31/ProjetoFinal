function Dinamica = calcularDinamica(Pcm)
    syms g

    
    n = Pcm{6}.n;
    syms Massa [1,n]
    syms Ixx   [n,1]
    syms Ixy   [n,1]
    syms Ixz   [n,1]
    syms Iyx   [n,1]
    syms Iyy   [n,1]
    syms Iyz   [n,1]
    syms Izx   [n,1]
    syms Izy   [n,1]
    syms Izz   [n,1]

    I = cell(1,n);
    for i=1:n
        I{i} = [ Ixx(i), -Ixy(i), -Ixz(i);
                -Iyx(i),  Iyy(i), -Iyz(i);
                -Izx(i), -Izy(i),  Izz(i)];
    end
    

    H1_aux             = sym(zeros(n,1));    
    Dinamica           = struct();
    Dinamica.G         = sym(zeros(n,1));
    Dinamica.M         = sym(zeros(n,n));
    Dinamica.H1        = sym(zeros(n,1));
    Dinamica.H2        = sym(zeros(n,1));
    Dinamica.Ec        = sym(zeros(1,n));
    Dinamica.Ep        = sym(zeros(1,n));
    Dinamica.EcT       = 0;
    Dinamica.EpT       = 0;
    Dinamica.gravidade = [0 0 g];

    
    for i = 1:n
        q              = Pcm{i}.juntas(1,:);
        dq             = Pcm{i}.juntas(2,:);
        Iaux           = Pcm{i}.T(1:3,1:3) * I{i} * Pcm{i}.T(1:3,1:3).';
        Dinamica.Ep(i) = (Massa(i)) * Dinamica.gravidade * (Pcm{i}.P);
        Dinamica.Ec(i) = (Massa(i)/2) * transpose(Pcm{i}.V) * Pcm{i}.V + (Massa(i)/2) * transpose(Pcm{i}.W) * Iaux * (Pcm{i}.W);
    end

    Dinamica.Lagrangeano = -Dinamica.EcT + Dinamica.EpT;

    for i = 1:n
        Dinamica.G(i,1)  =  diff(Dinamica.Ep(i),q(i));
        Dinamica.H2(i,1) = -diff(Dinamica.Ec(i),q(i));
        for j = 1:n
            Dinamica.M(i,j) = diff(diff(Dinamica.Ec(i),dq(i)),dq(j));
            H1_aux(j,1) = diff(diff(Dinamica.Ec(i),dq(i)),q(j))*dq(j);
        end

        Dinamica.H1(i,1) = sum(H1_aux);
    end

    Dinamica.H = Dinamica.H1 + Dinamica.H2;

end