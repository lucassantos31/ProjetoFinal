function Dinamica = calcularDinamica(Pcm)
    syms g

    
    n = Pcm{end}.n;
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

    q  = Pcm{end}.juntas(1,:);
    dq = Pcm{end}.juntas(2,:);

    for i = 1:n
        Iaux           = Pcm{i}.T(1:3,1:3) * I(i) * Pcm{i}.T(1:3,1:3);
        Dinamica.Ep(i) = (Massa(i)) * Dinamica.gravidade * (Pcm{i}.P);
        Dinamica.Ec(i) = (Massa(i)/2) * transpose(Pcm{i}.V) * Pcm{i}.V + (Massa(i)/2) * transpose(Pcm{i}.W) * Iaux * (Pcm{i}.W);
        Dinamica.EpT   = Dinamica.EpT + Dinamica.Ep(i);
        Dinamica.EcT   = Dinamica.EcT + Dinamica.Ec(i);
    end

    Dinamica.Lagrangeano = -Dinamica.EcT + Dinamica.EpT;

    for i = 1:n
        Dinamica.G(i,1) = 0;
        for k = 1:n
            Dinamica.G(i,1) = Dinamica.G(i,1) + diff(Dinamica.Ep(k), q(i));
        end

        Dinamica.H2(i,1) = -diff(Dinamica.EcT, q(i));
        H1_sum = sym(0);
        for j = 1:n
            Dinamica.M(i,j) = diff(diff(Dinamica.Ec(i),dq(i)),dq(j));
            H1_sum = H1_sum + diff(diff(Dinamica.EcT, dq(i)), q(j)) * dq(j);
        end

        Dinamica.H1(i,1) = H1_sum;
    end

    Dinamica.H = Dinamica.H1 + Dinamica.H2;

end