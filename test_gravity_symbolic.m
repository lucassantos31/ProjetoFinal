function test_gravity_symbolic()
%TEST_GRAVITY_SYMBOLIC Verifica se os torques gravitacionais aparecem nas
%juntas de base para um manipulador planar simples com dois elos.

    ordem = {'y','y'};
    excentricidades = [1 0 0; 1 0 0];

    robo = calcularTransformadas(ordem, excentricidades);
    cin = calcularCinematicaTotal(robo);

    Dinamica = calcularDinamica(cin);

    subsVars = {
        sym('Q1'),   sym('Q2'), ...
        sym('dQ1'),  sym('dQ2'), ...
        sym('d2Q1'), sym('d2Q2'), ...
        sym('Lx1'),  sym('Lx2'), ...
        sym('Lxcm1'), sym('Lxcm2'), ...
        sym('Lycm1'), sym('Lycm2'), ...
        sym('Lzcm1'), sym('Lzcm2'), ...
        sym('Massa1'), sym('Massa2'), ...
        sym('g')
    };

    vals = {
        0.3,  -0.2, ...   % posições articulares
        0.0,   0.0, ...   % velocidades articulares
        0.0,   0.0, ...   % acelerações articulares
        0.5,   0.4, ...   % comprimentos dos elos
        0.25,  0.2, ...   % centros de massa ao longo de X
        0.0,   0.0, ...   % centros de massa ao longo de Y
        0.0,   0.0, ...   % centros de massa ao longo de Z
        1.5,   1.0, ...   % massas
        9.81              % gravidade
    };

    Gval = double(subs(Dinamica.G, subsVars, vals));

    assert(abs(Gval(1)) > 1.0e-6, ...
        'O torque gravitacional na junta de base deve ser diferente de zero.');
    assert(abs(Gval(2)) > 1.0e-6, ...
        'O torque gravitacional na segunda junta deve ser diferente de zero.');
end
