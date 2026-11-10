function [p, result] = optimize_chainring(p)
% p: model_me의 매개변수 + rho0, NH, tau_iso, w_max, k_min

    th = transpose(((0:p.N-1) + 0.5) * 2*pi/p.N);

    B = zeros(p.N, 2*p.NH);
    D = B;
    C = B;

    for n = 1:p.NH
        B(:, 2*n-1) = cos(n*th);
        B(:, 2*n) = sin(n*th);
        D(:, 2*n-1) = -n*sin(n*th);
        D(:, 2*n) = n*cos(n*th);
        C(:, 2*n-1:2*n) = (1-n^2) * B(:, 2*n-1:2*n);
    end

    [k, J] = active_leg(th, p);
    h = 2*pi/p.N;
    hi = min(mod(th, pi) + h, pi-eps(pi));
    lo = max(mod(th, pi) - h, 0);
    [~, Jhi] = active_leg(hi, p);
    [~, Jlo] = active_leg(lo, p);
    Jp = (Jhi-Jlo) ./ (hi-lo);
    drive = abs(k) >= p.k_min;

    c0 = zeros(2*p.NH, 1);
    [result.E_circle, mean_circle] = objective(c0);
    [p.rho_c, result.E_opt] = fmincon(@objective, c0, [], [], [], [], [], [], @constraints);

    function [E, mean_u, min_capacity] = objective(c)
        rho = p.rho0 * (1 + B*c);
        rhop = p.rho0 * D*c;

        rho = max(rho(drive), p.rho0*1e-6);
        rhop = rhop(drive);
        w = (p.r*p.v0/p.R) ./ rho;
        a = -(p.r*p.v0/p.R)^2 * rhop ./ rho.^3;

        tp = rho/p.r * p.f*p.R + p.Ic*a;
        tau_h = (tp + J(drive).*a + 0.5*Jp(drive).*w.^2) ./ k(drive);
        capacity = max(p.tau_iso * (1-abs(k(drive).*w)/p.w_max), 0);
        u = abs(tau_h) ./ max(capacity, p.tau_iso*1e-6);

        weights = rho / sum(rho);  % dt = rho/V0 * dtheta
        mean_u = sum(weights .* u);
        E = sqrt(sum(weights .* (u-mean_u).^2));
        min_capacity = min(capacity) / p.tau_iso;
    end

    function [ineq, eq] = constraints(c)
        [~, mean_u, min_capacity] = objective(c);
        rho = 1 + B*c;
        curvature = 1 + C*c;

        n = transpose(1:p.NH);
        amplitude = hypot(c(1:2:end), c(2:2:end));
        mr = pi/p.N * sum(n .* amplitude);
        mc = pi/p.N * sum(n .* abs(1-n.^2) .* amplitude);

        ineq = [0.45+mr-rho; rho+mr-1.55; 0.4+mc-curvature;
                mean_u-mean_circle; 1e-6-min_capacity];
        eq = [];
    end
end

function [k, J] = active_leg(th, p)
% model_me의 fourbar, legJ와 같은 식
    th = mod(th, pi);

    P = 2*p.L1 * (p.d+p.l*cos(th));
    Q = 2*p.L1*p.l * sin(th);
    K = p.L1^2+p.l^2+p.d^2+2*p.d*p.l*cos(th)-p.L2^2;
    phi = atan2(Q, P) + acos(K ./ hypot(P, Q));

    Pp = -2*p.L1*p.l * sin(th);
    Qp = 2*p.L1*p.l * cos(th);
    Kp = -2*p.d*p.l * sin(th);
    k = (Kp-Pp.*cos(phi)-Qp.*sin(phi)) ./ (Q.*cos(phi)-P.*sin(phi));

    Mps = (p.L1^2*k.^2+p.l^2+2*p.L1*p.l*k.*cos(phi-th)) / 4;
    u = p.l*sin(th)-p.L1*sin(phi);
    v = -p.l*cos(th)-p.d+p.L1*cos(phi);
    up = p.l*cos(th)-p.L1*k.*cos(phi);
    vp = p.l*sin(th)-p.L1*k.*sin(phi);
    bps = (u.*vp-v.*up) / p.L2^2;

    J = p.m1*p.L1^2*k.^2/3 + p.m2*Mps + p.m2*p.L2^2*bps.^2/12;
end
