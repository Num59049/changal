clear; clc; close all;

//당연히 모든 단위는 SI unit이다.

// 파라미터
p.L1 = 0.44; p.L2 = 0.45; p.l = 0.17; p.d = 0.62;
//허벅지, 정강이, 페달의 길이, 페달 중심과 고관절 사이 거리

p.m1 = 8; p.m2 = 4; p.Ic = 0.05;
//허벅지, 정강이+발 질량, 체인링 회전관성모멘트

p.v0 = 0; p.f = 0; p.r = 0; p.M = 0; p.R = 0;//등속주행조건을 위해 필요
//v0: 속력, f: 마찰력, r: 스프로켓 반지름, M: 전체 질량, R: 바퀴반지름

p.N = 1440; //각도를 몇 번 쪼갤것인지
//p.rho: rho(theta)를 반환하는 체인링 형상 함수 핸들 (별도 지정)


// 원래 여기 검증이 있었는데 필요없어서 버림


// 페달의 회전각, 매개변수를 입력받아 고관절 각도, 고관절 각도를 페달의 회전각을 미분한것을 출력
function [phi, k] = fourbar(th,p)

    th = mod(th,pi); //반 바퀴마다 활성 다리를 교대: 0 <= th < pi
    P = 2*p.L1*(p.l*cos(th)+p.d);
    Q = 2*p.L1*p.l*sin(th);
    K = p.L1^2+p.l^2-p.L2^2+p.d^2+2*p.d*p.l*cos(th);
    //Pcos(phi)+Qsin(phi) = K가 성립하는 K를 찾는거임

    phi = acos(K./hypot(P,Q))+atan2(Q,P);
    //atan2(y,x) : 동경 (x,y)의 각을 [-pi,pi]에서 반환

    //-acos가 아닌 이유 : 기하상 그렇게 될 수밖에 없음. 사람의 무릎상.

    Pp = -2*p.L1*p.l*sin(th);
    Qp = 2*p.L1*p.l*cos(th);
    Kp = -2*p.d*p.l*sin(th);
    //prime의 의미로 p. dP/dth, dQ/dth, dK/dth를 나타냄.

    k=(Kp-Pp.*cos(phi)-Qp.*sin(phi))./(Q.*cos(phi)-P.*sin(phi));// 음함수 미분으로 k = dphi/dth

end

// 허벅지 + 종아리의 등가 회전관성모멘트. 정확히는 에너지 기준으로, 다리의 에너지가 1/2*J*w^2이도록 하는 J의함수임
function J = legJ(th,phi,k,p)
    th = mod(th,pi); //fourbar와 같은 활성 다리의 각도
    //Mps, bps는 안봄 뭔지아직모름
    Mps = (p.L1^2 * k.^2 + p.l^2 + 2 * p.L1 * p.l * k.*cos(phi-th))/4;
    //||dM/dth||^2인데 G prime square라 Mps(...), M = (B+C)/2


    u = p.l*sin(th)-p.L1*sin(phi);
    v = -p.l*cos(th)-p.d+p.L1*cos(phi);
    //(u,v) = C-B
    up = p.l*cos(th)-p.L1*k.*cos(phi);
    vp = p.l*sin(th)-p.L1*k.*sin(phi);

    bps = (u.*vp-v.*up)/p.L2^2;
    //dbeta/dth (제곱 전)

    J = p.m1*p.L1^2*k.^2/3 + p.m2*Mps + p.m2*p.L2^2*bps.^2/12;

end

function Jp = legJp(th, p)
    th = mod(th,pi);
    h = 2*pi/p.N;
    //교대 경계를 넘지 않게 차분점을 제한 (내부에서는 기존 중앙차분)
    thp = min(th+h,pi-eps(pi));
    thm = max(th-h,0);
    [phip, kp] = fourbar(thp,p);
    [phim, km] = fourbar(thm,p);
    Jp = (legJ(thp,phip,kp,p)-legJ(thm,phim,km,p))./(thp-thm);
    Jp(th==0) = NaN; //교대점에서는 주기 확장한 관성이 불연속이므로 미분 미정의
end

//모멘트팔 모양(rho+)
function rho1 = rho_function(th,p)

    rho1 = p.rho0*ones(size(th));
    for n = 1:numel(p.rho_c)/2
        rho1 = rho1 + p.rho0*(p.rho_c(2*n-1)*cos(n*th) + p.rho_c(2*n)*sin(n*th));
    end

end

function rhop = rho_prime(th,p)
    rhop = zeros(size(th));
    for n = 1:numel(p.rho_c)/2
        rhop = rhop + p.rho0*n*(-p.rho_c(2*n-1)*sin(n*th) + p.rho_c(2*n)*cos(n*th));
    end
end

// 시간에 따른 theta
function th = theta(t, rho, p)
    f = @(t, th) p.r*(p.v0)/(p.R)./rho(th);
    [~, th] = ode45(f, t, 0);
end

function w = theta_dot(th, rho, p)
    w = p.r*p.v0/p.R./rho;
end

function alpha = theta_dotdot(th, rho, rho_prime, p)
    alpha = -(p.r*p.v0/p.R)^2.*rho_prime./rho.^3;
end

//크랭크에 전달되는 페달 토크 (rho, rhop는 각도별 배열)
function tp = tau_p(rho,rhop,p)
    tp = (rho/p.r)*p.f*p.R - p.Ic*(p.r^2*p.v0^2)./(p.R^2*rho.^3).*rhop;
end


//끝판왕 함수, theta에 따른 rho와 시간을 받아, 그 시간에 따른 사람 고관절의 tau_h를 출력한다.
function tau_h = leg(t,p)

    rhof = @(th) rho_function(th,p);
    th = theta(t,rhof,p); //시간에 따른 체인링 회전각
    rho = rho_function(th, p);
    rhop = rho_prime(th,p);
    w = theta_dot(th,rho,p);
    a = theta_dotdot(th,rho,rhop,p);


    [phi, k] = fourbar(th, p);

    tau_h = (tau_p(rho,rhop,p) + legJ(th,phi,k,p).*a + 0.5*legJp(th,p).*w.^2)./k;

end
