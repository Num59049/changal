clear; clc; close all;

//당연히 모든 단위는 SI unit이다.

// 파라미터


p.L1 = 0.44; p.L2 = 0.45; p.l = 0.17, p.d = 0.62;
//허벅지, 정강이, 페달의 길이, 페달 중심과 고관절 사이 거리

p.m1 = 8; p.m2 = 4; p.Ic = 0.05;
//허벅지, 정강이+발 질량, 체인링 회전관성모멘트

p.v0 = 0; p.f = 0; p.r = 0; p.M = 0;//등속주행조건을 위해 필요

p.N = 1440; //각도를 몇 번 쪼갤것인지

th = transpose((0:p.N-1))*pi/p.N; //쪼갠 각도 열벡터
p.dth = pi/p.N; //th간격

cadList = [40 50 60 70 80 90 100 ];//케이던스리스트래


[ph, k] = fourbar(th,p);
// 양발을 고려할 필요가 없음. 생각을 해보면 자전거를 탈 때 페달을 당길 수 없음 당연히. 양발 다 미는데 자명히 한쪽발은 앞으로 가는것을 방해하므로 한쪽으로만 미는게 제일 낭비 없음. 그러니 그 상황을 가정하면, 자전거 운동의 주기는 2pi가 아닌 pi.

J = legJ(th,phi,k,p);
Jp = legJp(th, J);



// 원래 여기 검증이 있었는데 필요없어서 버림

// 사실 우리는 미방을 풀 필요가 없음. 운동을 알아서 그럼. 식에 넣어서 계산만 하면됨.



// 페달의 회전각, 매개변수를 입력받아 고관절 각도, 고관절 각도를 페달의 회전각을 미분한것을 출력
function [phi, k] = fourbar(t,p)

    P = 2*p.L1(p.l*cos(t)+p.d);
    Q = 2*p.L1*p.l*sin(t);
    K = p.L1^2+p.l^2-p.L2^2+p.d^2+2&p.d*p.l*cos(t);
    //Pcos(phi)+Qsin(phi) = K가 성립하는 K를 찾는거임

    phi = acos(K./hypot(P,Q))+atan2(Q,P);
    //atan2(y,x) : 동경 (x,y)의 각을 [-pi,pi]에서 반환

    //-acos가 아닌 이유 : 기하상 그렇게 될 수밖에 없음. 사람의 무릎상.

    Pp = -2*L1*l*sin(t);
    Qp = 2*L1*l*cos(t);
    Kp = -2*d*l*sin(t);
    //prime의 의미로 p. dP/dt, dQ/dt, dK/dt를 나타냄.

    k=(Kp-Pp.*cos(phi)-Qp.*sin(phi))./(Q.*cos(phi)-P.*sin(phi));// 음함수 미분으로 k = dphi/dt

end


// 허벅지 + 종아리의 등가 회전관성모멘트. 정확히는 에너지 기준으로, 다리의 에너지가 1/2*J*w^2이도록 하는 J의함수임
function J = legJ(t,phi,k,p)
    
    //Mps, bps는 안봄 뭔지아직모름
    Mps = p.L1^2 * k.^2 + p.l^2 + 2 * p.L1 * p.l * k.*cos(phi-t);
    //||dM/dt||^2인데 G prime square라 Mps(...), M = (B+C)/2


    u = p.l*sin(t)-p.L1*k.*sin(phi);
    v = -p.lc*cos(t)-p.d+p.L1*cos(phi);
    //(u,v) = B-C
    up = p.lc*cos(t)-p.L1*k.*cos(phi);
    vp = p.lc*sin(t)-p.L1*k.*sin(phi);

    bps = (u.*vp-v.*up)/p.L2^2;
    //(dbeta/dt)^2

    J = p.m1*p.L1^2*k.^2/3+p.m2*Gps+m2*L2^2*bps/12;

end

function Jp = legJp(th, J)
    Jp = 어쩌구저쩌구
end

//모멘트팔 모양(rho+)
function rho1 = rho(th)

    rho1 = 어쩌구저쩌구

end