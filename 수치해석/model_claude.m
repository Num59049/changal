%% =====================================================================
%  비원형 체인링 통합 모델 (중간성과 기준 버전)
%  A층 사람: 4절 링크 + 힘-속도 관계 -> 부담 최소화
%  B층 기하: 유효반지름(지지함수) -> 실제 윤곽, 유한 축간거리 보정
%  C층 구동계: 프리휠 포함 운동방정식 적분 -> 준정적 가정 검증
%  툴박스 없이 돌아감. (fmincon이 있으면 마지막 확인 단계도 실행)
%  =====================================================================
clear; clc; close all;

%% ---------------- 0. 파라미터 ----------------
p.L1 = 0.44;  p.L2 = 0.45;  p.lc = 0.17;  p.d = 0.62;  % 허벅지, 정강이, 크랭크, 고관절-크랭크축 [m]
p.rho0 = 0.083;                    % 기준(원형) 유효반지름 [m]
p.tau_iso = 200;  p.w_max = 15;    % 힘-속도 직선: 정지 토크 [N m], 한계 각속도 [rad/s]
p.P0 = 200;                        % 목표 동력 [W]
p.m1 = 8;  p.m2 = 4;  p.Ic = 0.05; % 허벅지, 정강이+발 질량 [kg], 크랭크 관성 [kg m^2]
p.N  = 1440;
% 각도 기준: theta = 페달이 맨 아래(하사점)에서 무릎 쪽(앞)으로 잰 각.
th = (0:p.N-1)''    *2*pi/p.N;  p.dth = 2*pi/p.N;
cadList = [40 50 60 70 80 90 100];

%% ---------------- 1. 4절 링크 ----------------
[phR,kR] = fourbar(th,p);   [phL,kL] = fourbar(th+pi,p);
K.th = th; K.kR = kR; K.kL = kL;

%% ---------------- 2. 다리 관성 ----------------
J = legJ(th,phR,kR,p) + legJ(th+pi,phL,kL,p) + p.Ic;
Jp = dper(J,p);
fprintf('[관성] 두 다리+크랭크 유효관성: %.3f ~ %.3f kg m^2\n\n', min(J), max(J));

%% ---------------- 3. 원형 기준선 ----------------
fprintf('--- 원형 체인링 기준선 (%d W) ---\n', p.P0);
E0 = zeros(size(cadList));
for i=1:numel(cadList)
    E0(i) = burden(p.rho0*ones(p.N,1),cadList(i),K,p);
    fprintf('%4d rpm : 부담 %.4f\n', cadList(i), E0(i));
end

%% ---------------- 4. 최적 형상 (2차 형상 격자 탐색) ----------------
// rho(th) = rho0 (1 + e cos(2 th - ph)),  최소 곡률반지름 = rho0 (1-3e)
kmin = 0.0;                    // 최소 곡률반지름 / rho0  (0 = 볼록 한계, 0.4 = 지름비 1.5)
emax = (1-kmin)/3 - 1e-4;
eG = linspace(0,emax,61);  phG = linspace(0,pi,181);
fprintf('\n--- 최적 형상 (최소 곡률반지름 >= %.2f rho0) ---\n', kmin);
fprintf('%5s | %7s %7s %7s | %6s %6s %7s | %8s %8s\n','rpm','원형','최적','개선','e','지름비','긴축deg','장력min원','장력min최');
opt = struct();
for i=1:numel(cadList)
    cad = cadList(i);  best = inf;
    for e = eG
        for ph = phG
            Eb = burden(shp(e,ph,th,p),cad,K,p);
            if Eb < best, best = Eb; be = e; bph = ph; end
        end
    end
    rho = shp(be,bph,th,p);  [~,im] = max(rho);
    Tc = tension(p.rho0*ones(p.N,1),cad,K,J,Jp,p);  To = tension(rho,cad,K,J,Jp,p);
    opt(i).e = be; opt(i).ph = bph; opt(i).E = best;
    fprintf('%5d | %7.4f %7.4f %6.2f%% | %6.3f %6.2f %7.1f | %8.1f %8.1f\n', cad, E0(i), best, ...
        100*(1-best/E0(i)), be, (1+be)/(1-be), rad2deg(th(im)), min(Tc), min(To));
end

%% ---------------- 5. 찌그러짐 정도별 개선폭 ----------------
fprintf('\n--- 찌그러짐 정도별 개선폭 ---\n');
eList = [0.05 0.0909 0.15 0.20 0.25 0.333];  cadT = [40 60 100];
imp = zeros(numel(eList),numel(cadT));
for a=1:numel(eList)
    for b=1:numel(cadT)
        bb = inf; for ph=phG, bb=min(bb,burden(shp(eList(a),ph,th,p),cadT(b),K,p)); end
        imp(a,b) = 100*(1-bb/burden(p.rho0*ones(p.N,1),cadT(b),K,p));
    end
    fprintf(' e=%.3f 지름비 %.2f : %5.2f%% %5.2f%% %5.2f%%  (40/60/100 rpm)\n', ...
        eList(a),(1+eList(a))/(1-eList(a)),imp(a,:));
end
figure; plot((1+eList)./(1-eList),imp,'o-','LineWidth',2); grid on;
xlabel('지름비 (긴 지름/짧은 지름)'); ylabel('개선폭 [%]'); legend('40 rpm','60 rpm','100 rpm','Location','northwest');
hold on; plot([1.2 1.2],ylim,'k--'); title('찌그러짐과 개선폭 (점선: 시판 수준)');

%% ---------------- 6. 실제 윤곽과 둘레 ----------------
i60 = find(cadList==60);  e=opt(i60).e; ph=opt(i60).ph;
s = th;                                           % 법선각
h  = p.rho0*(1+e*cos(pi-2*s-ph));                 % 지지함수 h(s) = rho(pi/2 - s)
hp = 2*p.rho0*e*sin(pi-2*s-ph);
x = h.*cos(s)-hp.*sin(s);  y = h.*sin(s)+hp.*cos(s);
per = sum(hypot(diff([x;x(1)]),diff([y;y(1)])));
fprintf('\n[둘레] 실제 %.5f m  vs  2*pi*rho0 %.5f m (평균반지름 고정 = 톱니 수 보존)\n', per, 2*pi*p.rho0);
figure; plot(1000*x,1000*y,'LineWidth',2); hold on;
plot(1000*p.rho0*cos(s),1000*p.rho0*sin(s),'--'); axis equal; grid on;
legend('최적 체인링','원형'); xlabel('mm'); ylabel('mm'); title('60 rpm 최적 형상 윤곽');

%% ---------------- 7. 축간거리 유한 보정 ----------------
D = 0.41;  rs = 0.025;
rhoInf = shp(e,ph,th,p);  rhoEx = zeros(p.N,1);
hf = @(s) p.rho0*(1+e*cos(pi-2*s-ph));
for j=1:p.N
    f = @(ps) hf(ps-th(j)) - D*cos(ps) - rs;
    psj = fzero(f,[pi/2-0.7, pi/2+0.7]);
    rhoEx(j) = hf(psj-th(j));
end
Eb0 = burden(p.rho0*ones(p.N,1),60,K,p);
fprintf('[축간거리 %.2f m] 개선폭 근사 %.2f%%  ->  실제 접선 %.2f%%\n', D, ...
    100*(1-burden(rhoInf,60,K,p)/Eb0), 100*(1-burden(rhoEx,60,K,p)/Eb0));

%% ---------------- 8. 구동계 운동방정식 (프리휠 포함) ----------------
q.M=85; q.Jrot=0.20; q.Rw=0.34; q.Crr=0.005; q.g=9.81; q.rhoA=1.225; q.CdA=0.40;
q.Meq = q.M + q.Jrot/q.Rw^2;
fprintf('\n--- 구동계 검증 (마지막 10초) ---\n');
cases = {60,0.030; 100,0.045};
for cI=1:size(cases,1)
    cad = cases{cI,1};  rsI = cases{cI,2};
    for kind=1:2
        if kind==1, rho = p.rho0*ones(p.N,1); lab='원형';
        else, ii=find(cadList==cad); rho=shp(opt(ii).e,opt(ii).ph,th,p); lab='최적'; end
        R = simDrive(rho,cad,rsI,K,J,Jp,p,q);
        fprintf('%4d rpm %s: 속도 %.3f (목표 %.3f) m/s, 체인속도 변동 %.2f%%, 프리휠 풀림 %.1f%%\n', ...
            cad,lab,R.v,R.vt,R.chain,R.free);
    end
end

%% ---------------- 9. (선택) 3차 성분까지 열어서 확인 ----------------
if exist('fmincon','file')==2
    fprintf('\n--- fmincon: 3차까지 열어도 2차만 남는지 (60 rpm, 최소곡률 0.4) ---\n');
    NH=3; obj=@(c) burden(shpF(c,th,p,NH),60,K,p);
    nlc=@(c) deal([0.4 - min(curvF(c,th,p,NH))/p.rho0; 0.45*p.rho0-min(shpF(c,th,p,NH)); ...
                   max(shpF(c,th,p,NH))-1.55*p.rho0],[]);
    o=optimoptions('fmincon','Display','off','Algorithm','sqp');
    c=fmincon(obj,[0 0 0.1 0.1 0 0]',[],[],[],[],[],[],nlc,o);
    fprintf(' 성분 크기 n1 %.4f  n2 %.4f  n3 %.4f\n', hypot(c(1),c(2)),hypot(c(3),c(4)),hypot(c(5),c(6)));
end

%% ===================== 함수 =====================
function [phi,k] = fourbar(t,p)
    P=2*p.L1*(p.d+p.lc*cos(t)); Q=2*p.L1*p.lc*sin(t);
    K=p.L1^2+p.lc^2+p.d^2+2*p.d*p.lc*cos(t)-p.L2^2;
    phi=atan2(Q,P)+acos(min(max(K./hypot(P,Q),-1),1));
    Pp=-2*p.L1*p.lc*sin(t); Qp=2*p.L1*p.lc*cos(t); Kp=-2*p.d*p.lc*sin(t);
    k=(Kp-Pp.*cos(phi)-Qp.*sin(phi))./(Q.*cos(phi)-P.*sin(phi));
end
function y = dper(f,p)                 % 주기 중앙차분
    y = (circshift(f,-1)-circshift(f,1))/(2*p.dth);
end
function J = legJ(t,phi,k,p)
    u=p.lc*sin(t)-p.L1*sin(phi); v=-p.lc*cos(t)-p.d+p.L1*cos(phi);
    ud=p.lc*cos(t)-p.L1*k.*cos(phi); vd=p.lc*sin(t)-p.L1*k.*sin(phi);
    lam=(u.*vd-v.*ud)/p.L2^2;
    J = p.m1*p.L1^2/3*k.^2 + p.m2/4*(p.L1^2*k.^2+p.lc^2+2*p.L1*p.lc*k.*cos(phi-t)) ...
        + p.m2*p.L2^2/12*lam.^2;
end
function y = cap(w,p)
    y = max(p.tau_iso*(1-abs(w)/p.w_max),0);
end
function [c,V0,td] = capTot(rho,cad,K,p)
    V0 = p.rho0*cad*2*pi/60;  td = V0./rho;
    c  = abs(K.kR).*cap(K.kR.*td,p) + abs(K.kL).*cap(K.kL.*td,p);
end
function E = burden(rho,cad,K,p)
    [c,V0] = capTot(rho,cad,K,p);
    E = 2*pi*p.rho0*(p.P0/V0) / (sum(c)*p.dth);
end
function rho = shp(e,ph,t,p)
    rho = p.rho0*(1+e*cos(2*t-ph));
end
function T = tension(rho,cad,K,J,Jp,p)
    [c,V0,td] = capTot(rho,cad,K,p);
    Q = 2*pi*p.rho0*(p.P0/V0)/(sum(c)*p.dth)*c;
    T = (Q - (J.*td.*dper(td,p) + 0.5*Jp.*td.^2))./rho;
end
function y = lookP(arr,x,p)            % 주기 선형보간
    s = mod(x,2*pi)/p.dth;  i0 = floor(s);  f = s-i0;
    y = arr(mod(i0,p.N)+1)*(1-f) + arr(mod(i0+1,p.N)+1)*f;
end
function R = simDrive(rho,cad,rs,K,J,Jp,p,q)
    [c,V0] = capTot(rho,cad,K,p);
    Qp = 2*pi*p.rho0*(p.P0/V0)/(sum(c)*p.dth)*c;   % 준정적 최적 배분 토크
    G = rho*q.Rw/rs;  Gp = dper(G,p);
    vt = V0*q.Rw/rs;
    sg = (p.P0/vt - q.Crr*q.M*q.g - 0.5*q.rhoA*q.CdA*vt^2)/(q.M*q.g);   % 목표 속도에서 P0가 되는 경사
    F = @(v) q.Crr*q.M*q.g + q.M*q.g*sg + 0.5*q.rhoA*q.CdA*v*abs(v);
    dt=2e-4; Tsim=40; n=round(Tsim/dt);
    th_=0; w=cad*2*pi/60; v=vt; eng=true;
    nrec = round(10/dt); rv=zeros(nrec,1); re=false(nrec,1); rc=zeros(nrec,1); r=0;
    for i=1:n
        Jf=lookP(J,th_,p); Jfp=lookP(Jp,th_,p); qq=lookP(Qp,th_,p);
        gg=lookP(G,th_,p); ggp=lookP(Gp,th_,p);
        if eng
            a=(qq-gg*F(v)-0.5*(Jfp+2*q.Meq*gg*ggp)*w^2)/(Jf+q.Meq*gg^2);
            Fd=q.Meq*(gg*a+ggp*w^2)+F(v);
            if Fd<0, eng=false; end
        end
        if ~eng
            a=(qq-0.5*Jfp*w^2)/Jf;  vd=-F(v)/q.Meq;
        end
        th_=th_+w*dt; w=w+a*dt;
        if eng
            v=gg*w;
        else
            v=v+vd*dt; g2=lookP(G,th_,p);
            if g2*w>v
                Jf2=lookP(J,th_,p);
                w=(Jf2*w+q.Meq*g2*v)/(Jf2+q.Meq*g2^2); v=g2*w; eng=true;
            end
        end
        if i>n-nrec
            r=r+1; rv(r)=v; re(r)=eng; rc(r)=lookP(rho,th_,p)*w;
        end
    end
    R.v=mean(rv); R.vt=vt; R.free=100*(1-mean(re));
    R.chain=100*(max(rc(re))-min(rc(re)))/V0;
end
function rho = shpF(c,t,p,NH)
    rho = p.rho0*ones(size(t));
    for n=1:NH, rho = rho + p.rho0*(c(2*n-1)*cos(n*t)+c(2*n)*sin(n*t)); end
end
function r = curvF(c,t,p,NH)
    r = p.rho0*ones(size(t));
    for n=1:NH, r = r + p.rho0*(1-n^2)*(c(2*n-1)*cos(n*t)+c(2*n)*sin(n*t)); end
end
