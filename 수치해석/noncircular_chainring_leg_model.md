# 비원형 체인링–다리 운동 모델 수식 정리

## 1. 모델의 가정

- 뒷바퀴는 일정한 각속도 \(\omega\)로 회전한다.
- 체인은 미끄러지지 않으며 뒷스프로킷은 뒷바퀴와 함께 회전한다.
- 앞 체인링의 유효 반지름은 크랭크 회전각 \(\theta\)에 따른 함수 \(\rho(\theta)\)이다.
- 다리는 허벅지와 종아리로 구성된 평면 2링크이며 고관절에서만 능동 토크가 발생한다.
- 한 시점에 페달을 누르는 다리 하나만 고려하고, 중력에 의한 위치에너지 변화 및 손실은 무시한다.

## 2. 체인링과 크랭크 운동

| 기호 | 의미 |
|---|---|
| \(t\) | 시간 |
| \(\theta(t)\) | 크랭크 회전각 |
| \(\omega\) | 뒷바퀴의 일정한 각속도 |
| \(r_s\) | 뒷스프로킷 반지름 |
| \(\rho(\theta)\) | 체인링의 유효 반지름(크랭크 중심에서 체인 접선까지의 수직거리) |

체인의 속도 일치 조건은

\[
\rho(\theta)\dot\theta=r_s\omega.
\]

따라서

\[
\boxed{\dot\theta(t)=\frac{r_s\omega}{\rho(\theta(t))}}
\tag{1}
\]

이고, \(\theta(0)=0\)일 때

\[
\boxed{t=\frac{1}{r_s\omega}\int_0^{\theta(t)}\rho(u)\,du}
\tag{2}
\]

이다(\(\omega\neq0\)). 시간에 대해 미분하면

\[
\boxed{\ddot\theta(t)=-\frac{(r_s\omega)^2\rho'(\theta(t))}{\rho(\theta(t))^3}}
\tag{3}
\]

을 얻는다. 여기서 \(\rho'=d\rho/d\theta\)이다.

## 3. 다리의 기하학적 구속조건

고관절 \(A\), 무릎 \(B\), 페달 \(C\), 크랭크 중심 \(O\)를 사용한다. 허벅지 길이는 \(L_1\), 종아리 길이는 \(L_2\), 크랭크 길이는 \(l_c\), 고관절 높이는 \(d\), 허벅지 각도는 \(\phi\)이다.

\[
O=(0,0),\qquad A=(0,d),
\]

\[
\boxed{B=(L_1\sin\phi,\ d-L_1\cos\phi)}
\tag{4}
\]

\[
\boxed{C=(l_c\sin\theta,\ -l_c\cos\theta)}
\tag{5}
\]

종아리 길이가 일정하므로

\[
\boxed{(L_1\sin\phi-l_c\sin\theta)^2+(d-L_1\cos\phi+l_c\cos\theta)^2=L_2^2}
\tag{6}
\]

이다. 다음을 정의한다.

\[
\begin{aligned}
P(\theta)&=2L_1(d+l_c\cos\theta),\\
Q(\theta)&=2L_1l_c\sin\theta,\\
K(\theta)&=L_1^2+l_c^2+d^2+2dl_c\cos\theta-L_2^2.
\end{aligned}
\]

그러면

\[
\boxed{P\cos\phi+Q\sin\phi=K}
\tag{7}
\]

이며, 기존 모델에서 선택한 기하학적 해는

\[
\boxed{\phi(\theta)=\operatorname{atan2}(Q,P)+\arccos\!\left(\frac{K}{\sqrt{P^2+Q^2}}\right)}
\tag{8}
\]

이다. 기하학적으로 가능한 범위와 해의 가지를 전제로 한다.

## 4. 각도비 \(k(\theta)\)

\[
\boxed{k(\theta)=\frac{d\phi}{d\theta}}
\tag{9}
\]

식 (7)을 미분하면

\[
\boxed{k(\theta)=\frac{K'-P'\cos\phi-Q'\sin\phi}{Q\cos\phi-P\sin\phi}}
\tag{10}
\]

이다. 여기서

\[
P'=-2L_1l_c\sin\theta,\quad Q'=2L_1l_c\cos\theta,\quad K'=-2dl_c\sin\theta.
\]

따라서

\[
\boxed{\dot\phi=k(\theta)\dot\theta}
\tag{11}
\]

이다.

## 5. 다리의 등가 관성모멘트

허벅지와 종아리는 각각 질량 \(m_1,m_2\)인 균일한 막대로 가정한다. 다리의 전체 운동에너지를

\[
\boxed{T_{\mathrm{leg}}=\frac12J_{\mathrm{leg}}(\theta)\dot\theta^2}
\tag{12}
\]

로 표현하여 \(J_{\mathrm{leg}}\)를 정의한다.

허벅지는 고관절을 중심으로 회전하므로

\[
\boxed{T_1=\frac16m_1L_1^2k(\theta)^2\dot\theta^2}
\tag{13}
\]

이다.

종아리의 방향각을 \(\psi\)라 하고 \(\dot\psi=\lambda(\theta)\dot\theta\)라 두면

\[
\boxed{\lambda(\theta)=\frac{uv'-vu'}{L_2^2}}
\tag{14}
\]

이다. 여기서

\[
\begin{aligned}
u&=l_c\sin\theta-L_1\sin\phi,\\
v&=-l_c\cos\theta-d+L_1\cos\phi,\\
u'&=l_c\cos\theta-L_1k\cos\phi,\\
v'&=l_c\sin\theta-L_1k\sin\phi.
\end{aligned}
\]

종아리 질량중심의 병진 운동과 질량중심 주위 회전을 합하면

\[
\boxed{\begin{aligned}
J_{\mathrm{leg}}(\theta)={}&\frac{m_1L_1^2}{3}k^2\\
&+\frac{m_2}{4}\left(L_1^2k^2+l_c^2+2L_1l_ck\cos(\phi-\theta)\right)\\
&+\frac{m_2L_2^2}{12}\lambda^2.
\end{aligned}}
\tag{15}
\]

\(\phi,k,\lambda\)는 모두 \(\theta\)의 함수이다.

## 6. 고관절 토크와 페달 토크

\(\tau_h(t)\)는 고관절 토크, \(\tau_p(\theta)\)는 크랭크에 전달되는 페달 토크이다. 다리의 운동에너지 변화율은

\[
\frac{dT_{\mathrm{leg}}}{dt}=\frac12J'_{\mathrm{leg}}(\theta)\dot\theta^3+J_{\mathrm{leg}}(\theta)\dot\theta\ddot\theta.
\]

고관절 입력 동력은 \(\tau_h\dot\phi\), 페달로 전달되는 동력은 \(\tau_p\dot\theta\)이므로, 손실과 중력을 무시하면

\[
\tau_h k(\theta)\dot\theta-\tau_p(\theta)\dot\theta=\frac{dT_{\mathrm{leg}}}{dt}.
\]

\(\dot\theta\neq0\)에서 정리하면

\[
\boxed{\tau_p(\theta)=k(\theta)\tau_h-J_{\mathrm{leg}}(\theta)\ddot\theta-\frac12J'_{\mathrm{leg}}(\theta)\dot\theta^2}
\tag{16}
\]

이다. 따라서 \(k(\theta)\neq0\)에서

\[
\boxed{\tau_h(t)=\frac{\tau_p(\theta(t))+J_{\mathrm{leg}}(\theta(t))\ddot\theta(t)+\frac12J'_{\mathrm{leg}}(\theta(t))\dot\theta(t)^2}{k(\theta(t))}}
\tag{17}
\]

이다. \(J'_{\mathrm{leg}}=dJ_{\mathrm{leg}}/d\theta\)이다.

## 7. 체인링 형상을 대입한 최종식

식 (1), (3)을 식 (17)에 대입하면

\[
\boxed{\begin{aligned}
\tau_h(t)=\frac{1}{k(\theta(t))}\Bigg[&\tau_p(\theta(t))\\
&-\frac{(r_s\omega)^2J_{\mathrm{leg}}(\theta(t))\rho'(\theta(t))}{\rho(\theta(t))^3}\\
&+\frac{(r_s\omega)^2J'_{\mathrm{leg}}(\theta(t))}{2\rho(\theta(t))^2}\Bigg].
\end{aligned}}
\tag{18}
\]

따라서 고관절 토크는 \(\tau_p\), \(k\), \(J_{\mathrm{leg}}\), \(J'_{\mathrm{leg}}\), \(\rho\), \(\rho'\), \(r_s\), \(\omega\)에 의해 결정된다.

## 8. 수치계산

\(J'_{\mathrm{leg}}\)는 중앙차분으로 근사할 수 있다.

\[
\boxed{J'_{\mathrm{leg}}(\theta)\approx\frac{J_{\mathrm{leg}}(\theta+h)-J_{\mathrm{leg}}(\theta-h)}{2h}}
\tag{19}
\]

\(\theta(t)\)는 다음 초기값 문제를 `ode45`로 계산할 수 있다.

\[
\boxed{\begin{cases}\displaystyle\frac{d\theta}{dt}=\frac{r_s\omega}{\rho(\theta)},\\[4pt]\theta(0)=0.\end{cases}}
\tag{20}
\]

수치적으로는 \(\theta(t)\)를 구한 후 \(\phi(\theta)\), \(k(\theta)\), \(J_{\mathrm{leg}}(\theta)\), \(J'_{\mathrm{leg}}(\theta)\)를 평가하여 식 (17) 또는 (18)에 대입한다.

## 9. 적용 범위 및 남은 과제

- 식 (17), (18)은 \(k(\theta)\neq0\)인 구간에서 직접 적용할 수 있다.
- 활성 다리 교대는 별도의 규칙으로 처리해야 하며, 교대 지점의 불연속성에 주의한다.
- 모델에서 구한 최적화 결과는 설정한 가정 아래의 결과이며, 실험으로 검증된 실제 자전거 성능과 동일시할 수 없다.
- 최적 체인링을 정하려면 최대 토크, 제곱평균 토크 등 최소화할 목적함수를 별도로 지정해야 한다.
