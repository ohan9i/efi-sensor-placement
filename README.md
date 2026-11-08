# EfI 센서 배치 시뮬레이터

Kim & Cho (2024), *Effective independence in optimal sensor placement associated with general Fisher information involving full error covariance matrix*, Mechanical Systems and Signal Processing 212, 111263 의 내용을 직접 만져 보며 이해하기 위한 인터랙티브 웹 시뮬레이션입니다.

**라이브:** https://ohan9i.github.io/efi-sensor-placement/

## 구성
1. **핵심 개념** — FIM $F=S^T\Sigma^{-1}S$, EfI 정의 $e_i = 1-\det F_{-i}/\det F$, 기존 근사식과 제안식
2. **3센서 예제 (논문 Table 1)** — S, Σ(분산·상관), 분해 방식 L을 바꿔 가며 정의식 / 제안식 Eq.(26) / 기존 근사 Eq.(18) 비교
3. **연속보 케이스 스터디 (논문 4장)** — 4경간 연속보 FE 모델(81개 후보, 10개 모드), Eq.(28)/(29) 상관 오차, 81→10 순차 제거를 4가지 방법(Alt.1–4)으로 실행해 엔트로피 −ln det F, 배치, 계산시간, 단계별 EfI 비교
4. **부록 A 증명 요약**

모든 계산은 브라우저에서 순수 JavaScript로 수행합니다 (`js/linalg.js`, `js/beam.js`, `js/efi.js`).
경간 길이(8–12–12–8)는 논문 그림의 총 길이 40, 최장 경간 12로부터 가정한 값입니다.

## 로컬 실행
```bash
python3 -m http.server 8000   # http://localhost:8000
node test/run.mjs             # 논문 Table 1 수치와 케이스 스터디 결과 확인
```
