# EfI 센서 배치 시뮬레이터

Kim & Cho (2024), *Effective independence in optimal sensor placement associated with general Fisher information involving full error covariance matrix*, Mechanical Systems and Signal Processing 212, 111263 의 내용을 직접 만져 보며 이해하기 위한 인터랙티브 웹 시뮬레이션입니다.

**라이브:** https://ohan9i.github.io/efi-sensor-placement/

## 배치 게임 (첫 화면)
논문 그림 1처럼 층마다 센서를 하나씩 둘 수 있는 전단 건물에서, 예산만큼 층을 골라 센서를 꽂습니다. det F 기반 정보 효율이 실시간 점수가 되고, 제출하면 논문의 정확한 EfI(Eq. 26)로 순차 제거한 봇, 촐레스키 근사(Eq. 18) 봇, 교환 탐색 최고 기록과 비교합니다.

- 스테이지 1: 9층 건물, 모드 3개, 센서 3개 (Σ = I)
- 스테이지 2: 20층 건물, 모드 4개, 센서 6개 (Eq. 28, λ = 3층)
- 스테이지 3: 30층 건물, 모드 5개, 센서 8개 (Eq. 28, λ = 4층)

## 원리 이해 (두 번째 탭)
왼쪽은 여섯 장의 설명, 오른쪽은 스크롤을 따라 바뀌는 하나의 시뮬레이션 스테이지입니다.

1. **81개 후보, 10개의 센서**: 4경간 연속보 FE 모델(Eq. 27)과 목표 모드 10개
2. **센서 하나의 몫, EfI**: 위치별 EfI와 다음에 제거될 센서
3. **오차가 서로 닮을 때**: Eq. 28/29 공분산, 상관 길이 λ, 센서 간 상관
4. **지름길의 함정**: 논문 Table 1의 3센서 예제, 연속보에서 네 방법이 고르는 센서 비교
5. **하나씩 지워 보기**: 81→10 순차 제거, 엔트로피 −ln det F, 배치, 계산 시간
6. **제안식이 정확한 이유**: 부록 A 증명 요약

디자인은 [unslop](https://github.com/heydamir/unslop) (Damir Sokolovsky, CC BY 4.0) 의 기준으로 다듬었습니다. 폰트는 IBM Plex Sans KR / IBM Plex Mono (OFL), 수식은 KaTeX를 함께 배포합니다.

모든 계산은 브라우저에서 순수 JavaScript로 수행합니다 (`js/linalg.js`, `js/beam.js`, `js/efi.js`).
경간 길이(8–12–12–8)는 논문 그림의 총 길이 40, 최장 경간 12로부터 가정한 값입니다.

## 로컬 실행
```bash
python3 -m http.server 8000   # http://localhost:8000
node test/run.mjs             # 논문 Table 1 수치와 케이스 스터디 결과 확인
```
