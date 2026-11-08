#!/bin/sh
# 배포 전 실행: 모든 JS/CSS 경로에 버전을 붙여 브라우저 캐시 때문에 옛 화면이 보이는 문제를 막음
V=${1:-$(date +%s)}
sed -i -E "s#(\./[a-z]+\.js)(\?v=[0-9]+)?'#\1?v=$V'#g" js/*.js
sed -i -E "s#(src|href)=\"((js|css)/[a-z]+\.(js|css))(\?v=[0-9]+)?\"#\1=\"\2?v=$V\"#g" index.html
