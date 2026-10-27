#!/bin/sh
# Build the exact MILP certificate checker viprchk (github.com/ambros-gleixner/VIPR at commit f221c3f) against
# GMP 6.3.0, with m4 1.4.19 built only to configure GMP. Sources are fetched from GNU and GitHub and their sha256
# checked. Everything goes under the directory given as $1 (default ./build); prints the path of viprchk.
#   needs: curl, git, tar, xz, make, a C and a C++ compiler
set -e
B=${1:-build}
mkdir -p "$B" && cd "$B" && B=$(pwd)
curl -fsSLO https://ftp.gnu.org/gnu/m4/m4-1.4.19.tar.xz
curl -fsSLO https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
printf '%s\n' '63aede5c6d33b6d9b13511cd0be2cac046f2e70fd0a07aa9573a04a82783af96  m4-1.4.19.tar.xz' \
  'a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898  gmp-6.3.0.tar.xz' | sha256sum -c
tar xf m4-1.4.19.tar.xz && tar xf gmp-6.3.0.tar.xz
(cd m4-1.4.19 && ./configure --prefix="$B/m4" && make && make install) > m4-build.log 2>&1
PATH="$B/m4/bin:$PATH"; export PATH
(cd gmp-6.3.0 && ./configure --prefix="$B/gmp" --enable-cxx && make && make install) > gmp-build.log 2>&1
rm -rf VIPR
git clone -q https://github.com/ambros-gleixner/VIPR && git -C VIPR checkout -q f221c3fd2274e8e859e8e24505278d6f7c1bb8ff
g++ -O2 -std=c++11 -I"$B/gmp/include" VIPR/code/viprchk.cpp -o viprchk -L"$B/gmp/lib" -lgmpxx -lgmp \
  -Wl,-rpath,"$B/gmp/lib"
echo "$B/viprchk"
