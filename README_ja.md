<a name= "japanese">

日本語 / [English](README.md)

<img src="doc/figs/komega.png" width="300">

## プロジェクト情報
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran-orange.svg)](https://fortran-lang.org/)
[![Platform: Linux](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)](https://github.com/k-yoshimi/Komega)
[![Version: 2.0.0](https://img.shields.io/badge/Version-2.0.0-green.svg)](https://github.com/k-yoshimi/Komega/releases)

## ビルドステータス
[![CI](https://github.com/k-yoshimi/Komega/workflows/CI/badge.svg)](https://github.com/k-yoshimi/Komega/actions)
[![Extended Test Matrix](https://github.com/k-yoshimi/Komega/workflows/Extended%20Test%20Matrix/badge.svg)](https://github.com/k-yoshimi/Komega/actions)
[![Quick Test](https://github.com/k-yoshimi/Komega/workflows/Quick%20Test/badge.svg)](https://github.com/k-yoshimi/Komega/actions)

# What is Kω ? 

Shifted-Krylov部分空間法に基づくソルバーライブラリと,
それを用いてHamiltonianと励起状態ベクトルから動的Green関数を計算するミニアプリである.

# Download

 * 最新版のパッケージはこちらから.
 
   https://github.com/k-yoshimi/Komega/releases/download/v2.0.0/komega-2.0.0.tar.gz

 * リリースノート
 
   https://github.com/k-yoshimi/Komega/releases/tag/v2.0.0
 * 過去のリリース
 
   https://github.com/k-yoshimi/Komega/releases

# Prerequisite

 * fortranコンパイラ
 * BLASライブラリ
 * LAPACKライブラリ(ミニアプリのみ使用)
 * MPIライブラリ(Optional)

# Documents

 * Manual for the Library
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/library/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/library/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/en/_build/latex/komega.pdf))
 * Manual for the sample program
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/software/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/ja/_build/latex/shiftk.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/software/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/en/_build/latex/shiftk.pdf))

# License

GNU LESSER GENERAL PUBLIC LICENSE Version 3

© 2016- The University of Tokyo. All rights reserved.

本ソフトウェアは2016年度 東京大学物性研究所 ソフトウェア高度化プロジェクトの支援を受け開発された。

(＊) 引用論文

[“Kω — Open-source library for the shifted Krylov subspace method of the form (zI−H)x=b”, Takeo Hoshi, Mitsuaki Kawamura, Kazuyoshi Yoshimi, Yuichi Motoyama, Takahiro Misawa, Youhei Yamaji, Synge Todo, Naoki Kawashima, Tomohiro Sogabe, Computer Physics Communications, Volume 258, January 2021, 107536.](https://www.sciencedirect.com/science/article/pii/S0010465520302551)

# Directory Tree

 * `app/`: ミニアプリ関連のディレクトリ
   * `src/`: ミニアプリソースコードのディレクトリ
   * `sample/`: ミニアプリサンプル用ディレクトリ
     * `Shiftk.nb`: テスト用Mathematicaノートブック(開発者向け)
     * `denovo/`: ハミルトニアンや右辺ベクトルの入力ファイルを用いない場合の例
     * `from_file/`: ハミルトニアンや右辺ベクトルの入力ファイルを用いる場合の例
 * `doc/`: ドキュメント用ディレクトリ
   * `index.html` : ドキュメントのトップページ
   * `library/`: ライブラリのドキュメントのディレクトリ
   * `software/`: ミニアプリのドキュメントのディレクトリ
 * `configure`: ビルド環境指定ファイル
 * `src/`: ライブラリ本体のディレクトリ
 * `test/`: ライブラリのテスト用ディレクトリ

# Install

もっともシンプルには次のとおりである.

 * コンソールで `$ ./configure --prefix=install_dir; make; make install` とタイプする.
   ただし`install_dir`はインストール先のディレクトリの絶対パスとする
   (以後各自のディレクトリ名で読み替えること).
 * `install_dir`で指定したディレクトリに以下のものが作られる.
   * `install_dir/lib/`内: 共有および静的ライブラリ
   * `install_dir/include/komega.h`: C/C++用ヘッダーファイル
   * `install_dir/bin/Shiftk.out`: ミニアプリ

詳しくはマニュアルを参照のこと.

# ミニアプリのテスト

 * `app/sample/denovo/`もしくは`app/sample/from_file/`ディレクトリに移動する.
 * `install_dir/bin/ShiftK.out namelist.def`とやる.
 * `dynamicalG.dat`などが作られると成功.
 * `namelist.def`の書式はマニュアル参照

# ライブラリの使用方法

## プログラム内での各ルーチンの呼び出し方

### fortran/C/C++の場合

マニュアル参照

## ライブラリのリンク方法

``` bash
$ gfortran myprog.f90 -L install_dir/src/shared -lshiftk -lblas -I install_dir/src/shared
$ gcc myprog.c -L install_dir/src/shared -lshiftk -lblas -I install_dir/src/shared
```
など.

動的リンクを行ったファイルを実行するときには,
環境変数`LD_LIBRARY_PATH`に`install_dir/lib`ディレクトリを追加しておく必要がある.

# ライブラリのテスト(Optional)

 * `test/`ディレクトリに移動する.
 * `$ ./solve_cc.x < complex_freq.in`とやる.
 * Residual vectorが十分小さくなっていれば成功
 * `solve_rc.x` も同様. `solve_cr.x`, `solve_rr.x` に関しては入力ファイル `real_freq.in` を使う.
 * `complex_freq.in`(名称は自由)のパラメーターは次の通り
   * `ndim`: 擬似ハミルトニアンの次元
   * `nl`: 射影のテスト用. 解ベクトルの先頭から`nl(<=ndim)`番目までを計算する.
   * `nz`: 振動数の点数
   * `itermax`: 最大反復回数
   * `threshold`: 収束判定のthreshold
   * `rnd_seed`: 擬似Hamiltonian生成のための乱数の種
   * このnamelistの後ろに振動数を1行に1個ずつ書く.
     
