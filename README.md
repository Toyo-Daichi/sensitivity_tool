## 基本情報(sensitivity_toolについて)
- アンサンブル手法に基づく簡易予報感度解析(Enomoto et al. 2015; 榎本ほか 2014)

### 基本となる考え方

場の状態ベクトル`x`の時間発展が非線形モデル`M(x)`で記述されたとする。初期時刻`t=0`における摂動`y(i=1:m)`を与えた、メンバー数`m`のアンサンブル予報において、時刻`t`における擾乱は次のように表せる。

```
z(i=1:m,t) = M(x+y(i=1:m)) - M(x)
```

**初期摂動の線形発展を仮定すると、"感度解析"は「検証領域において検証時刻に成長を表すメンバーの線形組み合わせの最適な係数`p(i=1:m)`を求めること」である。**

```
z(optimum) = p(i=1)z(i=1) + p(i=２)z(i=2) +　... + p(i=m)z(i=m)
p(i=1:m)   = (p(i=1), p(i=2), ..., p(i=m)).T
```
**同じ係数を用いて、発達する摂動に対応する初期摂動が求まる。**

```
y(optimum) = p(i=1)y(i=1) + p(i=２)y(i=2) +　... + p(i=m)y(i=m)
p(i=1:m)   = (p(i=1), p(i=2), ..., p(i=m)).T
```

## コードの紹介
- `anl_EnASA_rate.py` : アンサンブル随伴感度解析(Ensemble adjoint sensitivity analysis; EnASA)

`z(i=1:m)`を用いて**トータルエネルギーノルム`norm(i=1:m)`**(※1)を計算する。各メンバーから求めた`norm(i=1:m)`を用いて最適な係数`p(i=1:m)`を求める。

```python
for i in range(m):
  p(i) = norm(i)/np.sum(norm(1:m))
```

- `anl_EnSVSA_mode_svds.py` : アンサンブル特異ベクトル感度解析(Ensemble singular vector analysis; EnSVSA)

アンサンブル特異ベクトル法では、共分散`(y.T)G(y) = (p.T Y.T)G(Y p) = 1`の条件のもとで検証時刻における検証領域における擾乱`(p.T Z.T) H (Z)`を最大化する`p(i=1:m)`を求める。この問題ではラグランジュ関数からの微分から固有値問題を得る。

<img src="https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cleft%28%5Cmathrm%7BY%7D%5E%7B%5Cmathrm%7BT%7D%7D+%5Cmathrm%7BGY%7D%5Cright%29%5E%7B-1%7D+%5Cmathrm%7BZ%7D%5E%7B%5Cmathrm%7BT%7D%7D+%5Cmathrm%7BHZp%7D%3D%5CLambda+%5Cboldsymbol%7Bp%7D" 
alt="\left(\mathrm{Y}^{\mathrm{T}} \mathrm{GY}\right)^{-1} \mathrm{Z}^{\mathrm{T}} \mathrm{HZp}=\Lambda \boldsymbol{p}">

この時、`Y.T G Y`は互いに直交している行列とその転置行列を掛け合う(?)ので対角行列となる。したがって、`Z`の特異値問題に置き換えることができる。

<img src="https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+Z%3DU+%5CSigma+V%5E%7B%5Ctop%7D" 
alt="Z=U \Sigma V^{\top}">

この時、左特異値ベクトル`U`が共分散`Z.T H Z`の固有ベクトル、すなわち`p(i=1:m)`に相当する。

1. **検証領域**における検証時刻の摂動`z(i=1:m)`から抽出した`extraction_z（dims=1:ndims,i=1:m）`を作成する。

2. `extraction_z（dims=1:ndims,i=1:m）`で特異値分解して、`p(i=1:m)`を取得する。

### 計算時の注意点

- 摂動はアンサンブル平均からの差ではなく、コントロールランからの差で構成される。

- ノルムの計算は**検証領域**によって変化することが知られている。

## データ取得
`data_get.sh` : grib形式のデータを取得するコード(※1)

[京都大学生存圏研究所（RISH: Research Institute for Sustainable Humanosphere）](http://database.rish.kyoto-u.ac.jp/arch/jmadata/gpv-original.html)から取得することができる。詳しい格子情報や配信時間等は[気象行支援センター](http://database.rish.kyoto-u.ac.jp/arch/jmadata/gpv-original.html)に記載されている。

`data_grib2bin.sh, data_encode.py` : データ形式をGrads形式に変更するコード(※2)

### データにおける注意点
※1. データ取得サイト利用時は、下記の点の注意が必要。
>教育研究機関向けにデータを提供しています。企業活動等のためにデータを頻繁に必要とされる方は、気象業務支援センターからデータを直接購入し、データ提供スキーム全体の維持発展にご協力ください。（サイトより引用）

※2. 今回使用しているデータは週間予報アンサンブルGPVである。

2006030100以前はgrib形式、それ以降はgrib2形式で提供されている。これに伴いメンバー数が(25→27), 提供データの要素が一部変更されているので注意が必要である。
なお、 2020(令和2)年3月24日（火）以降のデータ提供は終了している。

>2020(令和2)年3月24日（火）を以て提供を終了しました。
>※最終提供プロダクトは、2020(令和2)年3月24日（火）12UTC初期値です。高分解能全球域GPV、高分解能日本域GPVをご利用ください。(サイトより引用)


## 作成情報
- 制作開始日　2020年7月5日

## 今後の計画
1. 参考文献の再現実験を行う。　2. 近年の大気場でも同様の実験で行ってみる。 

## 参考文献
Enomoto, T., S. Yamane, and W. Ohfuchi, 2015: Simple sensitivity analysis using ensemble forecasts. J. Meteor. Soc. Japan, 93, 199-	213.

榎本剛, 山根省三, 大淵済, 2014: アンサンブル手法に基づく簡易予報感度解析. 京都大学防災研究所年報, 57(B), 163-168. 
