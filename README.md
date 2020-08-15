## 基本情報(sensitivity_toolについて)
- アンサンブル手法に基づく簡易予報感度解析(Enomoto et al. 2015; 榎本ほか 2014; Matsueda et al. 2011)

#### 基本となる考え方

場の状態ベクトル`x`の時間発展が非線形モデル`M(x)`で記述されたとする。初期時刻`t=0`における摂動`y(i=1:m)`を与えた、メンバー数`m`のアンサンブル予報において、時刻`t`における擾乱は次のように表せる。

```
z(i=1:m,t) = M(x+y(i=1:m),t) - M(x,t)
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

<br>

## コードの紹介
### アンサンブル随伴感度解析(Ensemble adjoint sensitivity analysis; EnASA)  
- `anl_EnASA_rate.py` 

`z(i=1:m)`を用いて**トータルエネルギーノルム`norm(i=1:m)`** を計算する。各メンバーから求めた`norm(i=1:m)`を用いて最適な係数`p(i=1:m)`を求める。

```python
for i in range(m):
  p(i) = norm(i)/np.sum(norm(1:m))
```

<br>

### アンサンブル特異ベクトル感度解析(Ensemble singular vector analysis; EnSVSA)  
- `anl_EnSVSA_mode_svds.py`

アンサンブル特異ベクトル法では、共分散`(y.T)G(y) = (p.T Y.T)G(Y p) = 1`の条件のもとで検証時刻における検証領域における擾乱`(p.T Z.T) H (Z)`を最大化する`p(i=1:m)`を求める。この問題ではラグランジュ関数からの微分から固有値問題を得る。

<br>

ラグランジュ 関数を次のように定義する。  
![F(\boldsymbol{p}, \lambda)=\boldsymbol{p}^{\top} \boldsymbol{Z}^{\top} \mathbf{G}_{t} \boldsymbol{Z} \boldsymbol{p}+\lambda\left(1-\boldsymbol{p}^{\top} \mathbf{Y}^{\top} \mathbf{G}_{0} \mathbf{Y} \boldsymbol{p}\right)](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+F%28%5Cboldsymbol%7Bp%7D%2C+%5Clambda%29%3D%5Cboldsymbol%7Bp%7D%5E%7B%5Ctop%7D+%5Cboldsymbol%7BZ%7D%5E%7B%5Ctop%7D+%5Cmathbf%7BG%7D_%7Bt%7D+%5Cboldsymbol%7BZ%7D+%5Cboldsymbol%7Bp%7D%2B%5Clambda%5Cleft%281-%5Cboldsymbol%7Bp%7D%5E%7B%5Ctop%7D+%5Cmathbf%7BY%7D%5E%7B%5Ctop%7D+%5Cmathbf%7BG%7D_%7B0%7D+%5Cmathbf%7BY%7D+%5Cboldsymbol%7Bp%7D%5Cright%29)

<br>

pについて微分すると、

![\frac{\partial F(\boldsymbol{p}, \lambda)}{\partial \boldsymbol{p}}=2 \boldsymbol{p}^{\top} \mathbf{Z}^{\top} \mathbf{G}_{t} \mathbf{Z}-2 \lambda \boldsymbol{p}^{\top} \mathbf{Y}^{\top} \mathbf{G}_{0} \mathbf{Y}=\mathbf{0}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cfrac%7B%5Cpartial+F%28%5Cboldsymbol%7Bp%7D%2C+%5Clambda%29%7D%7B%5Cpartial+%5Cboldsymbol%7Bp%7D%7D%3D2+%5Cboldsymbol%7Bp%7D%5E%7B%5Ctop%7D+%5Cmathbf%7BZ%7D%5E%7B%5Ctop%7D+%5Cmathbf%7BG%7D_%7Bt%7D+%5Cmathbf%7BZ%7D-2+%5Clambda+%5Cboldsymbol%7Bp%7D%5E%7B%5Ctop%7D+%5Cmathbf%7BY%7D%5E%7B%5Ctop%7D+%5Cmathbf%7BG%7D_%7B0%7D+%5Cmathbf%7BY%7D%3D%5Cmathbf%7B0%7D)

したがって、下記の固有値問題に置き換えることができる。

![\left(\mathrm{Y}^{\mathrm{T}} \mathrm{GY}\right)^{-1} \mathrm{Z}^{\mathrm{T}} \mathrm{HZp}=\Lambda \boldsymbol{p}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cleft%28%5Cmathrm%7BY%7D%5E%7B%5Cmathrm%7BT%7D%7D+%5Cmathrm%7BGY%7D%5Cright%29%5E%7B-1%7D+%5Cmathrm%7BZ%7D%5E%7B%5Cmathrm%7BT%7D%7D+%5Cmathrm%7BHZp%7D%3D%5CLambda+%5Cboldsymbol%7Bp%7D)

はじめに、`Y.T G Y`は正規直交関数で対角行列になるので、`Y.T G Y`は考えなくて良い。`Z.T G Z`については2つの解法が考えられる。   

1. `Z.T G Z`の固有値問題として解く。 
`Z.T G Z`の行列サイズを確認してみると、`(m, dims) (dims, dims) (dims, m) = (m, m)`とメンバー数`m`[~O(10)]と等しくなり、簡単に固有値問題を解くことができる。

2. `Z`の特異値問題として解く。  
`Z.T G Z`の固有値問題を解く方法には、`Z`の特異値問題と置き換えることができる。Zは次のように分解することができる。
![Z=U \Sigma V^{\top}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+Z%3DU+%5CSigma+V%5E%7B%5Ctop%7D)

左特異値ベクトル`U`が共分散`Z Z.T`の固有ベクトル、右特異値ベクトル`V.T`が共分散`Z.T Z`の正規化された主成分の固有値ベクトルを表す。行列`sigma`の対角成分は特異値である。この時、共分散`Z.T Z`の正規化された主成分の固有値ベクトルが`p(i=1:m)`に相当する。  

1,2で作成した`p(i=1:m)`を下記のように初期場にかけて感度領域を作成する。固有値ベクトルの行列サイズが`(m, m)`であることに注意する。  
![\mathbf{y}=\mathbf{Y} \mathbf{p}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cmathbf%7Bx%7D%3D%5Cmathbf%7BY%7D+%5Cmathbf%7Bp%7D)

#### 実践的な手順

1. **検証領域**における検証時刻の摂動`z(i=1:m)`から抽出した`extraction_z（dims=1:ndims,i=1:m）`を作成する。

2. `extraction_z（dims=1:ndims,i=1:m）`で特異値分解して、`p(i=1:m)`を取得する。

<br>

### 計算時の注意点

- 摂動はアンサンブル平均からの差ではなく、コントロールランからの差で構成される。

- ノルムの計算は**検証領域**によって変化することが知られている。

<br>
<br>

## データ取得
`data_get.sh` : grib形式のデータを取得するコード(※1)

[京都大学生存圏研究所（RISH: Research Institute for Sustainable Humanosphere）](http://database.rish.kyoto-u.ac.jp/arch/jmadata/gpv-original.html)から取得することができる。詳しい格子情報や配信時間等は[気象行支援センター](http://database.rish.kyoto-u.ac.jp/arch/jmadata/gpv-original.html)に記載されている。

`data_grib2bin.sh, data_encode.py` : データ形式をGrads形式に変更するコード(※2)

<br>

#### データにおける注意点
※1. データ取得サイト利用時は、下記の点の注意が必要。
>教育研究機関向けにデータを提供しています。企業活動等のためにデータを頻繁に必要とされる方は、気象業務支援センターからデータを直接購入し、データ提供スキーム全体の維持発展にご協力ください。（サイトより引用）

※2. 今回使用しているデータは週間予報アンサンブルGPVである。

2006030100以前はgrib形式、それ以降はgrib2形式で提供されている。これに伴いメンバー数が(25→27), 提供データの要素が一部変更されているので注意が必要である。
なお、 2020(令和2)年3月24日（火）以降のデータ提供は終了している。

>2020(令和2)年3月24日（火）を以て提供を終了しました。
>※最終提供プロダクトは、2020(令和2)年3月24日（火）12UTC初期値です。高分解能全球域GPV、高分解能日本域GPVをご利用ください。(サイトより引用)

<br>
<br>

## 作成情報
- 制作開始日　2020年7月5日

## 今後の計画
1. 参考文献の再現実験を行う。　2. 近年の大気場でも同様の実験で行ってみる。 

## 参考文献
Enomoto, T., S. Yamane, and W. Ohfuchi, 2015: Simple sensitivity analysis using ensemble forecasts. J. Meteor. Soc. Japan, 93, 199-	213.  
榎本剛, 山根省三, 大淵済, 2014: アンサンブル手法に基づく簡易予報感度解析. 京都大学防災研究所年報, 57(B), 163-168. 
