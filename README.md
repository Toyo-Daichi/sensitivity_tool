## 基本情報(sensitivity_toolについて)
- アンサンブル手法に基づく簡易予報感度解析

## コードの紹介
`anl_EnASA_rate.py` : アンサンブル随伴感度解析

`anl_EnSVSA_mode_svds.py` : アンサンブル特異ベクトル感度解析

## データ取得
`data_get.sh` : grib形式のデータを取得するコード(※1)

[京都大学生存圏研究所（RISH: Research Institute for Sustainable Humanosphere）](http://database.rish.kyoto-u.ac.jp/arch/jmadata/gpv-original.html)から取得することができる。詳しい格子情報や配信時間等は[気象行支援センター](http://database.rish.kyoto-u.ac.jp/arch/jmadata/gpv-original.html)に記載されている。

`data_grib2bin.sh, data_encode.py` : データ形式をGrads形式に変更するコード(※2)

### 注意点
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
