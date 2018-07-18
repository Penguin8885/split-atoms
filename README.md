# split_atoms
原子集団(データ群)を慣性モーメントテンソルで軸を決めて修正再帰二分法で分割していくプログラム

### プログラムの内容
大規模データをMPIに分割するとき，空間局所性を維持してほぼ等分割するために開発された手法  
主に分子シュミレーションで用いられる．質量パラメータを利用することで各データに重み付けをすることができる．  

この手法では原子集団をある都合がいい軸に射影して，その軸上の値で原子集団に順序関係をつけて，2つのグループに分割する．
この軸に慣性モーメントテンソルの固有ベクトルを用いる(目的関数の最大化で導かれる)．  
グループ分けしたものをさらに修正再帰二分法を用いて設定された分割数まで分割していく．
修正再帰二分法は再帰二分法を拡張した方法で，再帰二分法が2のべき乗の数値のみ使用できるのを一般の数字で使用できるように拡張した方法．

一様分布2次元データに対する結果(1013原子, 11分割)
<img src="https://github.com/Penguin8885/split_atoms/blob/master/figure_0.png">

正規分布2次元データに対する結果(1013原子, 11分割)
<img src="https://github.com/Penguin8885/split_atoms/blob/master/figure_1.png">

### 元論文
A three-dimensional domain decomposition method for large-scale DFT electronic structure calculations  
https://arxiv.org/abs/1209.4506  
