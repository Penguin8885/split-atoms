import numpy as np
import matplotlib.pyplot as plt

def generate_data(num, type_n=0):
    # np.random.seed(0)

    if type_n == 0:
        # 一様分布でデータを発生
        x = np.random.uniform(-100, 100, num)
        y = np.random.uniform(-100, 100, num)
    else:
        # 正規分布でデータを発生
        mu = [10, 10]               # 平均
        cov = [[5, 1], [1, 5]]      # 共分散行列
        x, y = np.random.multivariate_normal(mu, cov, num).T

    z = np.zeros(num, dtype=np.int) # Z軸座標は全て0
    m = np.random.uniform(0,1,num)  # 質量パラメータは[0,1]の一様分布で生成
    data = list(zip(x, y, z, m))
    return data

def split_atoms(data, num_div):
    # 分割しなくていい(分割数1)の場合はそのままリストにして返す
    if num_div == 1:
        return [data]

    # 重心を計算
    center = np.zeros(3)
    total_mass = 0
    for datum in data:
        vec = np.array(datum[:3])   # 位置ベクトル
        mass = datum[3]             # 質量
        center += mass * vec
        total_mass += mass
    center /= total_mass            # 総和から重心を計算

    # 慣性モーメントテンソルを計算
    A = np.zeros([3,3])             # 慣性モーメントテンソル
    for x, y, z, m in data:
        x -= center[0]; y -= center[1]; z -= center[2]  # 重心を原点へシフト
        A[0][0] += m * (y**2 + z**2)
        A[0][1] -= m * x * y
        A[0][2] -= m * x * z
        A[1][0] -= m * y * x
        A[1][1] += m * (z**2 + x**2)
        A[1][2] -= m * y * z
        A[2][0] -= m * z * x
        A[2][1] -= m * z * y
        A[2][2] += m * (x**2 + y**2)

    # 慣性テンソルの固有値, 固有値ベクトルを計算して軸を決定
    eig_lam, eig_vec = np.linalg.eig(A)         # 固有値，固有ベクトル計算

    # 慣性テンソルから求めた3つの固有ベクトル方向の軸へ射影
    projection = np.empty([3,len(data)])
    objective = np.zeros(3)
    for i, ev in enumerate(eig_vec):            # 各固有ベクトルでループ
        for j, datum in enumerate(data):        # 各データでループ
            vec = np.array(datum[:3], np.float)           # 位置ベクトル
            mass = datum[3]                               # 質量
            vec -= center                                 # 重心を原点へシフト
            projection[i][j] = ev.dot(vec)                # 軸への射影値を計算
            objective[i] += mass * (projection[i][j])**2  # 目的関数を計算
    # print(objective)

    # 適切な軸を選んで射影値でソート
    max_index = np.argmax(objective)            # 目的関数を最大化する軸を選択
    ap = [(projection[max_index][i], datum) for i, datum in enumerate(data)]
                                                # (選択した軸への射影値, 元データ)の
                                                #             タップルのリストを作成
    ap = sorted(ap, key=lambda array:array[0])  # 射影値でソート

    # 修正再帰二分法の割合に合わせてデータ分割
    next_num_div1 = int(np.ceil(num_div/2))     # 次の分割数 1
    next_num_div2 = num_div - next_num_div1     # 次の分割数 2
    split_index = len(ap)*next_num_div1/num_div # 分割する位置を次の分割の割合から決定
    split_index = int(np.round(split_index))    # 四捨五入して整数型へ変換
    data_1 = [d for p, d in ap[:split_index]]   # 分割したデータ 1
    data_2 = [d for p, d in ap[split_index:]]   # 分割したデータ 2

    # 再帰してさらに分割させ，最後に統治
    s_1 = split_atoms(data_1, next_num_div1)    # 再帰分割
    s_2 = split_atoms(data_2, next_num_div2)    # 再帰分割
    split_list = s_1 + s_2                      # 統治

    return split_list

if __name__ == '__main__':
    # データ生成
    data = generate_data(num=1013, type_n=1)

    # 修正再帰二分法と慣性モーメントテンソルを用いた分割
    split_list = split_atoms(data, num_div=11)

    # 結果を表示
    for i, split_i in enumerate(split_list):
        x, y, z, m = zip(*split_i)  # データをx, y, zのリストに分割
        plt.scatter(x, y)           # x, y座標で散布図をプロット
        print(i, len(split_i))      # i番目のプロセスにわたすデータサイズを表示
    plt.show()
