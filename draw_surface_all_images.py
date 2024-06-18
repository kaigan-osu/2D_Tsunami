import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

mx = 750      # 東西方向の格子数
my = 495      # 南北方向の格子数
dx = 1620.0    # 格子間隔 [m]
dy = 1620.0    # 格子間隔 [m]

# 読込データのファイルリストを読み込む
case = np.loadtxt('flist.csv',delimiter=',', dtype='str')

for kt in range(len(case)):
    print(case[kt],' in process')
    eta = np.loadtxt('data/'+ case[kt], delimiter=',')
    etanew = np.flipud(eta[1:my+1, 1:mx+1])

    # 図形描画の準備
    fig = plt.figure(figsize=(11.2, 6))  # 図の枠の準備
    ax = fig.add_subplot(111)

    x = [dx * i / 1000.0 for i in range(mx)]    # x座標の値
    y = [dy * j / 1000.0 for j in range(my)]    # y座標の値
    X, Y = np.meshgrid(x, y)

    # 結果の出力
    ax.set_xlabel("X-axis [km]")
    ax.set_ylabel("Y-axis [km]")

    # カラーバーの範囲を設定
    vmin, vmax = -3.0, 3.0

    # 等高線の間隔を細かく設定
    levels = np.linspace(vmin, vmax, 100)

    # BoundaryNormを使用して、範囲外の色を指定
    cmap = plt.get_cmap('bwr')
    norm = BoundaryNorm(np.linspace(vmin, vmax, 101), cmap.N, clip=True)

    map = ax.contourf(X, Y, etanew, cmap=cmap, levels=levels, norm=norm,extend='both')   # 水位分布をプロット
    fig.colorbar(map, ax=ax, label='elevation [m]', ticks=np.arange(-3, 3.1, 0.5))

    # fig内でのaxes座標を取得，戻り値はBbox
    ax_pos = ax.get_position()
    # fig内座標でテキストを表示 Bboxは Bbox.x0, Bbox.x1, Bbox.y0, Bbox.y1で座標を取得できる
    fig.text(ax_pos.x1 - 0.11, ax_pos.y1 + 0.01, 'file= ' + case[kt])

    h = np.loadtxt('depth_1620-01.csv', delimiter=',')                # 水深データの読み込み
    hnew = np.flipud(h)
    ax.contour(X, Y, hnew, levels=[0.01], colors='black', linewidths=0.5)   # 地形の形状をプロット

    plt.savefig('images/'+ case[kt][:4] +'.jpg')
    plt.cla()

    #plt.show()