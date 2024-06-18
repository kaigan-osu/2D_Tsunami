import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

# 水深データと初期水位データを図示

eta = np.loadtxt('surface_1620.csv', delimiter=',')
h = np.loadtxt('depth_1620-01.csv', delimiter=',')       

eta = np.flipud(eta)
h = np.flipud(h)

mx = 750     # 東西方向の格子数
my = 495      # 南北方向の格子数
dx = 1620.0    # 格子間隔 [m]
dy = 1620.0    # 格子間隔 [m]

# 図形描画の準備
fig = plt.figure(figsize=(12, 12/mx*my))  # 図の枠の準備
ax = fig.add_subplot(111)

x = [dx * 2 * i / 1000.0 for i in range(mx)]    # x座標の値
y = [dy * 2 * j / 1000.0 for j in range(my)]    # y座標の値
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

map = ax.contourf(X, Y, eta, cmap=cmap, levels=levels, norm=norm,extend='both')   # 水位分布をプロット
fig.colorbar(map, ax=ax, label='elevation [m]', ticks=np.arange(-3, 3.1, 0.5))
ax.contour(X, Y, h, levels=[0.01], colors='black', linewidths=0.5)   # 地形の形状をプロット

plt.show()