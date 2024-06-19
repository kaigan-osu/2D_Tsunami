import numpy as np
import math 

# 208-211, 222, 224行目のバグを修正(2024.6.19)
# 2次元非線形長波方程式による津波伝播計算（2024.6.18）

#計算条件
ktend = 3600        # 総計算回数
output_step = 30    # 計算結果の出力間隔
mx = 750            # 東西方向の計算格子数
my = 495            # 南北方向の計算格子数
g = 9.80665         # 重力加速度 [m/s2]
dt = 2.0            # 1ステップの時間 [sec]
dx = 1620.0         # 東西方向の計算格子サイズ [m]
dy = 1620.0         # 南北方向の計算格子サイズ [m]
manning = 0.025     # マニングの粗度係数（摩擦の計算）
omega = 7.29 * 10 ** (-5)   # 地球の角速度 [rad/s]
f = 2.0 * omega * math.sin(math.pi / 6.0)  # コリオリ係数（北緯30°で代表）

#初期条件（計算領域の外側にゼロの境界をすべての方向に設定する）
etanew = np.array([[0.0] * (mx+2) for i in range(my+2)])
Mnew   = np.array([[0.0] * (mx+2) for i in range(my+2)])
M      = np.array([[0.0] * (mx+2) for i in range(my+2)])
Nnew   = np.array([[0.0] * (mx+2) for i in range(my+2)])
N      = np.array([[0.0] * (mx+2) for i in range(my+2)])

Mp1    = np.array([[0.0] * (mx+2) for i in range(my+2)])
Mp2    = np.array([[0.0] * (mx+2) for i in range(my+2)])
Mp3    = np.array([[0.0] * (mx+2) for i in range(my+2)])
Mp4    = np.array([[0.0] * (mx+2) for i in range(my+2)])
dm     = np.array([[0.0] * (mx+2) for i in range(my+2)])
N4     = np.array([[0.0] * (mx+2) for i in range(my+2)])

Np1    = np.array([[0.0] * (mx+2) for i in range(my+2)])
Np2    = np.array([[0.0] * (mx+2) for i in range(my+2)])
Np3    = np.array([[0.0] * (mx+2) for i in range(my+2)])
Np4    = np.array([[0.0] * (mx+2) for i in range(my+2)])
dn     = np.array([[0.0] * (mx+2) for i in range(my+2)])
M4     = np.array([[0.0] * (mx+2) for i in range(my+2)])
C      = np.array([[0.0] * (mx+2) for i in range(my+2)])

h = np.loadtxt("depth_1620-01.csv" ,  delimiter=',')   # 水深データ　＋：海下向き，－：陸上向き
h = np.vstack((h[[0],:], h, h[[my-1],:]))              # 上下に１行ずつコピーして追加(境界条件用)
h = np.hstack((h[:,[0]], h, h[:,[mx-1]]))              # 左右に１行ずつコピーして追加(境界条件用)

eta = np.loadtxt("surface_1620.csv" ,  delimiter=',')  # 初期水位データ　＋：上向き，－：下向き
eta = np.vstack((eta[[0],:], eta, eta[[my-1],:]))      # 上下に１行ずつコピーして追加(境界条件用)
eta = np.hstack((eta[:,[0]], eta, eta[:,[mx-1]]))      # 左右に１行ずつコピーして追加(境界条件用)

for j in range(my+2):
    for i in range (mx+2):
        if h[j,i] < 0.0:                        # 陸地の水深はゼロ
            h[j,i] = 0.0
        elif h[j,i] < 10.0 and h[j,i] > 0.0:    # 水深10mよりも浅い海はすべて水深10mにセット
            h[j,i] = 10.0

        C[j,i] = (g * h[j,i])**0.5              # 波速の計算

#最初の境界条件
M[:, 1] = 0.0
M[:, mx+1] = 0.0
M[0,:] = M[1,:]
M[my+1,:] = M[my,:]

N[1, :] = 0.0
N[my+1, :] = 0.0
N[:,0] = N[:,1]
N[:,mx+1] = N[:,mx]

#出力用ファイルのリストファイルを準備
flist = np.loadtxt('file_name.csv', dtype='str')

#時間の計算--------------------------------------------------------------------------
for kt in range (ktend):

    #水位の計算
    emax = -100.0
    emin = +100.0
    for i in range (1,mx+1):
        for j in range(1,my+1):
            if h[j,i] > 0.0:
                etanew[j,i] = eta[j,i] - dt * ((M[j,i+1] - M[j,i]) / dx + (N[j+1,i] - N[j,i]) / dy)
                if etanew[j,i] < -h[j,i]:                                   # 水位が水深よりも下回らないようにする
                    etanew[j,i] = -h[j,i]
            else:
                etanew[j,i] = 0.0
            
            if etanew[j,i] > emax:                                          # 水位の最大値を探す
                emax = etanew[j,i]
                jemax, iemax = j, i
            if etanew[j,i] < emin:                                          # 水位の最小値を探す
                emin = etanew[j,i]
                jemin, iemin = j, i

    print('step=',kt,'file=',flist[kt], '  etamax=','{:.3f}'.format(emax),'[',jemax,',',iemax,']',
          'etamin=','{:.3f}'.format(emin),'[',jemin,',',iemin,']')

    # 水位の境界条件の更新（ゾンマーフェルトの境界条件）
    etanew[:,0]    = eta[:,0]    - dt * C[:,1]  * (eta[:,0]    - eta[:,1]) / dx
    etanew[:,mx+1] = eta[:,mx+1] - dt * C[:,mx] * (eta[:,mx+1] - eta[:,mx]) / dx
    etanew[0,:]    = eta[0,:]    - dt * C[1,:]  * (eta[0,:]    - eta[1,:]) / dy
    etanew[my+1,:] = eta[my+1,:] - dt * C[my,:] * (eta[my+1,:] - eta[my,:]) / dy

#========================================================================================東西の流れの計算
    for i in range (2,mx+1):
        for j in range (1,my+1):
            dm[j,i] = max(etanew[j,i], etanew[j,i-1]) + max(h[j,i], h[j,i-1])  # 計算格子間の水位
            N4[j,i] = (N[j,i] + N[j+1,i] + N[j,i-1] + N[j+1,i-1]) / 4.0        # Mの計算点におけるNの値

    for i in range(2,mx+1):
        for j in range(1,my+1):

            if h[j,i] > 0.0 and h[j,i-1] > 0.0:                                # 計算点の両側が海であること
            
                if M[j,i] > 0.0 and dm[j,i] > 0.0 and dm[j,i-1] > 0.0:         # 流れの方向の確認と計算格子間に水があること
                    Mp1[j,i] = (M[j,i]**2 / dm[j,i] - M[j,i-1]**2 / dm[j,i-1]) / dx

                elif M[j,i] < 0.0 and dm[j,i+1] > 0.0 and dm[j,i] > 0.0:
                    Mp1[j,i] = (M[j,i+1]**2 / dm[j,i+1] - M[j,i]**2 / dm[j,i]) / dx

                else:
                    Mp1[j,i] = 0.0  # 流れの方向がない場合，計算格子間に水がない場合

                if N4[j,i] > 0.0 and dm[j,i] > 0.0 and dm[j-1,i] > 0.0:         # 流れの方向の確認と計算格子間に水があること
                    Mp2[j,i] = (M[j,i] * N4[j,i] / dm[j,i] - M[j-1,i] * N4[j-1,i] / dm[j-1,i]) / dy

                elif N4[j,i] < 0.0 and dm[j+1,i] > 0.0 and dm[j,i] > 0.0:
                    Mp2[j,i] = (M[j+1,i] * N4[j+1,i] / dm[j+1,i] - M[j,i] * N4[j,i] / dm[j,i]) / dy

                else:
                    Mp2[j,i] = 0.0  # 流れの方向がない場合，計算格子間に水がない場合

                # ------------------------------------------- 摩擦項の計算(浅すぎると摩擦項が発散する)
                if dm[j,i] >= 1.0:
                    Mp3[j,i] = g * manning**2 / dm[j,i]**(7.0/3.0) * M[j,i] * (M[j,i]**2 + N4[j,i]**2)**0.5
                elif dm[j,i] < 1.0:
                    Mp3_tmp = g * manning**2 / 1.0**(7.0/3.0) * M[j,i] * (M[j,i]**2 + N4[j,i]**2)**0.5
                    Mp3[j,i] = Mp3_tmp * dm[j,i]**2 
                #--------------------------------------------------------------------------ここまで

                Mp4[j,i] = - f * N4[j,i]  # コリオリ力の項

            else:   # 計算点の少なくともどちらかが陸の場合は岸を抜ける流れはない
                Mp1[j,i] = 0.0
                Mp2[j,i] = 0.0
                Mp3[j,i] = 0.0
                Mp4[j,i] = 0.0

    #----------------------------------------------------------------------- 東西方向の流れを更新
    for i in range (2,mx+1):
        for j in range (1,my+1):
            
            if h[j,i] > 0.0 and h[j,i-1] > 0.0:
                if dm[j,i] > 0.01:
                    Mnew[j,i] = M[j,i] - dt * (Mp1[j,i] + Mp2[j,i] + Mp3[j,i] + Mp4[j,i]
                                               + g * h[j,i] * (etanew[j,i] - etanew[j,i-1]) / dx)
                else:
                    Mnew[j,i] = 0.0
            else:
                Mnew[j,i] = 0.0  
    #-----------------------------------------------------------------------------------ここまで
                
    Mnew[0,:] = Mnew[1,:]      # 境界に沿った方向の流れは境界の摩擦はない境界
    Mnew[my+1,:] = Mnew[my,:]  # 境界に沿った方向の流れは境界の摩擦はない境界
    Mnew[:,1]    = M[:,1]    - dt * C[:,1]    * (M[:,1]    - M[:,2]) / dx  # 境界に垂直な流れ方向はゾンマーフェルト境界
    Mnew[:,mx+1] = M[:,mx+1] - dt * C[:,mx] * (M[:,mx+1] - M[:,mx]) / dx   # 境界に垂直な流れ方向はゾンマーフェルト境界

#========================================================================================南北の流れの計算
    for i in range (1,mx+1):
        for j in range (2,my+1):
            dn[j,i] = max(etanew[j,i], etanew[j-1,i]) + max(h[j,i], h[j-1,i])
            M4[j,i] = (M[j,i] + M[j-1,i] + M[j,i+1] + M[j-1,i+1]) / 4.0

    for i in range (1,mx+1):
        for j in range (2,my+1):

            if h[j,i] > 0.0 and h[j-1,i] > 0.0:

                if M4[j,i] > 0.0 and dn[j,i] > 0.0 and dn[j,i-1] > 0.0:
                    Np1[j,i] = (N[j,i] * M4[j,i] / dn[j,i] - N[j,i-1] * M4[j,i-1] / dn[j,i-1]) / dx

                elif M4[j,i] < 0.0 and dn[j,i+1] > 0.0 and dn[j,i] > 0.0:
                    Np1[j,i] = (N[j,i+1] * M4[j,i+1] / dn[j,i+1] - N[j,i] * M4[j,i] / dn[j,i]) / dx

                else:
                    Np1[j,i] = 0.0

                if N[j,i] > 0.0 and dn[j,i] > 0.0 and dn[j-1,i] > 0.0:
                    Np2[j,i] = (N[j,i]**2 / dn[j,i] - N[j-1,i]**2 / dn[j-1,i]) / dy

                elif N[j,i] < 0.0 and dn[j+1,i] > 0.0 and dn[j,i] > 0.0:
                    Np2[j,i] = (N[j+1,i]**2 / dn[j+1,i] - N[j,i]**2 / dn[j,i]) / dy

                else:
                    Np2[j,i] = 0.0

                # --------------------------------------------------------------------- 摩擦項の計算
                if dn[j,i] >= 1.0:    
                    Np3[j,i] = g * manning**2 /dn[j,i]**(7.0/3.0) * N[j,i] * (M4[j,i]**2+N[j,i]**2)**0.5
                elif dn[j,i] < 1.0:
                    Np3_tmp = g * manning**2 / 1.0**(7.0/3.0) * N[j,i] * (M4[j,i]**2 + N[j,i]**2)**0.5
                    Np3[j,i] = Np3_tmp * dn[j,i]**2 
                #--------------------------------------------------------------------------ここまで

                Np4[j,i] = f * M4[j,i]  # コリオリ力の項

            else:
                Np1[j,i] = 0.0
                Np2[j,i] = 0.0
                Np3[j,i] = 0.0
                Np4[j,i] = 0.0

    #--------------------------------------------------------------------------- 南北方向の流れを更新
    for i in range (1,mx+1):
        for j in range (2,my+1):

            if h[j,i] > 0.0 and h[j-1,i] > 0.0:
                if dn[j,i] > 0.01:
                    Nnew[j,i] = N[j,i] - dt * (Np1[j,i] + Np2[j,i] + Np3[j,i] + Np4[j,i]
                                                + g * dn[j,i] * (etanew[j,i] - etanew[j-1,i]) / dy)
                else:
                    Nnew[j,i] = 0.0
            else:
                Nnew[j,i] = 0.0
    #-----------------------------------------------------------------------------------ここまで

    Nnew[:,0] = Nnew[:,1]
    Nnew[:,mx+1] = Nnew[:,mx]
    Nnew[1,:]    = N[1,:]    - dt * C[1,:]  * (N[1,:] - N[2,:]) / dy
    Nnew[my+1,:] = N[my+1,:] - dt * C[my,:] * (N[my+1,:] - N[my,:]) / dy

#======================================================================================= 結果の出力
    if kt%output_step == 0:
        np.savetxt('data/' +flist[kt], etanew, delimiter=',')

#====================================================================================結果の入れ替え
    eta = etanew.copy()
    M = Mnew.copy()
    N = Nnew.copy()
