# ORB-SLAM3 `export_orb_preint_pack.cpp` 提炼说明

## 1) 导出逻辑（第一性原理）
- 目标：复用 ORB 预积分内部结果，转换到统一定义下比较，不重写 ORB 核心预积分。
- 关键做法：  
  1) 用 ORB 的 `Preintegrated` 计算 `dR/dP/dV/DT`。  
  2) 把 ORB 内部协方差顺序映射到 `z9=[dphi,dp,dv]` 与 `z6=[dba,dbg]`。  
  3) 构建 `Sigma_z15 = blockdiag(Sigma_z9, Sigma_bias_rw)`（按你当前要求不做 combined 15D 全交叉）。  
  4) 从 `dR,dP,dV,DT,JincBias` 构建 `J_e_preint/J_s_preint`。

### 核心代码
```cpp
ORB_SLAM3::IMU::Preintegrated pim(bias0, calib);
for (size_t k = 0; k + 1 < rows.size(); ++k) {
  const double dt = rows[k + 1].t - rows[k].t;
  if (!(dt > 0.0)) continue;
  const Eigen::Vector3d omega = 0.5 * (rows[k].w + rows[k + 1].w);
  const Eigen::Vector3d acc = 0.5 * (rows[k].a + rows[k + 1].a);
  pim.IntegrateNewMeasurement(acc.cast<float>(), omega.cast<float>(), static_cast<float>(dt));
}
```

```cpp
Eigen::Matrix<double,9,9> A9 = Eigen::Matrix<double,9,9>::Zero();
A9.block<3,3>(0,0) = Jr_inv.cast<double>();
A9.block<3,3>(3,6) = Eigen::Matrix3d::Identity();
A9.block<3,3>(6,3) = Eigen::Matrix3d::Identity();
Eigen::Matrix<double,9,9> Sigma_z9 = A9 * C9_orb * A9.transpose();
```

```cpp
const PreintFactorJacobians jac_preint =
    build_preint_factor_jacobians_local(dR_f.cast<double>(), dP_f.cast<double>(), dV_f.cast<double>(), DT, JincBias_ba_bg);
```

## 2) 生成文件内容（`orb_preint_pack.txt`）
- 核心块：  
  - `dR_orb`, `dP_orb`, `dV_orb`, `DT_orb`  
  - `Sigma_z9_orb`（9x9）  
  - `JincBias_ba_bg_orb`（9x6）  
  - `Sigma_bias_rw_orb`（6x6）  
  - `Sigma_z15_orb`（15x15，blockdiag）  
  - `J_e_preint_orb`（15x15）  
  - `J_s_preint_orb`（15x15）
- 说明：当前 `export_orb_preint_pack.cpp` 只输出上述块（不额外导出 ORB 的原始 `C`/`J*` 中间量）。

---

## 3) 详细版（步骤-代码-公式，一步一块）

### Step 1：把连续噪声密度转成 ORB 使用的离散参数
```cpp
const double sf = std::sqrt(1.0 / dt_nominal);
const float ng  = static_cast<float>(cfg.sigma_g_c * sf);
const float na  = static_cast<float>(cfg.sigma_a_c * sf);
const float ngw = static_cast<float>(cfg.sigma_gw_c / sf);
const float naw = static_cast<float>(cfg.sigma_aw_c / sf);
```
对应数学：
$$
n_g=\frac{\sigma_{g,c}}{\sqrt{\Delta t}},\quad
n_a=\frac{\sigma_{a,c}}{\sqrt{\Delta t}},\quad
n_{gw}=\sigma_{gw,c}\sqrt{\Delta t},\quad
n_{aw}=\sigma_{aw,c}\sqrt{\Delta t}.
$$

### Step 2：按 ORB Tracking 一致的 midpoint 方式积分
```cpp
const Eigen::Vector3d omega = 0.5 * (r0.w + r1.w);
const Eigen::Vector3d acc = 0.5 * (r0.a + r1.a);
pim.IntegrateNewMeasurement(acc.cast<float>(), omega.cast<float>(), static_cast<float>(dt));
```
对应数学（每小段）：
$$
\bar\omega_k=\frac{\omega_k+\omega_{k+1}}{2},\qquad
\bar a_k=\frac{a_k+a_{k+1}}{2}.
$$
再递推得到 \(\Delta R,\Delta p,\Delta v,\Delta t\)。

### Step 3：把 ORB 的 `C9` 顺序映射到目标 `z9=[dphi,dp,dv]`
```cpp
Eigen::Matrix<double, 9, 9> A9 = Eigen::Matrix<double, 9, 9>::Zero();
A9.block<3, 3>(0, 0) = Jr_inv.cast<double>();
A9.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
A9.block<3, 3>(6, 3) = Eigen::Matrix3d::Identity();
Sigma_z9 = A9 * C9_orb * A9.transpose();
```
ORB 原顺序是 \([d\theta,dv,dp]\)。映射关系：
$$
\begin{bmatrix}
d\phi\\dp\\dv
\end{bmatrix}

=
A_9
\begin{bmatrix}
d\theta\\dv\\dp
\end{bmatrix},\qquad
\Sigma_{z9}=A_9\,C9_{\text{orb}}\,A_9^\top.
$$

### Step 4：bias 协方差从 `[dbg,dba]` 交换到 `[dba,dbg]`
```cpp
P6.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
P6.block<3, 3>(3, 0) = Eigen::Matrix3d::Identity();
Sigma_bias_rw = P6 * Cb_orb * P6.transpose();
```
对应数学：
$$
\begin{bmatrix}dba\\dbg\end{bmatrix}
=
P_6\begin{bmatrix}dbg\\dba\end{bmatrix},\qquad
\Sigma_{bias}=P_6\,C_{b,\text{orb}}\,P_6^\top.
$$

### Step 5：构建 \(J_{\text{inc,bias}}\) 与 \(J_s,J_e\)
```cpp
const Eigen::Matrix3f dphi_dbg = Jr_inv * pim.JRg;
JincBias_ba_bg.block<3, 3>(0, 3) = dphi_dbg.cast<double>();
JincBias_ba_bg.block<3, 3>(3, 0) = pim.JPa.cast<double>();
...
const PreintFactorJacobians jac_preint =
    build_preint_factor_jacobians_local(dR_f.cast<double>(), dP_f.cast<double>(), dV_f.cast<double>(), DT, JincBias_ba_bg);
```
对应数学：
$$
\frac{\partial d\phi}{\partial b_g} = J_r^{-1}(\phi)\frac{\partial d\theta}{\partial b_g},
$$
$$
J_e=G^{-1},\qquad J_s=-G^{-1}\!\left(F+\begin{bmatrix}G_9J_{\text{inc,bias}}\\0\end{bmatrix}\right).
$$

### Step 6：输出 `z9 + z6 + z15(blockdiag)`（按当前任务要求）
```cpp
Sigma_z15.setZero();
Sigma_z15.block<9, 9>(0, 0) = Sigma_z9;
Sigma_z15.block<6, 6>(9, 9) = Sigma_bias_rw;
```
对应数学：
$$
\Sigma_{z15}=
\begin{bmatrix}
\Sigma_{z9} & 0\\
0 & \Sigma_{bias}
\end{bmatrix}.
$$

---

## 4) 通俗理解版本
- ORB 这条线最重要的现实点是：它原生就更像“9维惯性 + 6维bias随机游走”两个部分，不是一步到位的 combined 15D 全交叉模型。
- 所以 exporter 的策略是：尊重 ORB 原模型，把它先稳定导出，再统一到你的比较接口里。
- `Sigma_z15` 在这里是 `blockdiag(z9,z6)`，你可以把它理解成“把两块结果放到一个15维容器里”，不是重新发明一个 ORB 没有的全耦合协方差。
- `J_s/J_e` 的加入让比较更完整：不仅比较噪声协方差，也比较“残差对状态”的灵敏度是否一致。

---

## 5) 输入与配置（确保“跑出来就是你想要的”）

### 5.1 命令行参数
程序的 `usage()` 写死了 3 个参数（缺一不可）：
```bash
./export_orb_preint_pack \
  --imu_txt <imu_data_Tangent_0.txt> \
  --config_yaml <cpc_config_Tangent_0.yaml> \
  --out_txt <orb_preint_pack.txt>
```

### 5.2 `--imu_txt` 的文本格式
`read_imu_txt()` 期望每行是 7 列、空白分隔：
```
t  gx  gy  gz  ax  ay  az
```
- `t`：时间戳（double，单位秒）。
- `g*`：角速度（建议 rad/s；程序不做单位换算）。
- `a*`：加速度（建议 m/s^2；程序不做单位换算）。
- 空行或以 `#` 开头的行会被跳过。
- 会按 `t` 排序后再积分；若相邻两行 `dt<=0`，该小段会被跳过（但整体仍会继续）。

### 5.3 `--config_yaml` 的“极简 YAML”约定
`load_config_yaml()` 不是完整 YAML 解析器，而是基于字符串查找/逐行扫描：
- 必须包含 4 个标量（缺任意一个会直接报错）：
  - `sigma_g_c`, `sigma_a_c`, `sigma_gw_c`, `sigma_aw_c`
- 可选：`dt`（nominal dt，只用于把连续噪声密度离散化；如果没有，就用 IMU 文件前两帧的 `t` 差作为 `dt_nominal`）。
- 可选：bias（必须是“同一行的 inline list”形式；多行 list 不会被识别）：
```yaml
biases:
  gyro:  [bgx, bgy, bgz]
  accel: [bax, bay, baz]
```
- bias 的使用方式：构造 ORB 的 `IMU::Bias(bax,bay,baz,bwx,bwy,bwz)`；也就是说 YAML 里的 `accel` 对应 `ba`，`gyro` 对应 `bg`。

### 5.4 一个容易忽略但非常关键的点：`dt_nominal` 只影响噪声离散化
程序对“积分本身”用的是每一小段的真实 `dt=r1.t-r0.t`（midpoint）；
而 `dt_nominal` 仅用于把 `sigma_*_c` 转成 ORB `Calib` 里用的离散参数 `ng/na/ngw/naw`。

---

## 6) `build_preint_factor_jacobians_local()` 到底在构什么 Jacobian

### 6.1 两套 15 维量的顺序（这一步不弄清，后面全会看乱）
这个函数内部明确写了两种顺序：
- 状态增量（列顺序）`x`：`[dp, dtheta, dv, dba, dbg]`
- 残差/误差（行顺序）`z`：`[dphi, dp, dv, dba, dbg]`

其中 `dtheta` 与 `dphi` 的区别是：
- `dphi`：我们希望用来比较/统一的 SO(3) 李代数坐标（通过 `phi = Log(dR)` 定义在 `dR` 的右扰动附近）。
- `dtheta`：ORB 预积分内部更常用的旋转误差坐标。
- 两者通过右雅可比关联（代码里用 ORB 的 `InverseRightJacobianSO3`）：
  - `dtheta = Jr(phi) * dphi`
  - `dphi = Jr(phi)^{-1} * dtheta`

### 6.2 `F`：把“起点状态增量”推到“终点状态增量”的一阶近似
代码里给出了（省略 bias 行）：
```cpp
// x order: [dp, dtheta, dv, dba, dbg]
F.block<3,3>(0,3) = -skew(dP);
F.block<3,3>(0,6) = dt * I;
F.block<3,3>(6,3) = -skew(dV);
```
直观上可以读成：
- 终点 `dp` 里既包含起点 `dp`，也会受旋转小扰动影响（`-skew(dP)*dtheta`），并叠加 `dt*dv`。
- 终点 `dv` 同理会被 `-skew(dV)*dtheta` 影响。

### 6.3 `G`：把“残差坐标 z”映射回“状态坐标 x”
它做了两件事：
- 把 `dphi` 换成 `dtheta`（乘上 `Jr`）。
- 把 `dba/dbg` 的符号按残差定义取反（见 `-I`）。
所以最终用 `G_inv`（显式构的近似逆）把“状态增量”转回“残差坐标”，并得到：
- `J_e = G_inv`（残差对终点状态的雅可比）
- `J_s = -G_inv * J`（残差对起点状态的雅可比）

### 6.4 bias 的注入：`JincBias_ba_bg`（9x6）
`JincBias_ba_bg` 的行/列顺序在 exporter 里固定为：
- 行（9）：`[dphi, dp, dv]`
- 列（6）：`[dba, dbg]`
它来自 ORB `Preintegrated` 的 `JPa/JPg/JVa/JVg/JRg`，并额外做了 `dtheta->dphi` 的雅可比变换：
`dphi_dbg = Jr_inv * pim.JRg`。

---

## 7) `orb_preint_pack.txt` 的可解析格式（写 parser 时用）
这个文件是“按块输出矩阵”的纯文本：
- 注释/元信息行以 `#` 开头。
- 每个矩阵块格式：
  1) 第一行：`<name> (<rows>x<cols>)`
  2) 接下来 `rows` 行：每行 `cols` 个数，空格分隔（`fixed`，18 位小数）
  3) 一个空行作为块分隔

与顺序相关的最重要约定：
- `Sigma_z9_orb` 对应 `z9_order: [dphi, dp, dv]`
- `JincBias_ba_bg_orb`：行是 `[dphi, dp, dv]`，列是 `[dba, dbg]`
- `J_e_preint_orb` / `J_s_preint_orb`：行是 `z15=[dphi, dp, dv, dba, dbg]`；列是 `x15=[dp, dtheta, dv, dba, dbg]`
- `Sigma_z15_orb` 在当前 exporter 中是 `blockdiag(Sigma_z9_orb, Sigma_bias_rw_orb)`（惯性 9D 与 bias 6D 的交叉项被置零）
