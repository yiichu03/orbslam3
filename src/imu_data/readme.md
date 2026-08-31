ORB-SLAM3 的 exporter（orb_preint_pack.txt）里，compare_orbslam3_gtsam 实际只用这 8 个 block 名称：

dR_orb
dP_orb
dV_orb
DT_orb
Sigma_z9_orb
JincBias_ba_bg_orb
Sigma_bias_rw_orb
Sigma_z15_orb
其它像 C9_orb_raw、Cbias_orb_raw、JRg_orb_raw/JVa_orb_raw/... 也都是 debug，不参与 compare。

按它自身建模导出 9D 预积分测量 + 6D biasRW（以及 blockdiag 的 Sigma_z15），并用 compare_orbslam3_gtsam 对 GTSAM 9D+biasRW reference 全部通过，满足“按 GTSAM Tangent 误差状态定义拿到 jac/cov”的要求（只是它不是 Combined 15D 的那种带交叉项模型）。

## 交接复核

- ORB-SLAM3 核心预积分实现未修改；新增内容仅为独立 exporter、对比逻辑、输入和结果文件。
- exporter 使用 `InverseRightJacobianSO3(Log(dR))`，把 ORB-SLAM3 的局部旋转误差 `dtheta` 映射到 GTSAM Tangent 的 `dphi` 后，再导出协方差和 bias Jacobian。
- 2026-08-31 从当前源码重新编译并运行后，内置比较器和独立 `compare_orbslam3_gtsam` 均得到：

```text
[ OK ] Sigma_z9 (z9=[dphi,dp,dv])
[ OK ] JincBias_ba_bg (rows=[dphi,dp,dv])
```




文件夹有🔓时候 sudo chown -R $USER:$USER /路径/到/文件夹

docker start orbslam3_melodic
docker exec -it orbslam3_melodic bash



cd /workspace/ORB_SLAM3
chmod +x build.sh
./build.sh

/workspace/ORB_SLAM3/src/tools/export_orb_preint_pack \
   --imu_txt     /workspace/ORB_SLAM3/src/imu_data/imu_data_Tangent_0.txt \
   --config_yaml /workspace/ORB_SLAM3/src/imu_data/cpc_config_Tangent_0.yaml \
   --out_txt     /workspace/ORB_SLAM3/src/imu_data/orb_preint_pack.txt

# Export ORB-SLAM3 preintegration pack (9D preint + 6D bias RW)
# 生成的 orb_preint_pack.txt 需要拷贝到 swift_vio 的 imu_data 目录用于后续 gtsam reference + compare
#
# build 之后可执行文件在：
#   /workspace/ORB_SLAM3/src/tools/export_orb_preint_pack
#
# 用法示例：
/workspace/ORB_SLAM3/src/tools/export_orb_preint_pack \
   --imu_txt     /workspace/ORB_SLAM3/src/imu_data/imu_data_Tangent_0.txt \
   --config_yaml /workspace/ORB_SLAM3/src/imu_data/cpc_config_Tangent_0.yaml \
   --out_txt     /workspace/ORB_SLAM3/src/imu_data/orb_preint_pack.txt
