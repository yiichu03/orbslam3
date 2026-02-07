文件夹有🔓时候 sudo chown -R $USER:$USER /路径/到/文件夹

docker start orbslam3_melodic
docker exec -it orbslam3_melodic bash



cd /workspace/ORB_SLAM3
chmod +x build.sh
./build.sh

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
