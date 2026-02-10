#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ImuTypes.h"

namespace {

struct ImuRow {
  double t = 0.0;
  Eigen::Vector3d w = Eigen::Vector3d::Zero();
  Eigen::Vector3d a = Eigen::Vector3d::Zero();
};

struct Config {
  double sigma_g_c = 0.0;
  double sigma_a_c = 0.0;
  double sigma_gw_c = 0.0;
  double sigma_aw_c = 0.0;
  double dt_nominal = 0.0;
  Eigen::Vector3d bias_gyro = Eigen::Vector3d::Zero();
  Eigen::Vector3d bias_accel = Eigen::Vector3d::Zero();
};

static inline std::string trim(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
    b++;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
    e--;
  return s.substr(b, e - b);
}

static std::string slurp_file(const std::string& path) {
  std::ifstream ifs(path.c_str());
  if (!ifs.is_open()) {
    throw std::runtime_error("unable to open file: " + path);
  }
  std::stringstream buffer;
  buffer << ifs.rdbuf();
  return buffer.str();
}

static bool find_yaml_scalar_double(const std::string& content, const std::string& key, double& out) {
  const std::string needle = key + ":";
  size_t pos = content.find(needle);
  if (pos == std::string::npos)
    return false;
  pos += needle.size();
  while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos])))
    pos++;
  size_t end = pos;
  while (end < content.size()) {
    char c = content[end];
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '+' || c == '-' || c == '.' || c == 'e' || c == 'E') {
      end++;
      continue;
    }
    break;
  }
  if (end == pos)
    return false;
  out = std::stod(content.substr(pos, end - pos));
  return true;
}

static bool find_yaml_inline_list_under_section(const std::string& content, const std::string& section, const std::string& key,
                                                std::vector<double>& out) {
  std::istringstream iss(content);
  std::string line;
  bool in_section = false;
  const std::string section_hdr = section + ":";
  const std::string key_hdr = key + ":";

  while (std::getline(iss, line)) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();

    size_t first = 0;
    while (first < line.size() && (line[first] == ' ' || line[first] == '\t'))
      first++;
    if (first == line.size())
      continue;

    if (!in_section) {
      if (first == 0 && line.compare(0, section_hdr.size(), section_hdr) == 0) {
        in_section = true;
      }
      continue;
    }

    if (first == 0) {
      break;
    }

    if (line.compare(first, key_hdr.size(), key_hdr) != 0) {
      continue;
    }

    size_t lb = line.find('[', first + key_hdr.size());
    if (lb == std::string::npos)
      return false;
    size_t rb = line.find(']', lb + 1);
    if (rb == std::string::npos)
      return false;
    std::string inside = line.substr(lb + 1, rb - (lb + 1));

    std::vector<double> vals;
    std::stringstream ss(inside);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      std::string t = trim(tok);
      if (t.empty())
        continue;
      vals.push_back(std::stod(t));
    }
    out = std::move(vals);
    return !out.empty();
  }
  return false;
}

static Eigen::Vector3d to_vec3(const std::vector<double>& v, const std::string& what) {
  if (v.size() != 3) {
    throw std::runtime_error("expected 3 elements for " + what);
  }
  return Eigen::Vector3d(v[0], v[1], v[2]);
}

static Config load_config_yaml(const std::string& path) {
  const std::string content = slurp_file(path);
  Config cfg;

  if (!find_yaml_scalar_double(content, "sigma_g_c", cfg.sigma_g_c) || !find_yaml_scalar_double(content, "sigma_a_c", cfg.sigma_a_c) ||
      !find_yaml_scalar_double(content, "sigma_gw_c", cfg.sigma_gw_c) || !find_yaml_scalar_double(content, "sigma_aw_c", cfg.sigma_aw_c)) {
    throw std::runtime_error("config_yaml missing one of sigma_*_c: " + path);
  }

  // nominal dt (for noise discretization)
  if (!find_yaml_scalar_double(content, "dt", cfg.dt_nominal)) {
    cfg.dt_nominal = 0.0;
  }

  std::vector<double> bg, ba;
  if (find_yaml_inline_list_under_section(content, "biases", "gyro", bg) && bg.size() == 3) {
    cfg.bias_gyro = to_vec3(bg, "biases.gyro");
  }
  if (find_yaml_inline_list_under_section(content, "biases", "accel", ba) && ba.size() == 3) {
    cfg.bias_accel = to_vec3(ba, "biases.accel");
  }

  return cfg;
}

static std::vector<ImuRow> read_imu_txt(const std::string& path) {
  std::ifstream ifs(path.c_str());
  if (!ifs.is_open()) {
    throw std::runtime_error("unable to open imu_txt: " + path);
  }

  std::vector<ImuRow> rows;
  std::string line;
  while (std::getline(ifs, line)) {
    size_t i = 0;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i])))
      i++;
    if (i == line.size() || line[i] == '#')
      continue;

    std::istringstream iss(line);
    ImuRow r;
    double gx, gy, gz, ax, ay, az;
    if (!(iss >> r.t >> gx >> gy >> gz >> ax >> ay >> az))
      continue;
    r.w << gx, gy, gz;
    r.a << ax, ay, az;
    rows.push_back(r);
  }
  if (rows.size() < 2) {
    throw std::runtime_error("imu_txt has too few samples: " + path);
  }
  std::sort(rows.begin(), rows.end(), [](const ImuRow& a, const ImuRow& b) { return a.t < b.t; });
  return rows;
}

template <typename Derived>
static void appendMatrixBlock(std::ostream& os, const std::string& name, const Eigen::MatrixBase<Derived>& mat) {
  os << name << " (" << mat.rows() << "x" << mat.cols() << ")\n";
  os.setf(std::ios::fixed);
  os << std::setprecision(18);
  for (int r = 0; r < mat.rows(); ++r) {
    for (int c = 0; c < mat.cols(); ++c) {
      os << mat(r, c);
      if (c < mat.cols() - 1) {
        os << ' ';
      }
    }
    os << '\n';
  }
  os << '\n';
}

struct PreintFactorJacobians {
  Eigen::Matrix<double, 15, 15> J_s = Eigen::Matrix<double, 15, 15>::Zero();
  Eigen::Matrix<double, 15, 15> J_e = Eigen::Matrix<double, 15, 15>::Zero();
};

static inline Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
  Eigen::Matrix3d S;
  S << 0.0, -v.z(), v.y(), v.z(), 0.0, -v.x(), -v.y(), v.x(), 0.0;
  return S;
}

static PreintFactorJacobians build_preint_factor_jacobians_local(const Eigen::Matrix3d& dR, const Eigen::Vector3d& dP,
                                                                  const Eigen::Vector3d& dV, const double dt,
                                                                  const Eigen::Matrix<double, 9, 6>& JincBias_ba_bg) {
  const Eigen::Vector3f phi_f = Sophus::SO3f(dR.cast<float>()).log();
  const Eigen::Matrix3f Jr_inv_f = ORB_SLAM3::IMU::InverseRightJacobianSO3(phi_f);
  const Eigen::Matrix3d Jr_inv = Jr_inv_f.cast<double>();
  const Eigen::Matrix3d Jr = Jr_inv.inverse();
  const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

  // State x order: [dp, dtheta, dv, dba, dbg]
  Eigen::Matrix<double, 15, 15> F = Eigen::Matrix<double, 15, 15>::Zero();
  F.block<3, 3>(0, 0) = I;
  F.block<3, 3>(0, 3) = -skew(dP);
  F.block<3, 3>(0, 6) = dt * I;
  F.block<3, 3>(3, 3) = I;
  F.block<3, 3>(6, 3) = -skew(dV);
  F.block<3, 3>(6, 6) = I;
  F.block<3, 3>(9, 9) = I;
  F.block<3, 3>(12, 12) = I;

  // Residual z order: [dphi, dp, dv, dba, dbg]
  Eigen::Matrix<double, 15, 15> G = Eigen::Matrix<double, 15, 15>::Zero();
  G.block<3, 3>(0, 3) = I;     // dp_e <- dp
  G.block<3, 3>(3, 0) = Jr;    // dtheta_e <- Jr * dphi
  G.block<3, 3>(6, 6) = I;     // dv_e <- dv
  G.block<3, 3>(9, 9) = -I;    // dba_e <- -dba
  G.block<3, 3>(12, 12) = -I;  // dbg_e <- -dbg

  Eigen::Matrix<double, 15, 15> G_inv = G.transpose();
  G_inv.block<3, 3>(0, 3) = Jr_inv;

  Eigen::Matrix<double, 15, 15> J = F;
  J.topRightCorner<9, 6>() += G.topLeftCorner<9, 9>() * JincBias_ba_bg;

  PreintFactorJacobians out;
  out.J_e = G_inv;
  out.J_s = -G_inv * J;
  return out;
}

static void usage(const char* argv0) {
  std::cerr << "usage: " << argv0 << " --imu_txt <imu_data_Tangent_0.txt> --config_yaml <cpc_config_Tangent_0.yaml> --out_txt <orb_preint_pack.txt>\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string imu_txt;
    std::string config_yaml;
    std::string out_txt;

    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      auto next = [&]() -> std::string {
        if (i + 1 >= argc) {
          throw std::runtime_error("missing value after " + a);
        }
        return std::string(argv[++i]);
      };
      if (a == "--imu_txt") {
        imu_txt = next();
      } else if (a == "--config_yaml") {
        config_yaml = next();
      } else if (a == "--out_txt") {
        out_txt = next();
      } else if (a == "--help" || a == "-h") {
        usage(argv[0]);
        return 0;
      } else {
        throw std::runtime_error("unknown argument: " + a);
      }
    }

    if (imu_txt.empty() || config_yaml.empty() || out_txt.empty()) {
      usage(argv[0]);
      return 1;
    }

    const Config cfg = load_config_yaml(config_yaml);
    const std::vector<ImuRow> rows = read_imu_txt(imu_txt);

    const double dt_nominal = (cfg.dt_nominal > 0.0) ? cfg.dt_nominal : (rows[1].t - rows[0].t);
    if (!(dt_nominal > 0.0)) {
      throw std::runtime_error("invalid dt_nominal");
    }

    // ORB-SLAM3 Tracking.cc uses sf = sqrt(freq) for discretization.
    const double sf = std::sqrt(1.0 / dt_nominal);
    const float ng = static_cast<float>(cfg.sigma_g_c * sf);
    const float na = static_cast<float>(cfg.sigma_a_c * sf);
    const float ngw = static_cast<float>(cfg.sigma_gw_c / sf);
    const float naw = static_cast<float>(cfg.sigma_aw_c / sf);

    const ORB_SLAM3::IMU::Bias bias0(
        static_cast<float>(cfg.bias_accel.x()),
        static_cast<float>(cfg.bias_accel.y()),
        static_cast<float>(cfg.bias_accel.z()),
        static_cast<float>(cfg.bias_gyro.x()),
        static_cast<float>(cfg.bias_gyro.y()),
        static_cast<float>(cfg.bias_gyro.z()));

    const Sophus::SE3f Tbc(Sophus::SO3f(), Eigen::Vector3f::Zero());
    const ORB_SLAM3::IMU::Calib calib(Tbc, ng, na, ngw, naw);
    ORB_SLAM3::IMU::Preintegrated pim(bias0, calib);

    // Midpoint (trapezoidal) integration, consistent with ORB-SLAM3 Tracking.cc.
    for (size_t k = 0; k + 1 < rows.size(); ++k) {
      const ImuRow& r0 = rows[k];
      const ImuRow& r1 = rows[k + 1];
      const double dt = r1.t - r0.t;
      if (!(dt > 0.0))
        continue;
      const Eigen::Vector3d omega = 0.5 * (r0.w + r1.w);
      const Eigen::Vector3d acc = 0.5 * (r0.a + r1.a);
      pim.IntegrateNewMeasurement(acc.cast<float>(), omega.cast<float>(), static_cast<float>(dt));
    }

    const Eigen::Matrix3f dR_f = pim.GetOriginalDeltaRotation();
    const Eigen::Vector3f dP_f = pim.GetOriginalDeltaPosition();
    const Eigen::Vector3f dV_f = pim.GetOriginalDeltaVelocity();
    const double DT = static_cast<double>(pim.dT);

    const Eigen::Vector3f phi = Sophus::SO3f(dR_f).log();
    const Eigen::Matrix3f Jr_inv = ORB_SLAM3::IMU::InverseRightJacobianSO3(phi);

    // ORB C9 order: [dtheta, dv, dp]
    const Eigen::Matrix<double, 9, 9> C9_orb = pim.C.block<9, 9>(0, 0).cast<double>();

    // Target z9 order: [dphi, dp, dv]
    Eigen::Matrix<double, 9, 9> A9 = Eigen::Matrix<double, 9, 9>::Zero();
    A9.block<3, 3>(0, 0) = Jr_inv.cast<double>();
    A9.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity();
    A9.block<3, 3>(6, 3) = Eigen::Matrix3d::Identity();

    Eigen::Matrix<double, 9, 9> Sigma_z9 = A9 * C9_orb * A9.transpose();

    // ORB bias block order: [dbg, dba]
    const Eigen::Matrix<double, 6, 6> Cb_orb = pim.C.block<6, 6>(9, 9).cast<double>();
    Eigen::Matrix<double, 6, 6> P6 = Eigen::Matrix<double, 6, 6>::Zero();
    P6.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity(); // dba <- orb dba (at +3)
    P6.block<3, 3>(3, 0) = Eigen::Matrix3d::Identity(); // dbg <- orb dbg (at +0)
    Eigen::Matrix<double, 6, 6> Sigma_bias_rw = P6 * Cb_orb * P6.transpose(); // [dba, dbg]

    // Jacobians w.r.t biases in target z9 coordinates.
    const Eigen::Matrix3f dphi_dbg = Jr_inv * pim.JRg; // dphi / dbg
    Eigen::Matrix<double, 9, 6> JincBias_ba_bg = Eigen::Matrix<double, 9, 6>::Zero();
    // rows: [dphi, dp, dv], cols: [dba, dbg]
    JincBias_ba_bg.block<3, 3>(0, 3) = dphi_dbg.cast<double>();
    JincBias_ba_bg.block<3, 3>(3, 0) = pim.JPa.cast<double>();
    JincBias_ba_bg.block<3, 3>(3, 3) = pim.JPg.cast<double>();
    JincBias_ba_bg.block<3, 3>(6, 0) = pim.JVa.cast<double>();
    JincBias_ba_bg.block<3, 3>(6, 3) = pim.JVg.cast<double>();

    // Optional: 15D blockdiag covariance in z15=[dphi,dp,dv,dba,dbg] order.
    Eigen::Matrix<double, 15, 15> Sigma_z15 = Eigen::Matrix<double, 15, 15>::Zero();
    Sigma_z15.block<9, 9>(0, 0) = Sigma_z9;
    Sigma_z15.block<6, 6>(9, 9) = Sigma_bias_rw;
    const PreintFactorJacobians jac_preint =
        build_preint_factor_jacobians_local(dR_f.cast<double>(), dP_f.cast<double>(), dV_f.cast<double>(), DT, JincBias_ba_bg);

    // Write pack.
    std::ofstream ofs(out_txt.c_str(), std::ios::trunc);
    if (!ofs.is_open()) {
      throw std::runtime_error("failed to open for writing: " + out_txt);
    }
    ofs << "# export_orb_preint_pack\n";
    ofs << "# imu_txt: " << imu_txt << "\n";
    ofs << "# config_yaml: " << config_yaml << "\n";
    ofs << std::setprecision(17) << "# interval: ts=" << rows.front().t << " te=" << rows.back().t << " DT=" << DT << "\n";
    ofs << std::setprecision(17) << "# dt_nominal: " << dt_nominal << "\n";
    ofs << std::setprecision(17) << "# orb_noise_discrete: Ng=" << ng << " Na=" << na << " Ngw=" << ngw << " Naw=" << naw << "\n";
    ofs << "# z9_order: [dphi, dp, dv]\n";
    ofs << "# z6_order: [dba, dbg]\n";
    ofs << "# z15_order: [dphi, dp, dv, dba, dbg] (blockdiag)\n\n";

    appendMatrixBlock(ofs, "dR_orb", dR_f.cast<double>());
    appendMatrixBlock(ofs, "dP_orb", dP_f.cast<double>());
    appendMatrixBlock(ofs, "dV_orb", dV_f.cast<double>());
    Eigen::Matrix<double, 1, 1> DT_mat;
    DT_mat(0, 0) = DT;
    appendMatrixBlock(ofs, "DT_orb", DT_mat);

    appendMatrixBlock(ofs, "Sigma_z9_orb", Sigma_z9);
    appendMatrixBlock(ofs, "JincBias_ba_bg_orb", JincBias_ba_bg);
    appendMatrixBlock(ofs, "J_e_preint_orb", jac_preint.J_e);
    appendMatrixBlock(ofs, "J_s_preint_orb", jac_preint.J_s);
    appendMatrixBlock(ofs, "Sigma_bias_rw_orb", Sigma_bias_rw);
    appendMatrixBlock(ofs, "Sigma_z15_orb", Sigma_z15);

    // Raw internal blocks for debugging.
    appendMatrixBlock(ofs, "C9_orb_raw", C9_orb);
    appendMatrixBlock(ofs, "Cbias_orb_raw", Cb_orb);
    appendMatrixBlock(ofs, "JRg_orb_raw", pim.JRg.cast<double>());
    appendMatrixBlock(ofs, "JVa_orb_raw", pim.JVa.cast<double>());
    appendMatrixBlock(ofs, "JVg_orb_raw", pim.JVg.cast<double>());
    appendMatrixBlock(ofs, "JPa_orb_raw", pim.JPa.cast<double>());
    appendMatrixBlock(ofs, "JPg_orb_raw", pim.JPg.cast<double>());

    std::cout << std::setprecision(17) << "integrated interval: ts=" << rows.front().t << " te=" << rows.back().t << " DT=" << DT << "\n";
    std::cout << "wrote: " << out_txt << "\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "export_orb_preint_pack failed: " << e.what() << "\n";
    return 1;
  }
}
