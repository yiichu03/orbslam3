#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
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

static inline std::string trim(const std::string &s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
    b++;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
    e--;
  return s.substr(b, e - b);
}

static std::string slurp_file(const std::string &path) {
  std::ifstream ifs(path.c_str());
  if (!ifs.is_open()) {
    throw std::runtime_error("unable to open file: " + path);
  }
  std::stringstream buffer;
  buffer << ifs.rdbuf();
  return buffer.str();
}

static bool find_yaml_scalar_double(const std::string &content, const std::string &key, double &out) {
  const std::string needle = key + ":";
  size_t pos = content.find(needle);
  if (pos == std::string::npos)
    return false;
  pos += needle.size();
  while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos])))
    pos++;
  size_t end = pos;
  while (end < content.size()) {
    const char c = content[end];
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

static bool find_yaml_inline_list_under_section(const std::string &content, const std::string &section, const std::string &key,
                                                std::vector<double> &out) {
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

static Eigen::Vector3d to_vec3(const std::vector<double> &v, const std::string &what) {
  if (v.size() != 3) {
    throw std::runtime_error("expected 3 elements for " + what);
  }
  return Eigen::Vector3d(v[0], v[1], v[2]);
}

static Config load_config_yaml(const std::string &path) {
  const std::string content = slurp_file(path);
  Config cfg;

  if (!find_yaml_scalar_double(content, "sigma_g_c", cfg.sigma_g_c) || !find_yaml_scalar_double(content, "sigma_a_c", cfg.sigma_a_c) ||
      !find_yaml_scalar_double(content, "sigma_gw_c", cfg.sigma_gw_c) || !find_yaml_scalar_double(content, "sigma_aw_c", cfg.sigma_aw_c)) {
    throw std::runtime_error("config_yaml missing one of sigma_*_c: " + path);
  }

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

static std::vector<ImuRow> read_imu_txt(const std::string &path) {
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
    double gx = 0.0, gy = 0.0, gz = 0.0, ax = 0.0, ay = 0.0, az = 0.0;
    if (!(iss >> r.t >> gx >> gy >> gz >> ax >> ay >> az))
      continue;
    r.w << gx, gy, gz;
    r.a << ax, ay, az;
    rows.push_back(r);
  }
  if (rows.size() < 2) {
    throw std::runtime_error("imu_txt has too few samples: " + path);
  }
  std::sort(rows.begin(), rows.end(), [](const ImuRow &a, const ImuRow &b) { return a.t < b.t; });
  return rows;
}

template <typename Derived>
static void appendMatrixBlock(std::ostream &os, const std::string &name, const Eigen::MatrixBase<Derived> &mat) {
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

static bool read_next_data_line(std::istream &is, std::string &out) {
  std::string line;
  while (std::getline(is, line)) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    const std::string t = trim(line);
    if (t.empty() || t[0] == '#')
      continue;
    out = t;
    return true;
  }
  return false;
}

static bool parse_header_line(const std::string &line, std::string &name, int &rows, int &cols) {
  const size_t lb = line.find('(');
  const size_t rb = line.find(')', lb == std::string::npos ? 0 : lb + 1);
  if (lb == std::string::npos || rb == std::string::npos || rb <= lb + 1)
    return false;
  name = trim(line.substr(0, lb));
  const std::string inside = trim(line.substr(lb + 1, rb - (lb + 1)));
  const size_t x = inside.find('x');
  if (x == std::string::npos)
    return false;
  rows = std::stoi(trim(inside.substr(0, x)));
  cols = std::stoi(trim(inside.substr(x + 1)));
  return !name.empty() && rows > 0 && cols > 0;
}

static std::unordered_map<std::string, Eigen::MatrixXd> read_matrix_blocks(const std::string &path) {
  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    throw std::runtime_error("unable to open: " + path);
  }
  std::unordered_map<std::string, Eigen::MatrixXd> blocks;
  std::string line;
  while (read_next_data_line(ifs, line)) {
    std::string name;
    int rows = 0;
    int cols = 0;
    if (!parse_header_line(line, name, rows, cols)) {
      throw std::runtime_error("failed to parse matrix header line: '" + line + "' in " + path);
    }
    Eigen::MatrixXd mat(rows, cols);
    for (int r = 0; r < rows; ++r) {
      std::string row_line;
      if (!read_next_data_line(ifs, row_line)) {
        throw std::runtime_error("unexpected EOF reading matrix '" + name + "' from " + path);
      }
      std::istringstream iss(row_line);
      for (int c = 0; c < cols; ++c) {
        double v = 0.0;
        if (!(iss >> v)) {
          throw std::runtime_error("failed to parse matrix '" + name + "' row " + std::to_string(r) + " from " + path);
        }
        mat(r, c) = v;
      }
    }
    blocks[name] = std::move(mat);
  }
  return blocks;
}

static const Eigen::MatrixXd &get_block_or_throw(const std::unordered_map<std::string, Eigen::MatrixXd> &blocks, const std::string &name,
                                                 const std::string &path) {
  const auto it = blocks.find(name);
  if (it == blocks.end()) {
    throw std::runtime_error("missing block '" + name + "' in " + path);
  }
  return it->second;
}

static bool expect_near_abs_rel(const Eigen::MatrixXd &a, const Eigen::MatrixXd &b, double abs_tol, double rel_tol, const std::string &what) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    std::cerr << "[FAIL] " << what << ": shape mismatch: " << a.rows() << "x" << a.cols() << " vs " << b.rows() << "x" << b.cols() << "\n";
    return false;
  }

  double max_violation = -1.0;
  int max_r = 0;
  int max_c = 0;
  double max_diff = 0.0;
  double max_tol = 0.0;
  double max_a = 0.0;
  double max_b = 0.0;
  std::vector<std::pair<int, int>> fail_indices;
  std::vector<double> fail_a;
  std::vector<double> fail_b;
  std::vector<double> fail_diff;
  std::vector<double> fail_tol;

  double rel_needed_if_abs_fixed = 0.0;
  double abs_needed_if_rel_fixed = 0.0;

  for (int r = 0; r < a.rows(); ++r) {
    for (int c = 0; c < a.cols(); ++c) {
      const double va = a(r, c);
      const double vb = b(r, c);
      const double diff = std::abs(va - vb);
      const double scale = std::max(std::abs(va), std::abs(vb));
      const double tol = abs_tol + rel_tol * scale;
      const double violation = diff - tol;
      const double rel_needed_local = (diff > abs_tol && scale > 0.0) ? ((diff - abs_tol) / scale) : 0.0;
      const double abs_needed_local = std::max(0.0, diff - rel_tol * scale);
      rel_needed_if_abs_fixed = std::max(rel_needed_if_abs_fixed, rel_needed_local);
      abs_needed_if_rel_fixed = std::max(abs_needed_if_rel_fixed, abs_needed_local);
      if (violation > max_violation) {
        max_violation = violation;
        max_r = r;
        max_c = c;
        max_diff = diff;
        max_tol = tol;
        max_a = va;
        max_b = vb;
      }
      if (violation > 0.0) {
        fail_indices.emplace_back(r, c);
        fail_a.push_back(va);
        fail_b.push_back(vb);
        fail_diff.push_back(diff);
        fail_tol.push_back(tol);
      }
    }
  }

  if (max_violation > 0.0) {
    std::cerr << std::setprecision(18);
    std::cerr << "[FAIL] " << what << ": max violation at (" << max_r << "," << max_c << ")\n";
    std::cerr << "  a=" << max_a << " b=" << max_b << " |a-b|=" << max_diff << " tol=" << max_tol << " (abs=" << abs_tol
              << ", rel=" << rel_tol << ")\n";
    std::cerr << "  failing_entries_count=" << fail_indices.size() << "\n";
    for (size_t i = 0; i < fail_indices.size(); ++i) {
      std::cerr << "  (" << fail_indices[i].first << "," << fail_indices[i].second << ")"
                << " a=" << fail_a[i] << " b=" << fail_b[i] << " |a-b|=" << fail_diff[i] << " tol=" << fail_tol[i] << "\n";
    }
    std::cerr << "  to_pass_by_tuning_tolerance:\n";
    std::cerr << "    rel_needed_if_abs_fixed=" << rel_needed_if_abs_fixed << " (current_rel=" << rel_tol << ", abs_fixed=" << abs_tol
              << ")\n";
    std::cerr << "    abs_needed_if_rel_fixed=" << abs_needed_if_rel_fixed << " (current_abs=" << abs_tol << ", rel_fixed=" << rel_tol
              << ")\n";
    return false;
  }

  std::cout << "[ OK ] " << what << "\n";
  return true;
}

static std::string dirname_of(const std::string &path) {
  const size_t slash = path.find_last_of("/\\");
  if (slash == std::string::npos)
    return ".";
  if (slash == 0)
    return path.substr(0, 1);
  return path.substr(0, slash);
}

static std::string join_path(const std::string &dir, const std::string &name) {
  if (dir.empty())
    return name;
  const char last = dir.back();
  if (last == '/' || last == '\\')
    return dir + name;
  return dir + "/" + name;
}

static void compare_orbslam3_against_gtsam_ref(const std::string &config_yaml, const Eigen::Matrix<double, 9, 9> &Sigma_z9_orb,
                                               const Eigen::Matrix<double, 9, 6> &JincBias_ba_bg_orb) {
  const std::string imu_data_dir = dirname_of(config_yaml);
  const std::string gtsam_all = join_path(join_path(imu_data_dir, "gtsam_ref_out_orb"), "gtsam_ref_orb_preint_all.txt");
  const auto gtsam = read_matrix_blocks(gtsam_all);
  const Eigen::MatrixXd &Sigma_z9_gtsam = get_block_or_throw(gtsam, "Sigma_z9_gtsam", gtsam_all);
  const Eigen::MatrixXd &JincBias_gtsam = get_block_or_throw(gtsam, "JincBias_ba_bg_gtsam", gtsam_all);

  constexpr double abs_tol = 1e-4;
  constexpr double rel_tol = 1.5e-2;
  bool ok = true;
  ok &= expect_near_abs_rel(Eigen::MatrixXd(Sigma_z9_orb), Sigma_z9_gtsam, abs_tol, rel_tol, "Sigma_z9 (z9=[dphi,dp,dv])");
  ok &= expect_near_abs_rel(Eigen::MatrixXd(JincBias_ba_bg_orb), JincBias_gtsam, abs_tol, rel_tol, "JincBias_ba_bg (rows=[dphi,dp,dv])");
  if (!ok) {
    throw std::runtime_error("comparison failed");
  }
}

static Eigen::Matrix3d theta_to_phi_jacobian_from_dR(const Eigen::Matrix3f &dR) {
  const Eigen::Vector3f phi_hat = Sophus::SO3f(dR).log();
  return ORB_SLAM3::IMU::InverseRightJacobianSO3(phi_hat).cast<double>();
}

static void usage(const char *argv0) {
  std::cerr << "usage: " << argv0
            << " --imu_txt <imu_data_Tangent_0.txt> --config_yaml <cpc_config_Tangent_0.yaml> --out_txt <orb_preint_pack.txt>\n";
}

} // namespace

int main(int argc, char **argv) {
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

    const double sf = std::sqrt(1.0 / dt_nominal);
    const float ng = static_cast<float>(cfg.sigma_g_c * sf);
    const float na = static_cast<float>(cfg.sigma_a_c * sf);
    const float ngw = static_cast<float>(cfg.sigma_gw_c / sf);
    const float naw = static_cast<float>(cfg.sigma_aw_c / sf);

    const ORB_SLAM3::IMU::Bias bias0(static_cast<float>(cfg.bias_accel.x()), static_cast<float>(cfg.bias_accel.y()),
                                     static_cast<float>(cfg.bias_accel.z()), static_cast<float>(cfg.bias_gyro.x()),
                                     static_cast<float>(cfg.bias_gyro.y()), static_cast<float>(cfg.bias_gyro.z()));

    const Sophus::SE3f Tbc(Sophus::SO3f(), Eigen::Vector3f::Zero());
    const ORB_SLAM3::IMU::Calib calib(Tbc, ng, na, ngw, naw);
    ORB_SLAM3::IMU::Preintegrated pim(bias0, calib);

    for (size_t k = 0; k + 1 < rows.size(); ++k) {
      const ImuRow &r0 = rows[k];
      const ImuRow &r1 = rows[k + 1];
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

    const Eigen::Matrix<double, 9, 9> C9_orb = pim.C.block<9, 9>(0, 0).cast<double>(); // [dtheta,dv,dp]
    Eigen::Matrix3d T_theta_to_phi = Eigen::Matrix3d::Identity();
    // T_theta_to_phi = theta_to_phi_jacobian_from_dR(dR_f); // Comment this line to reproduce dphi == dtheta behavior.

    Eigen::Matrix<double, 9, 9> A9 = Eigen::Matrix<double, 9, 9>::Zero();
    A9.block<3, 3>(0, 0) = T_theta_to_phi;
    A9.block<3, 3>(3, 6) = Eigen::Matrix3d::Identity(); // dp <- orb dp
    A9.block<3, 3>(6, 3) = Eigen::Matrix3d::Identity(); // dv <- orb dv

    const Eigen::Matrix<double, 9, 9> Sigma_z9 = A9 * C9_orb * A9.transpose();

    const Eigen::Matrix<double, 6, 6> Cb_orb = pim.C.block<6, 6>(9, 9).cast<double>(); // [dbg,dba]
    Eigen::Matrix<double, 6, 6> P6 = Eigen::Matrix<double, 6, 6>::Zero();
    P6.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity(); // dba <- orb dba
    P6.block<3, 3>(3, 0) = Eigen::Matrix3d::Identity(); // dbg <- orb dbg
    const Eigen::Matrix<double, 6, 6> Sigma_bias_rw = P6 * Cb_orb * P6.transpose();

    Eigen::Matrix<double, 9, 6> JincBias_ba_bg = Eigen::Matrix<double, 9, 6>::Zero();
    JincBias_ba_bg.block<3, 3>(0, 3) = T_theta_to_phi * pim.JRg.cast<double>();
    JincBias_ba_bg.block<3, 3>(3, 0) = pim.JPa.cast<double>();
    JincBias_ba_bg.block<3, 3>(3, 3) = pim.JPg.cast<double>();
    JincBias_ba_bg.block<3, 3>(6, 0) = pim.JVa.cast<double>();
    JincBias_ba_bg.block<3, 3>(6, 3) = pim.JVg.cast<double>();
    std::cout << "[debug] T_theta_to_phi * pim.JRg:\n" << (T_theta_to_phi * pim.JRg.cast<double>()) << "\n[debug] pim.JRg:\n"
              << pim.JRg.cast<double>() << "\n";

    Eigen::Matrix<double, 15, 15> Sigma_z15 = Eigen::Matrix<double, 15, 15>::Zero();
    Sigma_z15.block<9, 9>(0, 0) = Sigma_z9;
    Sigma_z15.block<6, 6>(9, 9) = Sigma_bias_rw;

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
    appendMatrixBlock(ofs, "Sigma_bias_rw_orb", Sigma_bias_rw);
    appendMatrixBlock(ofs, "Sigma_z15_orb", Sigma_z15);
    compare_orbslam3_against_gtsam_ref(config_yaml, Sigma_z9, JincBias_ba_bg);

    std::cout << std::setprecision(17) << "integrated interval: ts=" << rows.front().t << " te=" << rows.back().t << " DT=" << DT << "\n";
    std::cout << "wrote: " << out_txt << "\n";
    return 0;

  } catch (const std::exception &e) {
    std::cerr << "export_orb_preint_pack failed: " << e.what() << "\n";
    return 1;
  }
}
