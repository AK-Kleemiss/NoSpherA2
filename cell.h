#pragma once

class cell
{
private:
  double as, bs, cs, a, b, c, alpha, beta, gamma, V;
  double ca, cb, cg, sa, sb, sg;
  double rcm[3][3];
  double cm[3][3];
  std::string crystal_system;
  std::vector < std::vector < std::vector <int> > > sym;
  const double PI = 3.14159265358979323846;
public:
  cell()
  {
    as = bs = cs = a = b = c = alpha = beta = gamma = V = ca = cb = cg = sa = sb = sg = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        rcm[i][j] = 0.0;
        cm[i][j] = 0.0;
      }
    sym.resize(3);
    crystal_system = "triclinic";
    for (int i = 0; i < 3; i++)
      sym[i].resize(3);
  };
  cell(const std::string filename)
  {
    sym.resize(3);
    for (int i = 0; i < 3; i++)
      sym[i].resize(3);
    read_CIF(filename);
    read_symm_CIF(filename);
  };
  cell(const std::string filename, std::ofstream& file)
  {
    file << "Reading: " << std::setw(44) << filename << std::flush;
    sym.resize(3);
    for (int i = 0; i < 3; i++)
      sym[i].resize(3);
    read_CIF(filename, file);
    read_symm_CIF(filename, file);
  };
  cell(const std::string filename, std::ofstream& file, const bool& debug)
  {
    if (debug)
      file << "starting to read cif!" << std::endl;
    file << "Reading: " << std::setw(44) << filename << std::flush;
    sym.resize(3);
    for (int i = 0; i < 3; i++)
      sym[i].resize(3);
    read_CIF(filename, file, debug);
    read_symm_CIF(filename, file, debug);
    if (debug) {
      file << "RCM done!" << std::endl;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
          file << std::setw(10) << std::fixed << get_rcm(i, j) / 2 / PI / 0.529177249 << ' ';
        file << std::endl;
      }
      file << "CM in 2*PI bohr:" << std::endl;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
          file << std::setw(10) << std::fixed << get_cm(i, j) << ' ';
        file << std::endl;
      }
    }
  };
  double get_as() { return as; };
  double get_bs() { return bs; };
  double get_cs() { return cs; };
  double get_ca() { return ca; };
  double get_cb() { return cb; };
  double get_cg() { return cg; };
  double get_sa() { return sa; };
  double get_sb() { return sb; };
  double get_sg() { return sg; };
  double get_rcm(int i, int j) { return rcm[i][j]; };
  double get_cm(int i, int j) { return cm[i][j]; };
  double get_sym(int i, int j, int k) { return sym[i][j][k]; };
  std::vector <std::vector < std::vector<int>>> get_sym() { return sym; };
  double get_a() { return a; };
  double get_b() { return b; };
  double get_c() { return c; };
  double get_V() { return V; };
  double get_angle(int i) const
  {
    switch (i) {
    case 0:
      return alpha;
      break;
    case 1:
      return beta;
      break;
    case 2:
      return gamma;
      break;
    default:
      std::cout << "Wrong angle!" << std::endl;
      return -400;
      break;
    }
  }
  double get_length(int i) const
  {
    switch (i) {
    case 0:
      return a;
      break;
    case 1:
      return b;
      break;
    case 2:
      return c;
      break;
    default:
      std::cout << "Wrong length!" << std::endl;
      return -1;
      break;
    }
  }

  void set_system()
  {
    if (beta != 90.0) {
      crystal_system = "triclinic";
      return;
    }
    else {
      if (alpha == beta && alpha == 90.0)
        crystal_system = "monoclinic";
    }
  }

  std::string get_crystal_system() const
  {
    return crystal_system;
  }

  double get_d_of_hkl(const std::vector<int> hkl) const 
  {
    double d, upper, lower;
    // d = sqrt( (1 - cos^2(alpha) - cos^2 (beta) - cos^2 (gamma) + 2ca*cb*cg) / (h^2 /a^2 *sin^2 (alpha) + k^2 / b^2 * sin^2 (beta) + l^2 /c^2 * sin^2(gamma) + 2 kl/bc (cos(beta)cos(gamma) - cos(alpha)) + 2 hl/ac (cos(alpha)cos(gamma) - cos(beta)) + 2 hk/ab (cos(beta)cos(alpha) - cos(gamma))) )
    upper = pow(hkl[0], 2) * pow(sa, 2) / pow(a, 2) + pow(hkl[1], 2) * pow(sb, 2) / pow(b, 2) + pow(hkl[2], 2) * pow(sg, 2) / pow(c, 2) + 2 * hkl[1] * hkl[2] / (b * c) * (cb * cg - ca) + 2 * hkl[0] * hkl[2] / (a * c) * (cg * ca - cb) + 2 * hkl[0] * hkl[1] / (a * b) * (ca * cb - cg);
    lower = 1 - pow(ca, 2) - pow(cb, 2) - pow(cg, 2) + 2 * ca * cb * cg;
    d = sqrt(lower / upper);
    return d;
  }

  double get_stl_of_hkl(const std::vector<int> hkl) const
  {
    return 1.0 / (2 * get_d_of_hkl(hkl));
  }

  void make_coords_cartesian(double* positions_cart, const double frac_x, const double frac_y, const double frac_z, const bool in_bohr = true)
  {
    positions_cart[0] = (a * frac_x + b * cg * frac_y + c * cb * frac_z);
    positions_cart[1] = (b * sg * frac_y + c * (ca - cb * cg) / sg * frac_z);
    positions_cart[2] = V / (a * b * sg) * frac_z;
    if (in_bohr)
      for (int i = 0; i < 3; i++)
        positions_cart[i] /= 0.529177249;
  }

  std::vector<double> get_coords_cartesian(const double frac_x, const double frac_y, const double frac_z, const bool in_bohr = true) const
  {
    std::vector<double> positions_cart{ 0., 0., 0. };
    positions_cart[0] = (a * frac_x + b * cg * frac_y + c * cb * frac_z);
    positions_cart[1] = (b * sg * frac_y + c * (ca - cb * cg) / sg * frac_z);
    positions_cart[2] = V / (a * b * sg) * frac_z;
    if (in_bohr)
      for (int i = 0; i < 3; i++)
        positions_cart[i] /= 0.529177249;
    return positions_cart;
  }

  bool read_CIF(std::string filename, std::ofstream& file, const bool& debug = false)
  {
    std::ifstream cif_input(filename.c_str(), std::ios::in);
    std::vector<bool> found;
    found.resize(7);
    for (int k = 0; k < 7; k++)
      found[k] = false;
    double v = 0.0;
    std::vector <std::string> cell_keywords;
    std::string line;
    cell_keywords.push_back("_cell_length_a");
    cell_keywords.push_back("_cell_length_b");
    cell_keywords.push_back("_cell_length_c");
    cell_keywords.push_back("_cell_angle_alpha");
    cell_keywords.push_back("_cell_angle_beta");
    cell_keywords.push_back("_cell_angle_gamma");
    cell_keywords.push_back("_cell_volume");
    if (debug)
      file << "\nStarting while !.eof()" << std::endl;
    while (!cif_input.eof()) {
      getline(cif_input, line);
      for (int k = 0; k < cell_keywords.size(); k++) {
        if (line.find(cell_keywords[k]) != std::string::npos) {
          switch (k) {
          case 0:
            a = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 1:
            b = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 2:
            c = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 3:
            alpha = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 4:
            beta = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 5:
            gamma = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 6:
            v = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          default:
            file << "This is weird... should never get here... aborting!" << std::endl;
            return false;
          }
          found[k] = true;
        }
      }
      if (found[0] == true && found[1] == true && found[2] == true && found[3] == true && found[4] == true && found[5] == true && found[6] == true)
        break;
    }
    ca = cos(0.0174532925199432944444444444444 * alpha);
    cb = cos(0.0174532925199432944444444444444 * beta);
    cg = cos(0.0174532925199432944444444444444 * gamma);
    sa = sin(0.0174532925199432944444444444444 * alpha);
    sb = sin(0.0174532925199432944444444444444 * beta);
    sg = sin(0.0174532925199432944444444444444 * gamma);
    V = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg);
    if (V / v > 1.1 || V / v < 0.9) {
      file << "Volume computed is more than 10% off, please check!" << std::endl;
      return false;
    }
    as = b * c * sa / V;
    bs = a * c * sb / V;
    cs = a * b * sg / V;

    if (debug)
      file << "Making cm and rcm" << std::endl
      << ca << " " << cb << " " << cg << " " << sa << " " << sb << " " << sg << " " << V << std::endl
      << as << " " << bs << " " << cs << std::endl;

    cm[0][0] = a / 0.529177249;
    cm[0][1] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
    cm[0][2] = sqrt(abs(a * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));

    cm[1][0] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
    cm[1][1] = b / 0.529177249;
    cm[1][2] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));

    cm[2][0] = sqrt(abs(a * c * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
    cm[2][1] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
    cm[2][2] = c / 0.529177249;

    rcm[0][0] = 2 * PI / a;
    rcm[0][1] = 0;
    rcm[0][2] = 0;

    rcm[1][0] = 2 * PI * -cg / (a * sg);
    rcm[1][1] = 2 * PI * 1 / (b * sg);
    rcm[1][2] = 0;

    rcm[2][0] = 2 * PI * b * c * (ca * cg - cb) / V / sg;
    rcm[2][1] = 2 * PI * a * c * (cb * cg - ca) / V / sg;
    rcm[2][2] = 2 * PI * a * b * sg / V;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (abs(rcm[i][j]) < 10e-10) {
          rcm[i][j] = 0.0;
          //cm[i][j] = 0.0;
        }
        else {
          rcm[i][j] = bohr2ang(rcm[i][j]);
          //cm[i][j] *= 0.529177249;
        }

    cif_input.close();
    return true;
  }
  bool read_CIF(std::string filename, const bool& debug = false)
  {
    std::ifstream cif_input(filename.c_str(), std::ios::in);
    std::vector<bool> found;
    found.resize(7);
    for (int k = 0; k < 7; k++)
      found[k] = false;
    double v = 0.0;
    std::vector <std::string> cell_keywords;
    std::string line;
    cell_keywords.push_back("_cell_length_a");
    cell_keywords.push_back("_cell_length_b");
    cell_keywords.push_back("_cell_length_c");
    cell_keywords.push_back("_cell_angle_alpha");
    cell_keywords.push_back("_cell_angle_beta");
    cell_keywords.push_back("_cell_angle_gamma");
    cell_keywords.push_back("_cell_volume");
    if (debug)
      std::cout << "\nStarting while !.eof()" << std::endl;
    while (!cif_input.eof()) {
      getline(cif_input, line);
      for (int k = 0; k < cell_keywords.size(); k++) {
        if (line.find(cell_keywords[k]) != std::string::npos) {
          switch (k) {
          case 0:
            a = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 1:
            b = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 2:
            c = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 3:
            alpha = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 4:
            beta = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 5:
            gamma = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          case 6:
            v = stod(line.substr(cell_keywords[k].length(), line.find("(")));
            break;
          default:
            std::cout << "This is weird... should never get here... aborting!" << std::endl;
            return false;
          }
          found[k] = true;
        }
      }
      if (found[0] == true && found[1] == true && found[2] == true && found[3] == true && found[4] == true && found[5] == true && found[6] == true)
        break;
    }
    ca = cos(0.0174532925199432944444444444444 * alpha);
    cb = cos(0.0174532925199432944444444444444 * beta);
    cg = cos(0.0174532925199432944444444444444 * gamma);
    sa = sin(0.0174532925199432944444444444444 * alpha);
    sb = sin(0.0174532925199432944444444444444 * beta);
    sg = sin(0.0174532925199432944444444444444 * gamma);
    V = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg);
    if (V / v > 1.1 || V / v < 0.9) {
      std::cout << "Volume computed is more than 10% off, please check!" << std::endl;
      return false;
    }
    as = b * c * sa / V;
    bs = a * c * sb / V;
    cs = a * b * sg / V;

    if (debug)
      std::cout << "Making cm and rcm" << std::endl
      << ca << " " << cb << " " << cg << " " << sa << " " << sb << " " << sg << " " << V << std::endl
      << as << " " << bs << " " << cs << std::endl;

    cm[0][0] = a / 0.529177249;
    cm[0][1] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
    cm[0][2] = sqrt(abs(a * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));

    cm[1][0] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
    cm[1][1] = b / 0.529177249;
    cm[1][2] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));

    cm[2][0] = sqrt(abs(a * c * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
    cm[2][1] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
    cm[2][2] = c / 0.529177249;

    rcm[0][0] = 2 * PI / a;
    rcm[0][1] = 0;
    rcm[0][2] = 0;

    rcm[1][0] = 2 * PI * -cg / (a * sg);
    rcm[1][1] = 2 * PI * 1 / (b * sg);
    rcm[1][2] = 0;

    rcm[2][0] = 2 * PI * b * c * (ca * cg - cb) / V / sg;
    rcm[2][1] = 2 * PI * a * c * (cb * cg - ca) / V / sg;
    rcm[2][2] = 2 * PI * a * b * sg / V;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (abs(rcm[i][j]) < 10e-10) {
          rcm[i][j] = 0.0;
          //cm[i][j] = 0.0;
        }
        else {
          rcm[i][j] = bohr2ang(rcm[i][j]);
          //cm[i][j] *= 0.529177249;
        }

    cif_input.close();
    return true;
  }

  void read_symm_CIF(std::string filename, std::ofstream& file, const bool& debug = false)
  {
    std::ifstream cif_input(filename.c_str(), std::ios::in);
    std::string line;
    cif_input.clear();
    cif_input.seekg(0, cif_input.beg);
    bool symm_found = false;
    int operation_field = 200;
    int count_fields = 0;
    while (!cif_input.eof() && !symm_found) {
      getline(cif_input, line);
      if (line.find("loop_") != std::string::npos) {
        //if(debug) file << "found loop!" << endl;
        while (line.find("_") != std::string::npos) {
          getline(cif_input, line);
          if (debug) file << "line in loop field definition: " << line << std::endl;
          if (line.find("space_group_symop_operation_xyz") != std::string::npos)
            operation_field = count_fields;
          else if (count_fields > 2 || (operation_field == 200 && count_fields != 0)) {
            if (debug) file << "I don't think this is the symmetry block.. moving on!" << std::endl;
            count_fields = 0;
            break;
          }
          count_fields++;
        }
        while (line.find("_") == std::string::npos && line.length() > 3 && count_fields != 0) {
          if (debug) file << "Reading operation!" << line << std::endl;
          symm_found = true;
          std::stringstream s(line);
          std::vector <std::string> fields;
          fields.resize(count_fields);
          int sym_from_cif[3][3];
          for (int x = 0; x < 3; x++)
            for (int y = 0; y < 3; y++)
              sym_from_cif[x][y] = 0;
          for (int i = 0; i < count_fields; i++)
            s >> fields[i];
          std::vector<std::string> vectors;
          vectors.resize(3);
          int column = 0;
          for (int c_ = 0; c_ < fields[operation_field].length(); c_++) {
            if (fields[operation_field][c_] != ',')
              vectors[column].push_back(fields[operation_field][c_]);
            else column++;
          }

          for (int x = 0; x < 3; x++) {
            if (vectors[x].find("X") != std::string::npos || vectors[x].find("x") != std::string::npos) {
              char sign = ' ';
              if (vectors[x].find("X") != std::string::npos && vectors[x].find("X") != 0)
                sign = vectors[x].at(vectors[x].find("X") - 1);
              else if (vectors[x].find("X") == 0)
                sign = '+';
              if (vectors[x].find("x") != std::string::npos && vectors[x].find("x") != 0)
                sign = vectors[x].at(vectors[x].find("x") - 1);
              else if (vectors[x].find("x") == 0)
                sign = '+';
              if (sign == '-')
                sym_from_cif[x][0] = -1;
              if (sign == '+')
                sym_from_cif[x][0] = 1;
            }
            if (vectors[x].find("Y") != std::string::npos || vectors[x].find("y") != std::string::npos) {
              char sign = ' ';
              if (vectors[x].find("Y") != std::string::npos && vectors[x].find("Y") != 0)
                sign = vectors[x].at(vectors[x].find("Y") - 1);
              else if (vectors[x].find("Y") == 0)
                sign = '+';
              if (vectors[x].find("y") != std::string::npos && vectors[x].find("y") != 0)
                sign = vectors[x].at(vectors[x].find("y") - 1);
              else if (vectors[x].find("y") == 0)
                sign = '+';
              if (sign == '-')
                sym_from_cif[x][1] = -1;
              if (sign == '+')
                sym_from_cif[x][1] = 1;
            }
            if (vectors[x].find("Z") != std::string::npos || vectors[x].find("z") != std::string::npos) {
              char sign = ' ';
              if (vectors[x].find("Z") != std::string::npos && vectors[x].find("Z") != 0)
                sign = vectors[x].at(vectors[x].find("Z") - 1);
              else if (vectors[x].find("Z") == 0)
                sign = '+';
              if (vectors[x].find("z") != std::string::npos && vectors[x].find("z") != 0)
                sign = vectors[x].at(vectors[x].find("z") - 1);
              else if (vectors[x].find("z") == 0)
                sign = '+';
              if (sign == '-')
                sym_from_cif[x][2] = -1;
              if (sign == '+')
                sym_from_cif[x][2] = 1;
            }
          }
          if (debug) {
            file << "Comparing ";
            for (int x = 0; x < 3; x++)
              for (int y = 0; y < 3; y++)
                file << sym_from_cif[x][y] << " ";
            file << std::endl;
          }
          bool already_known = false;
          for (int s_ = 0; s_ < sym[0][0].size(); s_++) {
            bool identical = true;
            bool inverse = true;
            for (int x = 0; x < 3; x++)
              for (int y = 0; y < 3; y++)
                if (sym[y][x][s_] != sym_from_cif[x][y])
                  identical = false;
            if (!identical)
              for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                  if (sym[y][x][s_] != sym_from_cif[x][y] * -1)
                    inverse = false;
            /*if (debug) {
              file << "Comparison with ";
              for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                  file << sym[x][y][s] << " ";
              file << "resulted in ";
              file << identical << " " << inverse << endl;
            }*/
            if (identical || inverse) {
              already_known = true;
              break;
            }
          }
          if (!already_known) {
            if (debug) file << "This is a new symmetry operation!" << std::endl;
            for (int x = 0; x < 3; x++)
              for (int y = 0; y < 3; y++)
                sym[y][x].push_back(sym_from_cif[x][y]);
          }
          getline(cif_input, line);
        }
      }
    }
  };

  void read_symm_CIF(std::string filename, const bool& debug = false)
  {
    std::ifstream cif_input(filename.c_str(), std::ios::in);
    std::string line;
    cif_input.clear();
    cif_input.seekg(0, cif_input.beg);
    bool symm_found = false;
    int operation_field = 200;
    int count_fields = 0;
    while (!cif_input.eof() && !symm_found) {
      getline(cif_input, line);
      if (line.find("loop_") != std::string::npos) {
        //if(debug) file << "found loop!" << endl;
        while (line.find("_") != std::string::npos) {
          getline(cif_input, line);
          if (debug) std::cout << "line in loop field definition: " << line << std::endl;
          if (line.find("space_group_symop_operation_xyz") != std::string::npos)
            operation_field = count_fields;
          else if (count_fields > 2 || (operation_field == 200 && count_fields != 0)) {
            if (debug) std::cout << "I don't think this is the symmetry block.. moving on!" << std::endl;
            count_fields = 0;
            break;
          }
          count_fields++;
        }
        while (line.find("_") == std::string::npos && line.length() > 3 && count_fields != 0) {
          if (debug) std::cout << "Reading operation!" << line << std::endl;
          symm_found = true;
          std::stringstream s(line);
          std::vector <std::string> fields;
          fields.resize(count_fields);
          int sym_from_cif[3][3];
          for (int x = 0; x < 3; x++)
            for (int y = 0; y < 3; y++)
              sym_from_cif[x][y] = 0;
          for (int i = 0; i < count_fields; i++)
            s >> fields[i];
          std::vector<std::string> vectors;
          vectors.resize(3);
          int column = 0;
          for (int c = 0; c < fields[operation_field].length(); c++) {
            if (fields[operation_field][c] != ',')
              vectors[column].push_back(fields[operation_field][c]);
            else column++;
          }

          for (int x = 0; x < 3; x++) {
            if (vectors[x].find("X") != std::string::npos || vectors[x].find("x") != std::string::npos) {
              char sign = ' ';
              if (vectors[x].find("X") != std::string::npos && vectors[x].find("X") != 0)
                sign = vectors[x].at(vectors[x].find("X") - 1);
              else if (vectors[x].find("X") == 0)
                sign = '+';
              if (vectors[x].find("x") != std::string::npos && vectors[x].find("x") != 0)
                sign = vectors[x].at(vectors[x].find("x") - 1);
              else if (vectors[x].find("x") == 0)
                sign = '+';
              if (sign == '-')
                sym_from_cif[x][0] = -1;
              else if (sign == '+')
                sym_from_cif[x][0] = 1;
            }
            if (vectors[x].find("Y") != std::string::npos || vectors[x].find("y") != std::string::npos) {
              char sign = ' ';
              if (vectors[x].find("Y") != std::string::npos && vectors[x].find("Y") != 0)
                sign = vectors[x].at(vectors[x].find("Y") - 1);
              else if (vectors[x].find("Y") == 0)
                sign = '+';
              if (vectors[x].find("y") != std::string::npos && vectors[x].find("y") != 0)
                sign = vectors[x].at(vectors[x].find("y") - 1);
              else if (vectors[x].find("y") == 0)
                sign = '+';
              if (sign == '-')
                sym_from_cif[x][1] = -1;
              else if (sign == '+')
                sym_from_cif[x][1] = 1;
            }
            if (vectors[x].find("Z") != std::string::npos || vectors[x].find("z") != std::string::npos) {
              char sign = ' ';
              if (vectors[x].find("Z") != std::string::npos && vectors[x].find("Z") != 0)
                sign = vectors[x].at(vectors[x].find("Z") - 1);
              else if (vectors[x].find("Z") == 0)
                sign = '+';
              if (vectors[x].find("z") != std::string::npos && vectors[x].find("z") != 0)
                sign = vectors[x].at(vectors[x].find("z") - 1);
              else if (vectors[x].find("z") == 0)
                sign = '+';
              if (sign == '-')
                sym_from_cif[x][2] = -1;
              else if (sign == '+')
                sym_from_cif[x][2] = 1;
            }
          }
          if (debug) {
            std::cout << "Comparing ";
            for (int x = 0; x < 3; x++)
              for (int y = 0; y < 3; y++)
                std::cout << sym_from_cif[x][y] << " ";
            std::cout << std::endl;
          }
          bool already_known = false;
          for (int s = 0; s < sym[0][0].size(); s++) {
            bool identical = true;
            bool inverse = true;
            for (int x = 0; x < 3; x++)
              for (int y = 0; y < 3; y++)
                if (sym[y][x][s] != sym_from_cif[x][y])
                  identical = false;
            if (!identical)
              for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                  if (sym[y][x][s] != sym_from_cif[x][y] * -1)
                    inverse = false;
            /*if (debug) {
              file << "Comparison with ";
              for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++)
                  file << sym[x][y][s] << " ";
              file << "resulted in ";
              file << identical << " " << inverse << endl;
            }*/
            if (identical || inverse) {
              already_known = true;
              break;
            }
          }
          if (!already_known) {
            if (debug) std::cout << "This is a new symmetry operation!" << std::endl;
            for (int x = 0; x < 3; x++)
              for (int y = 0; y < 3; y++)
                sym[y][x].push_back(sym_from_cif[x][y]);
          }
          getline(cif_input, line);
        }
      }
    }
  };
};