#pragma once

#include <vector>
#include <string>

class MO
{
private:
  int nr;
  double occ;
  double ener;
  std::vector<double> coefficients;
public:
  MO() {
    nr = 0;
    occ = 0.0;
    ener = 0.0;
  };
  MO(const int& number, const double& occupation, const double& energy) {
    nr = number;
    occ = occupation;
    ener = energy;
  };
  void push_back_coef(const double &val){
    coefficients.push_back(val);
  };
  bool erase_coef(const int &nr, const int &nex) {
    // Check, if Element (nr) is a valid Vector-Entry
    if (nr - 1 < 0)
      return false;
    if (size_t(nr) - 1 >= coefficients.size())
      return false;

    // consistencycheck with WFN
    if (nex != coefficients.size() - 1) 
      return false;

    // delete Vector-Entry and rearrange
    coefficients.erase(coefficients.begin() + (size_t(nr) - 1));

    return true;
  };
  double get_coefficient(const int& nr) const {
    if (nr >= coefficients.size()) {
      err_checkf(false, "Requested element outside of range!", std::cout);
    }
    if (nr < 0) {
      err_checkf(false, "Number below zero!", std::cout);
    }
    return coefficients[nr];
  }; 
  double get_coefficient_f(const int& nr) const {
    return coefficients[nr];
  };
  int get_nr() const {
    return coefficients[nr];
  };
  double* get_coefficient_ptr() {
    return coefficients.data();
  };
  void change_coefficient(int nr) {
    // Check, if Element (nr) is a valid Vector-Entry
    if (nr - 1 < 0)
      return;
    if (size_t(nr) - 1 >= coefficients.size())
      return;

    bool end = false;
    while (!end) {
      std::cout << "What is the value for the new exponent?\n";
      double temp_coef = 0.0;
      std::cin >> temp_coef;
      if (temp_coef > 1000 || temp_coef < -1000) {
        std::cout << "That's suspiciusly big for a MO-coefficient, please select something smaler...\n";
        continue;
      }
      coefficients[size_t(nr) - 1] = temp_coef;
      end = true;
    }
    cls();
  };
  bool change_coefficient(int nr, double value) {
    // Check, if Element (nr) is a valid Vector-Entry
    if (nr < 0) {
      err_checkf(false,"nr below 0!", std::cout);
      return false;
    }
    if (nr >= coefficients.size()) {
      err_checkf(false, "nr above size of MO!", std::cout);
      return false;
    }
    coefficients[nr] = value;
    return true;
  };
  void set_nr(const int inr) { nr = inr; };
  void set_occ(const int iocc) { occ = iocc; };
  double get_occ() const { return occ; };
  void set_ener(const double iener) { ener = iener; };
  int get_primitive_count() const { return coefficients.size(); };
  std::string hdr() {
    std::string temp;
    temp = "MO";
    if (nr < 10 && nr>0) temp.append("    ");
    else if (nr < 100) temp.append("   ");
    else temp.append("  ");
    temp.append(std::to_string(nr));
    temp.append("                  OCC NO =");
    if (occ == 2.0) temp.append("    2.00000000 ORB. ENERGY =");
    else if (occ == 0.0) temp.append("    0.00000000 ORB. ENERGY =");
    else if (occ < 2 && occ>0) {
      temp.append("    ");
      temp.append(std::to_string(occ));
      temp.append("  ORB. ENERGY =");
    }
    if (ener<0 && ener>-10) {
      temp.append("   ");
      temp.append(std::to_string(ener));
      temp.append("\n");
    }
    else if (ener<-10 && ener>-100) {
      temp.append("  ");
      temp.append(std::to_string(ener));
      temp.append("\n");
    }
    else if (ener < -100) {
      temp.append(" ");
      temp.append(std::to_string(ener));
      temp.append("\n");
    }
    else if (ener > 0) {
      temp.append("    ");
      temp.append(std::to_string(ener));
      temp.append("\n");
    }
    return temp;
  }
  double get_energy() const {
    return ener;
  };
  double* get_ptr_coefficients() { return &coefficients[0]; };
  const std::vector<double>& get_ptr_coef_vector() const { return coefficients; };
  void assign_coefs(const std::vector<double>& values) { coefficients.resize(values.size()); coefficients = values; }
};