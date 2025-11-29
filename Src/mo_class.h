#pragma once

#include <vector>
#include <string>

class MO
{
private:
  int nr; //acts as an index
  double occ;
  double ener;
  vec coefficients;
  int op; //0=alpha, 1=beta
public:
  MO() {
    nr = 0;
    op = 0;
    occ = 0.0;
    ener = 0.0;
  };
  MO(const int& number, const double& occupation, const double& energy, const int& oper = 0) {
    nr = number;
    occ = occupation;
    ener = energy;
    op = oper;
  };
  virtual ~MO() {};
  friend class MO_OCC;
  void push_back_coef(const double &val){
    coefficients.emplace_back(val);
  };
  const bool erase_coef(const int &_nr, const int &nex) {
    // Check, if Element (nr) is a valid Vector-Entry
    if (_nr - 1 < 0)
      return false;
    if (size_t(_nr) - 1 >= coefficients.size())
      return false;

    // consistencycheck with WFN
    if (nex != coefficients.size() - 1)
      return false;

    // delete Vector-Entry and rearrange
    coefficients.erase(coefficients.begin() + (size_t(_nr) - 1));

    return true;
  };
  const double& get_coefficient(const int& _nr) const {
    err_checkf(_nr < coefficients.size() && _nr >= 0, "Requested element outside of range! " + std::to_string(_nr), std::cout);
    return coefficients[_nr];
  };
  const vec& get_coefficients() const {
    return coefficients;
  };
  void assign_coefficients_size(int size)
  {
    coefficients.assign(size, 0.0);
  }
  void set_coefficients(vec coeff) {
    coefficients = coeff;
  };
  void insert_into_coefficients(std::ranges::input_range auto&& v) {
    coefficients.insert(coefficients.end(), std::ranges::begin(v), std::ranges::end(v));
  }


  const double& get_coefficient_f(const int& _nr) const {
    return coefficients[_nr];
  };
  double* get_coefficient_ptr() {
    return coefficients.data();
  };
  const bool set_coefficient(const int& _nr, const double& value) {
    // Check, if Element (nr) is a valid Vector-Entry
    if (_nr < 0) {
      err_checkf(false,"nr below 0!", std::cout);
      return false;
    }
    if (_nr >= (int) coefficients.size()) {
      err_checkf(false, "nr above size of MO!", std::cout);
      return false;
    }
    coefficients[_nr] = value;
    return true;
  };
  void set_nr(const int& inr) { nr = inr; };
  void set_occ(const int& iocc) { occ = iocc; };
  void set_occ(const double& iocc) { occ = iocc; };
  void set_op(const int& oper) { op = oper; };
  const double& get_occ() const { return occ; };
  const int& get_op() const { return op; };
  void set_ener(const double& iener) { ener = iener; };
  const int get_primitive_count() const { return (int) coefficients.size(); };
  const std::string hdr() {
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
  const double& get_energy() const {
    return ener;
  };
  const vec& get_ptr_coef_vector() const { return coefficients; };
  void assign_coefs(const vec& values) { coefficients.resize(values.size()); coefficients = values; }
};
