#pragma once

#ifndef __TSC_BLOCK_H__
#define __TSC_BLOCK_H__

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <complex>
#include "convenience.h"

class tsc_block
{
private:
  std::vector<std::vector<std::complex<double> > > sf; //[#scatterer] [#reflection]
  std::vector<std::string> scatterer; //Labels of reflections the correpsonding entry of sf belonds to
  std::vector<std::vector<int> > index; //[3] [index]
  std::string header;
  bool anomalous_dispersion;

public:
  template<typename numtype>
  tsc_block(
    std::vector<std::vector<numtype> >& given_sf,
    std::vector<std::string>& given_scatterer,
    std::vector<std::vector<int> >& given_index,
    std::string& given_header)
  {
    sf.resize(given_sf.size());
#pragma omp parallel for
    for (int i = 0; i < given_sf.size(); i++) {
      sf[i].resize(given_sf[i].size());
      for (int j = 0; j < given_sf[i].size(); j++) {
        sf[i][j] = given_sf[i][j];
      }
    }
    scatterer = given_scatterer;
    index = given_index;
    header = given_header;
    anomalous_dispersion = false;
  };
  template<typename numtype>
  tsc_block(
    std::vector<std::vector<numtype> >& given_sf,
    std::vector<std::string>& given_scatterer,
    std::vector<std::vector<int> >& given_index)
  {
    sf.resize(given_sf.size());
#pragma omp parallel for
    for (int i = 0; i < given_sf.size(); i++) {
      sf[i].resize(given_sf[i].size());
      for (int j = 0; j < given_sf[i].size(); j++) {
        sf[i][j] = given_sf[i][j];
      }
    }
    scatterer = given_scatterer;
    index = given_index;
    anomalous_dispersion = false;
  };
  template<typename numtype>
  tsc_block(
    std::vector<std::vector<numtype> >& given_sf,
    std::vector<std::string>& given_scatterer,
    hkl_list& given_index)
  {
    sf.resize(given_sf.size());
#pragma omp parallel for
    for (int i = 0; i < given_sf.size(); i++) {
      sf[i].resize(given_sf[i].size());
      for (int j = 0; j < given_sf[i].size(); j++) {
        sf[i][j] = given_sf[i][j];
      }
    }
    scatterer = given_scatterer;
    index.resize(3);
    for (const std::vector<int>& hkl : given_index) {
      for (int i = 0; i < 3; i++)
        index[i].push_back(hkl[i]);
    }
    anomalous_dispersion = false;
  };
  tsc_block() { anomalous_dispersion = false; };
  ~tsc_block()
  {
    for (int i = 0; i < sf.size(); i++) {
      sf[i].clear();
      std::vector<std::complex<double>>(sf[i]).swap(sf[i]);
    }
    sf.clear();
    shrink_vector<std::string>(scatterer);
    for (int i = 0; i < index.size(); i++) {
      shrink_vector<int>(index[i]);
    }
  };
  const std::vector<std::complex<double> > get_sf_for_scatterer(const unsigned int nr)
  {
    err_checkc(nr < scatterer.size(), "Wrong number in get SF");
    return sf[nr];
  };
  const std::vector<std::complex<double> > get_sf_for_scatterer(const unsigned int nr, std::ofstream& log)
  {
    err_checkf(nr < scatterer.size(), "Wrong number in get SF", log);
    return sf[nr];
  };
  const std::string get_scatterer(const unsigned int nr)
  {
    err_checkc(nr < scatterer.size(), "Invalid nr of scatterer");
    return scatterer[nr];
  };
  const std::string get_scatterer(const unsigned int nr, std::ofstream& log)
  {
    err_checkf(nr < scatterer.size(), "Invalid nr of scatterer", log);
    return scatterer[nr];
  };
  const std::vector<std::string> get_scatterers() { return scatterer; }
  void set_AD(const bool value) { anomalous_dispersion = value; };
  bool get_AD() { return anomalous_dispersion; };
  const std::vector<int> get_indices(const unsigned int nr)
  {
    err_checkc(nr < index[0].size(), "Invalid nr of index");
    return { index[0][nr], index[1][nr], index[2][nr] };
  };
  const std::vector<int> get_indices(const unsigned int nr, std::ofstream& file)
  {
    err_checkf(nr < index[0].size(), "Invalid nr of index", file);
    return { index[0][nr], index[1][nr], index[2][nr] };
  };
  const int get_index(const unsigned int dim, const unsigned int nr)
  {
    err_checkc(dim < 3, "invalid dimension for index");
    err_checkc(nr < index[dim].size(), "invalid nr for index");
    return index[dim][nr];
  };
  const int get_index(const unsigned int dim, const unsigned int nr, std::ofstream& file)
  {
    err_checkf(dim < 3, "invalid dimension for index", file);
    err_checkf(nr < index[dim].size(), "invalid nr for index", file);
    return index[dim][nr];
  };
  const bool is_empty() { return (sf.size() > 0 && scatterer.size() > 0 && index.size() > 0); };
  const unsigned int scatterer_size() { 
    if (sf.size() == scatterer.size()) 
      return sf.size(); 
    else 
      return 0; 
  };
  const unsigned int reflection_size() { 
    if (sf.size() == 0 || index.size() == 0) {
      return 0;
    }
    else if (sf[0].size() == index[0].size() && index[0].size() == index[1].size() && index[1].size() == index[2].size()) {
      return index[0].size();
    }
    else {
      return 0;
    }
  }
  void append(tsc_block& rhs, std::ofstream& log)
  {
    if (reflection_size() == 0) {
      *this = rhs;
      return;
    }
    //Appends the scatterers of rhs to the current set assuming same size of reflections.
    err_checkf(reflection_size() == rhs.reflection_size(), "Inconsistent number of reflections!", log);
    err_checkf(rhs.reflection_size() > 0, "Nothing to append or inconsinstency in given block detected, then please don't do it!", log);
#pragma omp parallel for
    for (int i = 0; i < rhs.reflection_size(); i++)
      for (int dim = 0; dim < 3; dim++)
        err_checkf(index[dim][i] == rhs.get_index(dim, i), "Mismatch in indices in append!", log);
    int new_scatterers = 0;
    std::vector<bool> is_new(rhs.scatterer_size(), true);
#pragma omp parallel for reduction(+:new_scatterers)
    for (int s = 0; s < rhs.scatterer_size(); s++) {
      for (unsigned int run = 0; run < scatterer_size(); run++)
        if (rhs.get_scatterer(s) == scatterer[run])
          is_new[s] = false;
      if (is_new[s] == false) continue;
      new_scatterers++;
    }
    const unsigned int old_size = sf.size();
    sf.resize(old_size + new_scatterers);
    scatterer.resize(old_size + new_scatterers);
#pragma omp parallel for
    for (int s = 0; s < rhs.scatterer_size(); s++) {
      if (is_new[s]) {
        unsigned int new_nr = old_size;
        for (int run = 0; run < s; run++)
          if (is_new[s]) new_nr++;
        sf[new_nr] = rhs.get_sf_for_scatterer(s, log);
        scatterer[new_nr] = rhs.get_scatterer(s, log);
      }
    }
  };
  void append(tsc_block rhs, std::ofstream& log)
  {
    if (reflection_size() == 0) {
      *this = rhs;
      return;
    }
    //Appends the scatterers of rhs to the current set assuming same size of reflections.
    err_checkf(reflection_size() == rhs.reflection_size(), "Inconsistent number of reflections!", log);
    err_checkf(rhs.reflection_size() > 0, "Nothing to append or inconsinstency in given block detected, then please don't do it!", log);
#pragma omp parallel for
    for (int i = 0; i < rhs.reflection_size(); i++)
      for (int dim = 0; dim < 3; dim++)
        err_checkf(index[dim][i] == rhs.get_index(dim, i), "Mismatch in indices in append!", log);
    int new_scatterers = 0;
    std::vector<bool> is_new(rhs.scatterer_size(), true);
#pragma omp parallel for reduction(+:new_scatterers)
    for (int s = 0; s < rhs.scatterer_size(); s++) {
      for (unsigned int run = 0; run < scatterer_size(); run++)
        if (rhs.get_scatterer(s) == scatterer[run])
          is_new[s] = false;
      if (is_new[s] == false) continue;
      new_scatterers++;
    }
    const unsigned int old_size = sf.size();
    sf.resize(old_size + new_scatterers);
    scatterer.resize(old_size + new_scatterers);
#pragma omp parallel for
    for (int s = 0; s < rhs.scatterer_size(); s++) {
      if (is_new[s]) {
        unsigned int new_nr = old_size;
        for (int run = 0; run < s; run++)
          if (is_new[s]) new_nr++;
        sf[new_nr] = rhs.get_sf_for_scatterer(s, log);
        scatterer[new_nr] = rhs.get_scatterer(s, log);
      }
    }
  };
  void write_tsc_file(std::string& cif)
  {
    std::ofstream tsc_file("experimental.tsc", std::ios::out);

    tsc_file << "TITLE: " << get_filename_from_path(cif).substr(0, cif.find(".cif")) << std::endl << "SYMM: ";
    tsc_file << "expanded";
    if (anomalous_dispersion) tsc_file << std::endl << "AD: TRUE";
    tsc_file << std::endl << "SCATTERERS:";
    for (int i = 0; i < scatterer.size(); i++)
      tsc_file << " " << scatterer[i];
    tsc_file << std::endl << "DATA:" << std::endl;

    for (int r = 0; r < index[0].size(); r++) {
      for (int h = 0; h < 3; h++)
        tsc_file << index[h][r] << " ";
      for (int i = 0; i < sf.size(); i++)
        tsc_file << std::scientific << std::setprecision(8) << real(sf[i][r]) << ","
        << std::scientific << std::setprecision(8) << imag(sf[i][r]) << " ";
      tsc_file << std::endl;
    }
    tsc_file.close();
    err_checkc(!tsc_file.bad(), "Error during writing of tsc file!");
  }
};


#endif