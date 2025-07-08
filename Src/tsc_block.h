#pragma once

#include "convenience.h"
#include <unordered_set>
template <typename numtype_index, typename numtype>
class tsc_block
{
private:
  cvec2 sf; //[#scatterer] [#reflection]
  svec scatterer;                // Labels of reflections the correpsonding entry of sf belonds to
  std::vector<std::vector<numtype_index>> index;     //[3] [miller_index]
  std::string header;
  bool anomalous_dispersion;

public:
  tsc_block(
      std::vector<std::vector<numtype>> &given_sf,
      svec &given_scatterer,
      std::vector<std::vector<numtype_index>> &given_index,
      std::string &given_header)
  {
    sf.resize(given_sf.size());
#pragma omp parallel for
    for (int i = 0; i < given_sf.size(); i++)
    {
      sf[i].resize(given_sf[i].size());
      for (int j = 0; j < given_sf[i].size(); j++)
      {
        err_checkf(is_nan(given_sf[i][j]) == false, "NaN in SF!", std::cout);
        sf[i][j] = given_sf[i][j];
      }
    }
    scatterer = given_scatterer;
    index = given_index;
    header = given_header;
    anomalous_dispersion = false;
  };
  tsc_block(
      std::vector<std::vector<numtype>> &given_sf,
      svec &given_scatterer,
      std::vector<std::vector<numtype_index>> &given_index)
  {
    sf.resize(given_sf.size());
#pragma omp parallel for
    for (int i = 0; i < given_sf.size(); i++)
    {
      sf[i].resize(given_sf[i].size());
      for (int j = 0; j < given_sf[i].size(); j++)
      {
        err_checkf(is_nan(given_sf[i][j]) == false, "NaN in SF!", std::cout);
        sf[i][j] = given_sf[i][j];
      }
    }
    scatterer = given_scatterer;
    index = given_index;
    anomalous_dispersion = false;
  };
  tsc_block(
      std::vector<std::vector<numtype>> &given_sf,
      svec &given_scatterer,
      hkl_list_d &given_index)
  {
    sf.resize(given_sf.size());
#pragma omp parallel for
    for (int i = 0; i < given_sf.size(); i++)
    {
      sf[i].resize(given_sf[i].size());
      for (int j = 0; j < given_sf[i].size(); j++)
      {
        err_checkf(is_nan(given_sf[i][j]) == false, "NaN in SF!", std::cout);
        sf[i][j] = given_sf[i][j];
      }
    }
    scatterer = given_scatterer;
    index.resize(3);
    for (const hkl_d &hkl : given_index)
    {
      for (int i = 0; i < 3; i++)
        index[i].push_back(static_cast<numtype_index>(hkl[i]));
    }
    anomalous_dispersion = false;
  };
  tsc_block(
      std::vector<std::vector<numtype>> &given_sf,
      svec &given_scatterer,
      hkl_list &given_index)
  {
    sf.resize(given_sf.size());
#pragma omp parallel for
    for (int i = 0; i < given_sf.size(); i++)
    {
      sf[i].resize(given_sf[i].size());
      for (int j = 0; j < given_sf[i].size(); j++)
      {
        err_checkf(is_nan(given_sf[i][j]) == false, "NaN in SF!", std::cout);
        sf[i][j] = given_sf[i][j];
      }
    }
    scatterer = given_scatterer;
    index.resize(3);
    for (const hkl_t &hkl : given_index)
    {
      for (int i = 0; i < 3; i++)
        index[i].push_back(hkl[i]);
    }
    anomalous_dispersion = false;
  };
  tsc_block(std::filesystem::path &file_name)
  {
    anomalous_dispersion = false;
    using namespace std;
    ifstream tsc_file(file_name.c_str(), ios::binary);

    auto charsize = sizeof(char);
    int head[1]{0};
    auto intsize = sizeof(head);
    auto indexsize = sizeof(numtype_index);
    tsc_file.read((char *)&head, intsize);
    char *cheader;
    string header_str;
    if (head[0] != 0)
    {
      cheader = new char[head[0]];
      tsc_file.read(cheader, head[0] * charsize);
      header_str = cheader;
      delete[] cheader;
    }
    header = header_str;
    // read scatterer labels and map onto scattterers list
    int sc_len[1]{0};
    tsc_file.read((char *)&sc_len, intsize);
    vector<char> scat_line(sc_len[0]);
    tsc_file.read((char *)scat_line.data(), sc_len[0] * charsize);
    string scat_str(scat_line.begin(), scat_line.end());
    // scat_str.resize(sc_len[0]);
    scatterer = split_string<string>(scat_str, string(" "));
    const int nr_scatterers = (int)scatterer.size();
    // read number of indices in tscb file
    int nr_hkl[1]{0};
    tsc_file.read((char *)&nr_hkl, intsize);
    // read indices and scattering factors row by row
    numtype_index rindex[3]{0, 0, 0};
    sf.resize(nr_scatterers);
    index.resize(3);
    for (int run = 0; run < *nr_hkl; run++)
    {
      tsc_file.read((char *)&rindex, 3 * indexsize);
      for (int i = 0; i < 3; i++)
      {
        index[i].push_back(rindex[i]);
      }
      vector<complex<double>> row(nr_scatterers);
      tsc_file.read((char *)row.data(), nr_scatterers * sizeof(complex<double>));
      for (int i = 0; i < nr_scatterers; i++)
      {
        sf[i].push_back(row[i]);
      }
    }
    tsc_file.close();
    err_checkf(!tsc_file.bad(), "TSCB file went bad!", std::cout);
  };
  tsc_block() { anomalous_dispersion = false; };
  ~tsc_block()
  {
    for (int i = 0; i < sf.size(); i++)
    {
      sf[i].clear();
      std::vector<std::complex<double>>(sf[i]).swap(sf[i]);
    }
    sf.clear();
    shrink_vector<std::string>(scatterer);
    for (int i = 0; i < index.size(); i++)
    {
      shrink_vector<numtype_index>(index[i]);
    }
  };
  const cvec get_sf_for_scatterer(const unsigned int nr)
  {
    err_checkf(nr < scatterer.size(), "Wrong number in get SF", std::cout);
    return sf[nr];
  };
  const cvec get_sf_for_scatterer(const unsigned int nr, std::ostream &log)
  {
    err_checkf(nr < scatterer.size(), "Wrong number in get SF", log);
    return sf[nr];
  };
  const std::string get_scatterer(const unsigned int nr)
  {
    err_checkf(nr < scatterer.size(), "Invalid nr of scatterer", std::cout);
    return scatterer[nr];
  };
  const std::string get_scatterer(const unsigned int nr, std::ostream &log)
  {
    err_checkf(nr < scatterer.size(), "Invalid nr of scatterer", log);
    return scatterer[nr];
  };
  const svec get_scatterers() const { return scatterer; }
  const std::string scatterers_string()
  {
    std::string result;
    for (int i = 0; i < scatterer.size(); i++)
    {
      result.append(scatterer[i]);
      if (i != scatterer.size() - 1)
        result.append(" ");
    }
    return result;
  }
  void set_AD(const bool value) { anomalous_dispersion = value; };
  bool get_AD() { return anomalous_dispersion; };
  const ivec get_indices(const unsigned int nr)
  {
    err_checkf(nr < index[0].size(), "Invalid nr of index", std::cout);
    return {index[0][nr], index[1][nr], index[2][nr]};
  };
  const ivec get_indices(const unsigned int nr, std::ofstream &file)
  {
    err_checkf(nr < index[0].size(), "Invalid nr of index", file);
    return {index[0][nr], index[1][nr], index[2][nr]};
  };
  const int get_index(const unsigned int dim, const unsigned int nr)
  {
    err_checkf(dim < 3, "invalid dimension for index", std::cout);
    err_checkf(nr < index[dim].size(), "invalid nr for index", std::cout);
    return index[dim][nr];
  };
  const int get_index(const unsigned int dim, const unsigned int nr, std::ofstream &file)
  {
    err_checkf(dim < 3, "invalid dimension for index", file);
    err_checkf(nr < index[dim].size(), "invalid nr for index", file);
    return index[dim][nr];
  };
  const bool is_empty() { return (sf.size() > 0 && scatterer.size() > 0 && index.size() > 0); };
  const int scatterer_size()
  {
    if (sf.size() == scatterer.size())
      return (int)sf.size();
    else
      return 0;
  };
  const int reflection_size()
  {
    if (sf.size() == 0 || index.size() == 0)
    {
      return 0;
    }
    else if (sf[0].size() == index[0].size() && index[0].size() == index[1].size() && index[1].size() == index[2].size())
    {
      return (int)index[0].size();
    }
    else
    {
      return 0;
    }
  }
  const std::vector<std::vector<numtype_index>> get_index_vector() const { return index; };
  void append(tsc_block &rhs, std::ostream &log)
  {
    if (reflection_size() == 0)
    {
      *this = rhs;
      return;
    }
    // Appends the scatterers of rhs to the current set assuming same size of reflections.
    err_checkf(reflection_size() == rhs.reflection_size(), "Inconsistent number of reflections!", log);
    err_checkf(rhs.reflection_size() > 0, "Nothing to append or inconsinstency in given block detected, then please don't do it!", log);
#pragma omp parallel for
    for (int i = 0; i < rhs.reflection_size(); i++)
      for (int dim = 0; dim < 3; dim++)
        err_checkf(index[dim][i] == rhs.get_index(dim, i), "Mismatch in indices in append!", log);
    int new_scatterers = 0;
    bvec is_new(rhs.scatterer_size(), true);
//#pragma omp parallel for reduction(+ : new_scatterers)
    for (int s = 0; s < (int)rhs.scatterer_size(); s++)
    {
      for (int run = 0; run < scatterer_size(); run++)
        if (rhs.get_scatterer(s) == scatterer[run])
          is_new[s] = false;
      if (is_new[s] == false)
        continue;
      new_scatterers++;
    }
    const int old_size = (int)sf.size();
    sf.resize(size_t(old_size + new_scatterers));
    scatterer.resize(size_t(old_size + new_scatterers));
#pragma omp parallel for
    for (int s = 0; s < rhs.scatterer_size(); s++)
    {
      if (is_new[s])
      {
        unsigned int new_nr = old_size;
        for (int run = 0; run < s; run++)
          if (is_new[s])
            new_nr++;
        sf[new_nr] = rhs.get_sf_for_scatterer(s, log);
        scatterer[new_nr] = rhs.get_scatterer(s, log);
      }
    }
    size_t sc_sf = sf.size();
    size_t nr_hkl_sf = sf[0].size();
    for(size_t i=0; i<sc_sf; i++) {
      if (sf[i].size() != nr_hkl_sf) {
        std::cerr << "Error: Inconsistent size in sf for scatterer " << i << "in append" << std::endl;
        return;
      }
    }
  };
  void append(tsc_block rhs, std::ostream &log)
  {
    if (reflection_size() == 0)
    {
      *this = rhs;
      return;
    }
    // Appends the scatterers of rhs to the current set assuming same size of reflections.
    err_checkf(reflection_size() == rhs.reflection_size(), "Inconsistent number of reflections!", log);
    err_checkf(rhs.reflection_size() > 0, "Nothing to append or inconsinstency in given block detected, then please don't do it!", log);
#pragma omp parallel for
    for (int i = 0; i < rhs.reflection_size(); i++)
      for (int dim = 0; dim < 3; dim++)
        err_checkf(index[dim][i] == rhs.get_index(dim, i), "Mismatch in indices in append!", log);

    std::unordered_set<std::string> existing_set(scatterer.begin(), scatterer.end());
    bvec is_new(rhs.scatterer_size(), false);
    std::atomic<int> new_scatterers(0);
//#pragma omp parallel for
    for (int s = 0; s < rhs.scatterer_size(); ++s)
    {
        const std::string& name = rhs.get_scatterer(s);

        if (existing_set.find(name) == existing_set.end())
        {
            is_new[s] = true;
            new_scatterers++;
        }
    }

    const unsigned int old_size = (int)sf.size();
    sf.resize(old_size + new_scatterers);
    scatterer.resize(old_size + new_scatterers);
#pragma omp parallel for
    for (int s = 0; s < rhs.scatterer_size(); s++)
    {
      if (is_new[s])
      {
        unsigned int new_nr = old_size;
        for (int run = 0; run < s; run++)
          if (is_new[s])
            new_nr++;
        sf[new_nr] = rhs.get_sf_for_scatterer(s, log);
        scatterer[new_nr] = rhs.get_scatterer(s, log);
      }
    }
    int64_t sc_sf = sf.size();
    int64_t nr_hkl_sf = sf[0].size();
    for(int64_t i=0; i<sc_sf; i++) {
      if (sf[i].size() != nr_hkl_sf) {
        std::cerr << "Error: Inconsistent size in sf for scatterer " << i << "in append2" << std::endl;
        return;
      }
    }
  };
  void write_tsc_file(const std::filesystem::path &cif, std::filesystem::path name = "experimental.tsc")
  {
    std::ofstream tsc_file(name, std::ios::out);

    tsc_file << "TITLE: " << cif.stem().string() << std::endl
             << "SYMM: ";
    tsc_file << "expanded";
    if (anomalous_dispersion)
      tsc_file << std::endl
               << "AD: TRUE";
    tsc_file << std::endl
             << "SCATTERERS:";
    for (int i = 0; i < scatterer.size(); i++)
      tsc_file << " " << scatterer[i];
    tsc_file << std::endl
             << "DATA:" << std::endl;

    for (int r = 0; r < index[0].size(); r++)
    {
      for (int h_loc = 0; h_loc < 3; h_loc++)
        tsc_file << index[h_loc][r] << " ";
      for (int i = 0; i < sf.size(); i++)
      {
        tsc_file << std::scientific << std::setprecision(8) << real(sf[i][r]) << ","
                 << std::scientific << std::setprecision(8) << imag(sf[i][r]) << " ";
      }
      tsc_file << std::endl;
    }
    tsc_file.close();
    err_checkf(!tsc_file.bad(), "Error during writing of tsc file!", std::cout);
  }
  void write_tsc_file_non_integer(const std::filesystem::path &cif, std::filesystem::path name = "experimental.tsc")
  {
    std::ofstream tsc_file(name, std::ios::out);

    tsc_file << "TITLE: " << cif.stem() << std::endl
             << "SYMM: ";
    tsc_file << "expanded";
    if (anomalous_dispersion)
      tsc_file << std::endl
               << "AD: TRUE";
    tsc_file << std::endl
             << "SCATTERERS:";
    for (int i = 0; i < scatterer.size(); i++)
      tsc_file << " " << scatterer[i];
    tsc_file << std::endl
             << "DATA:" << std::endl;

    for (int r = 0; r < index[0].size(); r++)
    {
      for (int h_loc = 0; h_loc < 3; h_loc++)
        tsc_file << std::fixed << std::setprecision(3) << index[h_loc][r] << " ";
      for (int i = 0; i < sf.size(); i++)
        tsc_file << std::scientific << std::setprecision(8) << real(sf[i][r]) << ","
                 << std::scientific << std::setprecision(8) << imag(sf[i][r]) << " ";
      tsc_file << std::endl;
    }
    tsc_file.close();
    err_checkf(!tsc_file.bad(), "Error during writing of tsc file!", std::cout);
  }
  void write_tscb_file(std::filesystem::path cif_name = "test.cif", std::filesystem::path name = "experimental.tscb")
  {
    try {  // Wrap in try-catch to help identify where segfaults occur
      //std::cerr << "Starting writing of tscb file!" << std::endl;
      
      // Remove the file if it exists
      if (std::filesystem::exists(name)) {
        std::filesystem::remove(name);
      }
      
      // Set to ensure buffering doesn't cause issues with I/O operations
      std::ios_base::sync_with_stdio(true);
      
      std::ofstream tsc_file(name.c_str(), std::ios::out | std::ios::binary);
      
      // Check if file was opened successfully
      err_checkf(tsc_file.is_open(), "Failed to open file for writing", std::cout);
      
      int head[1] = {static_cast<int>(header.size())};
      tsc_file.write((char *)&head, sizeof(head));
      tsc_file.flush(); 
      err_checkf(tsc_file.good(), "Problem with tsc file writing header size", std::cout);
      
      tsc_file.write(header.c_str(), head[0] * sizeof(char));
      tsc_file.flush();
      err_checkf(tsc_file.good(), "Problem with tsc file writing header content", std::cout);
      
      // Safely get scatterers string
      std::string sc;
      try {
        sc = scatterers_string();
        //std::cerr << "Got scatterers string, length: " << sc.size() << std::endl;
      } catch (const std::exception& e) {
        std::cerr << "Exception in scatterers_string(): " << e.what() << std::endl;
        tsc_file.close();
        return;
      }
      
      // Check if scatterers string is valid
      if (sc.empty()) {
        std::cerr << "Warning: Empty scatterers string" << std::endl;
      }
      
      head[0] = static_cast<int>(sc.size());
      tsc_file.write((char *)&head, sizeof(head));
      tsc_file.flush();
      err_checkf(tsc_file.good(), "Problem with tsc file writing scatterers size", std::cout);
      
      tsc_file.write(sc.c_str(), head[0] * sizeof(char));
      tsc_file.flush();
      err_checkf(tsc_file.good(), "Problem with tsc file writing scatterers content", std::cout);

      // Check if index vector is valid
      if (index.empty() || index[0].empty()) {
        std::cerr << "Error: Empty or invalid index vector" << std::endl;
        tsc_file.close();
        return;
      }
      
      // Validate index dimensions
      if (index.size() < 3) {
        std::cerr << "Error: Index vector has incorrect dimensions" << std::endl;
        tsc_file.close();
        return;
      }
      
      int nr_hkl[1] = {static_cast<int>(index[0].size())};
      //std::cerr << "Number of HKL indices: " << nr_hkl[0] << std::endl;
      
      tsc_file.write((char *)&nr_hkl, sizeof(nr_hkl));
      tsc_file.flush();
      err_checkf(tsc_file.good(), "Problem with tsc file writing nr_hkl", std::cout);
      
      // Safely get scatterer size
      const int scat_size = scatterer_size();
      //std::cerr << "Scatterer size: " << scat_size << std::endl;
      
      int64_t sc_sf = sf.size();
      int64_t nr_hkl_sf = sf[0].size();
      // Validate sf dimensions
      if (sf.empty() || sc_sf < scat_size || nr_hkl_sf < nr_hkl[0]) {
        std::cerr << "Error: Structure factor array has invalid dimensions" << std::endl;
        tsc_file.close();
        return;
      }
      for(int64_t i=0; i<sc_sf; i++) {
        if (sf[i].size() != nr_hkl_sf) {
          std::cerr << "Error: Inconsistent size in sf for scatterer " << i << std::endl;
          tsc_file.close();
          return;
        }
      }

      for (int run = 0; run < nr_hkl[0]; run++)
      {
        for (int i = 0; i < 3; i++)
        {
          if (i < index.size() && run < index[i].size()) {  // Bounds check
            tsc_file.write((char *)&(index[i][run]), sizeof(numtype_index));
          } else {
            std::cerr << "Error: Index array bounds exceeded at i=" << i << ", run=" << run << std::endl;
            tsc_file.close();
            return;
          }
        }
        
        for (int i = 0; i < scat_size; i++)
        {
          if (i < sf.size() && run < sf[i].size()) {  // Bounds check
            tsc_file.write((char *)&(sf[i][run]), sizeof(std::complex<double>));
          } else {
            std::cerr << "Error: SF array bounds exceeded at i=" << i << ", run=" << run << std::endl;
            tsc_file.close();
            return;
          }
        }
        
        // Add periodic flush to prevent buffer issues
        if (run % 1000 == 0) {
          tsc_file.flush();
        }
      }
      
      // Final flush and close
      tsc_file.flush();
      tsc_file.close();
      //std::cerr << "File closed successfully" << std::endl;
      
      // Verify the file was written correctly
      if (std::filesystem::exists(name)) {
        //std::cerr << "Successfully wrote " << std::filesystem::file_size(name) << " bytes to " << name << std::endl;
      } else {
        std::cerr << "Error: File not created successfully" << std::endl;
      }
      
    } catch (const std::exception& e) {
      std::cerr << "Exception in write_tscb_file: " << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Unknown exception in write_tscb_file" << std::endl;
    }
  }
};

inline bool merge_tscs(
    const std::string &mode,
    const pathvec &files,
    const bool old_tsc)
{
  // Currently only for Mode "pure merge"
  if (mode.empty() || false)
  {
    return false;
  }
  if (files.size() == 0)
    return false;
  svec labels;
  pathvec local_files(files);
  ivec2  indices;                    //[i] is h,k or l, [j] is the number of the reflection
  cvec2 form_fact; //[i] (len(labels)) scatterer, [j](len(indices[0])) reflection correpsonding to indices

  svec header;

  size_t offset = 0;
  indices.resize(3);
  for (int f = 0; f < files.size(); f++)
  {
    std::cout << "Reading file number: " << f + 1 << ": " << files[f] << std::endl;
    if (files[f].extension() == ".tscb")
    {
      std::filesystem::path name = files[f];
      std::filesystem::path new_name = name;
      new_name.replace_extension(".tsc");
      std::cout << "Converting to: " << new_name.string() << std::endl;
      tsc_block<int, cdouble> blocky(name);
      blocky.write_tsc_file(new_name, new_name);
      local_files[f] = new_name;
      std::cout << "Now reading converted: " << new_name << std::endl;
    }
    std::ifstream inf(local_files[f], std::ios::in);
    std::string line;
    bool data = false;
    bvec is_a_new_scatterer;
    int nr_scatterers = 0;
    int nr_new_scatterers;
    while (!data)
    {
      getline(inf, line);
      if (line.find("SCATTERERS:") != std::string::npos)
      {
        std::string temp_labels = line.substr(12, line.size()) + " ";
        const std::string delimiter = " ";
        size_t pos = 0;
        std::string new_label;
        while ((pos = temp_labels.find(delimiter)) != std::string::npos)
        {
          nr_scatterers++;
          new_label = temp_labels.substr(0, pos);
          bool is_new = true;
          for (int i = 0; i < labels.size(); i++)
            if (labels[i] == new_label)
              is_new = false;
          is_a_new_scatterer.push_back(is_new);
          if (is_new)
            labels.push_back(new_label);
          temp_labels.erase(0, pos + delimiter.length());
        }
        nr_new_scatterers = vec_sum(is_a_new_scatterer);
        std::cout << "Read " << nr_scatterers << " atoms, " << nr_new_scatterers << " are new." << std::endl;
        form_fact.resize(labels.size());
      }
      else if (line.find("DATA:") != std::string::npos)
        data = true;
      else if (f == 0)
      {
        if (line.find("AD: ") != std::string::npos)
          continue;
        header.push_back(line);
      }
    }
    std::cout << "Reading Data Block..." << std::endl;
    while (!inf.eof())
    {
      getline(inf, line);
      svec digest = split_string<std::string>(line, " ");
      if (digest.size() == 1 && digest[0] == "")
        continue;
      const int l_indices[3] = {stoi(digest[0]), stoi(digest[1]), stoi(digest[2])};
      if (l_indices[0] == 0 && l_indices[1] == 0 && l_indices[2] == 0)
        continue;
      if (f == 0)
      {
        // For the first file we only read what is not Freidel mates
        bool new_index = true;
#pragma omp parallel for reduction(&& : new_index)
        for (int run = 0; run < indices[0].size(); run++)
        {
          if (l_indices[0] == indices[0][run] && l_indices[1] == indices[1][run] && l_indices[2] == indices[2][run])
            new_index = false;
          else if (stoi(digest[0]) == -indices[0][run] && stoi(digest[1]) == -indices[1][run] && stoi(digest[2]) == -indices[2][run])
            new_index = false;
        }
        if (!new_index)
          continue;
        for (int i = 0; i < 3; i++)
          indices[i].push_back(l_indices[i]);
#pragma omp parallel for
        for (int i = 0; i < nr_scatterers; i++)
        {
          if (!is_a_new_scatterer[i])
            continue;
          size_t nr_in_new_scats = 0;
          for (int p = 0; p < i; p++)
            if (is_a_new_scatterer[p])
              nr_in_new_scats++;
          svec re_im = split_string<std::string>(digest[i + 3], ",");
          form_fact[offset + nr_in_new_scats].push_back(std::complex<double>(stod(re_im[0]), stod(re_im[1])));
        }
      }
      else
      {
        // Otherwise we put stuff where it belongs
#pragma omp parallel for
        for (int i = 0; i < nr_scatterers; i++)
        {
          if (!is_a_new_scatterer[i])
            continue;
          size_t nr_in_new_scats = 0;
          for (int p = 0; p < i; p++)
            if (is_a_new_scatterer[p])
              nr_in_new_scats++;
          form_fact[offset + nr_in_new_scats].resize(form_fact[0].size());
        }
        bool found = false;
        for (int run = 0; run < indices[0].size(); run++)
        {
          if (l_indices[0] == indices[0][run] && l_indices[1] == indices[1][run] && l_indices[2] == indices[2][run])
          {
#pragma omp parallel for
            for (int i = 0; i < nr_scatterers; i++)
            {
              if (!is_a_new_scatterer[i])
                continue;
              size_t nr_in_new_scats = 0;
              for (int p = 0; p < i; p++)
                if (is_a_new_scatterer[p])
                  nr_in_new_scats++;
              svec re_im = split_string<std::string>(digest[i + 3], ",");
              form_fact[offset + nr_in_new_scats][run] = std::complex<double>(stod(re_im[0]), stod(re_im[1]));
            }
            found = true;
          }
          else if (l_indices[0] == -indices[0][run] && l_indices[1] == -indices[1][run] && l_indices[2] == -indices[2][run])
          {
#pragma omp parallel for
            for (int i = 0; i < nr_scatterers; i++)
            {
              if (!is_a_new_scatterer[i])
                continue;
              size_t nr_in_new_scats = 0;
              for (int p = 0; p < i; p++)
                if (is_a_new_scatterer[p])
                  nr_in_new_scats++;
              svec re_im = split_string<std::string>(digest[i + 3], ",");
              form_fact[offset + nr_in_new_scats][run] = std::complex<double>(stod(re_im[0]), -stod(re_im[1]));
            }
            found = true;
          }
          if (found)
            break;
        }
      }
    }
    std::cout << "Data for " << form_fact[0].size() << " indices read." << std::endl;
    offset = labels.size();
  }
  std::cout << "Writing combined file..." << std::endl;
  std::string header_string("");
  for (size_t h_loc = 0; h_loc < header.size(); h_loc++)
    header_string += header[h_loc] + "\n";
  tsc_block<int, cdouble> combined(form_fact, labels, indices, header_string);
  if (!old_tsc)
  {
    combined.write_tscb_file("combined.tscb");
  }
  else
  {
    combined.write_tsc_file("combined", "combined.tsc");
  }
  std::cout << "Done!" << std::endl;

  return true;
}

inline bool merge_tscs_without_checks(
    const std::string &mode,
    const pathvec &files,
    const bool old_tsc)
{
  // Currently only for Mode "pure merge"
  if (mode.empty() || false)
  {
    return false;
  }
  if (files.size() == 0)
    return false;
  svec labels;
  pathvec local_files(files);
  ivec2 indices;                    //[i] is h,k or l, [j] is the number of the reflection
  cvec2 form_fact; //[i] (len(labels)) scatterer, [j](len(indices[0])) reflection correpsonding to indices
  svec header;

  size_t offset = 0;
  indices.resize(3);
  for (int f = 0; f < files.size(); f++)
  {
    std::cout << "Reading file number: " << f + 1 << ": " << files[f] << std::endl;
    if (files[f].extension() == ".tscb")
    {
      std::filesystem::path name = files[f];
      std::filesystem::path new_name = name.replace_extension(".tsc");
      std::cout << "Converting to: " << new_name << std::endl;
      tsc_block<int, cdouble> blocky(name);
      blocky.write_tsc_file(new_name, new_name);
      local_files[f] = new_name;
      std::cout << "Now reading converted: " << new_name << std::endl;
    }
    std::ifstream inf(local_files[f].c_str(), std::ios::in);
    std::string line;
    bool data = false;
    bvec is_a_new_scatterer;
    int nr_scatterers = 0;
    int nr_new_scatterers;
    while (!data)
    {
      getline(inf, line);
      if (line.find("SCATTERERS:") != std::string::npos)
      {
        std::string temp_labels = line.substr(12, line.size()) + " ";
        const std::string delimiter = " ";
        size_t pos = 0;
        std::string new_label;
        while ((pos = temp_labels.find(delimiter)) != std::string::npos)
        {
          nr_scatterers++;
          new_label = temp_labels.substr(0, pos);
          bool is_new = true;
          for (int i = 0; i < labels.size(); i++)
            if (labels[i] == new_label)
              is_new = false;
          is_a_new_scatterer.push_back(is_new);
          if (is_new)
            labels.push_back(new_label);
          temp_labels.erase(0, pos + delimiter.length());
        }
        nr_new_scatterers = vec_sum(is_a_new_scatterer);
        std::cout << "Read " << nr_scatterers << " atoms, " << nr_new_scatterers << " are new." << std::endl;
        form_fact.resize(labels.size());
      }
      else if (line.find("DATA:") != std::string::npos)
        data = true;
      else if (f == 0)
      {
        if (line.find("AD: ") != std::string::npos)
          continue;
        header.push_back(line);
      }
    }
    std::cout << "Reading Data Block..." << std::endl;
    while (!inf.eof())
    {
      getline(inf, line);
      svec digest = split_string<std::string>(line, " ");
      if (digest.size() == 1 && digest[0] == "")
        continue;
      if (f == 0)
        for (int i = 0; i < 3; i++)
          indices[i].push_back(stoi(digest[i]));
#pragma omp parallel for
      for (int i = 0; i < nr_scatterers; i++)
      {
        if (!is_a_new_scatterer[i])
          continue;
        size_t nr_in_new_scats = 0;
        for (int p = 0; p < i; p++)
          if (is_a_new_scatterer[p])
            nr_in_new_scats++;
        svec re_im = split_string<std::string>(digest[i + 3], ",");
        form_fact[offset + nr_in_new_scats].push_back(std::complex<double>(stod(re_im[0]), stod(re_im[1])));
      }
    }
    std::cout << "Data for " << form_fact.back().size() << " indices read." << std::endl;
    offset = labels.size();
    inf.close();
  }
  std::cout << "Writing combined file..." << std::endl;
  std::string header_string("");
  for (size_t h_loc = 0; h_loc < header.size(); h_loc++)
    header_string += header[h_loc] + "\n";
  tsc_block<int, cdouble> combined(form_fact, labels, indices, header_string);
  if (!old_tsc)
  {
    combined.write_tscb_file("combined.tscb");
  }
  else
  {
    combined.write_tsc_file("combined", "combined.tsc");
  }
  std::cout << "Done!" << std::endl;

  return true;
}
