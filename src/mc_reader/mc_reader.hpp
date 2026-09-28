// Blacklight simulation reader header

#ifndef MC_READER_H_
#define MC_READER_H_

// C++ headers
#include <fstream>  // ifstream
#include <iosfwd>   // streampos
#include <string>   // string

// Blacklight headers
#include "../blacklight.hpp"                 // enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../simulation_reader/simulation_reader.hpp"
#include "../utils/array.hpp"                // Array

//--------------------------------------------------------------------------------------------------

// OpacityTable reader
struct MCReader
{
  // Constructors and destructor
  MCReader(const InputReader *p_input_reader_, const SimulationReader *p_simulation_reader_);
  MCReader(const MCReader &source) = delete;
  MCReader &operator=(const MCReader &source) = delete;
  ~MCReader();

  // Pointers to other objects
  const InputReader *p_input_reader;
  const SimulationReader *p_simulation_reader;

  // Input data - general
  ModelType model_type;

  // Input data - opacity parameters
  std::string mc_file_name;
  std::string mc_freq_file_name;
  bool mc_input;
  Array<double> x1f, x2f, x3f;
  Array<double> x1v, x2v, x3v;
  Array<int> levels;
  
  FILE  *mc_file;
  FILE *mc_freq_file;

  // Data arrays
  int num_freqs;
  Array<float> *scattering_source_terms;
  Array<float> *grid_prim;
  Array<float> *scattering_error;

  Coordinates simulation_coord;
  SimulationFormat simulation_format;

  bool first_time = true;
  bool first_time_root_object_header = true;
  
  double simulation_m_msun;
  bool simulation_all_cgs;
  double simulation_r_rg;
  double simulation_rho_cgs;
  double simulation_v_c;
  bool simulation_hd_only;
  bool simulation_mc_temp;
  bool compton;
  bool stimulated_compton;
  bool mc_error;

  Array<double> freq_grid;
  //Array<double>* freq_grid_ptr = &freq_grid;
  Array<double> ln_freq_grid;

  double dlf;

  int ind_rho, ind_pgas, ind_kappa;
  int ind_uu1, ind_uu2, ind_uu3;
  int ind_bb1, ind_bb2, ind_bb3;
  PlasmaModel plasma_model;
  bool plasma_use_p;
  double plasma_thermal_frac;
  double plasma_kappa_frac;
  double plasma_power_frac;
  double plasma_mu;
  double plasma_ne_ni;
  double plasma_gamma;
  double plasma_gamma_i;
  double plasma_gamma_e;
  double plasma_rat_low;
  double plasma_rat_high;
  double bh_m, bh_a;

  std::ifstream data_stream;
  unsigned long int root_object_header_address;
  unsigned long int root_btree_address;
  unsigned long int root_name_heap_address;
  unsigned long int root_data_segment_address;

  int n_3_root;

  std::string *variable_names;
  int num_variable_names = 0;
  std::string *dataset_names;
  int num_dataset_names = 0;
  Array<int> num_variables;
   
  
  // Internal functions - mc_reader.cpp
  double Read(int snapshot);
  void ReadFreqFile();
  void Gradient(Array<float> &grad,Array<float> &f, Array<double> &x); 
  void CalculateSourceTerm(Array<float> &source_term,Array<float> &scattering,Array<float> &scattering_prime,Array<float> &scattering_prime_prime,Array<float> &scattering_error);

  
  // Internal functions - mc_hdf5_format_structure.cpp
  void ReadHDF5Superblock();
  void ReadHDF5RootGroupSymbolTableEntry();
  unsigned long int ReadHDF5Heap(unsigned long int heap_address);
  void ReadHDF5RootObjectHeader();
  void ReadHDF5FloatAttribute(const char *attribute_name, float *p_val);

  // Internal functions - mc_hdf5_format_metadata.cpp
  unsigned long int ReadHDF5DatasetHeaderAddress(const char *name, unsigned long int btree_address,
      unsigned long int data_segment_address);
  void ReadHDF5DataObjectHeader(unsigned long int data_object_header_address,
      unsigned char **p_datatype_raw, unsigned char **p_dataspace_raw, unsigned char **p_data_raw);
  static void ReadHDF5DataspaceDims(const unsigned char *dataspace_raw, unsigned long int **p_dims,
      int *p_num_dims);

  // Internal functions - mc_hdf5_format_arrays.cpp
  void ReadHDF5StringArray(const char *name, bool allocate, std::string **p_string_array,
      int *p_array_length);
  void ReadHDF5IntArray(const char *name, Array<int> &int_array);
  void ReadHDF5FloatArray(const char *name, Array<float> &float_array);
  void ReadHDF5FloatArray(const char *name, Array<double> &double_array);
  void ReadHDF5DoubleArray(const char *name, Array<double> &double_array);
  void ReadHDF5DoubleArray(const char *name, Array<float> &float_array);
  static void SetHDF5StringArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, bool allocate,
      std::string **string_array, int *p_array_length);
  static void SetHDF5IntArray(const unsigned char *datatype_raw, const unsigned char *dataspace_raw,
      const unsigned char *data_raw, Array<int> &int_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<float> &float_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      Array<double> &double_array);
  static void SetHDF5DoubleArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      Array<double> &double_array);
  static void SetHDF5DoubleArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      Array<float> &float_array);

  //Internal functions - simulation_geometry.cpp
  void ContravariantSimulationMetric(double x, double y, double z,double gcon[4][4]) const;
  void CovariantSimulationMetric(double x, double y, double z, double gcov[4][4]) const;
};

#endif
