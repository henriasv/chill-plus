#define NPY_NO_DEPRECATED_API

#include <vector>
#include "vec3.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/numeric/ublas/matrix.hpp>

class Chiller
{
public:
  Chiller();
  Chiller(double*, int, double*);
  ~Chiller();
  std::vector<int> getStatus();
  void constructCellLists(double);
  void cutoffNeighbors(double);
  void chillPlus();
  void movePositionsInsideBox();
  vec3 periodicDistance(vec3 & vec1, vec3 & vec2);
  vec3 periodicCell(vec3 cell);
  std::pair<float, float> polar_asimuthal(vec3 delta);
  std::complex<float> q_lm(int i, int m);
  std::complex<float> q_lm_from_matrix(int i, int m);
  std::vector<vec3> m_positions;
  std::vector<std::vector<int> > m_neighbors = {};
  std::vector<int>**** m_cells;
  std::vector<int> m_status;
  std::vector<vec3> m_simulationCell;
  boost::numeric::ublas::matrix<std::complex<float>> q_values;
  vec3 m_simulationCellSize;
  vec3 m_nBins;
  vec3 m_binSize;
  int m_numParticles = 0;
};
