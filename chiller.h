#define NPY_NO_DEPRECATED_API

#include <vector>
#include "vec3.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>

class Chiller
{
public:
  Chiller();
  Chiller(double*, int);
  std::vector<int> getStatus();
  void constructCellLists(double);
  void cutoffNeighbors(double);
  void chillPlus();
  vec3 periodicDistance(vec3 & vec1, vec3 & vec2, vec3 & cell1, vec3 & cell2);
  vec3 periodicCell(vec3 cell);
  std::vector<vec3> m_positions;
  std::vector<std::vector<int> > m_neighbors = {};
  std::vector<int>**** m_cells;
  std::vector<int> m_status;
  std::vector<vec3> m_simulationCell;
  vec3 m_simulationCellSize;
  vec3 m_nBins;
  vec3 m_binSize;
  int m_numParticles = 0;
};
