

#include "chiller.h"
#include <iostream>
#include <cmath>

Chiller::Chiller()
{

}

Chiller::~Chiller()
{
  for (int i = 0; i<m_nBins[0]; i++)
  {
    for (int j = 0; j<m_nBins[1]; j++)
    {
      for (int k = 0; k<m_nBins[2]; k++)
      {
        //m_cells[i][j][k]->clear();
        //std::cout << i << " " << j << " " << k << std::endl;
        if (m_cells[i][j][k] != NULL)
        {
          delete m_cells[i][j][k];
          m_cells[i][j][k] = NULL;
        }
        //std::cout << i << " " << j << " " << k << std::endl;
      }
      //delete m_cells[i][j];
    }
    //delete m_cells[i];
  }
  //delete m_cells;
}

Chiller::Chiller(double* positions, int numParticles, double* simulationCell)
{
  vec3 cellvec1(simulationCell[0], simulationCell[1], simulationCell[2]);
  vec3 cellvec2(simulationCell[3], simulationCell[4], simulationCell[5]);
  vec3 cellvec3(simulationCell[6], simulationCell[7], simulationCell[8]);
  vec3 cellorigin(simulationCell[9], simulationCell[10], simulationCell[11]);
  std::cout << cellvec1 << std::endl;
  std::cout << cellvec2 << std::endl;
  std::cout << cellvec3 << std::endl;
  std::cout << cellorigin << std::endl;
  m_simulationCell = {cellvec1, cellvec2, cellvec3, cellorigin};
  m_simulationCellSize = vec3(m_simulationCell[0].x(),
                              m_simulationCell[1].y(),
                              m_simulationCell[2].z());
  m_numParticles = numParticles;
  q_values = boost::numeric::ublas::matrix<std::complex<float>>(numParticles, 7);
  for (int i = 0; i<numParticles; i++)
  {
    m_positions.push_back(std::move(vec3(positions[3*i], positions[3*i+1], positions[3*i+2])));
    m_status.push_back(0);
    m_neighbors.push_back(std::vector<int>{});
  }
  movePositionsInsideBox();
  constructCellLists(4.0);
  cutoffNeighbors(4.0);
  chillPlus();
}

void Chiller::movePositionsInsideBox()
{
  for (int i = 0; i<m_numParticles; i++)
  {
    for (int j = 0; j<3; j++)
    {
      if (m_positions[i][j]>m_simulationCell[0][j]+m_simulationCell[3][j])
      {
        m_positions[i][j] = m_positions[i][j] - m_simulationCellSize[j];
      }
      if (m_positions[i][j] < m_simulationCell[3][j])
      {
        m_positions[i][j] = m_positions[i][j] + m_simulationCellSize[j];
      }
    }
  }
}

std::vector<int> Chiller::getStatus()
{
  return m_status;
}

void Chiller::constructCellLists(double cutoffLength)
{
  // First, construct the cell lists
  m_nBins = vec3( m_simulationCellSize[0]/cutoffLength,
                  m_simulationCellSize[1]/cutoffLength,
                  m_simulationCellSize[2]/cutoffLength);
  m_nBins.floor();
  m_binSize = m_simulationCellSize/m_nBins;
  std::cout << m_nBins << std::endl;

  m_cells = new std::vector<int>***[(int)m_nBins[0]];
  for (int i = 0; i<m_nBins[0]; i++)
  {
    m_cells[i] = new std::vector<int>**[(int)m_nBins[1]];
    for (int j = 0; j<m_nBins[1]; j++)
    {
      m_cells[i][j] = new std::vector<int>*[(int)m_nBins[2]];
      for (int k = 0; k<m_nBins[2]; k++)
      {
        m_cells[i][j][k] = new std::vector<int>;
      }
    }
  }



  for (int i = 0; i<m_numParticles; i++)
  {
    vec3 position = m_positions[i];
    vec3 binNumbers = (position-m_simulationCell[3])/m_binSize;
    m_cells[(int)binNumbers[0]][(int)binNumbers[1]][(int)binNumbers[2]]->push_back(i);
  }
}

void Chiller::cutoffNeighbors(double cutoffLength)
{
  double sqCutoff = cutoffLength*cutoffLength;
  for (int i = 0; i<m_nBins[0]; i++)
  {
    for (int j = 0; j<m_nBins[1]; j++)
    {
      for (int k = 0; k<m_nBins[2]; k++)
      {
        //std::cout << m_cells[i][j][k]->size() << " ";
        for (auto l : *m_cells[i][j][k])
        {
          vec3 thisCell = vec3((double)i, (double)j, (double)k);
          for (int m=-1; m<=1; m++)
          {
            for (int n=-1; n<=1; n++)
            {
              for (int o=-1; o<=1; o++)
                {
                  vec3 otherCell = periodicCell(thisCell+vec3((double)m, (double)n, (double)o));
                  int indX = (int) otherCell[0];
                  int indY = (int) otherCell[1];
                  int indZ = (int) otherCell[2];
                  for (auto p : *m_cells[indX][indY][indZ])
                  {
                    vec3 diff = periodicDistance(m_positions[p], m_positions[l]);
                    if (diff.sqnorm()<sqCutoff && l!=p)
                    {
                      m_neighbors[l].push_back(p);
                    }
                  }
                }
            }
          }
          while (m_neighbors[l].size() > 4)
          {
            double max_dist = 0;
            int max_index = 0;
            for (int p = 0; p<m_neighbors[l].size(); p++)
            {
              double tmp_dist = periodicDistance(m_positions[m_neighbors[l][p]], m_positions[l]).norm();
              if (tmp_dist > max_dist)
              {
                max_dist = tmp_dist;
                max_index = p;
              }
            }
            m_neighbors[l].erase(m_neighbors[l].begin()+max_index);
          }
          //std::cout << m_neighbors[l].size() << " ";
        }
      }
    }
  }
}

std::complex<float> Chiller::q_lm(int i, int m)
{
  std::complex<float> q = 0;
  for (int j = 0; j<4; j++)
  {
    vec3 posi = m_positions[i];
    vec3 posj = m_positions[m_neighbors[i][j]];
    vec3 delta = periodicDistance(posi, posj);

    std::pair<float, float> angles = polar_asimuthal(delta);
    q += boost::math::spherical_harmonic(3, m, angles.first, angles.second);
  }
  return q;
}

std::complex<float> Chiller::q_lm_from_matrix(int i, int m)
{
  //std::cout << "Getting q(" << i << "," << m << ") from matrix"<< std::endl;
  return q_values(i, m+3);
}

void Chiller::chillPlus()
{
  for (int i = 0; i<m_numParticles; i++)
  {
    if (m_neighbors[i].size() == 4)
    {
      for (int m = -3; m<=3; m++)
      {
        //std::cout << "Setting q(" << i << "," << m << ") to  matrix"<< std::endl;
        q_values(i,m+3) = q_lm(i,m);
      }
    }
  }
  for (int i = 0; i<m_numParticles; i++)
  {
    if (m_neighbors[i].size() == 4)
    {
      int num_eclipsed = 0;
      bool eclipsed[] = {false, false, false, false};
      for (int j = 0; j<m_neighbors[i].size(); j++)
      {
        std::complex<float> c1 = 0;
        std::complex<float> c2 = 0;
        std::complex<float> c3 = 0;
        std::complex<float> q_i = 0;
        std::complex<float> q_j = 0;
        if (m_neighbors[m_neighbors[i][j]].size() == 4)
        {
          for (int m = -3; m<=3; m++)
          {
            q_i = q_lm_from_matrix(i, m);
            q_j = q_lm_from_matrix(m_neighbors[i][j], m);
            c1 += q_i*std::conj(q_j);
            c2 += q_i*std::conj(q_i);
            c3 += q_j*std::conj(q_j);
          }
          std::complex<float> c_ij = c1/(std::sqrt(c2)*std::sqrt(c3));
          if (std::real(c_ij) > -0.35 && std::real(c_ij) < 0.25)
          {
            eclipsed[j] = true;
            num_eclipsed ++;
          }
        }
      }

      if (m_status[i]<num_eclipsed)
      {
        m_status[i] = num_eclipsed;
      }

/*
      for (int j = 0; j<m_neighbors[i].size(); j++)
      {
        if (m_status[m_neighbors[i][j]] < num_eclipsed-1)
        {
          m_status[m_neighbors[i][j]] = num_eclipsed-1;
        }
      }
*/
      //std::cout << num_eclipsed << std::endl;
    }
  }
}


std::pair<float, float> Chiller::polar_asimuthal(vec3 delta)
{
  float asimuthal = std::atan2(delta.y(), delta.x());
  float xy_distance = std::sqrt(delta.x()*delta.x()+delta.y()*delta.y());
  float polar = std::atan2(xy_distance, delta.z());
  return std::pair<float, float>(polar, asimuthal);
}

vec3 Chiller::periodicDistance(vec3 & vec1, vec3 & vec2)
{
    vec3 diff = vec2-vec1;
    diff[0] = diff[0] - m_simulationCellSize[0]*(diff[0]>(m_simulationCellSize[0]/2)) + m_simulationCellSize[0]*(diff[0] < -(m_simulationCellSize[0]/2));
    diff[1] = diff[1] - m_simulationCellSize[1]*(diff[1]>(m_simulationCellSize[1]/2)) + m_simulationCellSize[1]*(diff[1] < -(m_simulationCellSize[1]/2));
    diff[2] = diff[2] - m_simulationCellSize[2]*(diff[2]>(m_simulationCellSize[2]/2)) + m_simulationCellSize[2]*(diff[2] < -(m_simulationCellSize[2]/2));
    return diff;
}

vec3 Chiller::periodicCell(vec3 cell)
{
    vec3 ret = cell;
    for (int i = 0; i < 3; i++)
    {
        if (ret[i] == -1)
            ret[i] = m_nBins[i]-1;
        if (ret[i] == m_nBins[i])
            ret[i] = 0;
    }
    return ret;
}
