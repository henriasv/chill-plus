

#include "chiller.h"
#include <iostream>

Chiller::Chiller()
{

}

Chiller::Chiller(double* positions, int numParticles)
{
  m_simulationCell = {vec3(288.72, 0, 0), vec3(0, 288.72, 0), vec3(0, 0, 288.72), vec3(0, 0, 0)};
  m_simulationCellSize = vec3(m_simulationCell[0].x()-m_simulationCell[3].x(),
                              m_simulationCell[1].y()-m_simulationCell[3].y(),
                              m_simulationCell[2].z()-m_simulationCell[3].z());
  m_numParticles = numParticles;
  for (int i = 0; i<numParticles; i++)
  {
    m_positions.push_back(std::move(vec3(positions[3*i], positions[3*i+1], positions[3*i+2])));
    m_status.push_back(0);
    m_neighbors.push_back(std::vector<int>{});
  }
  constructCellLists(3.5);
  cutoffNeighbors(3.5);
}

std::vector<int> Chiller::getStatus()
{
  std::vector<int> ret;
  for (int i = 0; i<m_numParticles; i++)
  {
    ret.push_back(i);
  }
  return ret;
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
    //binNumbers.floor();
    if (binNumbers[0] >= m_nBins[0] || binNumbers[1]>=m_nBins[1] || binNumbers[2] >= m_nBins[2])
    {std::cout << "too large bin number " << binNumbers << std::endl << position << std::endl;}
    //std::cout << binNumbers << std::endl;
    m_cells[(int)binNumbers[0]][(int)binNumbers[1]][(int)binNumbers[2]]->push_back(i);
  }
  std::cout << "hei" << std::endl;
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
                    vec3 diff = periodicDistance(m_positions[p], m_positions[l], thisCell, otherCell);
                    if (diff.sqnorm()<sqCutoff && l!=p)
                    {
                      m_neighbors[l].push_back(p);
                    }
                  }
                }
            }
          }
          if (m_neighbors[l].size() > 4)
          {
            std::cout << m_positions[l] << std::endl;
          }
          //std::cout << m_neighbors[l].size() << " ";
        }
      }
    }
  }
}

void Chiller::chillPlus()
{
  for (int i = 0; i<m_numParticles; i++)
  {
    for (int j = 0; j<m_neighbors[i].size(); j++)
    {

    }
  }
}

vec3 Chiller::periodicDistance(vec3 & vec1, vec3 & vec2, vec3 & cell1, vec3 & cell2)
{
    vec3 diff = vec2-vec1;
    vec3 cellDiff = cell2-cell1;
    diff[0] = diff[0] - m_simulationCellSize[0]*(cellDiff[0]==(m_nBins[0]-1)) + m_simulationCellSize[0]*(cellDiff[0]==(-m_nBins[0]+1));
    diff[1] = diff[1] - m_simulationCellSize[1]*(cellDiff[1]==(m_nBins[1]-1)) + m_simulationCellSize[1]*(cellDiff[1]==(-m_nBins[1]+1));
    diff[2] = diff[2] - m_simulationCellSize[2]*(cellDiff[2]==(m_nBins[2]-1)) + m_simulationCellSize[2]*(cellDiff[2]==(-m_nBins[2]+1));
    return diff;
}

vec3 Chiller::periodicCell(vec3 cell)
{
    vec3 ret = cell;
    for (int index = 0; index < 3; index++)
    {
        if (ret[index] == -1)
            ret[index] = m_nBins[index];
        if (ret[index] == m_nBins[index])
            ret[index] = 0;
    }
    return ret;
}
