#ifndef __VESTECSEEDINGALGORITHM_H
#define __VESTECSEEDINGALGORITHM_H

#include <vtkAlgorithm.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkPolyData.h>
#include <vtkCell.h>

#include <map>
#include <string>
#include <vector>
#include <random>

#include <Eigen/Dense>

//Matrix to compute the determinant
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 1, 4, 4> DynamicMatrix;

//Forward declaration of vtkDataSet
class vtkDataSet;

/// This class extends the algorithm described in:
/// "Detection and classification of critical points in piecewise linear vector fields" by Wang et al.
/// (that extends the following method) "Robust Detection of Singularities in Vector Fields" by Bhatia et al./// 
/// -----
/// other references here:
/// (possible extension using a different index for classification) "Morse set classification and hierarchical refinement using conley index" by Chen G et al.
/// (survey) "A survey of topology‐based methods in visualization" by Heine et al.
/// (survey) "From numerics to combinatorics: a survey of topological methods for vector field visualization" by Wang et al.
class CriticalPointExtractor {
  public:
    /**
     * Store vector and points in internal data structure 
     */
    CriticalPointExtractor(vtkDataSet* input, double* currentSingularity, int mpiRank, bool pertubate = true);

    /**
     * Identify the critical cells 
     */
    void ComputeCriticalCells();

    /*
     * Save output into a vtk file
     */ 
    void writeCriticalCells(vtkSmartPointer<vtkDataSet> output);

    /**
     * Point types used for classification (2D and 3D)
     */
    enum PointType {      
      UNCLASSIFIED_SINGULARITY=-1,
      ATTRACTING_NODE=0, 
      ATTRACTING_FOCUS=1,     
      ATTRACTING_NODE_SADDLE=2, 
      ATTRACTING_FOCUS_SADDLE=3, 
      REPELLING_NODE_SADDLE=4,
      REPELLING_FOCUS_SADDLE=5,       
      REPELLING_NODE=6, 
      REPELLING_FOCUS=7,
      CENTER_NODE=8,
      NODE_SADDLE_2D=9,
      REGULAR_POINT=10
    };

    /**
     * Meta data about a critical point
     * id = simplex id in the global array
     */
    struct CriticalPoint {
      vtkIdType id;
      PointType type;

      CriticalPoint(vtkIdType i, PointType t) { 
        id = i; 
        type = t;
      }
    };

    /**
     *  Meta data about the vtk input dataset
     */
    struct DataSetMetadata {
      double local_bounds[6] = {-1,-1,-1,-1,-1,-1};
      vtkIdType max_global_id=LLONG_MIN, min_global_id=LLONG_MAX, max_local_id=LLONG_MIN, min_local_id=LLONG_MAX;
    };

    ~CriticalPointExtractor() {     
      delete []position;
	    delete []vector;
	    delete []vecCellIds;
    }
    
  private:
    /**
     * Sort the vector integers but returns the needed swap operations
     */
    int  Sort(vtkIdType* ids, int n);

    /**
     * Sort the vector with 3 integers but returns the needed swap operations
     */
    int	 Sort3(vtkIdType* ids);

    /**
     * Sort the vector with 4 integers but returns the needed swap operations
     */
	  int  Sort4(vtkIdType* ids);
    
    /**
     * Compute the determinant
     * tmpIds: Contains the vertex indices needed to fill the matrix. Is a copy since they will be sorted
     * grid: Gives access to vector field
     * currentSingularity: The singularity e.g. zero vector (0,0,0)
     * vecMatrix: The matrix used to compute the determinant
     */
    double ComputeDeterminant(std::array<vtkIdType, 4> &tmpIds, DynamicMatrix &vecMatrix, bool usePoints, int perturbationID = -1);
   
    /**
     * Check if the singularity is in cell
     * ids: The vertex ids spanning the cell
     * vecMatrix: The matrix used to compute the determinant
     */ 
    PointType PointInCell(const vtkIdType* ids, DynamicMatrix &vecMatrix);

    /**
     * Check if direction of the determinant is positive (counter-clockwise) 
     */
    bool DeterminantCounterClockWise(double& det);

    /**
     * Perturbation based on a unique hash which is normalized
     */
    void Perturbate(double* values, vtkIdType &hash, DataSetMetadata &dm);

    /**
     * Compute a hash for a given position. Two same positions in a distributed setup should produce the same hash
     */
    vtkIdType ComputeHash(double* pos);
    
    /**
     * Advanced Critical Simplex Classification
     * Implementing method from: 
     *    A. Globus, C. Levit and T. Lasinski, 
     *    "A tool for visualizing the topology of three-dimensional vector fields," 
     *    Proceeding Visualization '91, 1991, pp. 33-40, doi: 10.1109/VISUAL.1991.175773.
     */
    CriticalPointExtractor::PointType ClassifyCriticalSimplex(const vtkIdType* ids);

    /**
     * Initialize the matrices for critical simplex classification (either 3x3 or 2x2)
     */
    template<class T> void InitializeMatrices(const vtkIdType* ids, T& coordsMatrix, T& vectorsMatrix);

    /**
     * Compute statistic of the real and imaginary parts of the eigenvalues
     */
    template<class T> void CheckEigenvalues(T& ev, int &posReal, int &negReal, int &zeroImag, int &complexImag);

    /**
     * Classify the critical simplex
     */
    template<class T> CriticalPointExtractor::PointType GetCriticalSimplexType(T& ev, int &posReal, int &negReal, int &zeroImag, int &complexImag);

    /**
     * Computes the barycentric coordinates of a simplex
     * The method implemented is nicely explained at the following link: https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_triangles
     */
    double* ComputeBarycentricCoordinates(const vtkIdType* ids);
    
private:
    vtkIdType ZERO_ID;  //!< Vertex ID of the singularity: always number of vertices + 1 
    int iExchangeIndex; //!< The row id of the matrix to exchange with the singularity    
    double* position;   //!< Store point coordinates
	  double* vector;     //!< Store vector field

    vtkIdType* vecCellIds;              //!< The point ids for each cell
    int numCellIds;                     //!< The number of vertices per simplex
    double singularity[3];              //!< The singularity to identify normally (0,0,0)
    int numThreads;                     //!< Number of OpenMP threads
    double eps = 1 / std::pow(10,14);   //!< A very small constant used for perturbation
	  double delta = 4;                   //!< A constant used for perturbatino (must be greater than >= numCellIds)   
    vtkIdType numSimplicesPerCell;      //!< The algorithms just works with simplices. In case of a regular grid we have to create e.g. 5 tetra per voxel
	  std::vector<CriticalPoint> vecCriticalCellIDs;  //!< vector for storing the critical simplices
    vtkIdType numSimplices;                         //!< Total number of simplices
    
    std::map<vtkIdType,int> global_id_uniqueness_map; // for debug only  
public: // for debug only
    int local_deg_cases  = 0;
    int global_deg_cases = 0;
};


template <class T>
inline void hash_combine(vtkIdType& seed, const T& v)
{
    std::hash<T> constexpr hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <> struct std::hash<double*> {
    vtkIdType operator()(double* const& p) const {
        vtkIdType v = 0x778abe;
        hash_combine(v, p[0]);
        hash_combine(v, p[1]);
        hash_combine(v, p[2]);
		return v;
    }
};

template<class T> void CriticalPointExtractor::InitializeMatrices(const vtkIdType* ids, T& coordsMatrix, T& vectorsMatrix) {
  int mSize = coordsMatrix.rows();
  for (int i = 0; i < mSize; i++)
	{
		for (int j = 0; j < mSize; j++)
		{
			coordsMatrix(i,j) = position[ids[j]*3+i] - position[ids[mSize]*3+i];
			vectorsMatrix(i,j) = vector[ids[j]*3+i] - vector[ids[mSize]*3+i];
		}				
	}
}

template<class T> void CriticalPointExtractor::CheckEigenvalues(T& ev, int &posReal, int &negReal, int &zeroImag, int &complexImag) {
  auto real = ev.real();
	auto imag = ev.imag();		

	for (int i = 0; i < real.size(); i++)
	{
		if (imag[i] == 0.0)
			zeroImag++;
		else
			complexImag++;

		if (real[i] < 0)
			negReal++;
		else if (real[i] > 0)
			posReal++;
	}
}

template<class T> CriticalPointExtractor::PointType CriticalPointExtractor::GetCriticalSimplexType(T& ev, int &posReal, int &negReal, int &zeroImag, int &complexImag) {
  int numIds = numCellIds;

  if(numIds == 4) { //3D case - tetrahedra
    if (posReal + negReal == 3)
		{
			switch (posReal)
			{
			case 0:
				if (complexImag == 0)
					return ATTRACTING_NODE;
				else
					return ATTRACTING_FOCUS;
			case 1:
				if (complexImag == 0)
					return ATTRACTING_NODE_SADDLE;
				else
					return ATTRACTING_FOCUS_SADDLE;
			case 2:
				if (complexImag == 0)
					return REPELLING_NODE_SADDLE;
				else
					return REPELLING_FOCUS_SADDLE;
			case 3:
				if (complexImag == 0)
					return REPELLING_NODE;
				else
					return REPELLING_FOCUS;
			default:
				{
					std::cout << "[ERROR] Cannot classify the critical simplex with the following eigenvalues: "<<std::endl
			  	  <<" real: "<<ev.real()[0]<<" "<<ev.real()[1]<<" "<<ev.real()[2]<<std::endl
			  	  <<" imag: "<<ev.imag()[0]<<" "<<ev.imag()[1]<<" "<<ev.imag()[2]<<std::endl;
					return UNCLASSIFIED_SINGULARITY;
				}
			}
		} else if (complexImag > 0) {
			return CENTER_NODE;
		} else {
			std::cout << "[ERROR] Cannot classify the critical simplex with the following eigenvalues: "<<std::endl
				  <<" real: "<<ev.real()[0]<<" "<<ev.real()[1]<<" "<<ev.real()[2]<<std::endl
			 	  <<" imag: "<<ev.imag()[0]<<" "<<ev.imag()[1]<<" "<<ev.imag()[2]<<std::endl;	
      return UNCLASSIFIED_SINGULARITY;
    }     
  } else { //2D case - triangles
    if (posReal + negReal == 2)
		{
			switch (posReal)
			{
			case 0:
				if (complexImag == 0)
          return ATTRACTING_NODE;
				else 
          return ATTRACTING_FOCUS;
			case 1:
				return NODE_SADDLE_2D;
			case 2:
				if (complexImag == 0)
          return REPELLING_NODE;
				else
          return REPELLING_FOCUS;
			default:
				{
					std::cout << "[ERROR] Cannot classify the critical simplex with the following eigenvalues: "<<std::endl
				    <<" real: "<<ev.real()[0]<<" "<<ev.real()[1]<<std::endl
				    <<" imag: "<<ev.imag()[0]<<" "<<ev.imag()[1]<<std::endl;
					return UNCLASSIFIED_SINGULARITY;
				}
			}
		}
		else if (complexImag == 2) {
			return CENTER_NODE;
		} else {
			std::cout << "[ERROR] Cannot classify the critical simplex with the following eigenvalues: "<<std::endl
				    <<" real: "<<ev.real()[0]<<" "<<ev.real()[1]<<std::endl
				    <<" imag: "<<ev.imag()[0]<<" "<<ev.imag()[1]<<std::endl;
      return UNCLASSIFIED_SINGULARITY;
    }  
  }
}

/**
 * VTK Wrapper to execute the above algorithm  
 */
class VTK_EXPORT VestecCriticalPointExtractionAlgorithm : public vtkDataSetAlgorithm
{
public:
  static VestecCriticalPointExtractionAlgorithm *New();
  vtkTypeMacro(VestecCriticalPointExtractionAlgorithm, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  VestecCriticalPointExtractionAlgorithm();
  ~VestecCriticalPointExtractionAlgorithm();

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector );

  virtual int FillOutputPortInformation(int port, vtkInformation* info) override;

private:
  VestecCriticalPointExtractionAlgorithm( const VestecCriticalPointExtractionAlgorithm& ); // Not implemented.
  void operator = ( const VestecCriticalPointExtractionAlgorithm& );  // Not implemented.
};


#endif