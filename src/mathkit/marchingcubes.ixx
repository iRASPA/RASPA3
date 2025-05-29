/**
 * @file    MarchingCubes.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MarchingCubes Algorithm
 */
//________________________________________________

module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <vector>
#endif

export module marching_cubes;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <cstddef>;
import <vector>;
#endif

import int3;

//_____________________________________________________________________________
// types
/** unsigned char alias */
typedef unsigned char uchar ;
/** signed char alias */
typedef   signed char schar ;
/** isovalue alias */
typedef        double real  ;

//-----------------------------------------------------------------------------
// Vertex structure
/** \struct Vertex "MarchingCubes.h" MarchingCubes
 * Position and normal of a vertex
 * \brief vertex structure
 * \param x X coordinate
 * \param y Y coordinate
 * \param z Z coordinate
 * \param nx X component of the normal
 * \param ny Y component of the normal
 * \param nz Z component of the normal
 */
export typedef struct
{
  real  x,  y,  z ;  /**< Vertex coordinates */
  real nx, ny, nz ;  /**< Vertex normal */
} Vertex ;

//-----------------------------------------------------------------------------
// Triangle structure
/** \struct Triangle "MarchingCubes.h" MarchingCubes
 * Indices of the oriented triange vertices
 * \brief triangle structure
 * \param v1 First vertex index
 * \param v2 Second vertex index
 * \param v3 Third vertex index
 */
export typedef struct
{
  int v1,v2,v3 ;  /**< Triangle vertices */
} Triangle ;
//_____________________________________________________________________________



//_____________________________________________________________________________
/** Marching Cubes algorithm wrapper */
/** \class MarchingCubes
  * \brief Marching Cubes algorithm.
  */
export class MarchingCubes
//-----------------------------------------------------------------------------
{
// Constructors
public :
  /**
   * Main and default constructor
   * \brief constructor
   * \param size_x width  of the grid
   * \param size_y depth  of the grid
   * \param size_z height of the grid
   */
  MarchingCubes ( const int size_x = -1, const int size_y = -1, const int size_z = -1 ) ;

//-----------------------------------------------------------------------------
// Accessors
public :
  /** accesses the number of vertices of the generated mesh */
  inline size_t nverts() const { return _vertices.size() ; }
  /** accesses the number of triangles of the generated mesh */
  inline size_t ntrigs() const { return _triangles.size() ; }
  /** accesses a specific vertex of the generated mesh */
  inline const Vertex   * vert( const ptrdiff_t i ) const { if( i < 0  || i >= static_cast<ptrdiff_t>(nverts()) ) return static_cast<Vertex *>(nullptr) ; return _vertices.data()  + i ; }
  /** accesses a specific triangle of the generated mesh */
  inline const Triangle * trig( const ptrdiff_t i ) const { if( i < 0  || i >= static_cast<ptrdiff_t>(ntrigs()) ) return static_cast<Triangle*>(nullptr) ; return _triangles.data() + i ; }

  /** accesses the vertex buffer of the generated mesh */
  inline Vertex   *vertices () { return _vertices.data()  ; }
  /** accesses the triangle buffer of the generated mesh */
  inline Triangle *triangles() { return _triangles.data() ; }

  /**  accesses the width  of the grid */
  inline int size_x() const { return _size_x ; }
  /**  accesses the depth  of the grid */
  inline int size_y() const { return _size_y ; }
  /**  accesses the height of the grid */
  inline int size_z() const { return _size_z ; }

  /**
   * changes the size of the grid
   * \param size_x width  of the grid
   * \param size_y depth  of the grid
   * \param size_z height of the grid
   */
  inline void set_resolution( const int size_x, const int size_y, const int size_z ) { _size_x = size_x ;  _size_y = size_y ;  _size_z = size_z ; }
  /**
   * selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes
   * \param originalMC true for the original Marching Cubes
   */
  inline void set_method    ( const bool originalMC = false ) { _originalMC = originalMC ; }

  // Data access
  /**
   * accesses a specific cube of the grid
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
    inline double get_data(const int3 &coord) const {
        return _data[ static_cast<size_t>(coord.x + coord.y*_size_x + coord.z*_size_x*_size_y)];
    }

  /**
   * sets a specific cube of the grid
   * \param val new value for the cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_data  ( const real val, const size_t i, const size_t j, const size_t k ) { _data[ i + j*static_cast<size_t>(_size_x) + k*static_cast<size_t>(_size_x)*static_cast<size_t>(_size_y)] = val ; }

  // Data initialization
  /** inits temporary structures (must set sizes before call) : the grid and the vertex index per cube */
  void init_temps () ;
  /** inits all structures (must set sizes before call) : the temporary structures and the mesh buffers */
  void init_all   () ;


//-----------------------------------------------------------------------------
// Algorithm
public :
  /**
   * Main algorithm : must be called after init_all
   * \param iso isovalue
   */
  void run( real iso = 0.0 ) ;

protected :
  /** tesselates one cube */
  void process_cube (double *cube);
  /** tests if the components of the tesselation of the cube should be connected through the interior of the cube */
  bool test_interior( schar s, double *cube )    ;


//-----------------------------------------------------------------------------
// Operations
protected :
  /**
   * computes almost all the vertices of the mesh by interpolation along the cubes edges
   * \param iso isovalue
   */
  void compute_intersection_points( real iso ) ;

  /**
   * routine to add a triangle to the mesh
   * \param trig the code for the triangle as a sequence of edges index
   * \param n    the number of triangles to produce
   * \param v12  the index of the interior vertex to use, if necessary
   */
  void add_triangle ( const char* trig, char n, int v12 = -1 ) ;

  int add_vertex(const int3 &grid_coord, const int3&dir, int corner, double *cube);
  /** adds a vertex inside the current cube */
  int add_c_vertex() ;

  /**
   * interpolates the horizontal gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_x_grad( const int i, const int j, const int k ) const ;
  /**
   * interpolates the longitudinal gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_y_grad( const int i, const int j, const int k ) const ;
  /**
   * interpolates the vertical gradient of the implicit function at the lower vertex of the specified cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  real get_z_grad( const int i, const int j, const int k ) const ;

  /**
   * accesses the pre-computed vertex index on the lower horizontal edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_x_vert( const int i, const int j, const int k ) const { return _x_verts[static_cast<size_t>( i + j*_size_x + k*_size_x*_size_y)] ; }
  /**
   * accesses the pre-computed vertex index on the lower longitudinal edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_y_vert( const int i, const int j, const int k ) const { return _y_verts[static_cast<size_t>( i + j*_size_x + k*_size_x*_size_y)] ; }
  /**
   * accesses the pre-computed vertex index on the lower vertical edge of a specific cube
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline int   get_z_vert( const int i, const int j, const int k ) const { return _z_verts[static_cast<size_t>( i + j*_size_x + k*_size_x*_size_y)] ; }

  /**
   * sets the pre-computed vertex index on the lower horizontal edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_x_vert( const int val, const int i, const int j, const int k ) { _x_verts[static_cast<size_t>( i + j*_size_x + k*_size_x*_size_y)] = val ; }
  /**
   * sets the pre-computed vertex index on the lower longitudinal edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_y_vert( const int val, const int i, const int j, const int k ) { _y_verts[static_cast<size_t>( i + j*_size_x + k*_size_x*_size_y)] = val ; }
  /**
   * sets the pre-computed vertex index on the lower vertical edge of a specific cube
   * \param val the index of the new vertex
   * \param i abscisse of the cube
   * \param j ordinate of the cube
   * \param k height of the cube
   */
  inline void  set_z_vert( const int val, const int i, const int j, const int k ) { _z_verts[static_cast<size_t>( i + j*_size_x + k*_size_x*_size_y)] = val ; }

//-----------------------------------------------------------------------------
// Elements
protected :
  bool      _originalMC ;   /**< selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes */

  int       _size_x     ;  /**< width  of the grid */
  int       _size_y     ;  /**< depth  of the grid */
  int       _size_z     ;  /**< height of the grid */
  std::vector<double> _data;

    std::vector<int> _x_verts    ;  /**< pre-computed vertex indices on the lower horizontal   edge of each cube */
    std::vector<int> _y_verts    ;  /**< pre-computed vertex indices on the lower longitudinal edge of each cube */
    std::vector<int> _z_verts    ;  /**< pre-computed vertex indices on the lower vertical     edge of each cube */

  std::vector<Vertex> _vertices   ;  /**< vertex   buffer */
    std::vector<Triangle> _triangles  ;  /**< triangle buffer */

  int       _i          ;  /**< abscisse of the active cube */
  int       _j          ;  /**< height of the active cube */
  int       _k          ;  /**< ordinate of the active cube */

  uchar     _lut_entry  ;  /**< cube sign representation in [0..255] */
  uchar     _case       ;  /**< case of the active cube in [0..15] */
  uchar     _config     ;  /**< configuration of the active cube */
  uchar     _subconfig  ;  /**< subconfiguration of the active cube */

public:
  static const char cases[256][2];
static const char tiling1[16][3];
static const char tiling2[24][6];
static const char test3[24];
static const char tiling3_1[24][6];
static const char tiling3_2[24][12];
static const char test4[8];
static const char tiling4_1[8][6];
static const char tiling4_2[8][18];
static const char tiling5[48][9];
static const char test6[48][3];
static const char tiling6_1_1[48][9];
static const char tiling6_1_2[48][27];
static const char tiling6_2[48][15];
static const char test7[16][5];
static const char tiling7_1[16][9];
static const char tiling7_2[16][3][15];
static const char tiling7_3[16][3][27];
static const char tiling7_4_1[16][15];
static const char tiling7_4_2[16][27];
static const char tiling8[6][6];
static const char tiling9[8][12];
static const char test10[6][3];
static const char tiling10_1_1[6][12];
static const char tiling10_1_1_[6][12];
static const char tiling10_1_2[6][24];
static const char tiling10_2[6][24];
static const char tiling10_2_[6][24];
static const char tiling11[12][12];
static const char test12[24][4];
static const char tiling12_1_1[24][12];
static const char tiling12_1_1_[24][12];
static const char tiling12_1_2[24][24];
static const char tiling12_2[24][24];
static const char tiling12_2_[24][24];
static const char test13[2][7];
static const char subconfig13[64];
static const char tiling13_1[2][12];
static const char tiling13_1_[2][12];
static const char tiling13_2[2][6][18];
static const char tiling13_2_[2][6][18];
static const char tiling13_3[2][12][30];
static const char tiling13_3_[2][12][30];
static const char tiling13_4[2][4][36];
static const char tiling13_5_1[2][4][18];
static const char tiling13_5_2[2][4][30];
static const char tiling14[12][12];
static const char casesClassic[256][16];
bool test_face( schar face, double *cube );
};
//_____________________________________________________________________________

