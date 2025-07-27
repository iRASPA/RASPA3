module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <string>
#include <vector>
#endif

module skseitzintegermatrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import double3x3;
import skonethirdseitzmatrix;
import skrotationmatrix;

SKSeitzIntegerMatrix::SKSeitzIntegerMatrix() {}

SKSeitzIntegerMatrix::SKSeitzIntegerMatrix(SKRotationMatrix rotation, int3 translation)
{
  this->rotation = rotation;
  this->translation = translation;  // denominator = 24
}

SKSeitzIntegerMatrix::SKSeitzIntegerMatrix(char xvalue, char yvalue, char zvalue)
{
  char referenceValue = '0';

  std::size_t x = std::size_t(xvalue - referenceValue);
  std::size_t y = std::size_t(yvalue - referenceValue);
  std::size_t z = std::size_t(zvalue - referenceValue);

  // qDebug() << "xyz: " << x <<", " << y <<"," << z;

  // assert(x < SKSeitzIntegerMatrix::SeitzData.size());
  // assert(y < SKSeitzIntegerMatrix::SeitzData.size());
  // assert(z < SKSeitzIntegerMatrix::SeitzData.size());

  int3 r1 = int3(SKSeitzIntegerMatrix::SeitzData[x].r1(), SKSeitzIntegerMatrix::SeitzData[y].r1(),
                 SKSeitzIntegerMatrix::SeitzData[z].r1());
  int3 r2 = int3(SKSeitzIntegerMatrix::SeitzData[x].r2(), SKSeitzIntegerMatrix::SeitzData[y].r2(),
                 SKSeitzIntegerMatrix::SeitzData[z].r2());
  int3 r3 = int3(SKSeitzIntegerMatrix::SeitzData[x].r3(), SKSeitzIntegerMatrix::SeitzData[y].r3(),
                 SKSeitzIntegerMatrix::SeitzData[z].r3());
  this->rotation = SKRotationMatrix(r1, r2, r3);
  this->translation = int3(SKSeitzIntegerMatrix::SeitzData[x].t(), SKSeitzIntegerMatrix::SeitzData[y].t(),
                           SKSeitzIntegerMatrix::SeitzData[z].t());
}

std::vector<SKSeitzIntegerMatrix> SKSeitzIntegerMatrix::SeitzMatrices(std::string encoding)
{
  std::size_t m = encoding.size() / 3;

  // qDebug() << QString::fromStdString(encoding) << ", " << m;

  std::vector<SKSeitzIntegerMatrix> matrices = std::vector<SKSeitzIntegerMatrix>();
  matrices.resize(m);

  for (std::size_t i = 0; i < m; i++)
  {
    char x = encoding[3 * i];
    char y = encoding[3 * i + 1];
    char z = encoding[3 * i + 2];

    matrices[i] = SKSeitzIntegerMatrix(x, y, z);
  }

  return matrices;
}

std::vector<SKOneThirdSeitzMatrix> SKSeitzIntegerMatrix::SeitzData = std::vector<SKOneThirdSeitzMatrix>{
    SKOneThirdSeitzMatrix(" x", '0', 1, 0, 0, 0),         SKOneThirdSeitzMatrix(" y", '1', 0, 1, 0, 0),
    SKOneThirdSeitzMatrix(" z", '2', 0, 0, 1, 0),         SKOneThirdSeitzMatrix("-x", '3', -1, 0, 0, 0),
    SKOneThirdSeitzMatrix("-y", '4', 0, -1, 0, 0),        SKOneThirdSeitzMatrix("-z", '5', 0, 0, -1, 0),
    SKOneThirdSeitzMatrix("x-y", '6', 1, -1, 0, 0),       SKOneThirdSeitzMatrix("-x+y", '7', -1, 1, 0, 0),
    SKOneThirdSeitzMatrix("x+1/2", '8', 1, 0, 0, 12),     SKOneThirdSeitzMatrix("x+1/3", '9', 1, 0, 0, 8),
    SKOneThirdSeitzMatrix("x+1/4", ':', 1, 0, 0, 6),      SKOneThirdSeitzMatrix("x+2/3", ';', 1, 0, 0, 16),
    SKOneThirdSeitzMatrix("x+3/4", '<', 1, 0, 0, 18),     SKOneThirdSeitzMatrix("y+1/2", '=', 0, 1, 0, 12),
    SKOneThirdSeitzMatrix("y+1/3", '>', 0, 1, 0, 8),      SKOneThirdSeitzMatrix("y+1/4", '?', 0, 1, 0, 6),
    SKOneThirdSeitzMatrix("y+2/3", '@', 0, 1, 0, 16),     SKOneThirdSeitzMatrix("y+3/4", 'A', 0, 1, 0, 18),
    SKOneThirdSeitzMatrix("z+1/2", 'B', 0, 0, 1, 12),     SKOneThirdSeitzMatrix("z+1/3", 'C', 0, 0, 1, 8),
    SKOneThirdSeitzMatrix("z+1/4", 'D', 0, 0, 1, 6),      SKOneThirdSeitzMatrix("z+1/6", 'E', 0, 0, 1, 4),
    SKOneThirdSeitzMatrix("z+2/3", 'F', 0, 0, 1, 16),     SKOneThirdSeitzMatrix("z+3/4", 'G', 0, 0, 1, 18),
    SKOneThirdSeitzMatrix("z+5/6", 'H', 0, 0, 1, 20),     SKOneThirdSeitzMatrix("-x+1/2", 'I', -1, 0, 0, 12),
    SKOneThirdSeitzMatrix("-x+1/3", 'J', -1, 0, 0, 8),    SKOneThirdSeitzMatrix("-x+1/4", 'K', -1, 0, 0, 6),
    SKOneThirdSeitzMatrix("-x+2/3", 'L', -1, 0, 0, 16),   SKOneThirdSeitzMatrix("-x+3/4", 'M', -1, 0, 0, 18),
    SKOneThirdSeitzMatrix("-y+1/2", 'N', 0, -1, 0, 12),   SKOneThirdSeitzMatrix("-y+1/3", 'O', 0, -1, 0, 8),
    SKOneThirdSeitzMatrix("-y+1/4", 'P', 0, -1, 0, 6),    SKOneThirdSeitzMatrix("-y+2/3", 'Q', 0, -1, 0, 16),
    SKOneThirdSeitzMatrix("-y+3/4", 'R', 0, -1, 0, 18),   SKOneThirdSeitzMatrix("-z+1/2", 'S', 0, 0, -1, 12),
    SKOneThirdSeitzMatrix("-z+1/3", 'T', 0, 0, -1, 8),    SKOneThirdSeitzMatrix("-z+1/4", 'U', 0, 0, -1, 6),
    SKOneThirdSeitzMatrix("-z+1/6", 'V', 0, 0, -1, 4),    SKOneThirdSeitzMatrix("-z+2/3", 'W', 0, 0, -1, 16),
    SKOneThirdSeitzMatrix("-z+3/4", 'X', 0, 0, -1, 18),   SKOneThirdSeitzMatrix("-z+5/6", 'Y', 0, 0, -1, 20),
    SKOneThirdSeitzMatrix("x-y+1/3", 'Z', 1, -1, 0, 8),   SKOneThirdSeitzMatrix("x-y+2/3", '[', 1, -1, 0, 16),
    SKOneThirdSeitzMatrix("-x+y+1/3", '\\', -1, 1, 0, 8), SKOneThirdSeitzMatrix("-x+y+2/3", ']', -1, 1, 0, 16)};
