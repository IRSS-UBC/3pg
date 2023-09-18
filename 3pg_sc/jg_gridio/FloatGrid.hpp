#ifndef FLOATGRID_H
#define FLOATGRID_H

//#ident "@(#)FloatGrid.hpp     1.3 98/12/23 14:08:31 JCG"

#include <windows.h>
#include "Grid.hpp"


class FloatGrid : public Grid {
private:
  static float defaultnodata;
  char filename[100];
  
  
public:
  float nodata;
  float* z;
  FILE* fp;
 
  HANDLE hFile; 
  HANDLE hMapFile; 
  LPVOID lpMapAddress; 
  bool bFlag;
  bool runPeriodDone;
  int currRow;
  int lastrow;
  int lastcell;

  FloatGrid() : Grid() { z = 0; nodata = defaultnodata; }
  FloatGrid(Grid &grid) : Grid(grid) { z = 0; nodata = defaultnodata; }
  FloatGrid(FloatGrid &grid) : Grid(grid) { z = 0; nodata = defaultnodata; }

  //New constructor for the purpose of opening a file for writing
  FloatGrid(char* fname);

  ~FloatGrid();

  FloatGrid& operator=(FloatGrid &from);
  
  //void Allocate();
  void Initialise();

  
  int InLimits(int i, int j);
  int InLimits(double x, double y);

  bool Exists(char *fname);  //ANL
  void Read(char *fname); // throw (Exception);
  float GetVal(int cellIndex);
  void Write(char *fname); //throw (Exception);
  void Allocate();
  void ResetGrid();
  void CloseGrid();
  void PrintGridState(void);
};


#endif
