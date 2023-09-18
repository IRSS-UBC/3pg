#ifndef BYTEGRID_H
#define BYTEGRID_H

//#ident "@(#)ByteGrid.hpp 1.1 98/03/31 13:27:15 JCG"

#include "Grid.hpp"


typedef unsigned char Byte;

class ByteGrid : public Grid {
public:
  Byte nodata;
  Byte *b;

  void Initialise();
  
  ByteGrid() : Grid() { b = 0; }
  ByteGrid(Grid &grid) : Grid(grid) { b = 0; }
  ByteGrid(ByteGrid &grid) : Grid(grid) { b = 0; nodata = grid.nodata; }
  
  ~ByteGrid();

  ByteGrid& operator=(ByteGrid &from);

  bool Exists(char *fname);
  void Allocate();
  void ResetGrid();
  void CloseGrid();
  void PrintGridState(void);
  void Read(char *fname);// throw (Exception);
  void Write(char *fname);// throw (Exception);
};


#endif
