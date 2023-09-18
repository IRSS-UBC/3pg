#ifndef GRID_H
#define GRID_H

//#ident "%W% %D% %T% JCG"

// Turn of 'C++ Exception Specification ignored' warnings on VC++6. 
#ifdef WIN32
#pragma warning( disable : 4290 ) 
#endif

class Exception;
class BasicPoint;
/*
  Defines a Grid data structure without any actual data storage.
  Derived classes can declare their own type of storage and
  provide appropriate Read and Write functions.

  The Grid is designed to have a single cell border that is
  outside the area specified by the xmin/xmax/ymin/ymax rectangle.
  These values specify the centres of cells inside the border
  */


class Grid {
public:
  double cellsize;
  double xmin, xmax, ymin, ymax;
  int nrows, ncols;
  long npoints;
  int offset[8];

  Grid();                     // constructor
  Grid(Grid &oldgrid);        // assignment constructor
  Grid(int ncols, int nrows, double xmin, double ymin, double cellsize);
    // creates from specifications; adds border
  Grid(double xmin, double xmax, double ymin, double ymax, double cellsize);
    // creates from specifications; adds border
  virtual ~Grid() {}

  Grid& operator=(Grid& from);
  
  virtual bool Exists(char *fname);// throw (Exception) = 0; // ANL
  virtual void Allocate();// throw (Exception) = 0;  // ANL
  virtual void Read(char *fname);// throw (Exception) = 0;
  virtual void Write(char *fname);// throw (Exception) = 0;
  virtual void ResetGrid(void);// throw (Exception) = 0;
  virtual void CloseGrid(void);// throw (Exception) = 0;
  virtual void PrintGridState(void);// throw (Exception) = 0;


  

  /* InLimits functions return 0 if not within grid limits,
     1 if within grid limits and 2 if within grid limits and
     a valid value is available at that point.  The Grid implementation
     does not involve nodata, so will only check bounds and return
     0 or 1.  Note that the one cell border is EXCLUDED from the
     region considered inside the limits.
     */
  virtual int InLimits(int i, int j);
  virtual int InLimits(double x, double y);

  /* the following functions throw exceptions if the supplied
     values are outside the grid limits */
  void XYtoIJ(double x, double y, int &i, int &j);// throw (Exception);
  void IJtoXY(int i, int j, double &x, double &y);// throw (Exception);
  long IJtoIndex(int i, int j);// throw (Exception);
  void IndextoIJ(long k, int &i, int &j);// throw (Exception);
  long XYtoIndex(double x, double y);// throw (Exception);
  /*  BasicPoint* CellEdge(double x0, double y0, double x1, double y1,
                       int *newi = 0, int *newj = 0, int trace = 0); */

};



#endif
