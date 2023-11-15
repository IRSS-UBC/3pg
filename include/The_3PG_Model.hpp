// 3PG model run routine.
#include "Data_io.hpp"
#include "Params.hpp"

void runTreeModel( MYDate spMinMY, MYDate spMaxMY, bool spatial, long cellIndex, std::vector<PPPG_PARAM>& params, std::vector<PPPG_SERIES_PARAM>& series, std::vector<PPPG_OP_VAR>& outputs);
