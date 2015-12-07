#pragma once
#include <string>
#include "vertex.h"
#include "cellToVert.h"

using namespace std;

void readSettingData(string, vertex *&v, cellToVert *&CV, int *V_N, int *CELL_NUM);
