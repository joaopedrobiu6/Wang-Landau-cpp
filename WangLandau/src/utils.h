#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <map>


void save_data(std::vector<double> &data, std::string filename, float E_min);

std::string GetFilename(int L, float J, int q);
