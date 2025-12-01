#include "utils.h"


void save_data(std::vector<double> &data, std::string filename, float E_min)
{
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < data.size(); ++i)
    {
        file << E_min + i << "\t" << data[i] << std::endl;
    }
    file.close();
};

std::string GetFilename(int L, float J, int q)
{
    std::string J_str = std::to_string(J);
    J_str = J_str.substr(0, J_str.find(".") + 3);
    std::string filename = "L" + std::to_string(L) + "_J" + J_str + "_q" + std::to_string(q) + ".txt";
    return filename;
};
