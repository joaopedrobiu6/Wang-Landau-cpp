#include "WangLandau.h"
#include <omp.h>
#include <filesystem>

bool isFlat(const std::vector<int> &hist, double h_tol)
{
    // Remove entries with zero counts
    std::vector<int> values;
    for (const auto &entry : hist)
    {
        if (entry > 0)
        {
            values.push_back(entry);
        }
    }

    if (values.empty())
    {
        return false;
    }

    double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double min_value = *std::min_element(values.begin(), values.end());

    if (min_value < h_tol * mean)
    {
        // std::cout << "it's not flat yet" << std::endl;
        return false;
    }
    return true;
}

void info_dump(std::string folder, int MC_N, int L, int q, double f_tol, double h_tol, bool NoLog, int sampleInterval, int f_factor)
{
    std::ofstream info_file(folder + "/info.txt");
    info_file << "MC samples: " << MC_N << std::endl;
    info_file << "Lattice size: " << L << std::endl;
    info_file << "q: " << q << std::endl;
    info_file << "f tolerance: " << f_tol << std::endl;
    info_file << "Histogram tolerance: " << h_tol << std::endl;
    if (NoLog)
    {
        info_file << "Exp form" << std::endl;
    }
    else
    {
        info_file << "Log form" << std::endl;
    }
    info_file << "Sample Check Interval: " << sampleInterval << std::endl;
    info_file << "f update factor: " << f_factor << std::endl;
    info_file.close();
}

void info_print(int MC_N, int L, int q, double f_tol, double h_tol, bool NoLog, int sampleInterval, int f_factor)
{
    std::cout << "\n*******************************************" << std::endl;
    std::cout << "\nStarting Wang-Landau simulation" << std::endl;
    std::cout << "MC samples: " << MC_N << std::endl;
    std::cout << "Lattice size: " << L << std::endl;
    std::cout << "q: " << q << std::endl;
    std::cout << "f tolerance: " << f_tol << std::endl;
    std::cout << "Histogram tolerance: " << h_tol << std::endl;
    if (NoLog)
    {
        std::cout << "Exp form" << std::endl;
    }
    else
    {
        std::cout << "Log form" << std::endl;
    }
    std::cout << "Sample Check Interval: " << sampleInterval << std::endl;
    std::cout << "f update factor: " << f_factor << std::endl;
    std::cout << "\n*******************************************\n"
              << std::endl;
}

std::map<int, double> WangLandauPotts(PottsLattice lat, int MC_N, int q, double f_tol, double h_tol, bool NoLog, int sampleInterval, int f_factor, std::string filename, std::string save_path)
{
    int L = lat.lattice.size();

    info_print(MC_N, L, q, f_tol, h_tol, NoLog, sampleInterval, f_factor);

    // create output folder if it doesn't exist with name equal to filename without extension and in result folder
    std::string folder_name = save_path + "/" + filename.substr(0, filename.find_last_of("."));
    std::filesystem::create_directories(folder_name);

    // create info file in output folder with parameters
    info_dump(folder_name, MC_N, L, q, f_tol, h_tol, NoLog, sampleInterval, f_factor);
    srand(time(NULL));
    std::pair<float, float> E_limit = lat.Energy_Limit();
    float E_min = E_limit.first;
    float E_max = E_limit.second;
    int num_bins = (int)(E_max - E_min) + 1;

    std::vector<double> lng(num_bins, 0.0);
    std::vector<double> g(num_bins, 1.0);
    std::vector<int> hist(num_bins, 0);

    double f = std::exp(1.0);
    double lnf = std::log(f);

    // Initialize RNG
    std::uniform_int_distribution<int> dist_L(0, L - 1);
    std::uniform_int_distribution<int> dist_q(1, q);
    std::uniform_real_distribution<double> dist_real(0.0, 1.0);

    // Outer loop: check flatness only between parallel regions
    if (NoLog)
    {
        while (f - 1 > f_tol)
        {
            std::fill(hist.begin(), hist.end(), 0);

// Parallel Region
#pragma omp parallel
            {
                int tid = omp_get_thread_num();

                // Each thread has its own walker (lattice)
                PottsLattice local_lat = lat; // Copy constructor

                // Seed local RNG differently for each thread
                std::random_device rd;
                local_lat.rng.seed(rd() ^ tid);

                // Thread-local copy of g to reduce atomic reads
                std::vector<double> g_local = g;

                double Current_E = local_lat.Potts_Energy();

                for (int i = 0; i < MC_N; i++)
                {
                    int x = dist_L(local_lat.rng);
                    int y = dist_L(local_lat.rng);
                    int s0 = local_lat.lattice[x][y];

                    int s1 = dist_q(local_lat.rng);

                    // Calculate energy change without modifying lattice
                    double delta_E = local_lat.Energy_Change(x, y, s1);
                    double New_E = Current_E + delta_E;

                    int idx_current = (int)(Current_E - E_min);
                    int idx_new = (int)(New_E - E_min);

                    // Use thread-local g copy for reads
                    double g_current = g_local[idx_current];
                    double g_new = g_local[idx_new];

                    if (dist_real(local_lat.rng) < std::min(1.0, g_current / g_new))
                    {
                        local_lat.lattice[x][y] = s1;
                        Current_E = New_E;
                        idx_current = idx_new;
                    }

#pragma omp atomic
                    lng[idx_current] += lnf;
#pragma omp atomic
                    hist[idx_current] += 1;
#pragma omp atomic
                    g[idx_current] *= f;
                    g_local[idx_current] *= f;
                }
            } // End parallel region

            // Check flatness and update f (outside parallel region)
            std::cout << "f: " << f << std::endl;
            bool flat = isFlat(hist, h_tol);

            if (flat)
            {
                f = pow(f, 1.0 / f_factor);
            }
        }
    }
    else
    {
        while (lnf > f_tol)
        {
            std::fill(hist.begin(), hist.end(), 0);

// Parallel Region
#pragma omp parallel
            {
                int tid = omp_get_thread_num();

                // Each thread has its own walker (lattice)
                PottsLattice local_lat = lat; // Copy constructor

                // Seed local RNG differently for each thread
                std::random_device rd;
                local_lat.rng.seed(rd() ^ tid);

                double Current_E = local_lat.Potts_Energy();

                for (int i = 0; i < MC_N; i++)
                {
                    int x = dist_L(local_lat.rng);
                    int y = dist_L(local_lat.rng);
                    int s0 = local_lat.lattice[x][y];

                    int s1 = dist_q(local_lat.rng);

                    double delta_E = local_lat.Energy_Change(x, y, s1);
                    double New_E = Current_E + delta_E;

                    int idx_current = (int)(Current_E - E_min);
                    int idx_new = (int)(New_E - E_min);

                    // Use global lng for reads
                    double lng_current = lng[idx_current];
                    double lng_new = lng[idx_new];

                    if (dist_real(local_lat.rng) < std::exp(lng_current - lng_new))
                    {
                        local_lat.lattice[x][y] = s1;
                        Current_E = New_E;
                        idx_current = idx_new;
                    }

#pragma omp atomic
                    lng[idx_current] += lnf;
#pragma omp atomic
                    hist[idx_current] += 1;
                }
            } // End parallel region

            // Check flatness and update lnf (outside parallel region)
            bool flat = isFlat(hist, h_tol);

            if (flat)
            {
                lnf = lnf / f_factor;
            }
        }
    }

    filename = save_path + "/" + filename;
    save_data(lng, filename, E_min);

    // Convert back to map for return type compatibility (or change return type)
    // Changing return type would require changing header and main.
    // For now, let's construct the map to return.
    std::map<int, double> result_map;
    for (int i = 0; i < lng.size(); ++i)
    {
        if (lng[i] != 0.0)
        { // Only store visited states? Or all? Original stored all in range.
            result_map[(int)(E_min + i)] = lng[i];
        }
    }
    return result_map;
}