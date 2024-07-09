
#include "main.h" 
#include <boost/program_options.hpp> 
#include <Eigen/Dense> 
#include <cstring>
#include <iostream> 
#include <fstream>
#include <vector> 
#include <unistd.h>


namespace po = boost::program_options;
SimConfig sim;

Eigen::VectorXd quants(0);   // quantities 
Eigen::VectorXd prices(0); 
Eigen::VectorXd quants_0(0); 
Eigen::VectorXd prices_0(0); 
Eigen::MatrixXd Elas(0, 0);  // elasticity matrix 
Eigen::MatrixXd Elas_Inv(0, 0); 
Eigen::Matrix<double,Eigen::Dynamic,4> Prod_Costs(0, 4); 





int main(int argc, char* argv[]) {

    parse_options(argc, argv);
    std::cout << sim.to_string() << std::endl; 

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < 8; i++) {
        add_good(gen); 
        std::cout << Elas << "\n\n"; 
    }
    


    return 0;
}


void add_good(std::mt19937 gen) {

    // a + rand*(b-a) transforms the uniform distribution to produce rand outputs between a and b 
    std::uniform_real_distribution<> dist(0, 1); 

    // generate new entry for quants_0 and prices_0 vectors 
    int N = Elas.rows(); 
    double MPS = std::pow(10, LOG_MPS_MIN + dist(gen) * (LOG_MPS_MAX - LOG_MPS_MIN)); 
    double r = std::pow(10, LOG_R_MIN + dist(gen) * (LOG_R_MAX - LOG_R_MIN)); 
    double quant_0 = 2 * r; 
    double price_0 = MPS / r; 
    quants_0.conservativeResize(N + 1); 
    prices_0.conservativeResize(N + 1); 
    quants_0(N) = quant_0; 
    prices_0(N) = price_0; 

    // add new row/col to Elasticity matrix & recalculate its inverse 
    double ela_ii = quant_0 / price_0;  // diagonal entry of elasticity matrix for this good's row 
    Elas.conservativeResize(N+1, N+1); 
    Elas(N, N) = ela_ii; 
    for (int i = 0; i < N; i++) {
        Elas(N, i) =  (sim.compl_min + dist(gen) * (sim.compl_max - sim.compl_min)) * ela_ii; 
        Elas(i, N) = ((sim.compl_min + dist(gen) * (sim.compl_max - sim.compl_min)) * Elas(i, i)); 
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(Elas); 
    Elas_Inv = cod.pseudoInverse(); 

    quants.conservativeResize(N + 1); 
    prices.conservativeResize(N + 1); 
    quants(N) = 0; 
    prices(N) = Elas_Inv.row(N).dot(quants_0 - quants);

    // generate new production cost parameters (c0->fixed costs; c1->linear costs; c2->quadratic; c3->cubic) 
    double c0 = dist(gen) * MPS / 5; 
    double c3 = dist(gen) * 16 / 3 * price_0 / (quant_0*quant_0); 
    std::vector<double> intervals = {0.0, quant_0}; 
    std::vector<double> weights = {1.0, 0.0};  // max probability at 0, decreasing to zero at quant_0 
    std::piecewise_linear_distribution<> dist2(intervals.begin(), intervals.end(), weights.begin()); 

    double vertex_x = dist2(gen); 
    Eigen::VectorXd quants_temp = quants; 
    quants_temp(N) = vertex_x; 
    double vertex_y_max = Elas_Inv.row(N).dot(quants_0 - quants_temp); 
    double vertex_y = dist(gen) * vertex_y_max; 

    double c2 = -3 * c3 * vertex_x; 
    double c1 = vertex_y + 3 * c3 * vertex_x  * vertex_x; 
    Prod_Costs.conservativeResize(Prod_Costs.rows() + 1, Eigen::NoChange); 
    Prod_Costs.row(Prod_Costs.rows() - 1) << c0, c1, c2, c3; 
}



void parse_options(int argc, char* argv[]) {
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message") 
        ("config,c", po::value<std::string>(), "(str) config file to use") 
        ("timestep,t", po::value<std::string>(), "(str) timestep (day, week, month, or year)") 
        ("runtime,T", po::value<int>(), "(int) how long the sim should run, in years") 
        ("record,r", "write sim data to output files") 
        ("rec_step,R", po::value<std::string>(), "(int) output data every ___ timesteps") 
        ("log,l", "output debug and other log data") 
        ("warmup,w", po::value<int>(), "(int) number of warmup timesteps (no data is collected for these)")

        ("num_firms,f", po::value<int>(), "(int) number of firms in the simulation")
        ("ease_of_entry,e", po::value<float>(), "(float) ease of entry for new firms")
        ("compl_sign,s", po::value<std::string>(), "(string) possible signs for complementarity values (pos, neg, zero, or any)")
        ("compl_min,m", po::value<float>(), "(float) min complementarity magnitude")
        ("compl_max,M", po::value<float>(), "(float) max complementarity magnitude")
        ("flat_demand,d", "forces elasticity matrix to be zero, resulting in constant demand curve")
        ("uniform_demand,u", "q_0, p_0 entries for all firms will be same")
        ("uniform_costs,U", "all firms have same production cost parameters")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(1);
    }

    // Read options from the config file if provided.
    if (vm.count("config")) {
        std::string config_file = "config/" + vm["config"].as<std::string>() + ".cfg"; 
        std::ifstream ifs(config_file.c_str()); 
        if (ifs) {
            po::store(po::parse_config_file(ifs, desc), vm);
        } else {
            std::cerr << "Cannot open config file: " << config_file << "\n";
            exit(1);
        }
    }

    po::notify(vm);

    if (vm.count("timestep")) {
        std::string timestep_str = vm["timestep"].as<std::string>(); 
        if       (timestep_str == "day"  )  { sim.timestep = Timestep::DAY;   } 
        else if  (timestep_str == "week" )  { sim.timestep = Timestep::WEEK;  } 
        else if  (timestep_str == "month")  { sim.timestep = Timestep::MONTH; } 
        else if  (timestep_str == "year" )  { sim.timestep = Timestep::YEAR;  } 
        else  { 
            std::cerr << "Invalid timestep: " << timestep_str << "\n"; 
            exit(1); 
        }
    }

    if  (vm.count("runtime"))  { sim.runtime = vm["runtime"].as<int>();         }
    if  (vm.count("record"))   { sim.record = true;                             }
    if  (vm.count("rec_step")) { sim.recording_step = vm["rec_step"].as<int>(); }
    if  (vm.count("log"))      { sim.log = true;                                }
    if  (vm.count("warmup"))   { sim.warmup = vm["warmup"].as<int>();           }

    if  (vm.count("num_firms"))      { sim.num_firms = vm["num_firms"].as<int>();           }
    if  (vm.count("ease_of_entry"))  { sim.ease_of_entry = vm["ease_of_entry"].as<float>(); }
    if  (vm.count("compl_sign"))     { sim.compl_sign = vm["compl_sign"].as<std::string>(); } 
    if  (vm.count("compl_min"))      { sim.compl_min = vm["compl_min"].as<float>();         }
    if  (vm.count("compl_max"))      { sim.compl_max = vm["compl_max"].as<float>();         }
    if  (vm.count("flat_demand"))    { sim.flat_demand = true;                              }
    if  (vm.count("uniform_demand")) { sim.uniform_demand = true;                           }
    if  (vm.count("uniform_costs"))  { sim.uniform_costs = true;                            }
}


std::string SimConfig::to_string() const {
    return "Simulation Config Parameters: \n=============================" 
            + std::string("\nnum_firms: ") + std::to_string(num_firms) 
            + std::string("\nease_of_entry: ") + std::to_string(ease_of_entry) 
            + std::string("\ncompl_sign: ") + compl_sign
            + std::string("\ncompl_min: ") + std::to_string(compl_min) 
            + std::string("\ncompl_max: ") + std::to_string(compl_max) 
            + std::string("\nflat_demand: ") + std::to_string(flat_demand) 
            + std::string("\nuniform_demand: ") + std::to_string(uniform_demand) 
            + std::string("\nuniform_costs: ") + std::to_string(uniform_costs) 
            + std::string("\ntimestep: ") + std::to_string(timestep) 
            + std::string("\nruntime: ") + std::to_string(runtime) 
            + std::string("\nrecord: ") + std::to_string(record) 
            + std::string("\nrecording_step: ") + std::to_string(recording_step) 
            + std::string("\nlog: ") + std::to_string(log) 
            + std::string("\nwarmup: ") + std::to_string(warmup); 
}



