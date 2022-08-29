export module isotherm_fitting_simulation;

import <vector>;
import <span>;
import <tuple>;
import <string>;
import <fstream>;

import input_reader;
import component;
import system;
import multi_site_isotherm;
import isotherm_fitting;


export struct IsothermFittingSimulation
{
  public:
    IsothermFittingSimulation(InputReader &inputreader);

    void run();

  private:
    std::vector<System> systems;
};
