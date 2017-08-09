/*
  Copyright (C) 2017 Sven Willner <sven.willner@gmail.com>

  This file is part of DICE++.

  DICE++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  DICE++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with DICE++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "DICE.h"
#include "settingsnode.h"

using Time = size_t;
using Value = double;

static void print_usage(const char* program_name) {
    std::cerr << "DICE++\n"
                 "   C++-Implementation of the DICE Dynamic Integrated Climate-Economy\n"
                 "   Model of the Economics of Global Warming by W. Nordhaus as described in:\n\n"
                 "   Nordhaus, William D. The Climate Casino: Risk, Uncertainty, and Economics\n"
                 "   for a Warming World. Yale University Press (2013).\n\n"
                 "Version:  " DICEPP_VERSION
                 "\n"
                 "Author:   Sven Willner <sven.willner@pik-potsdam.de>\n\n"
                 "Source:   https://github.com/swillner/dicepp\n"
                 "License:  AGPL, (c) 2017 Sven Willner (see LICENSE file)\n\n"
                 "Usage:    "
              << program_name
              << " (<option> | <settingsfile>)\n"
                 "Options:\n"
                 "   -h, --help     Print this help text\n"
                 "   -v, --version  Print version"
              << std::endl;
}

int main(int argc, char* argv[]) {
#ifndef DEBUG
    try {
#endif
        if (argc != 2) {
            print_usage(argv[0]);
            return 1;
        }
        const std::string arg = argv[1];
        if (arg.length() > 1 && arg[0] == '-') {
            if (arg == "--version" || arg == "-v") {
                std::cout << DICEPP_VERSION << std::endl;
            } else if (arg == "--help" || arg == "-h") {
                print_usage(argv[0]);
            } else {
                print_usage(argv[0]);
                return 1;
            }
        } else {
            settings::SettingsNode settings;
            if (arg == "-") {
                std::cin >> std::noskipws;
                settings = settings::SettingsNode(std::cin);
            } else {
                std::ifstream settings_file(arg);
                settings = settings::SettingsNode(settings_file);
            }
            dice::DICE<Value, Time> dice(settings);
            dice.initialize();
            dice.run();
            dice.output();
        }
#ifndef DEBUG
    } catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 255;
    }
#endif
    return 0;
}
