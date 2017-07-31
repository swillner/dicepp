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
                 "Version:  " DICEPP_VERSION
                 "\n"
                 "Author:   Sven Willner <sven.willner@pik-potsdam.de>\n\n"
                 "Source:   https://github.com/swillner/dicepp\n"
                 "License:  AGPL, (c) 2017 Sven Willner (see LICENSE file)\n\n"
                 "Usage: "
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
