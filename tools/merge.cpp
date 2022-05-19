#include <cstdlib>
#include <iostream>
#include <vector>

#include <range/v3/action/push_back.hpp>

#include "arguments.hpp"

[[noreturn]] static void Usage() {

    std::cout << "Usage:\n";
    std::cout << "merge [-z,--gzip] -i,--input file1 file2 ... [-o,--output filename]\n";
    std::cout << "  -i,--input     List of input files\n";
    std::cout << "  -z,--gzip      Input files are gzip compressed\n";
    std::cout << "  -o,--output    Save the output to filename (default is merged.pb)\n";

    std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {

    std::cerr << "Run with -h or --help to see usage.\n";

    std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {

    Arguments args = GetArguments(argc, argv);

    bool gzip = false;
    std::vector<std::string_view> input_filenames;
    std::string result_filename = "merged.pb";

    for (auto [name, params] : args) {
        if (name == "-h" or name == "--help") {
            Usage();
        } else if (name == "-i" or name == "--input") {
            ranges::action::push_back(input_filenames, params);
        } else if (name == "-z" or name == "--gzip") {
            gzip = true;
        } else if (name == "-o" or name == "--output") {
            if (params.empty()) {
                std::cerr << "Specify result file name.\n";
                Fail();
            }
            result_filename = *params.begin();
        }
    }

    if (input_filenames.size() < 2) {
        std::cerr << "Specify at least two input file names.\n";
        Fail();
    }

    std::cout << "Merging ";
    if (gzip) {
        std::cout << "gzip compressed ";
    }
    for (auto& i : input_filenames) {
        std::cout << i << " ";
    }
    std::cout << "into " << result_filename << "\n";

    return EXIT_SUCCESS;
}
