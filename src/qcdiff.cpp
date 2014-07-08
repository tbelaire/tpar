
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
int main(int argc, char *argv[]) {

    bool disp_log;

    po::options_description desc("Allowed Options");
    desc.add_options()
        ("help", "displays help")
        ("verbose,v", "Display additional logging")
        ;
    po::options_description hidden_desc("Never see this");
    hidden_desc.add_options()
        ("left-file", po::value<string>(), "")
        ("right-file", po::value<string>(), "")
        ;
    po::positional_options_description pos_desc;
    pos_desc.add("left-file", 1).add("right-file", 2);

    po::options_description all_options;
    all_options.add(desc).add(hidden_desc);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(pos_desc)
                .run(), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << desc << endl;
            return 0;
        }
    }
    catch(std::exception& e)
    {
        cout << "Error: " << e.what() << endl;
        cout << desc << endl;
        return 1;
    }

    disp_log = vm.count("verbose");

    cout << "Logging: " << (disp_log ? "Yes" : "No") << endl;

    if (vm.count("left-file") && vm.count("right-file")) {
        string left_file(vm["left-file"].as<string>());
        string right_file(vm["right-file"].as<string>());


        cout << "Comparing " << left_file << " with " << right_file << endl;
    }
    else
    {
        cout << "Nope, fail" << endl;
    }
}
