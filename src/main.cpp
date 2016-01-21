#include <iostream>
#include <array>

#include "include/weather_data.hpp"
#include "include/color_tables.hpp"

void run_visualization(std::string input_file) {
	weather_data quantity;
	quantity.load_data(input_file);
	quantity.create_grid({1025, 1025});
	quantity.shepard_interpolation(2.0);
	//quantity.hardy_interpolation(2.0);
	quantity.compute_colormap(jet);
	quantity.compute_isolines(6);
	quantity.export_colormap();
	quantity.export_isolines();
	return;
}

int main(int argc, char** argv) {
	if (argc != 2)
		std::cout << "Usage: ./application_name.exe scattered_data_filename.csv" << std::endl;
	else {
		run_visualization(argv[1]);
	}
	return 0;
}