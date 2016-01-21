#ifndef _WEATHER_
#define _WEATHER_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>

#include "cell.hpp"
#include "isoline.hpp"

typedef struct {
	std::array<unsigned int, 2> size;
	std::array<double, 2> spacing;
	std::array<double, 4> boundary;
	std::vector<std::vector<cell>> cells;
	std::vector<std::vector<node>> nodes;
} uniform_data;

typedef struct {
	unsigned int size;
	std::array<double, 4> boundary;
	std::vector<std::array<double, 3>> values;
} scattered_data;

class weather_data {
public:
	weather_data();
	void load_data(std::string input_file);
	void create_grid(std::array<unsigned int, 2> size);
	weather_data& shepard_interpolation(double p);
	weather_data& hardy_interpolation(double R);
	weather_data& compute_colormap(std::vector<std::array<int, 4>> color_table);
	weather_data& compute_isolines(unsigned int n);
	void export_colormap();
	void export_isolines();
	~weather_data();
private:
	std::string name;
	uniform_data uniform;
	scattered_data scattered;
	std::vector<unsigned char> colormap;
	std::vector<isoline> isolines;
};

#endif