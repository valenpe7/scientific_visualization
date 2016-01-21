#ifndef _CELL_
#define _CELL_

#include <iostream>
#include <array>

#include "node.hpp"
#include "segment.hpp"

class cell {
	friend class weather_data;
public:
	cell() = default;
	segment top_bottom();
	segment top_left();
	segment top_right();
	segment bottom_left();
	segment bottom_right();
	segment left_right();
	cell& compute_average();
	~cell() = default;
private:
	int id;
	double average;
	std::array<int, 4> RGBA;
	std::array<node*, 4> nodes;
};

#endif
