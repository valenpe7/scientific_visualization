#ifndef _ISOLINE_
#define _ISOLINE_

#include <iostream>
#include <vector>

#include "segment.hpp"

class isoline {
	friend class weather_data;
public:
	isoline() = default;
	~isoline() = default;
private:
	double value;
	std::vector<segment> segments;
};

#endif